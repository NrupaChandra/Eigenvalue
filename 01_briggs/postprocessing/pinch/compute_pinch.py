
import json
import os
import sys
import numpy as np
from scipy.linalg import eig
from couette_eigenvaluesspectrum import alpha_spectrum

# ---------------------------------------------------------------------------
# SETTINGS -- must match the Briggs run
# ---------------------------------------------------------------------------
RE = 2000.0
BETA = 0.0
NUM_MODES = 100

RADIUS = 0.2
N_POINTS = 24
DEGREE = 6
LEVELS = 3


# ---------------------------------------------------------------------------
# TEMPORAL SOLVER: given alpha, calculate omega eigenvalues
# ---------------------------------------------------------------------------
class Temporal:
    def __init__(self, num_modes=NUM_MODES, Re=RE, beta=BETA):
        self.Re, self.beta, self.N = Re, complex(beta), num_modes
        n = num_modes

        yc = np.cos(np.arange(n) * np.pi / (n - 1))
        y = 0.5 - 0.5 * yc

        D0 = np.zeros((n, n))
        for j in range(n):
            D0[:, j] = np.cos(j * np.arccos(yc))

        D1 = np.zeros((n, n))
        D2 = np.zeros((n, n))
        D3 = np.zeros((n, n))
        D4 = np.zeros((n, n))

        D1[:, 1] = D0[:, 0]
        D1[:, 2] = 4 * D0[:, 1]
        D2[:, 2] = 4 * D0[:, 0]

        for j in range(3, n):
            D1[:, j] = 2*j*D0[:, j-1] + j*D1[:, j-2]/(j-2)
            D2[:, j] = 2*j*D1[:, j-1] + j*D2[:, j-2]/(j-2)
            D3[:, j] = 2*j*D2[:, j-1] + j*D3[:, j-2]/(j-2)
            D4[:, j] = 2*j*D3[:, j-1] + j*D4[:, j-2]/(j-2)

        self.D0 = D0
        self.D1 = D1 / (-0.5)
        self.D2 = D2 / 0.25
        self.D4 = D4 / 0.0625
        self.U = y[:, None]

    def spectrum(self, alpha):
        a = complex(alpha)
        k2 = a*a + self.beta*self.beta
        n, Re = self.N, self.Re
        D0, D1, D2, D4 = self.D0, self.D1, self.D2, self.D4

        A = (-1j*a*self.U*D2 + 1j*a*self.U*k2*D0
             + D4/Re - 2/Re*k2*D2 + k2*k2/Re*D0)
        B = -1j*D2 + 1j*k2*D0

        top = np.vstack([D0[0:1, :], D1[0:1, :]])
        bot = np.vstack([D1[n-1:n, :], D0[n-1:n, :]])

        A = np.vstack([-200j*top, A[2:n-2, :], -200j*bot])
        B = np.vstack([top,       B[2:n-2, :],       bot])

        w = eig(A, B, right=False)
        return w[np.isfinite(w.real) & np.isfinite(w.imag)]

    def nearest(self, alpha, omega_ref):
        """Follow one omega branch by choosing the eigenvalue nearest omega_ref."""
        w = self.spectrum(alpha)
        return w[np.argmin(np.abs(w - omega_ref))]


def cx(values):
    return np.array([complex(z["re"], z["im"]) for z in values])


# ---------------------------------------------------------------------------
# 1. CALCULATE INITIAL SEED
# ---------------------------------------------------------------------------
def calculate_seed(last):
    L = cx(last["L"])
    omega_F = cx(last["omega_F"])
    F = cx(last["F"])
    alpha_u = cx(last["alpha_L_u"])
    alpha_l = cx(last["alpha_L_l"])

    # F point whose omega has the largest imaginary part
    i = int(np.argmax(omega_F.imag))

    # nearest L point in the omega-plane
    j = int(np.argmin(np.abs(L - omega_F[i])))

    # midpoint of upper and lower alpha branches
    branch_mid = 0.5 * (alpha_u[j] + alpha_l[j])

    # seed = midpoint between F point and branch midpoint
    seed = 0.5 * (F[i] + branch_mid)

    return seed, omega_F[i], i, j


# ---------------------------------------------------------------------------
# 2. CIRCLE 
# ---------------------------------------------------------------------------
def pinch_from_circle(solver, centre, omega_ref, radius):
    # 24 equally spaced alpha points on the circle
    theta = 2*np.pi*np.arange(N_POINTS) / N_POINTS
    alpha = centre + radius*np.exp(1j*theta)

    # calculate omega(alpha) while staying on one branch
    omega = np.empty(N_POINTS, complex)
    ref = omega_ref

    for k, a in enumerate(alpha):
        ref = solver.nearest(a, ref)
        omega[k] = ref

    # normalized coordinate: circle becomes |s| = 1
    s = (alpha - centre) / radius

    # fit omega(s) with a degree-DEGREE polynomial
    V = np.vander(s, DEGREE + 1, increasing=True)
    c, *_ = np.linalg.lstsq(V, omega, rcond=None)

    fit_error = np.max(np.abs(V @ c - omega))
    dc = np.array([k*c[k] for k in range(1, DEGREE + 1)])#  derivative of fitted polynomial

    # solve d(omega)/ds = 0 -> also d(omega)/d(alpha) = 0
    roots = np.roots(dc[::-1])
    roots = roots[np.abs(roots) < 1.2]

    if len(roots) == 0:
        raise RuntimeError(
            "No pinch candidate found near the circle. Try a larger RADIUS."
        )

    # choose stationary point closest to circle centre
    s_p = roots[np.argmin(np.abs(roots))]

    alpha_p = centre + radius*s_p
    omega_p = sum(c[k]*s_p**k for k in range(len(c)))

    # Return the alpha circle points and their corresponding omega values.
    return alpha_p, omega_p, fit_error, alpha, omega


# ---------------------------------------------------------------------------
# 3. REPEAT WITH SMALLER RADII
# ---------------------------------------------------------------------------
def find_pinch(solver, seed, omega_seed):
    centre = seed
    omega_ref = omega_seed
    radius = RADIUS

    alpha_estimates = []
    level_results = []

    print("\nlevel   radius        alpha_p                         fit error      change")
    print("----------------------------------------------------------------------------")

    previous = None

    for level in range(1, LEVELS + 1):
        alpha_p, omega_p, fit_error, circle_points, omega_points = pinch_from_circle(
            solver, centre, omega_ref, radius
        )

        change = np.nan if previous is None else abs(alpha_p - previous)
        change_text = "---" if previous is None else f"{change:.3e}"

        print(
            f"{level:3d}   {radius:8.5f}   "
            f"{alpha_p.real:+.9f} {alpha_p.imag:+.9f}i   "
            f"{fit_error:.3e}   {change_text}"
        )

        level_results.append({
            "level": level,
            "centre": centre,
            "radius": float(radius),
            "circle_points": circle_points.copy(),
            "omega_points": omega_points.copy(),
            "alpha_p": alpha_p,
            "omega_p": omega_p,
            "fit_error": float(fit_error),
            "change_from_previous": None if previous is None else float(change),
        })

        alpha_estimates.append(alpha_p)

        # next search: centre on current pinch estimate and halve the radius
        centre = alpha_p
        omega_ref = omega_p
        previous = alpha_p
        radius /= 2

    # maximum disagreement between the estimates from all radii
    error = max(abs(a-b) for a in alpha_estimates for b in alpha_estimates)

    return alpha_p, omega_p, error, level_results


# ---------------------------------------------------------------------------
# 4. WRITE JSON OUTPUT
# ---------------------------------------------------------------------------
def write_output_json(output_path, source_path, last, seed, alpha_p, omega_p,
                      error, level_results):
    out = {
        "source_file": os.path.basename(source_path),
        "iteration": int(last["iteration"]),
        "settings": {
            "Re": RE,
            "beta": BETA,
            "num_modes": NUM_MODES,
            "initial_radius": RADIUS,
            "n_points": N_POINTS,
            "degree": DEGREE,
            "levels": LEVELS,
        },
        "seed": {
            "re": float(seed.real),
            "im": float(seed.imag),
        },
        "alpha_p": {
            "re": float(alpha_p.real),
            "im": float(alpha_p.imag),
        },
        "omega_p": {
            "re": float(omega_p.real),
            "im": float(omega_p.imag),
        },
        "radius_check_error": float(error),
        "level_results": [],
    }

    for result in level_results:
        out["level_results"].append({
            "level": result["level"],
            "centre": {
                "re": float(result["centre"].real),
                "im": float(result["centre"].imag),
            },
            "radius": result["radius"],
            "circle_points": [
                {
                    "re": float(z.real),
                    "im": float(z.imag)
                }
                for z in result["circle_points"]
            ],
            "alpha_p": {
                "re": float(result["alpha_p"].real),
                "im": float(result["alpha_p"].imag),
            },
            "omega_p": {
                "re": float(result["omega_p"].real),
                "im": float(result["omega_p"].imag),
            },
            "fit_error": result["fit_error"],
            "change_from_previous": result["change_from_previous"],
        })

    with open(output_path, "w") as f:
        json.dump(out, f, indent=2)


# ---------------------------------------------------------------------------
# 5. MAKE PLOT
# ---------------------------------------------------------------------------
def make_plot(output_path, last, seed, alpha_p, omega_p, level_results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    L = cx(last["L"])
    omega_F = cx(last["omega_F"])
    F = cx(last["F"])
    alpha_u = cx(last["alpha_L_u"])
    alpha_l = cx(last["alpha_L_l"])

    fig, ax = plt.subplots(1, 2, figsize=(16, 7.5))

    # omega-plane
    ax[0].plot(L.real, L.imag, "b", lw=1.6, label="L")
    ax[0].plot(omega_F.real, omega_F.imag, "r", lw=1.6, label=r"$\omega_F$")
    ax[0].plot(
        omega_p.real, omega_p.imag,"*", ms=5, label=r"$\omega_p$")
    ax[0].set_xlabel(r"$\omega_r$")
    ax[0].set_ylabel(r"$\omega_i$")
    ax[0].set_title(r"$\omega$-plane")
    ax[0].grid(True, alpha=0.4)
    ax[0].legend(loc="best")

    # alpha-plane
    ax[1].plot(F.real, F.imag, "b", lw=1.6, label="F")
    ax[1].plot(alpha_u.real, alpha_u.imag, "r", lw=1.6, label=r"$\alpha_L^u$")
    ax[1].plot(alpha_l.real, alpha_l.imag, "g", lw=1.6, label=r"$\alpha_L^l$")

    # plot every circle used in the refinement
    for result in level_results:
        centre = result["centre"]
        radius = result["radius"]
        points = result["circle_points"]
        level = result["level"]

        # Plot points first and use the same automatically chosen colour
        # for the corresponding circle.
        point_line, = ax[1].plot(
            points.real,
            points.imag,
            ".",
            ms=6,
            label=f"circle {level}: R={radius:g}",
        )

        circle = plt.Circle(
            (centre.real, centre.imag),
            radius,
            fill=False,
            lw=1.6,
            ls="--",
            edgecolor=point_line.get_color(),
        )
        ax[1].add_patch(circle)

    # seed and final pinch
    ax[1].plot(seed.real, seed.imag, "o", ms=5, label="seed")
    ax[1].plot(alpha_p.real, alpha_p.imag, "*", ms=5, label=r"$\alpha_p$")

    ax[1].set_xlabel(r"$\alpha_r$")
    ax[1].set_ylabel(r"$\alpha_i$")
    ax[1].set_title(r"$\alpha$-plane")
    ax[1].grid(True, alpha=0.4)
    ax[1].legend(loc="best")

    fig.suptitle(
        rf"$\alpha_p$ = {alpha_p.real:+.7f} {alpha_p.imag:+.7f}i"
        rf"      $\omega_p$ = {omega_p.real:+.7f} {omega_p.imag:+.7f}i",
        fontsize=12,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output_path, dpi=160)
    plt.close(fig)



# ---------------------------------------------------------------------------
# 6. MAKE SECOND PLOT 
# ---------------------------------------------------------------------------
def make_plot_with_omega_samples(output_path, last, seed, alpha_p, omega_p,
                                 level_results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    L = cx(last["L"])
    omega_F = cx(last["omega_F"])
    F = cx(last["F"])
    alpha_u = cx(last["alpha_L_u"])
    alpha_l = cx(last["alpha_L_l"])

    fig, ax = plt.subplots(1, 2, figsize=(16, 7.5))

    # omega-plane: same as the first plot
    ax[0].plot(L.real, L.imag, "b", lw=1.6, label="L")
    ax[0].plot(omega_F.real, omega_F.imag, "r", lw=1.6, label=r"$\omega_F$")
    ax[0].plot(
        omega_p.real, omega_p.imag, "*", ms=5, label=r"$\omega_p$"
    )

    # Add the omega values produced by the alpha sample points.
    # Each level/circle gets one colour.
    level_colours = {}

    for result in level_results:
        level = result["level"]
        radius = result["radius"]
        omega_points = result["omega_points"]

        omega_line, = ax[0].plot(
            omega_points.real,
            omega_points.imag,
            ".",
            ms=6,
            label=rf"$\omega(\alpha_k)$ circle {level}: R={radius:g}",
        )
        level_colours[level] = omega_line.get_color()

    ax[0].set_xlabel(r"$\omega_r$")
    ax[0].set_ylabel(r"$\omega_i$")
    ax[0].set_title(r"$\omega$-plane with circle-sample $\omega$ values")
    ax[0].grid(True, alpha=0.4)
    ax[0].legend(loc="best")

    # alpha-plane: same circles and sample points as the first plot
    ax[1].plot(F.real, F.imag, "b", lw=1.6, label="F")
    ax[1].plot(alpha_u.real, alpha_u.imag, "r", lw=1.6, label=r"$\alpha_L^u$")
    ax[1].plot(alpha_l.real, alpha_l.imag, "g", lw=1.6, label=r"$\alpha_L^l$")

    for result in level_results:
        centre = result["centre"]
        radius = result["radius"]
        points = result["circle_points"]
        level = result["level"]
        colour = level_colours[level]

        ax[1].plot(
            points.real,
            points.imag,
            ".",
            ms=6,
            color=colour,
            label=f"circle {level}: R={radius:g}",
        )

        circle = plt.Circle(
            (centre.real, centre.imag),
            radius,
            fill=False,
            lw=1.6,
            ls="--",
            edgecolor=colour,
        )
        ax[1].add_patch(circle)

    ax[1].plot(seed.real, seed.imag, "o", ms=5, label="seed")
    ax[1].plot(alpha_p.real, alpha_p.imag, "*", ms=5, label=r"$\alpha_p$")

    ax[1].set_xlabel(r"$\alpha_r$")
    ax[1].set_ylabel(r"$\alpha_i$")
    ax[1].set_title(r"$\alpha$-plane")
    ax[1].grid(True, alpha=0.4)
    ax[1].legend(loc="best")

    fig.suptitle(
        rf"$\alpha_p$ = {alpha_p.real:+.7f} {alpha_p.imag:+.7f}i"
        rf"      $\omega_p$ = {omega_p.real:+.7f} {omega_p.imag:+.7f}i",
        fontsize=12,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output_path, dpi=160)
    plt.close(fig)

# ---------------------------------------------------------------------------
# 7. MAKE THIRD PLOT: ALPHA SPECTRUM AT OMEGA_P
# ---------------------------------------------------------------------------
def make_plot_with_alpha_spectrum(output_path, last, alpha_p, omega_p):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # -------------------------------------------------
    # Briggs data from final iteration
    # -------------------------------------------------
    F = cx(last["F"])
    alpha_u = cx(last["alpha_L_u"])
    alpha_l = cx(last["alpha_L_l"])

    # -------------------------------------------------
    # Calculate complete alpha spectrum at omega_p
    # -------------------------------------------------
    alpha_spec = alpha_spectrum(omega_p)

    print(
        f"calculated {len(alpha_spec)} finite alpha eigenvalues "
        f"at omega_p = {omega_p}" , alpha_spec
    )

    # =================================================
    # THREE ALPHA-PLANE WINDOWS, WIDEST TO TIGHTEST
    # =================================================

    # -------------------------------------------------
    # 1. SPECTRUM window.
    #
    # The full spectrum has no natural outer edge: it runs out to
    # |alpha| ~ 1e8, which is the unresolved tail of the N = NUM_MODES
    # Chebyshev discretisation, not physics.  So the panel shows every
    # eigenvalue with |alpha| <= SPECTRUM_MAX_ABS and reports how many
    # were left out.  Raise or lower this one number to see more or less.
    # -------------------------------------------------
    SPECTRUM_MAX_ABS = 20.0

    spec_kept = alpha_spec[np.abs(alpha_spec) <= SPECTRUM_MAX_ABS]

    print(
        f"{len(spec_kept)} of {len(alpha_spec)} alpha eigenvalues "
        f"have |alpha| <= {SPECTRUM_MAX_ABS:g}  "
        f"(largest |alpha| in the spectrum: {np.max(np.abs(alpha_spec)):.3e})"
    )

    spec_view = np.concatenate([spec_kept, np.array([alpha_p])])

    spec_r_pad = 0.05 * (np.max(spec_view.real) - np.min(spec_view.real))
    spec_i_pad = 0.05 * (np.max(spec_view.imag) - np.min(spec_view.imag))

    spec_x = (np.min(spec_view.real) - spec_r_pad,
              np.max(spec_view.real) + spec_r_pad)
    spec_y = (np.min(spec_view.imag) - spec_i_pad,
              np.max(spec_view.imag) + spec_i_pad)

    # -------------------------------------------------
    # 2. BRIGGS window: the contours of the final iteration.
    #
    # 0.05 = matplotlib's default autoscale margin, so this window comes
    # out identical to the autoscaled axes of plots 1 and 2.
    # -------------------------------------------------
    alpha_view = np.concatenate([
        F,
        alpha_u,
        alpha_l,
        np.array([alpha_p])
    ])

    alpha_r_pad = 0.05 * (np.max(alpha_view.real) - np.min(alpha_view.real))
    alpha_i_pad = 0.05 * (np.max(alpha_view.imag) - np.min(alpha_view.imag))

    briggs_x = (np.min(alpha_view.real) - alpha_r_pad,
                np.max(alpha_view.real) + alpha_r_pad)
    briggs_y = (np.min(alpha_view.imag) - alpha_i_pad,
                np.max(alpha_view.imag) + alpha_i_pad)

    # -------------------------------------------------
    # 3. ZOOM window: centred on alpha_p, half-width ZOOM_HALF_WIDTH
    # -------------------------------------------------
    ZOOM_HALF_WIDTH = 0.05

    zoom_x = (alpha_p.real - ZOOM_HALF_WIDTH, alpha_p.real + ZOOM_HALF_WIDTH)
    zoom_y = (alpha_p.imag - ZOOM_HALF_WIDTH, alpha_p.imag + ZOOM_HALF_WIDTH)

    # =================================================
    # CREATE FIGURE
    # =================================================
    fig, ax = plt.subplots(1, 3, figsize=(24, 7.5))

    panels = [
        (0, spec_x, spec_y,
         rf"$\alpha$-plane spectrum: $|\alpha| \leq$ {SPECTRUM_MAX_ABS:g}"
         rf"  ({len(spec_kept)} of {len(alpha_spec)} modes)",
         (briggs_x, briggs_y), "Briggs window"),
        (1, briggs_x, briggs_y,
         r"$\alpha$-plane",
         (zoom_x, zoom_y), "zoom window"),
        (2, zoom_x, zoom_y,
         rf"$\alpha$-plane zoom: $\alpha_p \pm$ {ZOOM_HALF_WIDTH:g}",
         None, None),
    ]

    for k, xlim, ylim, title, box, box_label in panels:

        # points inside this window only
        visible = (
            (alpha_spec.real >= xlim[0])
            & (alpha_spec.real <= xlim[1])
            & (alpha_spec.imag >= ylim[0])
            & (alpha_spec.imag <= ylim[1])
        )

        # F contour
        ax[k].plot(
            F.real,
            F.imag,
            "b",
            lw=1.6,
            label="F"
        )

        # Upper branch
        ax[k].plot(
            alpha_u.real,
            alpha_u.imag,
            "r",
            lw=1.6,
            label=r"$\alpha_L^u$"
        )

        # Lower branch
        ax[k].plot(
            alpha_l.real,
            alpha_l.imag,
            "g",
            lw=1.6,
            label=r"$\alpha_L^l$"
        )

        # -------------------------------------------------
        # Alpha spectrum at omega_p
        # -------------------------------------------------
        ax[k].plot(
            alpha_spec[visible].real,
            alpha_spec[visible].imag,
            "k.",
            ms=7,
            label=r"$\alpha(\omega_p)$"
        )

        # -------------------------------------------------
        # Final pinch
        # -------------------------------------------------
        ax[k].plot(
            alpha_p.real,
            alpha_p.imag,
            "*",
            ms=5,
            label=r"$\alpha_p$"
        )

        # -------------------------------------------------
        # Outline of the window shown in the next panel
        # -------------------------------------------------
        if box is not None:
            (bx0, bx1), (by0, by1) = box
            ax[k].add_patch(
                plt.Rectangle(
                    (bx0, by0),
                    bx1 - bx0,
                    by1 - by0,
                    fill=False,
                    lw=1.0,
                    ls="--",
                    edgecolor="k",
                    label=box_label
                )
            )

        ax[k].set_xlim(xlim[0], xlim[1])
        ax[k].set_ylim(ylim[0], ylim[1])

        ax[k].set_xlabel(r"$\alpha_r$")
        ax[k].set_ylabel(r"$\alpha_i$")
        ax[k].set_title(title)

        ax[k].grid(True, alpha=0.4)
        ax[k].legend(loc="best")


    # =================================================
    # COMMON TITLE
    # =================================================

    fig.suptitle(
        rf"$\alpha_p$ = {alpha_p.real:+.7f} {alpha_p.imag:+.7f}i"
        rf"      $\omega_p$ = {omega_p.real:+.7f} {omega_p.imag:+.7f}i",
        fontsize=12
    )

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    fig.savefig(
        output_path,
        dpi=160
    )

    plt.close(fig)
# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        print("Usage: python find_pinch_simple.py <big_history.json>")
        sys.exit(1)

    path = sys.argv[1]
    stem = os.path.splitext(os.path.abspath(path))[0]

    print(f"reading {os.path.basename(path)} ...")

    with open(path) as f:
        frames = json.load(f)

    if not isinstance(frames, list):
        raise ValueError(
            "Use the BIG Briggs history JSON, not a *_pinch*.json output file."
        )

    last = frames[-1]
    seed, omega_seed, i, j = calculate_seed(last)

    print(f"using last frame: iteration {last['iteration']}")
    print(f"F node: {i+1}, L node: {j+1}")
    print(f"seed = {seed.real:+.9f} {seed.imag:+.9f}i")

    solver = Temporal()

    alpha_p, omega_p, error, level_results = find_pinch(
        solver, seed, omega_seed
    )

    print("\nFINAL")
    print(f"alpha_p = {alpha_p.real:+.9f} {alpha_p.imag:+.9f}i")
    print(f"omega_p = {omega_p.real:+.9f} {omega_p.imag:+.9f}i")
    print(f"radius-check error = {error:.3e}")

    # write JSON
    json_path = stem + "_pinch_simple.json"
    write_output_json(
        json_path, path, last, seed, alpha_p, omega_p,
        error, level_results
    )
    print(f"wrote {json_path}")

    # make plot
    png_path = stem + "_pinch_simple.png"
    make_plot(
        png_path, last, seed, alpha_p, omega_p, level_results
    )
    print(f"wrote {png_path}")


    # make second plot with omega values from all alpha-circle sample points
    omega_samples_png_path = stem + "_pinch_simple_omega_samples.png"
    make_plot_with_omega_samples(
        omega_samples_png_path, last, seed, alpha_p, omega_p, level_results
    )
    print(f"wrote {omega_samples_png_path}")

    # make third plot with complete alpha spectrum at omega_p
    alpha_spectrum_png_path = stem + "_pinch_simple_alpha_spectrum.png"

    make_plot_with_alpha_spectrum(
        alpha_spectrum_png_path,
        last,
        alpha_p,
        omega_p
    )

    print(f"wrote {alpha_spectrum_png_path}")

    print(f"wrote {alpha_spectrum_png_path}")

if __name__ == "__main__":
    main()