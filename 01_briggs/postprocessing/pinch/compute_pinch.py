
import json
import os
import sys
import numpy as np
import mpmath as mp
import couette_eigenvaluesspectrum as spatial
from couette_eigenvaluesspectrum import alpha_spectrum


# ---------------------------------------------------------------------------
# PRECISION
#
# Every calculation in this script is done with mpmath at HIGH_PREC_DPS
# decimal digits.  There is no double-precision arithmetic left anywhere in
# the pinch calculation: the temporal eigenvalue problem, the circle sweep,
# the polynomial fit, the root finding and the final Newton refinement all
# run at the same working precision.
#
# numpy is still imported, but only as a container for values that are on
# their way into matplotlib or into the JSON file.  Those conversions are
# marked where they happen and never feed back into a calculation.
# ---------------------------------------------------------------------------
HIGH_PREC_DPS = 60
mp.mp.dps = HIGH_PREC_DPS
spatial.set_precision(HIGH_PREC_DPS)

I_UNIT = mp.mpc(0, 1)          # exact imaginary unit, replaces the literal 1j
BC_SHIFT = mp.mpc(0, -200)     # replaces the literal -200j of the boundary rows


# ---------------------------------------------------------------------------
# SETTINGS -- must match the Briggs run
# ---------------------------------------------------------------------------
RE = 2000.0
BETA = 0.0
NUM_MODES = 100

RADIUS = 0.1
N_POINTS = 24
DEGREE = 6
LEVELS = 3

# A single arbitrary-precision eigenvalue solve at NUM_MODES = 100 takes about
# a minute, and the circle sweep does N_POINTS * LEVELS of them, so a stock run
# is measured in hours.  PROGRESS writes a counter to stderr during the sweep
# so the run is not silent; stdout is unaffected.
PROGRESS = True


# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------
def to_np(values):
    """
    Convert a sequence of mpmath numbers into a numpy complex array.

    Used only at the matplotlib boundary.  matplotlib is a double-precision
    library, so the high-precision values have to be rounded before they can
    be drawn; nothing that comes out of here is ever used in a calculation.
    """
    return np.array([complex(z) for z in values], dtype=complex)


def cx(values):
    """
    Read a list of {"re": ..., "im": ...} entries from the Briggs JSON.

    json.load gives back Python floats, and mp.mpf(float) is exact, so this
    carries the stored values across without adding or losing a single bit.
    """
    return [mp.mpc(mp.mpf(z["re"]), mp.mpf(z["im"])) for z in values]


# ---------------------------------------------------------------------------
# TEMPORAL SOLVER: given alpha, calculate omega eigenvalues
#
# One arbitrary-precision solver, used everywhere.  It builds the Chebyshev
# operators once, assembles A(alpha) and B(alpha) on demand, and provides
# both the eigenvalue spectrum and the determinant form of the dispersion
# relation that the final Newton refinement needs.
# ---------------------------------------------------------------------------
class Temporal:
    def __init__(self, num_modes=NUM_MODES, Re=RE, beta=BETA):
        mp.mp.dps = HIGH_PREC_DPS

        self.Re = mp.mpf(str(Re))
        self.beta = mp.mpc(str(beta))
        self.N = num_modes
        n = num_modes

        yc = [
            mp.cos(mp.mpf(k) * mp.pi / (n - 1))
            for k in range(n)
        ]

        y = [
            mp.mpf("0.5") - mp.mpf("0.5") * yy
            for yy in yc
        ]

        D0 = mp.matrix(n, n)
        D1 = mp.matrix(n, n)
        D2 = mp.matrix(n, n)
        D3 = mp.matrix(n, n)
        D4 = mp.matrix(n, n)

        # Chebyshev basis
        for r in range(n):
            for j in range(n):
                D0[r, j] = mp.cos(j * mp.acos(yc[r]))

        # Initial derivative columns
        for r in range(n):
            D1[r, 1] = D0[r, 0]

            if n > 2:
                D1[r, 2] = 4 * D0[r, 1]
                D2[r, 2] = 4 * D0[r, 0]

        # Higher derivative columns
        for j in range(3, n):
            factor = mp.mpf(j) / (j - 2)

            for r in range(n):
                D1[r, j] = 2*j*D0[r, j-1] + factor*D1[r, j-2]
                D2[r, j] = 2*j*D1[r, j-1] + factor*D2[r, j-2]
                D3[r, j] = 2*j*D2[r, j-1] + factor*D3[r, j-2]
                D4[r, j] = 2*j*D3[r, j-1] + factor*D4[r, j-2]

        # Map [-1,1] -> [0,1]
        for r in range(n):
            for j in range(n):
                D1[r, j] /= mp.mpf("-0.5")
                D2[r, j] /= mp.mpf("0.25")
                D4[r, j] /= mp.mpf("0.0625")

        self.D0 = D0
        self.D1 = D1
        self.D2 = D2
        self.D4 = D4
        self.U = y

    def matrices(self, alpha):
        """
        Construct the generalized eigenvalue matrices A(alpha), B(alpha).
        """

        a = mp.mpc(alpha)
        beta = self.beta
        Re = self.Re

        n = self.N
        k2 = a*a + beta*beta

        A = mp.matrix(n, n)
        B = mp.matrix(n, n)

        # Interior Orr-Sommerfeld operator
        for r in range(n):
            U = self.U[r]

            for j in range(n):

                A[r, j] = (
                    -I_UNIT*a*U*self.D2[r, j]
                    + I_UNIT*a*U*k2*self.D0[r, j]
                    + self.D4[r, j]/Re
                    - 2*k2*self.D2[r, j]/Re
                    + k2*k2*self.D0[r, j]/Re
                )

                B[r, j] = (
                    -I_UNIT*self.D2[r, j]
                    + I_UNIT*k2*self.D0[r, j]
                )

        # ---------------------------------------------------------------
        # Boundary conditions
        #
        # Rows 0, 1 carry D0 and D1 at the first wall, rows n-2, n-1 carry
        # D1 and D0 at the last wall.  The -200j shift parks the four
        # constraint modes at omega = -200j instead of at infinity.
        # ---------------------------------------------------------------

        for j in range(n):

            # first boundary
            A[0, j] = BC_SHIFT*self.D0[0, j]
            B[0, j] = self.D0[0, j]

            A[1, j] = BC_SHIFT*self.D1[0, j]
            B[1, j] = self.D1[0, j]

            # last boundary
            A[n-2, j] = BC_SHIFT*self.D1[n-1, j]
            B[n-2, j] = self.D1[n-1, j]

            A[n-1, j] = BC_SHIFT*self.D0[n-1, j]
            B[n-1, j] = self.D0[n-1, j]

        return A, B

    def spectrum(self, alpha):
        """
        Eigenvalues of A(alpha) x = omega B(alpha) x.

        B is non-singular for this problem (all n eigenvalues are finite,
        the four boundary modes sitting at omega = -200j), so the pencil is
        reduced to the standard problem B^-1 A and solved with mp.eig.
        """

        A, B = self.matrices(alpha)

        C = B**-1 * A

        w = mp.eig(C, left=False, right=False)

        return [z for z in w if mp.isfinite(z)]

    def nearest(self, alpha, omega_ref):
        """Follow one omega branch by choosing the eigenvalue nearest omega_ref."""
        w = self.spectrum(alpha)
        return min(w, key=lambda z: abs(z - omega_ref))

    def dispersion(self, alpha, omega):
        """
        Numerical dispersion relation:

              D(alpha, omega) = det(A(alpha) - omega B(alpha))

        An eigenvalue satisfies D = 0.
        """

        A, B = self.matrices(alpha)

        M = A - mp.mpc(omega)*B

        return mp.det(M)


# ---------------------------------------------------------------------------
# PINCH REFINEMENT
# ---------------------------------------------------------------------------
def refine_pinch_high_precision(alpha_guess, omega_guess, solver=None):
        """
        Refine the circle-fit result by solving

            D(alpha, omega) = 0
            dD/dalpha = 0

        The circle fit is already at full working precision, so this is a
        change of equation rather than a change of precision: it replaces
        the polynomial model of omega(alpha) by the determinant itself.
        """

        mp.mp.dps = HIGH_PREC_DPS

        print("\nPINCH REFINEMENT")
        print("--------------------------------------------")
        print(f"working precision = {HIGH_PREC_DPS} decimal digits")

        hp = Temporal() if solver is None else solver

        alpha0 = mp.mpc(alpha_guess)
        omega0 = mp.mpc(omega_guess)

        def D(alpha, omega):
            return hp.dispersion(alpha, omega)

        def D_alpha(alpha, omega):
            return mp.diff(
                lambda aa: D(aa, omega),
                alpha
            )

        # Two complex equations, two complex unknowns
        alpha_hp, omega_hp = mp.findroot(
            (D, D_alpha),
            (alpha0, omega0),
            tol=mp.mpf("1e-50"),
            maxsteps=30
        )

        residual_D = abs(D(alpha_hp, omega_hp))
        residual_Dalpha = abs(D_alpha(alpha_hp, omega_hp))

        print("\nRefined result:")
        print("alpha_p =")
        print(mp.nstr(alpha_hp, 55))

        print("\nomega_p =")
        print(mp.nstr(omega_hp, 55))

        print("\nChecks:")
        print("|D|       =", mp.nstr(residual_D, 10))
        print("|D_alpha| =", mp.nstr(residual_Dalpha, 10))

        return alpha_hp, omega_hp, residual_D, residual_Dalpha


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
    i = max(range(len(omega_F)), key=lambda k: omega_F[k].imag)

    # nearest L point in the omega-plane
    j = min(range(len(L)), key=lambda k: abs(L[k] - omega_F[i]))

    # midpoint of upper and lower alpha branches
    branch_mid = mp.mpf("0.5") * (alpha_u[j] + alpha_l[j])

    # seed = midpoint between F point and branch midpoint
    seed = mp.mpf("0.5") * (F[i] + branch_mid)

    return seed, omega_F[i], i, j


# ---------------------------------------------------------------------------
# 2. CIRCLE
# ---------------------------------------------------------------------------
def pinch_from_circle(solver, centre, omega_ref, radius):
    # N_POINTS equally spaced alpha points on the circle.
    # expjpi(2k/N) is exp(2*pi*i*k/N) evaluated without ever rounding pi.
    alpha = [
        centre + radius*mp.expjpi(2*mp.mpf(k)/N_POINTS)
        for k in range(N_POINTS)
    ]

    # calculate omega(alpha) while staying on one branch
    omega = []
    ref = omega_ref

    for k, a in enumerate(alpha):
        ref = solver.nearest(a, ref)
        omega.append(ref)

        if PROGRESS:
            print(
                f"    alpha sample {k+1:4d} / {N_POINTS}",
                end="\r",
                file=sys.stderr,
                flush=True,
            )

    if PROGRESS:
        print(" " * 40, end="\r", file=sys.stderr, flush=True)

    # normalized coordinate: circle becomes |s| = 1
    s = [(a - centre) / radius for a in alpha]

    # fit omega(s) with a degree-DEGREE polynomial
    V = mp.matrix(N_POINTS, DEGREE + 1)
    for r in range(N_POINTS):
        power = mp.mpc(1)
        for k in range(DEGREE + 1):
            V[r, k] = power
            power *= s[r]

    # least squares in arbitrary precision (QR); the numpy lstsq it replaces
    # solved the same full-rank overdetermined system
    c_solution, _ = mp.qr_solve(V, mp.matrix(omega))
    c = [c_solution[k] for k in range(DEGREE + 1)]

    fit_error = max(
        abs(sum((V[r, k]*c[k] for k in range(DEGREE + 1)), mp.mpc(0)) - omega[r])
        for r in range(N_POINTS)
    )

    dc = [k*c[k] for k in range(1, DEGREE + 1)]  # derivative of fitted polynomial

    # solve d(omega)/ds = 0 -> also d(omega)/d(alpha) = 0
    dc_descending = list(reversed(dc))

    # mirror numpy.roots, which drops leading zero coefficients
    while len(dc_descending) > 1 and dc_descending[0] == 0:
        dc_descending.pop(0)

    roots = mp.polyroots(
        dc_descending,
        maxsteps=200,
        extraprec=2*mp.mp.prec
    )

    roots = [r for r in roots if abs(r) < mp.mpf("1.2")]

    if len(roots) == 0:
        raise RuntimeError(
            "No pinch candidate found near the circle. Try a larger RADIUS."
        )

    # choose stationary point closest to circle centre
    s_p = min(roots, key=abs)

    alpha_p = centre + radius*s_p
    omega_p = sum((c[k]*s_p**k for k in range(len(c))), mp.mpc(0))

    # Return the alpha circle points and their corresponding omega values.
    return alpha_p, omega_p, fit_error, alpha, omega


# ---------------------------------------------------------------------------
# 3. REPEAT WITH SMALLER RADII
# ---------------------------------------------------------------------------
def find_pinch(solver, seed, omega_seed):
    centre = seed
    omega_ref = omega_seed
    radius = mp.mpf(str(RADIUS))

    alpha_estimates = []
    level_results = []

    print("\nlevel   radius        alpha_p                         fit error      change")
    print("----------------------------------------------------------------------------")

    previous = None

    for level in range(1, LEVELS + 1):
        alpha_p, omega_p, fit_error, circle_points, omega_points = pinch_from_circle(
            solver, centre, omega_ref, radius
        )

        change = mp.nan if previous is None else abs(alpha_p - previous)
        change_text = "---" if previous is None else f"{change:.3e}"

        print(
            f"{level:3d}   {radius:8.5f}   "
            f"{alpha_p.real:+.9f} {alpha_p.imag:+.9f}i   "
            f"{fit_error:.3e}   {change_text}"
        )

        level_results.append({
            "level": level,
            "centre": centre,
            "radius": radius,
            "circle_points": list(circle_points),
            "omega_points": list(omega_points),
            "alpha_p": alpha_p,
            "omega_p": omega_p,
            "fit_error": fit_error,
            "change_from_previous": None if previous is None else change,
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
#
# Every numeric field the original file had is still written, unchanged and
# still as a JSON number, so anything that already reads these files keeps
# working.  Alongside each one there is now a "_str" field holding the same
# quantity as a 50-digit decimal string, which is the only way to get the
# full working precision through JSON.
# ---------------------------------------------------------------------------
def write_output_json(output_path, source_path, last, seed, alpha_p, omega_p,
                      error, level_results,
                      alpha_hp, omega_hp,
                      residual_D, residual_Dalpha):

    def complex_entry(z):
        return {
            "re": float(mp.re(z)),
            "im": float(mp.im(z)),
            "re_str": mp.nstr(mp.re(z), 50),
            "im_str": mp.nstr(mp.im(z), 50),
        }

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
            "working_precision_dps": HIGH_PREC_DPS,
        },
        "seed": complex_entry(seed),
        "alpha_p_high_precision": {
            "re": mp.nstr(mp.re(alpha_hp), 50),
            "im": mp.nstr(mp.im(alpha_hp), 50),
        },
        "omega_p_high_precision": {
            "re": mp.nstr(mp.re(omega_hp), 50),
            "im": mp.nstr(mp.im(omega_hp), 50),
        },
        "high_precision_checks": {
            "D_residual": mp.nstr(residual_D, 20),
            "D_alpha_residual": mp.nstr(residual_Dalpha, 20),
        },
        "radius_check_error": float(error),
        "radius_check_error_str": mp.nstr(error, 50),
        "level_results": [],
    }

    for result in level_results:
        change = result["change_from_previous"]

        out["level_results"].append({
            "level": result["level"],
            "centre": complex_entry(result["centre"]),
            "radius": float(result["radius"]),
            "radius_str": mp.nstr(result["radius"], 50),
            "circle_points": [
                complex_entry(z)
                for z in result["circle_points"]
            ],
            "alpha_p": complex_entry(result["alpha_p"]),
            "omega_p": complex_entry(result["omega_p"]),
            "fit_error": float(result["fit_error"]),
            "fit_error_str": mp.nstr(result["fit_error"], 50),
            "change_from_previous": None if change is None else float(change),
            "change_from_previous_str": None if change is None else mp.nstr(change, 50),
        })

    with open(output_path, "w") as f:
        json.dump(out, f, indent=2)


# ---------------------------------------------------------------------------
# 5. MAKE PLOT
#
# From here on the values are handed to matplotlib, which is a
# double-precision library.  to_np() and complex() round the high-precision
# results once, at the point where they are drawn; the axis limits computed
# from them are plot geometry, not physics.
# ---------------------------------------------------------------------------
def make_plot(output_path, last, seed, alpha_p, omega_p, level_results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    L = to_np(cx(last["L"]))
    omega_F = to_np(cx(last["omega_F"]))
    F = to_np(cx(last["F"]))
    alpha_u = to_np(cx(last["alpha_L_u"]))
    alpha_l = to_np(cx(last["alpha_L_l"]))

    seed = complex(seed)
    alpha_p = complex(alpha_p)
    omega_p = complex(omega_p)

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
        centre = complex(result["centre"])
        radius = float(result["radius"])
        points = to_np(result["circle_points"])
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

    L = to_np(cx(last["L"]))
    omega_F = to_np(cx(last["omega_F"]))
    F = to_np(cx(last["F"]))
    alpha_u = to_np(cx(last["alpha_L_u"]))
    alpha_l = to_np(cx(last["alpha_L_l"]))

    seed = complex(seed)
    alpha_p = complex(alpha_p)
    omega_p = complex(omega_p)

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
        radius = float(result["radius"])
        omega_points = to_np(result["omega_points"])

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
        centre = complex(result["centre"])
        radius = float(result["radius"])
        points = to_np(result["circle_points"])
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
    # Calculate complete alpha spectrum at omega_p
    #
    # Done before anything is rounded: alpha_spectrum() takes the
    # high-precision omega_p and returns high-precision alpha.
    # -------------------------------------------------
    alpha_spec_hp = alpha_spectrum(omega_p)

    print(
        f"calculated {len(alpha_spec_hp)} finite alpha eigenvalues "
        f"at omega_p = {omega_p}" , alpha_spec_hp
    )

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

    spec_max_abs_hp = mp.mpf(str(SPECTRUM_MAX_ABS))

    spec_kept_hp = [z for z in alpha_spec_hp if abs(z) <= spec_max_abs_hp]

    largest_abs_hp = max(abs(z) for z in alpha_spec_hp)

    print(
        f"{len(spec_kept_hp)} of {len(alpha_spec_hp)} alpha eigenvalues "
        f"have |alpha| <= {SPECTRUM_MAX_ABS:g}  "
        f"(largest |alpha| in the spectrum: {largest_abs_hp:.3e})"
    )

    # -------------------------------------------------
    # Everything below here is drawing: values are rounded to double once,
    # here, and the window limits are plot geometry.
    # -------------------------------------------------
    F = to_np(cx(last["F"]))
    alpha_u = to_np(cx(last["alpha_L_u"]))
    alpha_l = to_np(cx(last["alpha_L_l"]))

    alpha_p = complex(alpha_p)
    omega_p_plot = complex(omega_p)

    alpha_spec = to_np(alpha_spec_hp)
    spec_kept = to_np(spec_kept_hp)

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
        rf"      $\omega_p$ = {omega_p_plot.real:+.7f} {omega_p_plot.imag:+.7f}i",
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
    # ---------------------------------------------------------
    # FINAL REFINEMENT
    # ---------------------------------------------------------
    alpha_hp, omega_hp, residual_D, residual_Dalpha = (
        refine_pinch_high_precision(alpha_p, omega_p, solver)
    )

    print("\nFINAL")
    print(f"alpha_p = {alpha_p.real:+.9f} {alpha_p.imag:+.9f}i")
    print(f"omega_p = {omega_p.real:+.9f} {omega_p.imag:+.9f}i")
    print(f"radius-check error = {error:.3e}")

    # write JSON
    json_path = stem + "_pinch_simple.json"
    write_output_json(
        json_path,
        path,
        last,
        seed,
        alpha_p,
        omega_p,
        error,
        level_results,
        alpha_hp,
        omega_hp,
        residual_D,
        residual_Dalpha
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
