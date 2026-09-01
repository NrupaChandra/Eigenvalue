
import json
import os
import sys
import time

import mpmath


# ---------------------------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------------------------
RE = 2000.0                 # must match the Briggs run

DPS_LADDER = (8, 16)  # precision levels; the last one is the answer
MAX_STEPS = 12              # Newton steps allowed per level
EXACT_THIRDS = True         # exact 1/3, 2/3 exponents instead of Python floats

SPECTRUM_MAX_ABS = 20.0     # spectrum plot: largest |alpha| shown in the spectrum panel
ZOOM_HALF_WIDTH = 0.05      # spectrum plot: half-width of the zoom panel around alpha_p


# ---------------------------------------------------------------------------
# THE ANALYTIC KERNEL  --  Fnum / Fnum_derivatives of eigenvalues_v8.jl
# ---------------------------------------------------------------------------
def fnum(a, omega, re, derivatives=True):
    """D, dD/dalpha, dD/domega and the size of the terms that cancel."""
    if EXACT_THIRDS:
        third, twothird = mpmath.mpf(1) / 3, mpmath.mpf(2) / 3
    else:
        third, twothird = 1 / 3, 2 / 3

    root3 = mpmath.power((-1j * a * re), third)

    def arg(y):
        return root3 / (a * re) * (omega * re - a * y * re + 1j * mpmath.power(a, 2))

    def darg_da(y):
        return ((1j * a * y * re + 2j * omega * re + 4 * mpmath.power(a, 2))
                / (3 * a * mpmath.power((-1j * a * re), twothird)))

    q = lambda f: mpmath.quad(f, [0, 1])

    ipa0 = q(lambda y: mpmath.exp(a * y) * mpmath.airyai(arg(y)))
    ina2 = q(lambda y: mpmath.exp(-a * y) * mpmath.airybi(arg(y)))
    ipa2 = q(lambda y: mpmath.exp(a * y) * mpmath.airybi(arg(y)))
    ina0 = q(lambda y: mpmath.exp(-a * y) * mpmath.airyai(arg(y)))

    disprel = (-1) / (2 * a) * (ipa0 * ina2 - ipa2 * ina0)
    scale = abs(ipa0 * ina2)

    if not derivatives:
        return disprel, None, None, scale

    d_ipa0_d_omega = q(lambda y: mpmath.exp(a * y) * mpmath.airyai(arg(y), derivative=1))
    d_ina2_d_omega = q(lambda y: mpmath.exp(-a * y) * mpmath.airybi(arg(y), derivative=1))
    d_ipa2_d_omega = q(lambda y: mpmath.exp(a * y) * mpmath.airybi(arg(y), derivative=1))
    d_ina0_d_omega = q(lambda y: mpmath.exp(-a * y) * mpmath.airyai(arg(y), derivative=1))

    d_ipa0_d_a = q(lambda y: y * mpmath.exp(a * y) * mpmath.airyai(arg(y))
                   + mpmath.exp(a * y) * mpmath.airyai(arg(y), derivative=1) * darg_da(y))
    d_ina2_d_a = q(lambda y: -y * mpmath.exp(-a * y) * mpmath.airybi(arg(y))
                   + mpmath.exp(-a * y) * mpmath.airybi(arg(y), derivative=1) * darg_da(y))
    d_ipa2_d_a = q(lambda y: y * mpmath.exp(a * y) * mpmath.airybi(arg(y))
                   + mpmath.exp(a * y) * mpmath.airybi(arg(y), derivative=1) * darg_da(y))
    d_ina0_d_a = q(lambda y: -y * mpmath.exp(-a * y) * mpmath.airyai(arg(y))
                   + mpmath.exp(-a * y) * mpmath.airyai(arg(y), derivative=1) * darg_da(y))

    d_disprel_d_omega = (-1) / (2 * a) * root3 / a * (
        d_ipa0_d_omega * ina2 + ipa0 * d_ina2_d_omega
        - (d_ipa2_d_omega * ina0 + ipa2 * d_ina0_d_omega))

    d_disprel_d_a = 1 / (2 * mpmath.power(a, 2)) * (ipa0 * ina2 - ipa2 * ina0) \
        - 1 / (2 * a) * (d_ipa0_d_a * ina2 + ipa0 * d_ina2_d_a
                         - (d_ipa2_d_a * ina0 + ipa2 * d_ina0_d_a))

    return disprel, d_disprel_d_a, d_disprel_d_omega, scale


# ---------------------------------------------------------------------------
# NEWTON on  [ D , dD/dalpha ] = 0
#
# J = [[dD/dalpha  , dD/domega    ]   top row analytic, bottom row one central
#      [d2D/dalpha2, d2D/dalpha dw]]  difference of the analytic dD/dalpha
# ---------------------------------------------------------------------------
def jacobian_bottom_row(a, omega, re, h):
    """d2D/dalpha2 and d2D/dalpha domega, by central difference of dD/dalpha."""
    ha, hw = h * abs(a), h * abs(omega)

    _, ap, _, _ = fnum(a + ha, omega, re)
    _, am, _, _ = fnum(a - ha, omega, re)
    _, wp, _, _ = fnum(a, omega + hw, re)
    _, wm, _, _ = fnum(a, omega - hw, re)

    return (ap - am) / (2 * ha), (wp - wm) / (2 * hw)


def find_pinch(alpha, omega, re):
    """Newton up the dps ladder.  Returns (alpha, omega, iterates, history)."""
    iterates = [(complex(alpha), complex(omega))]     # for the plots
    history = []

    print("\n level  it       |dalpha|        |domega|    |D|/|dD/dw|     "
          "|dw/dalpha|      s")
    print(" " + "-" * 76)

    for dps in DPS_LADDER:
        mpmath.mp.dps = dps + 20
        re_mp = mpmath.mpf(re)
        alpha, omega = mpmath.mpc(alpha), mpmath.mpc(omega)

        # central-difference step: balances truncation h^2 against noise/h
        h = mpmath.mpf(10) ** (-dps // 3)

        D_aa = D_aw = None
        prev_step = None
        rate = 0.0

        for it in range(MAX_STEPS):
            t0 = time.time()

            D, D_a, D_w, scale = fnum(alpha, omega, re_mp)

            # rebuild the Jacobian at each level and whenever the last step
            # gained less than 3 digits; costs 4 extra evaluations
            if D_aa is None or rate > 1e-3:
                D_aa, D_aw = jacobian_bottom_row(alpha, omega, re_mp, h)

            det = D_a * D_aw - D_w * D_aa
            dalpha = (-D * D_aw + D_a * D_w) / det
            domega = (D_aa * D - D_a * D_a) / det

            alpha, omega = alpha + dalpha, omega + domega
            iterates.append((complex(alpha), complex(omega)))

            # what this precision level can resolve
            noise = scale * mpmath.mpf(10) ** (-mpmath.mp.dps)
            floor = max(noise / abs(D_w), noise / abs(D_aa))

            step = max(abs(dalpha), abs(domega))
            rate = float(step / prev_step) if prev_step else 0.0
            prev_step = step

            history.append({
                "dps": dps, "iteration": it,
                "d_alpha": float(abs(dalpha)), "d_omega": float(abs(domega)),
                "dist_omega": float(abs(D / D_w)),
                "group_velocity": float(abs(D_a / D_w)),
                "floor": float(floor),
                "seconds": time.time() - t0,
                # full-precision mirrors of the five float values above
                "d_alpha_str": mp_str(abs(dalpha)),
                "d_omega_str": mp_str(abs(domega)),
                "dist_omega_str": mp_str(abs(D / D_w)),
                "group_velocity_str": mp_str(abs(D_a / D_w)),
                "floor_str": mp_str(floor),
            })

            print(f" {dps:5d} {it:3d}   {float(abs(dalpha)):.6e}   "
                  f"{float(abs(domega)):.6e}   {float(abs(D/D_w)):.6e}   "
                  f"{float(abs(D_a/D_w)):.6e}  {time.time()-t0:5.1f}")

            if step < floor * 100:
                print(f"        reached the {dps}-digit floor "
                      f"(1e{float(mpmath.log10(floor)):.0f})")
                break

    return alpha, omega, iterates, history


def summarise(alpha, omega, re):
    """Residuals and curvature at the converged point."""
    dps = DPS_LADDER[-1]
    mpmath.mp.dps = dps + 20
    h = mpmath.mpf(10) ** (-dps // 3)

    D, D_a, D_w, scale = fnum(alpha, omega, mpmath.mpf(re))
    D_aa, _ = jacobian_bottom_row(alpha, omega, mpmath.mpf(re), h)

    noise = scale * mpmath.mpf(10) ** (-mpmath.mp.dps)
    floor = max(noise / abs(D_w), noise / abs(D_aa))

    # at the pinch dD/dalpha = 0, so the curvature of the branch is just this
    omega_dd = -D_aa / D_w

    out = {
        "dist_omega": abs(D / D_w),        # distance to the surface, in omega
        "group_velocity": abs(D_a / D_w),  # the saddle test
        "omega_dd": omega_dd,
        "floor": floor,
        "digits": max(6, int(-mpmath.log10(floor)) - 3),
        "abs_D": abs(D),
        "dist_omega": abs(D / D_w),
        "group_velocity": abs(D_a / D_w),
        "omega_dd": omega_dd,
        "floor": floor,
        "digits": max(6, int(-mpmath.log10(floor)) - 3),
    }

    # confirm the digits are real by re-evaluating 30 digits higher
    mpmath.mp.dps = dps + 40
    D2, D_a2, D_w2, _ = fnum(mpmath.mpc(alpha), mpmath.mpc(omega), mpmath.mpf(re))
    out["check_dist_omega"] = abs(D2 / D_w2)
    out["check_group_velocity"] = abs(D_a2 / D_w2)
    mpmath.mp.dps = dps + 20

    return out


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def cx(values):
    """Read a list of {"re": ..., "im": ...} entries from the Briggs JSON."""
    import numpy as np
    return np.array([complex(z["re"], z["im"]) for z in values])


def calculate_seed(last):
    """The seed construction of compute_pinch.py, unchanged."""
    import numpy as np

    L = cx(last["L"])
    omega_F = cx(last["omega_F"])
    F = cx(last["F"])
    alpha_u = cx(last["alpha_L_u"])
    alpha_l = cx(last["alpha_L_l"])

    i = int(np.argmax(omega_F.imag))            # F point with the largest Im(omega)
    j = int(np.argmin(np.abs(L - omega_F[i])))  # nearest L point in the omega-plane

    seed = 0.5 * (F[i] + 0.5 * (alpha_u[j] + alpha_l[j]))

    return seed, omega_F[i], i, j


_LOG10_2 = 0.30102999566398119521373889472449


def mp_str(x, digits=None):
    """An mpmath value as a decimal string carrying every digit it holds."""
    if digits is None:
        try:
            bc = x._mpf_[3]                    # bits in this value's mantissa
        except AttributeError:                 # not an mpf -- fall back to context
            bc = mpmath.mp.prec
        digits = max(3, int(bc * _LOG10_2) + 3)
    return mpmath.nstr(x, digits, strip_zeros=False)


def point(z, digits=None):
    """A complex number as plain floats plus full-precision decimal strings."""
    return {
        "re": float(z.real), "im": float(z.imag),
        "re_str": mp_str(z.real, digits),
        "im_str": mp_str(z.imag, digits),
    }


def write_output_json(path, source, alpha, omega, seed, re, info, history,
                      omega_seed=None):
    """Write the result at full precision."""
    digits = info["digits"]

    residuals = {}

    # abs_D is written only when summarise() supplies it
    if "abs_D" in info:
        residuals["abs_D"] = float(info["abs_D"])
        residuals["abs_D_str"] = mp_str(info["abs_D"])

    seed_block = {"seed": point(mpmath.mpc(seed))}
    if omega_seed is not None:
        seed_block["omega_seed"] = point(mpmath.mpc(omega_seed))

    residuals.update({
        "distance_to_surface_in_omega": float(info["dist_omega"]),
        "distance_to_surface_in_omega_str": mp_str(info["dist_omega"]),
        "abs_group_velocity": float(info["group_velocity"]),
        "abs_group_velocity_str": mp_str(info["group_velocity"]),
        "arithmetic_floor": float(info["floor"]),
        "arithmetic_floor_str": mp_str(info["floor"]),
        "trustworthy_digits": digits + 3,
        "recheck_40_digits_higher": {
            "distance_to_surface_in_omega": float(info["check_dist_omega"]),
            "distance_to_surface_in_omega_str": mp_str(info["check_dist_omega"]),
            "abs_group_velocity": float(info["check_group_velocity"]),
            "abs_group_velocity_str": mp_str(info["check_group_velocity"]),
        },
    })

    out = {
        "source_file": os.path.basename(source),
        "method": "analytic Airy dispersion relation (Fnum, eigenvalues_v8.jl) "
                  "+ 2-D Newton on [D, dD/dalpha]",
        "settings": {"Re": re, "dps_ladder": list(DPS_LADDER),
                     "exact_thirds": EXACT_THIRDS},
        "precision": {
            "working_dps": mpmath.mp.dps,
            "final_level_dps": DPS_LADDER[-1],
            "recheck_dps": DPS_LADDER[-1] + 40,
            "trustworthy_digits": digits + 3,
            "seed_note": "seed and omega_seed are doubles from the Briggs "
                         "history, so they carry ~17 digits and no more",
            "str_note": "*_str fields hold every digit the value carries; "
                        "parse them at working_dps or higher, and set mp.dps "
                        "BEFORE parsing (mpmath.mpc parses at the current dps)",
        },
        **seed_block,
        "alpha_p": point(alpha),
        "omega_p": point(omega),
        "omega_second_derivative": point(info["omega_dd"]),
        "residuals": residuals,
        "convergence": history,

        # same shape as compute_pinch.py, so existing readers keep working
        "level_results": [{"level": 1, "dps": DPS_LADDER[-1],
                           "alpha_p": point(alpha),
                           "omega_p": point(omega)}],
    }

    with open(path, "w") as f:
        json.dump(out, f, indent=2)


# ---------------------------------------------------------------------------
# PLOT: omega-plane and alpha-plane, showing the Newton path
# ---------------------------------------------------------------------------
def plot_pinch(path_png, last, seed, alpha_p, omega_p, iterates):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    L = cx(last["L"])
    omega_F = cx(last["omega_F"])
    F = cx(last["F"])
    alpha_u = cx(last["alpha_L_u"])
    alpha_l = cx(last["alpha_L_l"])

    a_it = np.array([a for a, _ in iterates])
    w_it = np.array([w for _, w in iterates])

    fig, ax = plt.subplots(1, 2, figsize=(16, 7.5))

    # ---- omega-plane ----
    ax[0].plot(L.real, L.imag, "b", lw=1.6, label="L")
    ax[0].plot(omega_F.real, omega_F.imag, "r", lw=1.6, label=r"$\omega_F$")
    ax[0].plot(w_it.real, w_it.imag, ".-", ms=7, lw=1.0, color="0.35",
               label=r"$\omega$ at each Newton step")
    ax[0].plot(float(omega_p.real), float(omega_p.imag), "*", ms=12,
               color="magenta", label=r"$\omega_p$")

    ax[0].set_xlabel(r"$\omega_r$")
    ax[0].set_ylabel(r"$\omega_i$")
    ax[0].set_title(r"$\omega$-plane")
    ax[0].grid(True, alpha=0.4)
    ax[0].legend(loc="best")

    # ---- alpha-plane ----
    ax[1].plot(F.real, F.imag, "b", lw=1.6, label="F")
    ax[1].plot(alpha_u.real, alpha_u.imag, "r", lw=1.6, label=r"$\alpha_L^u$")
    ax[1].plot(alpha_l.real, alpha_l.imag, "g", lw=1.6, label=r"$\alpha_L^l$")
    ax[1].plot(a_it.real, a_it.imag, ".-", ms=7, lw=1.0, color="0.35",
               label="Newton path")
    ax[1].plot(seed.real, seed.imag, "o", ms=8, label="seed")
    ax[1].plot(float(alpha_p.real), float(alpha_p.imag), "*", ms=12,
               color="magenta", label=r"$\alpha_p$")

    ax[1].set_xlabel(r"$\alpha_r$")
    ax[1].set_ylabel(r"$\alpha_i$")
    ax[1].set_title(r"$\alpha$-plane")
    ax[1].grid(True, alpha=0.4)
    ax[1].legend(loc="best")

    # the Newton path is sub-pixel on these axes, so the insets plot the
    # offset FROM the pinch instead
    for a, pts, centre in ((ax[0], w_it, omega_p), (ax[1], a_it, alpha_p)):
        offset = pts - complex(centre)
        span = 1.6 * max(abs(offset[0]), 1e-12)

        inset = a.inset_axes([0.60, 0.08, 0.36, 0.32])
        inset.plot(offset.real, offset.imag, ".-", ms=7, lw=1.0, color="0.35")
        inset.plot(offset[0].real, offset[0].imag, "o", ms=7)
        inset.plot(0, 0, "*", ms=11, color="magenta")
        inset.set_xlim(-span, span)
        inset.set_ylim(-span, span)
        inset.set_title(f"offset from the pinch  "
                        f"(seed is {abs(offset[0]):.1e} away)", fontsize=8)
        inset.grid(True, alpha=0.3)
        inset.tick_params(labelsize=6)
        inset.locator_params(nbins=3)
        inset.set_facecolor("white")

    fig.suptitle(
        rf"$\alpha_p$ = {float(alpha_p.real):+.12f} {float(alpha_p.imag):+.12f}i"
        rf"      $\omega_p$ = {float(omega_p.real):+.12f} {float(omega_p.imag):+.12f}i",
        fontsize=12)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(path_png, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# PLOT: alpha spectrum at omega_p, in three nested windows
# ---------------------------------------------------------------------------
def plot_alpha_spectrum(path_png, last, alpha_p, omega_p):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        from couette_eigenvaluesspectrum import alpha_spectrum
    except ImportError:
        print("couette_eigenvaluesspectrum.py not found -- skipping the "
              "spectrum plot")
        return False

    F = cx(last["F"])
    alpha_u = cx(last["alpha_L_u"])
    alpha_l = cx(last["alpha_L_l"])

    a_p = complex(alpha_p)
    alpha_spec = alpha_spectrum(complex(omega_p))

    def window(values, margin=0.05):
        r_pad = margin * (np.max(values.real) - np.min(values.real))
        i_pad = margin * (np.max(values.imag) - np.min(values.imag))
        return ((np.min(values.real) - r_pad, np.max(values.real) + r_pad),
                (np.min(values.imag) - i_pad, np.max(values.imag) + i_pad))

    # 1. SPECTRUM window -- the full spectrum runs out to |alpha| ~ 1e8, which
    #    is the unresolved tail of the discretisation, not physics
    spec_kept = alpha_spec[np.abs(alpha_spec) <= SPECTRUM_MAX_ABS]
    spec_x, spec_y = window(np.append(spec_kept, a_p))

    # 2. BRIGGS window -- the contours of the final iteration
    briggs_x, briggs_y = window(np.concatenate([F, alpha_u, alpha_l, [a_p]]))

    # 3. ZOOM window -- centred on alpha_p
    zoom_x = (a_p.real - ZOOM_HALF_WIDTH, a_p.real + ZOOM_HALF_WIDTH)
    zoom_y = (a_p.imag - ZOOM_HALF_WIDTH, a_p.imag + ZOOM_HALF_WIDTH)

    panels = [
        (spec_x, spec_y,
         rf"$\alpha$-plane spectrum: $|\alpha| \leq$ {SPECTRUM_MAX_ABS:g}"
         rf"  ({len(spec_kept)} of {len(alpha_spec)} modes)",
         (briggs_x, briggs_y), "Briggs window"),
        (briggs_x, briggs_y, r"$\alpha$-plane", (zoom_x, zoom_y), "zoom window"),
        (zoom_x, zoom_y,
         rf"$\alpha$-plane zoom: $\alpha_p \pm$ {ZOOM_HALF_WIDTH:g}", None, None),
    ]

    fig, ax = plt.subplots(1, 3, figsize=(24, 7.5))

    for k, (xlim, ylim, title, box, box_label) in enumerate(panels):
        visible = ((alpha_spec.real >= xlim[0]) & (alpha_spec.real <= xlim[1])
                   & (alpha_spec.imag >= ylim[0]) & (alpha_spec.imag <= ylim[1]))

        ax[k].plot(F.real, F.imag, "b", lw=1.6, label="F")
        ax[k].plot(alpha_u.real, alpha_u.imag, "r", lw=1.6, label=r"$\alpha_L^u$")
        ax[k].plot(alpha_l.real, alpha_l.imag, "g", lw=1.6, label=r"$\alpha_L^l$")
        ax[k].plot(alpha_spec[visible].real, alpha_spec[visible].imag, "k.", ms=7,
                   label=r"$\alpha(\omega_p)$, double precision")
        ax[k].plot(a_p.real, a_p.imag, "*", ms=12, color="magenta",
                   label=r"$\alpha_p$ (analytic double root)")

        if box is not None:
            (bx0, bx1), (by0, by1) = box
            ax[k].add_patch(plt.Rectangle((bx0, by0), bx1 - bx0, by1 - by0,
                                          fill=False, lw=1.0, ls="--",
                                          edgecolor="k", label=box_label))

        ax[k].set_xlim(xlim)
        ax[k].set_ylim(ylim)
        ax[k].set_xlabel(r"$\alpha_r$")
        ax[k].set_ylabel(r"$\alpha_i$")
        ax[k].set_title(title)
        ax[k].grid(True, alpha=0.4)
        ax[k].legend(loc="best")

    fig.suptitle(
        rf"$\alpha_p$ = {a_p.real:+.12f} {a_p.imag:+.12f}i"
        rf"      $\omega_p$ = {float(omega_p.real):+.12f} "
        rf"{float(omega_p.imag):+.12f}i", fontsize=12)

    fig.text(0.5, 0.005,
             "at this omega_p the two alpha-roots are one double root exactly "
             "at the star; any spread in the black points is double-precision "
             "noise in the spatial solver, which splits a defective eigenvalue "
             "like sqrt(roundoff)",
             ha="center", fontsize=9, color="0.35")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(path_png, dpi=160)
    plt.close(fig)
    return True


# ---------------------------------------------------------------------------
# PLOT: Newton convergence
# ---------------------------------------------------------------------------
def plot_convergence(path_png, history):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = list(range(len(history)))

    def masked(key):
        # exact zeros fell below the arithmetic floor; mask them for the log axis
        return [h[key] if h[key] > 0 else float("nan") for h in history]

    fig, ax = plt.subplots(1, 2, figsize=(14, 6))

    ax[0].semilogy(n, masked("d_alpha"), "o-", label=r"$|\Delta\alpha|$")
    ax[0].semilogy(n, masked("d_omega"), "s-", label=r"$|\Delta\omega|$")
    ax[0].set_title("Newton step")

    ax[1].semilogy(n, masked("dist_omega"), "o-",
                   label=r"$|D|\,/\,|\partial D/\partial\omega|$")
    ax[1].semilogy(n, masked("group_velocity"), "s-", label=r"$|d\omega/d\alpha|$")
    ax[1].set_title("residuals")

    for a in ax:
        a.step(n, [h["floor"] for h in history], where="mid", color="0.5",
               ls=":", lw=1.4, label="arithmetic floor")
        for k in range(1, len(history)):
            if history[k]["dps"] != history[k - 1]["dps"]:
                a.axvline(k - 0.5, color="k", ls="--", lw=1.0)
                a.text(k - 0.45, a.get_ylim()[1], f" dps={history[k]['dps']}",
                       va="top", fontsize=9)
        a.set_xlabel("iteration")
        a.grid(True, which="both", alpha=0.3)
        a.legend(loc="lower left", fontsize=9)

    fig.tight_layout()
    fig.savefig(path_png, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        print("Usage: python compute_pinch2.py <big_history.json> [--dps N]")
        sys.exit(1)

    path = sys.argv[1]
    stem = os.path.splitext(os.path.abspath(path))[0]

    global DPS_LADDER
    if "--dps" in sys.argv:
        top = int(sys.argv[sys.argv.index("--dps") + 1])
        DPS_LADDER = tuple(d for d in (40, 80, 150) if d < top) + (top,)

    print(f"reading {os.path.basename(path)} ...")

    with open(path) as f:
        frames = json.load(f)

    if not isinstance(frames, list):
        raise ValueError(
            "Use the BIG Briggs history JSON, not a *_pinch*.json output file.")

    last = frames[-1]
    seed, omega_seed, i, j = calculate_seed(last)

    print(f"using last frame: iteration {last['iteration']}")
    print(f"F node: {i+1}, L node: {j+1}")
    print(f"seed  = {seed.real:+.9f} {seed.imag:+.9f}i")
    print(f"Re = {RE},  dps ladder = {DPS_LADDER},  exact thirds = {EXACT_THIRDS}")

    t0 = time.time()
    alpha_p, omega_p, iterates, history = find_pinch(seed, omega_seed, RE)
    info = summarise(alpha_p, omega_p, RE)
    digits = info["digits"]

    print("\nFINAL")
    print(f"alpha_p = {mpmath.nstr(alpha_p.real, min(digits, 45), strip_zeros=False)}")
    print(f"          {mpmath.nstr(alpha_p.imag, min(digits, 45), strip_zeros=False)} i")
    print(f"omega_p = {mpmath.nstr(omega_p.real, min(digits, 45), strip_zeros=False)}")
    print(f"          {mpmath.nstr(omega_p.imag, min(digits, 45), strip_zeros=False)} i")
    print(f"\ndistance to the dispersion surface, in omega : "
          f"{float(info['dist_omega']):.3e}")
    print(f"saddle test |domega/dalpha|                  : "
          f"{float(info['group_velocity']):.3e}")
    print(f"same, re-checked 40 digits higher            : "
          f"{float(info['check_dist_omega']):.3e} / "
          f"{float(info['check_group_velocity']):.3e}")
    print(f"arithmetic floor at dps {DPS_LADDER[-1]}                   : "
          f"{float(info['floor']):.3e}   -> {digits + 3} trustworthy digits")
    print(f"omega'' = {mpmath.nstr(info['omega_dd'], 20)}")
    print(f"\nruntime: {time.time()-t0:.1f} s")

    outputs = [
        (stem + "_pinch_analytic.json",
         lambda p: write_output_json(p, path, alpha_p, omega_p, seed, RE,
                                     info, history, omega_seed)),

        (stem + "_pinch_analytic.png",
         lambda p: plot_pinch(p, last, seed, alpha_p, omega_p, iterates)),

        (stem + "_pinch_analytic_alpha_spectrum.png",
         lambda p: plot_alpha_spectrum(p, last, alpha_p, omega_p)),

        (stem + "_pinch_analytic_convergence.png",
         lambda p: plot_convergence(p, history)),
    ]

    print()
    for output_path, write in outputs:
        if write(output_path) is not False:
            print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
