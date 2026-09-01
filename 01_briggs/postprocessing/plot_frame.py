
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------- settings --
HERE = Path(__file__).resolve().parent
JSON_FILE = HERE.parent / "results" / "json" / "contour_iteration_v4.4.json"
ITERATION = 353       # None = last entry in the file; or a number
POINTS = True         # draw the discrete sample points
MARKEVERY = 1         # show every n-th point
ZOOM_PAD = 6          # close-up half-width, in units of the gap
ALPHA_REF = None      # e.g. 0.5725194852012776 - 3.038860744115944j
SHOW = True           # also open interactive windows
# -----------------------------------------------------------------------------

if len(sys.argv) > 1:
    ITERATION = int(sys.argv[1])


def to_complex(arr):
    return np.array([z["re"] for z in arr]) + 1j * np.array([z["im"] for z in arr])


# --- load one iteration ------------------------------------------------------
data = json.load(open(JSON_FILE))
if ITERATION is None:
    entry = data[-1]
else:
    hits = [d for d in data if d["iteration"] == ITERATION]
    if not hits:
        raise SystemExit(f"iteration {ITERATION} not in {JSON_FILE.name} "
                         f"(1 .. {data[-1]['iteration']})")
    entry = hits[0]

it = entry["iteration"]
L = to_complex(entry["L"])
wF = to_complex(entry["omega_F"])
F = to_complex(entry["F"])
au = to_complex(entry["alpha_L_u"])
al = to_complex(entry["alpha_L_l"])

# --- branch gap --------------------------------------------------------------
gap = np.abs(au - al)
i = int(gap.argmin())
d = gap[i]
mid = 0.5 * (au[i] + al[i])

mk = dict(marker="o", markersize=2.5, markevery=MARKEVERY) if POINTS else {}


def draw_alpha(ax, msize):

    """alpha-plane curves + dashed segment across the gap; returns its midpoint."""
    m = dict(mk, markersize=msize) if POINTS else {}
    ax.plot(F.real, F.imag, "b-", label=r"$F$", **m)
    ax.plot(au.real, au.imag, "r-", label=r"$\alpha_L^u$", **m)
    ax.plot(al.real, al.imag, "g-", label=r"$\alpha_L^l$", **m)
    ax.plot([au.real[i], al.real[i]], [au.imag[i], al.imag[i]], "k.--", lw=1.2, zorder=5)
    if ALPHA_REF is not None:
        ax.plot(ALPHA_REF.real, ALPHA_REF.imag, "mx", ms=10, mew=2, zorder=6,
                label=r"$\alpha_p$")
    ax.set_xlabel(r"$\alpha_r$")
    ax.set_ylabel(r"$\alpha_i$")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    return mid.real, mid.imag


# --- figure 1: the two panels ------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

ax1.plot(L.real, L.imag, "b-", label=r"$L$", **mk)
ax1.plot(wF.real, wF.imag, "r-", label=r"$\omega_F$", **mk)
ax1.plot(L.real[i], L.imag[i], "k*", ms=12, zorder=5, label=r"$\omega$ of min. gap")
ax1.set_xlabel(r"$\omega_r$")
ax1.set_ylabel(r"$\omega_i$")
ax1.set_title(r"$\omega$-plane")
ax1.grid(True, alpha=0.3)
ax1.legend(loc="best")

xm, ym = draw_alpha(ax2, 2.5)
ax2.annotate(rf"$d_{{\min}} = {d:.4e}$", xy=(xm, ym), xytext=(8, 8),
             textcoords="offset points", zorder=7,
             bbox=dict(boxstyle="round,pad=0.3", fc="w", ec="k", alpha=0.85))
ax2.set_title(r"$\alpha$-plane")

fig.suptitle(rf"iteration {it}   $\omega_i = {L.imag[0]:+.6e}$   "
             rf"$d_{{\min}} = {d:.4e}$  at  $\omega_r = {L.real[i]:.5f}$  (index {i})")
fig.tight_layout()
out = JSON_FILE.with_name(f"{JSON_FILE.stem}_iter{it}.png")
fig.savefig(out, dpi=150)

# --- figure 2: alpha-plane close-up ------------------------------------------
figz, axz = plt.subplots(figsize=(7.5, 7))
xm, ym = draw_alpha(axz, 4)
axz.annotate(rf"$d_{{\min}} = {d:.6e}$", xy=(xm, ym), xytext=(14, 14),
             textcoords="offset points", fontsize=11, zorder=7,
             arrowprops=dict(arrowstyle="-", lw=0.8),
             bbox=dict(boxstyle="round,pad=0.35", fc="w", ec="k", alpha=0.9))
half = ZOOM_PAD * d
axz.set_xlim(xm - half, xm + half)
axz.set_ylim(ym - half, ym + half)
axz.set_aspect("equal", adjustable="box")
axz.text(0.02, 0.02, f"a_u[{i}]    = {au[i]:.6f}\na_l[{i}]    = {al[i]:.6f}\n"
                     f"midpoint = {mid:.6f}",
         transform=axz.transAxes, fontsize=9, va="bottom", family="monospace",
         bbox=dict(boxstyle="round,pad=0.3", fc="w", ec="0.6", alpha=0.85))
axz.set_title(rf"$\alpha$-plane close-up  --  iteration {it}, "
              rf"$\omega_r = {L.real[i]:.5f}$")
figz.tight_layout()
outz = out.with_name(f"{out.stem}_alpha_zoom.png")
figz.savefig(outz, dpi=150)

# --- numbers -----------------------------------------------------------------
print(f"iteration {it} of {JSON_FILE.name}  ({len(data)} entries)")
print(f"  omega_i   = {L.imag[0]:+.9e}")
print(f"  d_min     = {d:.6e}   at index {i}, omega_r = {L.real[i]:.6f}")
print(f"  alpha_u   = {au[i]:.6f}")
print(f"  alpha_l   = {al[i]:.6f}")
print(f"  midpoint  = {mid:.6f}"
      + (f"   |mid - alpha_p| = {abs(mid - ALPHA_REF):.3e}" if ALPHA_REF else ""))
print(f"  saved     -> {out}")
print(f"            -> {outz}")

if SHOW:
    plt.show()
plt.close("all")
