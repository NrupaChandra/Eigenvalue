import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load("summary_v46.npz")
it = d["it"]
OUT = "/sessions/magical-upbeat-hamilton/mnt/Eigenvalue/01_briggs/vib_ribbon/results/"
HL1, RL0, RL1, AR0, AR1, RR0, RR1, HR0 = 94, 95, 254, 255, 314, 315, 474, 475
ARC = (AR0 + AR1) // 2

plt.rcParams.update({"font.size": 9, "figure.dpi": 150, "axes.grid": True,
                     "grid.alpha": .3, "savefig.bbox": "tight"})

# ---------------------------------------------------------------- figure 1
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
for a, its, ttl in [(ax[0], [11, 12, 13], "iterations 11-13  (before the hop)"),
                    (ax[1], [14, 15, 16], "iterations 14-16  (after the hop)")]:
    for k, i in enumerate(its):
        au, al, F = d[f"au_{i}"], d[f"al_{i}"], d[f"F_{i}"]
        c = plt.cm.viridis(k / 3)
        a.plot(au.real, au.imag, "-", color=c, lw=1.4, label=f"$\\alpha_u$, it {i}")
        a.plot(F.real, F.imag, "--", color="0.5", lw=.9)
    a.plot(d["au_13"][ARC].real, d["au_13"][ARC].imag, "ko", ms=6, mfc="none")
    a.plot(d["au_14"][ARC].real, d["au_14"][ARC].imag, "rs", ms=6, mfc="none")
    a.axhline(0, color="k", lw=.6)
    a.set_xlim(0.15, 0.8); a.set_ylim(-0.15, 0.30)
    a.set_xlabel(r"$\Re\,\alpha$"); a.set_ylabel(r"$\Im\,\alpha$")
    a.set_title(ttl); a.legend(fontsize=7, loc="upper left")
fig.suptitle(r"v4.6 upper branch $\alpha_u$ along the ribbon contour; grey dashed = $F$."
             "\n" r"circle / square = value at the arc ($\omega=\omega_f$)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig(OUT + "diag_v46_1_branch_hop.png")

# ---------------------------------------------------------------- figure 2
fig, ax = plt.subplots(2, 2, figsize=(10, 6.4))
a = ax[0, 0]
a.plot(it, np.abs(d["au_botL"]), lw=1)
a.axvline(14, color="r", ls=":"); a.set_yscale("linear")
a.set_xscale("log"); a.set_xlabel("iteration")
a.set_ylabel(r"$|\alpha_u|$ at $\omega=0.24+i\omega_i$")
a.set_title("end of the left horizontal: step change at it 14")

a = ax[0, 1]
a.plot(it, d["au_arc"].real, label=r"$\Re\,\alpha_u(\omega_f)$")
a.plot(it, d["au_arc"].imag, label=r"$\Im\,\alpha_u(\omega_f)$")
a.axvline(14, color="r", ls=":"); a.set_xscale("log")
a.set_xlabel("iteration"); a.legend(fontsize=7)
a.set_title(r"value at the arc ($\omega_f=0.25$): two distinct values only")

a = ax[1, 0]
a.plot(it, d["dpair"], label=r"$\min_{i,j}|\alpha_u^i-\alpha_l^j|$")
a.axhline(1e-4, color="r", ls="--", lw=.8, label="pinch_tol $=10^{-4}$")
a.set_yscale("log"); a.set_xscale("log"); a.set_xlabel("iteration")
a.legend(fontsize=7); a.set_title("pinch criterion is never approached")

a = ax[1, 1]
a.plot(it, d["om_i"], label=r"$\omega_i$ of L")
a.plot(it, d["omF_max"], "--", label=r"$\max\Im\,\omega_F$")
a.plot(it, [f.imag.mean() if False else 0 for f in it], alpha=0)
a.set_xscale("log"); a.set_xlabel("iteration"); a.legend(fontsize=7)
a.set_title(r"L rides directly on $\omega_F$ (the stall)")
fig.tight_layout()
fig.savefig(OUT + "diag_v46_2_timeseries.png")

# ---------------------------------------------------------------- figure 3
fig, ax = plt.subplots(figsize=(6.4, 4.2))
lv13 = {48: 0.40636 + 0.00911j, 96: 0.40636 + 0.00911j, 192: 0.40636 + 0.00911j}
lv14 = {48: 0.65360 + 0.08051j, 96: 0.40647 + 0.00375j, 192: 0.40647 + 0.00375j}
n = [48, 96, 192]
ax.plot(n, [abs(lv13[k] - lv13[192]) for k in n], "o-", label=r"$\omega_i=-0.04730$ (it 13)")
ax.plot(n, [abs(lv14[k] - lv14[192]) for k in n], "s-", label=r"$\omega_i=-0.05130$ (it 14)")
ax.axvline(48, color="r", ls=":", label="resolution used by v4.6")
ax.set_yscale("symlog", linthresh=1e-6)
ax.set_xscale("log"); ax.set_xticks(n); ax.set_xticklabels(n)
ax.set_xlabel("nodes on the left horizontal (seed $\\to\\ \\omega_r=0.24$)")
ax.set_ylabel(r"$|\alpha(0.24) - \alpha_{\rm converged}(0.24)|$")
ax.set_title("the it-14 branch is a sampling artefact:\nit disappears as soon as the "
             "horizontal is refined")
ax.legend(fontsize=8)
fig.savefig(OUT + "diag_v46_3_resolution.png")

# ---------------------------------------------------------------- figure 4
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
for a, i in [(ax[0], 13), (ax[1], 1000)]:
    au, al, F = d[f"au_{i}"], d[f"al_{i}"], d[f"F_{i}"]
    a.plot(F.real, F.imag, "k--", lw=1.2, label="$F$")
    a.plot(au.real, au.imag, "r-", lw=1.3, label=r"$\alpha_u$")
    a.plot(al.real, al.imag, "g-", lw=1.3, label=r"$\alpha_l$")
    a.axhline(0, color="0.3", lw=.6)
    a.set_xlabel(r"$\Re\,\alpha$"); a.set_ylabel(r"$\Im\,\alpha$")
    a.set_title(f"iteration {i}"); a.legend(fontsize=8)
ax[0].set_ylim(-3.6, 0.4); ax[1].set_ylim(-3.6, 0.4)
fig.suptitle("F is dragged bodily down to $\\Im\\,\\alpha\\approx-1.15$ while "
             "$\\alpha_l$ stays on the wall-mode ladder near $-3.3i$", fontsize=10)
fig.tight_layout()
fig.savefig(OUT + "diag_v46_4_F_dragged.png")
print("wrote 4 figures to", OUT)
