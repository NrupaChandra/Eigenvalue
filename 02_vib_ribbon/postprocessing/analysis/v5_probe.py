"""Why does v5 still hop, at N = 300?

v5 frame 14 (clean)  is at omega_i = -0.0498
v5 frame 15 (hopped) is at omega_i = -0.0538
Neither level was covered by the earlier tests (which used -0.04730 / -0.05130).

Track alpha along the v5 path -- seed node -> right along the horizontal ->
up the left riser -> onto the arc -- at BOTH levels, at v5's own sampling, and
find the first node where the two curves part company.  That says whether the
hop is on the horizontal (my original diagnosis) or on the riser/arc.
"""
import numpy as np, sys
from multiprocessing import Pool
from os_spatial import spectrum

WR_SEED, WR_ARC = 0.12374581939799331, 0.24     # N=300 seed node, riser foot
OMEGA_F, R_ARC = 0.25, 0.01
N_HOR, N_VER, N_ARC = 140, 160, 15              # v5 spacing on the horizontal
DA_MAX = 0.012


def path(wi):
    hor = np.linspace(WR_SEED, WR_ARC, N_HOR) + 1j * wi
    ris = (OMEGA_F - R_ARC) + 1j * np.linspace(wi, 0.0, N_VER)
    arc = OMEGA_F + R_ARC * np.exp(1j * np.linspace(np.pi, np.pi / 2, N_ARC))
    return np.concatenate([hor, ris[1:], arc])


def track(nodes, spectra, seed, da_max=DA_MAX):
    a1 = a2 = complex(seed); w1 = w2 = nodes[0]
    out = []; have2 = False; worst = 0.0
    for j, w in enumerate(nodes):
        pred = a1
        if have2 and abs(w1 - w2) > 1e-14:
            d = (w - w1) / (w1 - w2) * (a1 - a2)
            if abs(d) > da_max:
                d *= da_max / abs(d)
            pred = a1 + d
        ev = spectra[j]
        cand = ev[np.abs(ev - pred).argmin()]
        worst = max(worst, abs(cand - a1))
        out.append(cand)
        a2, w2 = a1, w1
        a1, w1 = cand, w
        have2 = True
    return np.array(out), worst


if __name__ == "__main__":
    levels = [float(x) for x in sys.argv[1:]] or [-0.0498, -0.0538]
    res = {}
    for wi in levels:
        nodes = path(wi)
        with Pool(2) as p:
            spectra = p.map(spectrum, list(nodes))
        # seed: the branch sitting next to F near alpha ~ 0.26 - 0.05i
        s0 = spectra[0][np.abs(spectra[0] - (0.26 - 0.05j)).argmin()]
        a, worst = track(nodes, spectra, s0)
        res[wi] = (nodes, a, worst)
        print(f"omega_i = {wi:+.5f}  seed alpha = {s0:.5f}  max step = {worst:.4f}")
        print(f"   end of horizontal (idx {N_HOR-1}, omega_r=0.24): {a[N_HOR-1]:.5f}")
        print(f"   riser top         (omega=0.24+0i)             : {a[N_HOR+N_VER-2]:.5f}")
        print(f"   arc apex          (omega=0.25+0.01i)          : {a[-1]:.5f}")

    if len(levels) == 2:
        n0, a0, _ = res[levels[0]]
        n1, a1, _ = res[levels[1]]
        d = np.abs(a1 - a0)
        print("\nwhere do the two levels part company?  |alpha(-0.0538) - alpha(-0.0498)|")
        seg = lambda j: ("horizontal" if j < N_HOR else
                         "riser" if j < N_HOR + N_VER - 1 else "arc")
        for j in range(len(d)):
            if j < 4 or (j % 10 == 0 and j < N_HOR + 40) or d[j] > 0.02 and d[j - 1] <= 0.02:
                print(f"   idx {j:4d} [{seg(j):10s}] omega={n0[j]:.5f} "
                      f"a0={a0[j]:+.5f} a1={a1[j]:+.5f}  |diff|={d[j]:.4f}")
        k = np.argmax(d > 0.02) if (d > 0.02).any() else -1
        print(f"\nfirst node with |diff| > 0.02: idx {k} in the {seg(k)} "
              f"at omega = {n0[k]:.5f}" if k >= 0 else "\nnever exceeds 0.02")
        print(f"max |diff| over the whole path = {d.max():.4f} at idx {d.argmax()} ({seg(d.argmax())})")
