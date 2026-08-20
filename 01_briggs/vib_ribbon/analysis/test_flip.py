"""Decisive test for the iteration-14 upper-branch flip in v4.6.

Walks the left horizontal of the ribbon contour outward from the seed node
(index 48, omega_r = 0.12121) at omega_i = -0.05130 (iteration 14) using
(a) PURE continuation  -- nearest eigenvalue to the secant predictor, no filter
(b) v4.6 continuation  -- side-of-F filter applied first, then nearest
and compares both with what the run actually recorded.

If (a) and (b) diverge at the node where the branch crosses the real alpha
axis, the flip is caused by the side-of-F filter, not by the detour.
"""
import numpy as np
from multiprocessing import Pool
from os_spatial import spectrum

d = np.load("summary_v46.npz")
IT = 14
au_run = d[f"au_{IT}"]
L = d[f"L_{IT}"]
F13 = d["F_13"]
F14 = d[f"F_{IT}"]

I0, I1 = 48, 95          # seed .. riser bottom
nodes = L[I0:I1 + 1]
print("omega_i =", L[0].imag, " nodes", nodes[0], "->", nodes[-1])


def normals(F):
    t = np.empty_like(F)
    t[1:-1] = F[2:] - F[:-2]
    t[0] = t[1]; t[-1] = t[-2]
    n = 1j * t / np.abs(t)
    return n


def sides(evs, F, nF):
    k = np.abs(evs[:, None] - F[None, :]).argmin(axis=1)
    return np.real(np.conj(nF[k]) * (evs - F[k]))


def continue_branch(nodes, spectra, sd_list, seed, use_side, da_max=0.03, sep_frac=0.6):
    vals = np.empty(len(nodes), complex)
    flags = np.zeros(len(nodes), bool)
    a1 = a2 = complex(seed); w1 = w2 = nodes[0]
    have2 = False
    nfb = 0
    for j, w in enumerate(nodes):
        pred = a1
        if have2 and abs(w1 - w2) > 1e-14:
            disp = (w - w1) / (w1 - w2) * (a1 - a2)
            if abs(disp) > da_max:
                disp *= da_max / abs(disp)
            pred = a1 + disp
        evs = spectra[j]; sd = sd_list[j]
        keep = np.arange(len(evs))
        if use_side:
            k = np.where(sd > 0.0)[0]
            if len(k) == 0:
                nfb += 1
            else:
                keep = k
        dist = np.abs(evs[keep] - pred)
        k1 = dist.argmin()
        d1 = dist[k1]
        d2 = np.sort(dist)[1] if len(dist) > 1 else np.inf
        cand = evs[keep[k1]]
        vals[j] = cand
        flags[j] = (abs(cand - a1) > da_max) or (d1 > sep_frac * d2)
        a2, w2 = a1, w1
        a1, w1 = cand, w
        have2 = True
    return vals, flags, nfb


if __name__ == "__main__":
    with Pool(2) as p:
        spectra = p.map(spectrum, list(nodes))
    for name, F in [("F(it13)", F13), ("F(it14)", F14)]:
        nF = normals(F)
        sd = [sides(s, F, nF) for s in spectra]
        seed = au_run[I0]
        v_pure, f_pure, _ = continue_branch(nodes, spectra, sd, seed, use_side=False)
        v_side, f_side, nfb = continue_branch(nodes, spectra, sd, seed, use_side=True)
        print(f"\n=== side projections taken against {name}; side fallbacks = {nfb} ===")
        print(f"{'idx':>4}{'om_r':>9} {'PURE continuation':>21} {'SIDE-filtered':>21} "
              f"{'v4.6 recorded':>21} {'|side-pure|':>11} {'sd(pure)':>10} {'flag':>5}")
        for j in range(len(nodes)):
            i = I0 + j
            if not (i < 62 or i % 8 == 0 or i >= 93):
                continue
            k = np.abs(spectra[j] - v_pure[j]).argmin()
            print(f"{i:4d}{nodes[j].real:9.5f} {v_pure[j].real:+.4f}{v_pure[j].imag:+.4f}j"
                  f"      {v_side[j].real:+.4f}{v_side[j].imag:+.4f}j"
                  f"      {au_run[i].real:+.4f}{au_run[i].imag:+.4f}j"
                  f"   {abs(v_side[j]-v_pure[j]):11.4f} {sd[j][k]:+10.5f}"
                  f" {'F' if f_side[j] else '.':>5}")
        print(f"end of left horizontal: pure {v_pure[-1]:.5f} | side {v_side[-1]:.5f} "
              f"| run {au_run[I1]:.5f}")
