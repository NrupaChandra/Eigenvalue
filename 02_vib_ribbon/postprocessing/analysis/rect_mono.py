"""Closed-rectangle monodromy test in the omega-plane + inspection of what the
lower branch actually is.

Rectangle: omega_r in [0.12121, R], omega_i in [-0.05130, -0.04730]
           (the two horizontal levels of v4.6 iterations 14 and 13).
Track alpha by pure continuation all the way round.  If alpha does not come
back to itself, a branch point of alpha(omega) lies inside -- which is what
makes the it-13 -> it-14 upper branch jump.
Shrinking R localises its real part.
"""
import os
import numpy as np
from multiprocessing import Pool
from os_spatial import spectrum

WI_TOP, WI_BOT = -0.04730, -0.05130      # iterations 13 and 14
WR_L = 0.12121212121212122               # the seed node
SEED_ALPHA = 0.2560 - 0.0381j            # recorded alpha_u at the seed, it 13


def leg(a, b, n):
    return list(a + (b - a) * np.linspace(0, 1, n)[1:])


def loop_nodes(R, nh=48, nv=10):
    c0 = complex(WR_L, WI_TOP)
    p = [c0]
    p += leg(c0, complex(R, WI_TOP), nh)          # right along it-13 level
    p += leg(complex(R, WI_TOP), complex(R, WI_BOT), nv)   # down
    p += leg(complex(R, WI_BOT), complex(WR_L, WI_BOT), nh)  # left along it-14
    p += leg(complex(WR_L, WI_BOT), c0, nv)       # up, back to start
    return p


def track(nodes, spectra, seed, da_max=0.03):
    a1 = a2 = complex(seed); w1 = w2 = nodes[0]
    out = []
    have2 = False
    worst = 0.0
    for j, w in enumerate(nodes):
        pred = a1
        if have2 and abs(w1 - w2) > 1e-14:
            disp = (w - w1) / (w1 - w2) * (a1 - a2)
            if abs(disp) > da_max:
                disp *= da_max / abs(disp)
            pred = a1 + disp
        ev = spectra[j]
        cand = ev[np.abs(ev - pred).argmin()]
        worst = max(worst, abs(cand - a1))
        out.append(cand)
        a2, w2 = a1, w1
        a1, w1 = cand, w
        have2 = True
    return np.array(out), worst


if __name__ == "__main__":
    print("=== closed-rectangle monodromy, omega_i from -0.04730 (it13) to -0.05130 (it14) ===")
    print(f"{'R (right edge)':>15} {'alpha start':>20} {'alpha after loop':>20} "
          f"{'|monodromy|':>12} {'max step':>9}")
    import sys
    for R in [float(x) for x in sys.argv[1:]]:
        nodes = loop_nodes(R, nh=int(os.environ.get('NH','48')), nv=int(os.environ.get('NV','10')))
        with Pool(2) as p:
            spectra = p.map(spectrum, nodes)
        a, worst = track(nodes, spectra, SEED_ALPHA)
        print(f"{R:15.4f} {a[0].real:+.5f}{a[0].imag:+.5f}j  {a[-1].real:+.5f}{a[-1].imag:+.5f}j "
              f"{abs(a[-1]-a[0]):12.5f} {worst:9.4f}")

    if len(sys.argv)>2: raise SystemExit
    print("\n=== what is the LOWER branch?  spectrum at omega = 0.12121 - 0.04730i ===")
    ev = spectrum(complex(WR_L, WI_TOP))
    below = ev[ev.imag < 0]
    above = ev[ev.imag > 0]
    print(f"{len(ev)} finite eigenvalues: {len(above)} with Im a > 0, {len(below)} with Im a < 0")
    print("  eigenvalues with Im alpha < 0, sorted by |alpha|:")
    for z in below[np.argsort(np.abs(below))][:10]:
        print(f"     {z:+.6f}   |a| = {abs(z):.4f}")
    print("  eigenvalues with Im alpha < 0 and 0 <= Re alpha <= 1 (i.e. genuinely 'under F'):")
    m = below[(below.real >= 0) & (below.real <= 1)]
    for z in m[np.argsort(np.abs(m.imag))][:10]:
        print(f"     {z:+.6f}")
