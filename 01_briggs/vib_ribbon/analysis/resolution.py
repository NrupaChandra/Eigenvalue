"""Is the iteration-13 -> 14 upper-branch flip a sampling artefact?

Continue alpha_u from the seed node (omega_r = 0.12121) rightwards to
omega_r = 0.24 along each of the two horizontal levels, at the run's own
resolution and at 4x finer, and compare the endpoint.

If the endpoint depends on the sampling, the tracker is under-resolving a
near-degeneracy and neither curve can be trusted.
"""
import numpy as np, sys
from multiprocessing import Pool
from os_spatial import spectrum
from rect_mono import track

WR_L, WR_R = 0.12121212121212122, 0.24
SEEDS = {-0.04730: 0.2560 - 0.0381j, -0.05130: 0.2593 - 0.0470j}

if __name__ == "__main__":
    wi = float(sys.argv[1])
    print(f"level omega_i = {wi}   seed alpha = {SEEDS[wi]}")
    print(f"{'nodes':>7} {'d(omega_r)':>11} {'alpha at omega_r=0.24':>24} {'max step':>10}")
    for n in [48, 96, 192]:
        nodes = list(np.linspace(WR_L, WR_R, n) + 1j * wi)
        with Pool(2) as p:
            spectra = p.map(spectrum, nodes)
        a, worst = track(nodes, spectra, SEEDS[wi])
        print(f"{n:7d} {(WR_R-WR_L)/(n-1):11.5f}   {a[-1].real:+.5f}{a[-1].imag:+.5f}j"
              f"        {worst:10.4f}")
        # where does Im(alpha) cross zero?
        s = np.sign(a.imag)
        cr = np.where(np.diff(s) != 0)[0]
        for c in cr:
            print(f"          Im(alpha) crosses 0 at omega_r ~ {nodes[c].real:.5f}, "
                  f"alpha ~ {a[c]:.5f}")
