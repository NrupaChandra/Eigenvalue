"""How accurate are the eigenvalues themselves?  The Chebyshev D4 recursion at
num_modes=150 is badly conditioned; if the eigenvalue noise floor is larger than
the local branch spacing, no tracker can resolve the branches there."""
import numpy as np, os_spatial as M
from numpy.linalg import cond

def spec_nm(nm, w):
    M.NM = nm
    M.Y, M.D0, M.D1, M.D2, M.D3, M.D4 = M._build_operators(nm)
    M.U = M.Y[:, None] * np.ones((1, nm))
    return M.spectrum(w, cap=1e3), M.D4

for w in [0.1465 - 0.0513j, 0.25 + 0.01j]:
    print(f"\nomega = {w}")
    print(f"{'nm':>5} {'cond(D4)':>11}   nearest eigenvalues to 0.32 and to 0.55")
    for nm in [100, 130, 150, 180, 220]:
        ev, D4 = spec_nm(nm, w)
        a = ev[np.argsort(np.abs(ev - 0.32))[:3]]
        b = ev[np.argsort(np.abs(ev - 0.553 - 0.255j))[:1]]
        print(f"{nm:5d} {cond(D4):11.3e}   " + "  ".join(f"{z:+.6f}" for z in a)
              + f"   |   {b[0]:+.6f}")
