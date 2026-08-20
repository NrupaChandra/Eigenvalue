"""Independent re-implementation of the spatial Orr-Sommerfeld eigenvalue
problem used by briggsv4.6_ribbon.jl (mode :alpha_collocation).

Plane Couette flow, u(y) = y on y in [0,1], Re = 2000, beta = 0, v_g = 0,
num_modes = 150 Chebyshev collocation points, clamped BCs phi = phi' = 0 at
both walls, quadratic-in-alpha problem linearised by companion form.

Written from the Julia source only (not translated line-by-line from any
existing Python port), so agreement is a real cross-check.
"""
import numpy as np
from scipy.linalg import eig

Re = 2000.0
beta = 0.0 + 0.0j
NM = 150
v_g = 0.0 + 0.0j
start, terminate = 0.0, 1.0


def _build_operators(nm=NM):
    yc = np.cos(np.arange(nm) * np.pi / (nm - 1))          # 1 -> -1
    y = (start + terminate) / 2 - yc * (terminate - start) / 2   # 0 -> 1

    D0 = np.zeros((nm, nm))
    for j in range(nm):
        D0[:, j] = np.cos(j * np.arccos(yc))               # T_j(yc)

    D1 = np.zeros((nm, nm)); D2 = np.zeros((nm, nm))
    D3 = np.zeros((nm, nm)); D4 = np.zeros((nm, nm))
    # T_0' = 0, T_1' = T_0, T_2' = 4 T_1 ; T_2'' = 4 T_0
    D1[:, 1] = D0[:, 0]
    D1[:, 2] = 4 * D0[:, 1]
    D2[:, 2] = 4 * D0[:, 0]
    for j in range(3, nm):          # Julia index j+1, so (j_jl - 1) = j
        n = j
        D1[:, j] = 2 * n * D0[:, j - 1] + n * D1[:, j - 2] / (n - 2)
        D2[:, j] = 2 * n * D1[:, j - 1] + n * D2[:, j - 2] / (n - 2)
        D3[:, j] = 2 * n * D2[:, j - 1] + n * D3[:, j - 2] / (n - 2)
        D4[:, j] = 2 * n * D3[:, j - 1] + n * D4[:, j - 2] / (n - 2)

    s = -(terminate - start) / 2
    return y, D0, D1 / s, D2 / s**2, D3 / s**3, D4 / s**4


Y, D0, D1, D2, D3, D4 = _build_operators()
U = Y[:, None] * np.ones((1, NM))       # u(y) broadcast over columns
D2U = 0.0
NB = 200.0                              # the -200i boundary-row scaling


def spatial_pencil(omega, nm=None):
    nm = NM if nm is None else nm
    """Return (A, B) of the linearised quadratic alpha-EVP, A x = alpha B x."""
    om = complex(omega)
    A11 = (-2j * om * D1 - 4 / Re * D3 + 4 / Re * beta**2 * D1
           - 1j * U * D2 + 1j * beta**2 * U * D0 + 1j * D2U * D0
           - 1j * v_g * D2 + 1j * v_g * beta**2 * D0).astype(complex)
    A12 = (1j * om * D2 - 1j * om * beta**2 * D0 + 1 / Re * D4
           - 2 / Re * beta**2 * D2 + 1 / Re * beta**4 * D0).astype(complex)
    B11 = (-4 / Re * D2 - 2j * U * D1 + 2j * v_g * D1).astype(complex)

    A11[:2, :] = 0.0; A11[-2:, :] = 0.0
    B11[:2, :] = 0.0; B11[-2:, :] = 0.0
    A12[0, :] = -NB * 1j * D0[0, :]
    A12[1, :] = -NB * 1j * D1[0, :]
    A12[-2, :] = -NB * 1j * D1[-1, :]
    A12[-1, :] = -NB * 1j * D0[-1, :]

    I = np.eye(nm, dtype=complex); Z = np.zeros((nm, nm), dtype=complex)
    A = np.block([[A11, A12], [I, Z]])
    B = np.block([[B11, Z], [Z, I]])
    return A, B


def spectrum(omega, cap=1.0e3):
    A, B = spatial_pencil(omega, NM)
    ev = eig(A, B, right=False)
    ev = ev[np.isfinite(ev.real) & np.isfinite(ev.imag)]
    return ev[np.abs(ev) < cap]


if __name__ == "__main__":
    import sys
    for w in [0.0, 0.25 + 0.01j, 0.1212121212]:
        ev = spectrum(w)
        print(f"omega = {w}  -> {len(ev)} finite eigenvalues")
        idx = np.argsort(np.abs(ev))
        for k in idx[:8]:
            print(f"   {ev[k]:.8f}")
