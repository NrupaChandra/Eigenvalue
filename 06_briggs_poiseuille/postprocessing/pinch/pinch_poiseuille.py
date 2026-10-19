"""
pinch_poiseuille.py -- independent check of the pinch points of the Poiseuille_v1.x.jl problem.

Same discretisation as the Julia file (Chebyshev coefficients, D-matrix recursion,
-200i boundary rows, diagonal equilibration, exp(-alpha*y) quadratic for the spatial
problem):  U = 1 - y^2 on [-1, 1], Re = 5000, beta = 0.

What it does
  1. Zero-frequency saddle.  On alpha = -i*kappa the temporal eigenvalues are purely
     imaginary (omega = i*sigma) or come in pairs.  The top sigma(kappa) has a minimum;
     there d omega / d alpha = 0, i.e. a double root of the spatial problem at omega = i*sigma_min.
  2. TS saddle.  Quadratic-fit iteration on the TS temporal branch (continued from alpha = 1).
  3. For both: check the double root in the spatial solver, the square-root law
     d = 2 sqrt(2 |d omega| / |omega''|), and the pinch test (follow the two roots up to
     omega_i = +1e-3 at fixed omega_r; one must end above, one below the real alpha axis).

Usage:  python pinch_poiseuille.py            (150 modes, ~2 min)
        python pinch_poiseuille.py 60         (60 modes, same digits shown, ~15 s)
"""
import sys
import numpy as np
import scipy.linalg as sla

RE = 5000.0
NM = int(sys.argv[1]) if len(sys.argv) > 1 else 150


def build(nm=NM, re=RE):
    x = np.cos(np.arange(nm) * np.pi / (nm - 1))
    y = -x                                   # start = -1, terminate = 1
    D0 = np.cos(np.outer(np.arccos(np.clip(x, -1, 1)), np.arange(nm)))
    D = [D0] + [np.zeros((nm, nm)) for _ in range(4)]
    D[1][:, 1] = D0[:, 0]; D[1][:, 2] = 4 * D0[:, 1]; D[2][:, 2] = 4 * D0[:, 0]
    for jj in range(4, nm + 1):              # same recursion as the Julia file
        j = jj - 1
        for k in range(1, 5):
            D[k][:, j] = 2 * (jj - 1) * D[k - 1][:, j - 1] + (jj - 1) * D[k][:, j - 2] / (jj - 3)
    for k in range(1, 5):
        D[k] = D[k] / (-1.0) ** k            # (-(terminate-start)/2)^k
    return dict(nm=nm, re=re, D=D, u=1 - y ** 2, d2u=-2 * np.ones(nm))


def equilibrate(*mats, iters=12):
    W = sum(np.abs(M) for M in mats); n = W.shape[0]
    R = np.ones(n); S = np.ones(n)
    for _ in range(iters):
        r = np.sqrt(np.maximum(W.max(1), 1e-300)); R /= r; W = W / r[:, None]
        c = np.sqrt(np.maximum(W.max(0), 1e-300)); S /= c; W = W / c[None, :]
    return R, S


def bc_rows(P, scale):
    D0, D1 = P['D'][0], P['D'][1]; n = P['nm']
    return scale * np.vstack([D0[:1], D1[:1]]), scale * np.vstack([D1[n - 1:n], D0[n - 1:n]])


def temporal(alpha, P):
    D0, D1, D2, D3, D4 = P['D']; n = P['nm']; re = P['re']
    U = P['u'][:, None]; U2 = P['d2u'][:, None]; k2 = alpha ** 2
    A = (-1j * alpha * U * D2 + 1j * alpha * U * k2 * D0 + 1j * alpha * U2 * D0
         + D4 / re - 2 / re * k2 * D2 + k2 ** 2 / re * D0)
    B = -1j * D2 + 1j * k2 * D0
    ta, ba = bc_rows(P, -200j); tb, bb = bc_rows(P, 1.0)
    A = np.vstack([ta, A[2:n - 2], ba]); B = np.vstack([tb, B[2:n - 2], bb]).astype(complex)
    R, S = equilibrate(A, B)
    ev = sla.eigvals(R[:, None] * A * S, R[:, None] * B * S)
    return ev[np.isfinite(ev)]


def spatial(omega, P):
    D0, D1, D2, D3, D4 = P['D']; n = P['nm']; re = P['re']
    U = P['u'][:, None]; U2 = P['d2u'][:, None]; Z2 = np.zeros((2, n))
    A11 = -2j * omega * D1 - 4 / re * D3 - 1j * U * D2 + 1j * U2 * D0
    A12 = 1j * omega * D2 + D4 / re
    B11 = -4 / re * D2 - 2j * U * D1
    ta, ba = bc_rows(P, -200j)
    A11 = np.vstack([Z2, A11[2:n - 2], Z2]).astype(complex)
    A12 = np.vstack([ta, A12[2:n - 2], ba]).astype(complex)
    B11 = np.vstack([Z2, B11[2:n - 2], Z2]).astype(complex)
    R, S = equilibrate(B11, A11, A12)
    A11, A12, B11 = (R[:, None] * M * S for M in (A11, A12, B11))
    I = np.eye(n); Z = np.zeros((n, n))
    ev = sla.eigvals(np.block([[A11, A12], [I, Z]]), np.block([[B11, Z], [Z, I]]))
    return ev[np.isfinite(ev)]


def follow_temporal(P, a0, w0, a1, nstep):
    w = w0
    for t in np.linspace(0, 1, nstep + 1)[1:]:
        ev = temporal(a0 + t * (a1 - a0), P); w = ev[np.argmin(abs(ev - w))]
    return w


def pinch_test(P, alpha_s, omega_s, label):
    ev = spatial(omega_s, P); o = np.argsort(abs(ev - alpha_s))
    print(f'  {label}: spatial roots at omega_s nearest alpha_s: {ev[o[0]]:.8f}, {ev[o[1]]:.8f} '
          f'(separation {abs(ev[o[0]] - ev[o[1]]):.2e})')
    ev = spatial(omega_s + 1e-3j, P); o = np.argsort(abs(ev - alpha_s)); r = ev[o[:2]].copy()
    print(f'    at omega_s + 1e-3 i the pair is {abs(r[0] - r[1]):.5f} apart')
    for wi in np.linspace(omega_s.imag + 1e-3, 1e-3, 400):
        ev = spatial(omega_s.real + 1j * wi, P)
        r = np.array([ev[np.argmin(abs(ev - x))] for x in r])
    sides = ['above' if x.imag > 0 else 'below' for x in r]
    print(f'    followed up to omega_i = +1e-3 at omega_r = {omega_s.real:.6f}: '
          f'{r[0]:.5f} ({sides[0]}), {r[1]:.5f} ({sides[1]})  ->  '
          + ('GENUINE alpha+/alpha- pinch' if sides[0] != sides[1] else 'NOT a pinch (same side)'))


def main():
    P = build()
    print(f'Poiseuille, Re = {RE:g}, beta = 0, {NM} Chebyshev modes')

    # ---------------- 1. zero-frequency saddle on the imaginary alpha axis
    def top_sigma(k):
        ev = temporal(-1j * k, P)
        return np.max(ev[np.abs(ev.real) < 1e-6].imag)
    k0 = 0.30
    for h in (0.02, 0.005, 0.001, 0.0002):
        kk = np.array([k0 - h, k0, k0 + h]); ss = np.array([top_sigma(x) for x in kk])
        c = np.polyfit(kk, ss, 2); k0 = -c[1] / (2 * c[0])
    a0 = -1j * k0; w0 = 1j * top_sigma(k0); w2 = -1j * 2 * c[0]
    print(f'\n1. zero-frequency saddle:  alpha_0 = {a0.imag:+.8f} i   omega_0 = {w0.imag:+.8f} i   '
          f"omega'' = {w2.imag:+.4f} i")
    pinch_test(P, a0, w0, 'zero-frequency')

    # ---------------- 2. TS saddle (quadratic fit on the TS branch)
    ev = temporal(1.0, P); c = ev / 1.0
    m = (c.real > 0.15) & (c.real < 0.45); w = ev[m][np.argmax(ev[m].imag)]
    a = 0.62 - 1.18j; w = follow_temporal(P, 1.0 + 0j, w, a, 70); h = 0.02
    for _ in range(40):
        wp = follow_temporal(P, a, w, a + h, 2); wm = follow_temporal(P, a, w, a - h, 2)
        d1 = (wp - wm) / (2 * h); d2 = (wp - 2 * w + wm) / h ** 2
        da = -d1 / d2
        if abs(da) > 0.15:
            da *= 0.15 / abs(da)
        w = follow_temporal(P, a, w, a + da, max(2, int(abs(da) / 0.02))); a += da
        if abs(da) < 1e-11:
            break
    print(f"\n2. TS saddle:  alpha_s = {a:.8f}   omega_s = {w:.8f}   omega'' = {d2:.4f}")
    pinch_test(P, a, w, 'TS')
    print(f'\n   Im(omega_0) - Im(omega_s) = {w0.imag - w.imag:+.4f}  (the larger Im omega is the absolute frequency)')


if __name__ == '__main__':
    main()
