import mpmath as mp

# -------------------------------------------------
# WORKING PRECISION
#
# Every number in this module is an mpmath arbitrary-precision value.
# There is no double-precision (numpy/scipy) arithmetic anywhere here.
#
# HIGH_PREC_DPS must agree with the value used by compute_pinch.py; that
# script calls set_precision() below at import time so that the two files
# can never drift apart.
# -------------------------------------------------
HIGH_PREC_DPS = 60
mp.mp.dps = HIGH_PREC_DPS

I_UNIT = mp.mpc(0, 1)          # exact imaginary unit, replaces the literal 1j
BC_SHIFT = mp.mpc(0, -200)     # replaces the literal -200j of the boundary rows


# -------------------------------------------------
# INFINITE EIGENVALUES
#
# The pencil (A, B) of the spatial problem is singular: B has four
# identically zero rows (the boundary-condition rows), so some alpha
# eigenvalues are infinite.  scipy's QZ reports those as inf/nan, which the
# original code removed with np.isfinite.
#
# mpmath has no QZ, so the problem is solved in the reversed form
#
#       B v = mu A v ,      mu = 1 / alpha
#
# where the infinite alpha become mu = 0.  In exact arithmetic those mu are
# exactly zero; at finite working precision they come out at the rounding
# level instead.  An eigenvalue is therefore treated as mu = 0, i.e.
# alpha = infinity, when
#
#       |mu| <= INFINITE_ALPHA_TOL * max|mu|
#
# The default (1e-30 at 60 digits) sits in a very wide empty gap: the
# rounding-level mu are around 1e-48, while the smallest physically
# meaningful |mu| is around 1e-9.  Rebuilt by set_precision().
# -------------------------------------------------
INFINITE_ALPHA_TOL = mp.mpf(10) ** (-(HIGH_PREC_DPS // 2))


# -------------------------------------------------
# SETUP
# -------------------------------------------------
Re = mp.mpf("2000")
beta = mp.mpc(0, 0)
num_modes = 100
start = mp.mpf(0)
terminate = mp.mpf(1)
v_g = mp.mpc(0, 0)


def build_operators():
    """
    Build the Chebyshev collocation operators at the current working
    precision and publish them as module-level names, exactly as the
    original script did at import time.
    """

    global y_colloc_points, y_colloc_points_new
    global D0_static, D0_static_V2
    global D1_static_V2, D2_static_V2, D3_static_V2, D4_static_V2
    global D1_static_V3, D2_static_V3, D3_static_V3, D4_static_V3
    global D0, D1, D2, D3, D4
    global u, d2u, u_ones, d2u_ones

    n = num_modes

    y_colloc_points = [
        mp.cos(mp.mpf(j - 1) * mp.pi / (n - 1))
        for j in range(1, n + 1)
    ]

    half_span = (terminate - start) / 2

    y_colloc_points_new = [
        ((start + terminate) / 2) - yy * half_span
        for yy in y_colloc_points
    ]

    D0_static = mp.matrix(n, n)
    for r in range(n):
        for j in range(1, n + 1):
            D0_static[r, j - 1] = mp.cos((j - 1) * mp.acos(y_colloc_points[r]))

    D1_static_V2 = mp.matrix(n, n)
    D2_static_V2 = mp.matrix(n, n)
    D3_static_V2 = mp.matrix(n, n)
    D4_static_V2 = mp.matrix(n, n)
    D0_static_V2 = D0_static

    for r in range(n):
        D1_static_V2[r, 1] = D0_static[r, 0]
        D1_static_V2[r, 2] = 4 * D0_static[r, 1]
        D2_static_V2[r, 2] = 4 * D0_static[r, 0]

    for j in range(4, n + 1):
        factor = mp.mpf(j - 1) / (j - 3)

        for r in range(n):
            D1_static_V2[r, j-1] = 2*(j-1)*D0_static_V2[r, j-2] + factor*D1_static_V2[r, j-3]
            D2_static_V2[r, j-1] = 2*(j-1)*D1_static_V2[r, j-2] + factor*D2_static_V2[r, j-3]
            D3_static_V2[r, j-1] = 2*(j-1)*D2_static_V2[r, j-2] + factor*D3_static_V2[r, j-3]
            D4_static_V2[r, j-1] = 2*(j-1)*D3_static_V2[r, j-2] + factor*D4_static_V2[r, j-3]

    # Map [-1, 1] -> [start, terminate]
    scale = -(terminate - start) / 2

    D1_static_V3 = mp.matrix(n, n)
    D2_static_V3 = mp.matrix(n, n)
    D3_static_V3 = mp.matrix(n, n)
    D4_static_V3 = mp.matrix(n, n)

    for r in range(n):
        for j in range(n):
            D1_static_V3[r, j] = D1_static_V2[r, j] / scale**1
            D2_static_V3[r, j] = D2_static_V2[r, j] / scale**2
            D3_static_V3[r, j] = D3_static_V2[r, j] / scale**3
            D4_static_V3[r, j] = D4_static_V2[r, j] / scale**4

    D0 = D0_static
    D1 = D1_static_V3
    D2 = D2_static_V3
    D3 = D3_static_V3
    D4 = D4_static_V3

    u = y_colloc_points_new
    d2u = mp.mpf(0)

    # u * ones(1, length(u))  and  d2u * ones(1, length(u))
    u_ones = mp.matrix(n, n)
    d2u_ones = mp.matrix(n, n)

    for r in range(n):
        for j in range(n):
            u_ones[r, j] = u[r]
            d2u_ones[r, j] = d2u


def set_precision(dps):
    """
    Set the working precision of this module and rebuild the operators.

    compute_pinch.py calls this so that both files always run at the same
    number of digits.  A repeated call with the same value is a no-op.
    """

    global HIGH_PREC_DPS, INFINITE_ALPHA_TOL

    if dps == HIGH_PREC_DPS and mp.mp.dps == dps:
        return

    HIGH_PREC_DPS = dps
    mp.mp.dps = dps
    INFINITE_ALPHA_TOL = mp.mpf(10) ** (-(HIGH_PREC_DPS // 2))

    build_operators()


build_operators()


# -------------------------------------------------
# SPATIAL EIGENVALUE PROBLEM
# -------------------------------------------------

def couetteflow_alpha(omega, eigenvectors=True):
    """
    Return the alpha eigenvalues (and, unless eigenvectors=False, the right
    eigenvectors) of the spatial problem at the given omega.

    Infinite eigenvalues are returned as mpc(inf, inf) so that mp.isfinite()
    removes exactly the modes that np.isfinite() removed before.

    eigenvectors=False skips the eigenvector accumulation, which is by far
    the most expensive part of the solve.  alpha_spectrum() discards the
    eigenvectors, so it uses that path.
    """

    n = num_modes

    w = mp.mpc(omega)
    b2 = beta**2
    b4 = beta**4
    inv_Re = 1 / Re

    A = mp.matrix(2*n, 2*n)
    B = mp.matrix(2*n, 2*n)

    for r in range(n):

        # A21 = eye(num_modes) , B22 = eye(num_modes)
        A[n + r, r] = mp.mpf(1)
        B[n + r, n + r] = mp.mpf(1)

        # Rows 0, 1, n-2, n-1 of A11 and B11 are zero (they carry the
        # boundary conditions through A12 instead), so only the interior
        # rows are assembled.
        interior = (2 <= r < n - 2)

        for j in range(n):

            if interior:
                A[r, j] = (
                    -2 * I_UNIT * w * D1[r, j]
                    - 4 * inv_Re * D3[r, j]
                    + 4 * inv_Re * b2 * D1[r, j]
                    - I_UNIT * u_ones[r, j] * D2[r, j]
                    + I_UNIT * b2 * u_ones[r, j] * D0[r, j]
                    + I_UNIT * d2u_ones[r, j] * D0[r, j]
                    - I_UNIT * v_g * D2[r, j]
                    + I_UNIT * v_g * b2 * D0[r, j]
                )

                B[r, j] = (
                    -4 * inv_Re * D2[r, j]
                    - 2 * I_UNIT * u_ones[r, j] * D1[r, j]
                    + 2 * I_UNIT * v_g * D1[r, j]
                )

            # A12 block
            A[r, n + j] = (
                I_UNIT * w * D2[r, j]
                - I_UNIT * w * b2 * D0[r, j]
                + inv_Re * D4[r, j]
                - 2 * inv_Re * b2 * D2[r, j]
                + inv_Re * b4 * D0[r, j]
            )

    # Boundary rows of the A12 block
    for j in range(n):
        A[0, n + j] = BC_SHIFT * D0[0, j]
        A[1, n + j] = BC_SHIFT * D1[0, j]
        A[n - 2, n + j] = BC_SHIFT * D1[n - 1, j]
        A[n - 1, n + j] = BC_SHIFT * D0[n - 1, j]

    # ---------------------------------------------
    # Reversed pencil:  B v = mu A v ,  mu = 1/alpha.
    # A is non-singular, B is not, so this is the orientation that keeps
    # the solve well posed and sends the infinite alpha to mu = 0.
    # ---------------------------------------------
    C = A**-1 * B

    if eigenvectors:
        mu, eigvecs = mp.eig(C, left=False, right=True)
    else:
        mu = mp.eig(C, left=False, right=False)
        eigvecs = None

    mu_scale = max(abs(m) for m in mu)
    cutoff = INFINITE_ALPHA_TOL * mu_scale

    eigvals = [
        mp.mpc(mp.inf, mp.inf) if abs(m) <= cutoff else 1 / m
        for m in mu
    ]

    return eigvals, eigvecs


def alpha_spectrum(omega):

    eigvals, _ = couetteflow_alpha(omega, eigenvectors=False)

    return [z for z in eigvals if mp.isfinite(z)]
