import numpy as np
from scipy.linalg import eig

# -------------------------------------------------
# SETUP
# -------------------------------------------------
Re = 2000.0
beta = 0.0 + 0.0j
num_modes = 100      
start = 0
terminate = 1
v_g = 0.0 + 0.0j

y_colloc_points = np.array([np.cos((j - 1) * np.pi / (num_modes - 1)) for j in range(1, num_modes + 1)])
y_colloc_points_new = ((start + terminate) / 2) - y_colloc_points * ((terminate - start) / 2)

D0_static = np.zeros((num_modes, num_modes))
for j in range(1, num_modes + 1):
    D0_static[:, j - 1] = np.cos((j - 1) * np.arccos(y_colloc_points))

D1_static_V2 = np.zeros((num_modes, num_modes))
D2_static_V2 = np.zeros((num_modes, num_modes))
D3_static_V2 = np.zeros((num_modes, num_modes))
D4_static_V2 = np.zeros((num_modes, num_modes))
D0_static_V2 = D0_static

D1_static_V2[:, 1] = D0_static[:, 0]
D1_static_V2[:, 2] = 4 * D0_static[:, 1]
D2_static_V2[:, 2] = 4 * D0_static[:, 0]

for j in range(4, num_modes + 1):
    D1_static_V2[:, j-1] = 2*(j-1)*D0_static_V2[:, j-2] + (j-1)*D1_static_V2[:, j-3]/(j-3)
    D2_static_V2[:, j-1] = 2*(j-1)*D1_static_V2[:, j-2] + (j-1)*D2_static_V2[:, j-3]/(j-3)
    D3_static_V2[:, j-1] = 2*(j-1)*D2_static_V2[:, j-2] + (j-1)*D3_static_V2[:, j-3]/(j-3)
    D4_static_V2[:, j-1] = 2*(j-1)*D3_static_V2[:, j-2] + (j-1)*D4_static_V2[:, j-3]/(j-3)

D1_static_V3 = D1_static_V2 / (-(terminate - start) / 2)**1
D2_static_V3 = D2_static_V2 / (-(terminate - start) / 2)**2
D3_static_V3 = D3_static_V2 / (-(terminate - start) / 2)**3
D4_static_V3 = D4_static_V2 / (-(terminate - start) / 2)**4

D0 = D0_static
D1 = D1_static_V3
D2 = D2_static_V3
D3 = D3_static_V3
D4 = D4_static_V3

u = y_colloc_points_new
d2u = 0.0

# u * ones(1, length(u))  and  d2u * ones(1, length(u))
u_ones = u[:, None] * np.ones((1, len(u)))
d2u_ones = d2u * np.ones((1, len(u)))


# -------------------------------------------------
# SPATIAL EIGENVALUE PROBLEM
# -------------------------------------------------

def couetteflow_alpha(omega):

        A11 = -2 * 1j * omega * D1 - 4 / Re * D3 + 4 / Re * beta**2 * D1 - 1j * u_ones * D2 + 1j * beta**2 * u_ones * D0 + 1j * d2u_ones * D0 - 1j * v_g * D2 + 1j * v_g * beta**2 * D0
        A12 = 1j * omega * D2 - 1j * omega * beta**2 * D0 + 1 / Re * D4 - 2 / Re * beta**2 * D2 + 1 / Re * beta**4 * D0
        A11 = np.vstack([np.zeros((2, num_modes), complex),                    A11[2:num_modes - 2, :],    np.zeros((2, num_modes), complex)])
        A12 = np.vstack([-200j * np.vstack([D0[0:1, :], D1[0:1, :]]),   A12[2:num_modes - 2, :],    -200j * np.vstack([D1[num_modes-1:num_modes, :], D0[num_modes-1:num_modes, :]])])
        A21 = 1 * np.eye(num_modes, dtype=complex)
        A22 = np.zeros((num_modes, num_modes), complex)
        A = np.block([[A11, A12], [A21, A22]])
        B11 = - 4 / Re * D2 - 2 * 1j * u_ones * D1 + 2 * 1j * v_g * D1
        B11 = np.vstack([np.zeros((2, num_modes), complex),  B11[2:num_modes - 2, :],    np.zeros((2, num_modes), complex)])
        B12 = np.zeros((num_modes, num_modes), complex)
        B21 = np.zeros((num_modes, num_modes), complex)
        B22 = 1 * np.eye(num_modes, dtype=complex)
        B = np.block([[B11, B12], [B21, B22]])

        eigvals, eigvecs = eig(A, B)

        return eigvals, eigvecs

def alpha_spectrum(omega):

    eigvals, _ = couetteflow_alpha(omega)

    mask = (
        np.isfinite(eigvals.real)
        & np.isfinite(eigvals.imag)
    )

    return eigvals[mask]