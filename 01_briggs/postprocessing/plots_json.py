import json
import matplotlib.pyplot as plt


json_file = "contour_iteration_algebraic_testv6.json"
with open(json_file, "r") as f:
    data = json.load(f)

# Chose which frame i.e which iteration -1 means the last frame :D 
frame = data[-150]
iteration = frame.get("iteration", len(data))


def get_xy(points):
    x = [p["re"] for p in points]
    y = [p["im"] for p in points]
    return x, y

# Exact pinch point for the algebraic test case
alpha_exact = (-0.2, -0.4)   # alpha_r, alpha_i
omega_exact = (-0.88, 0.8)   # omega_r, omega_i

L_x, L_y = get_xy(frame["L"])
omegaF_x, omegaF_y = get_xy(frame["omega_F"])
F_x, F_y = get_xy(frame["F"])
alpha_u_x, alpha_u_y = get_xy(frame["alpha_L_u"])
alpha_l_x, alpha_l_y = get_xy(frame["alpha_L_l"])

fig, axes = plt.subplots(1, 2, figsize=(16, 7))


# Left plot: omega-plane

ax = axes[0]
ax.plot(L_x, L_y, "b.-", markersize=3, linewidth=1.5, label="L")
ax.plot(omegaF_x, omegaF_y, "r.-", markersize=3, linewidth=1.5, label=r"$\omega_F$")
ax.plot(
    omega_exact[0],
    omega_exact[1],
    "ko",
    markersize=7,
    label=r"exact $\omega$ pinch"
)

ax.set_title(r"$\omega$-plane, iteration {}".format(iteration))
ax.set_xlabel(r"$\omega_r$")
ax.set_ylabel(r"$\omega_i$")
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_aspect("equal", adjustable="box")

# Optional: set limits
ax.set_xlim(-1.25, -0.35)
ax.set_ylim(0.10, 1.55)


# Right plot: alpha-plane

ax = axes[1]
ax.plot(F_x, F_y, "b.-", markersize=3, linewidth=1.5, label="F")
ax.plot(alpha_u_x, alpha_u_y, "r.-", markersize=3, linewidth=1.5, label=r"$\alpha_L^u$")
ax.plot(alpha_l_x, alpha_l_y, "g.-", markersize=3, linewidth=1.5, label=r"$\alpha_L^l$")
ax.plot(
    alpha_exact[0],
    alpha_exact[1],
    "ko",
    markersize=7,
    label=r"exact $\alpha$ pinch"
)

ax.set_title(r"$\alpha$-plane, iteration {}".format(iteration))
ax.set_xlabel(r"$\alpha_r$")
ax.set_ylabel(r"$\alpha_i$")
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_aspect("equal", adjustable="box")

# Optional: set limits 
ax.set_xlim(-1.1, 0.65)
ax.set_ylim(-1.15, 0.40)


plt.tight_layout()

# Save figure if needed
plt.savefig("Intermediate2", dpi=300, bbox_inches="tight")

plt.show()