# Plot one iteration of contour_iteration_v8.jsonl: omega plane (left), alpha plane (right).
# Usage:  python plot_v8_frames.py 25      (no number -> last iteration)
import json, sys
import matplotlib.pyplot as plt

LOG = "contour_iteration_v8.jsonl"

entries = []
for line in open(LOG):
    try:
        entries.append(json.loads(line))
    except json.JSONDecodeError:      # partial last line from a running job
        pass

it = int(sys.argv[1]) if len(sys.argv) > 1 else entries[-1]["iteration"]
e = next(x for x in entries if x["iteration"] == it)

def cx(key):                          # [{"re":..,"im":..}, ...] -> (re, im)
    arr = e.get(key) or []
    return [p["re"] for p in arr], [p["im"] for p in arr]

fig, (ax_w, ax_a) = plt.subplots(1, 2, figsize=(14, 6))

# omega plane (left)
ax_w.plot(*cx("L"), "o-", ms=3, lw=0.8, label="L")
ax_w.plot(*cx("omega_F"), "s-", ms=3, lw=0.8, label="omega_F")
ax_w.set_xlabel(r"$\omega_r$"); ax_w.set_ylabel(r"$\omega_i$")
ax_w.set_title(f"omega plane, iteration {it}")
ax_w.legend(fontsize=8); ax_w.grid(alpha=0.3)

# alpha plane (right)
if "cloud_u" not in e:
    print("this log has no clouds (written before cloud logging was added)")
ax_a.plot(*cx("cloud_u"), ".", ms=2, alpha=0.4, color="tab:blue", label="cloud upper")
ax_a.plot(*cx("cloud_l"), ".", ms=2, alpha=0.4, color="tab:red", label="cloud lower")
ax_a.plot(*cx("F"), "k.-", ms=3, lw=0.8, label="F")
ax_a.plot(*cx("alpha_L_u"), "-", lw=1.5, color="tab:blue", label="alpha_L_u (frontier)")
ax_a.plot(*cx("alpha_L_l"), "-", lw=1.5, color="tab:red", label="alpha_L_l (frontier)")
ax_a.plot(*cx("pair_u"), "x", ms=4, color="navy", label="pair_u")
ax_a.plot(*cx("pair_l"), "x", ms=4, color="darkred", label="pair_l")
ax_a.set_xlabel(r"$\alpha_r$"); ax_a.set_ylabel(r"$\alpha_i$")
ax_a.set_title(f"alpha plane, iteration {it}, d_branch = {e['d_branch']}")
ax_a.legend(fontsize=8); ax_a.grid(alpha=0.3)

fig.tight_layout()
plt.show()
