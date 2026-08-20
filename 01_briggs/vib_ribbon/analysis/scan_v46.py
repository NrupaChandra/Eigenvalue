#!/usr/bin/env python3
"""Streaming scan of contour_iteration_v4.6_ribbon.json (94 MB, 1001 frames).

Memory-safe: one frame decoded at a time via JSONDecoder.raw_decode over a
sliding buffer.  Reduces each frame to a small set of scalars + a few stored
full curves, written to summary_v46.npz.

Ribbon geometry (from briggsv4.6_ribbon.jl, N=100, omega_r=LinRange(0,0.5,100),
omega_f=0.25, r_arc=0.01, n_sub=2, n_ver=160, n_arc=60):
    idx   0.. 94  left horizontal   (omega_r 0 .. 0.23737, level omega_i)
    idx  95..254  left riser        (0.24 + i*linspace(omega_i, 0))   UP
    idx 255..314  arc               (0.25 + 0.01*exp(i*theta), theta pi->0)
    idx 315..474  right riser       (0.26 + i*linspace(0, omega_i))   DOWN
    idx 475..569  right horizontal  (omega_r 0.26263 .. 0.5, level omega_i)
"""
import json, numpy as np, sys

PATH = "/sessions/magical-upbeat-hamilton/mnt/Eigenvalue/01_briggs/vib_ribbon/results/contour_iteration_v4.6_ribbon.json"
OUT = "/sessions/magical-upbeat-hamilton/mnt/outputs/summary_v46.npz"

HL0, HL1 = 0, 94        # left horizontal
RL0, RL1 = 95, 254      # left riser (bottom -> top)
AR0, AR1 = 255, 314     # arc
RR0, RR1 = 315, 474     # right riser (top -> bottom)
HR0, HR1 = 475, 569     # right horizontal
I_ARC_TOP = (AR0 + AR1) // 2   # omega ~ 0.25 + 0.01i, the forcing frequency

KEEP_FRAMES = set([1, 2, 5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 25, 30,
                   50, 100, 200, 300, 500, 700, 1000])


def cx(lst):
    return np.array([c["re"] for c in lst]) + 1j * np.array([c["im"] for c in lst])


def frames(path):
    dec = json.JSONDecoder()
    buf = ""
    pos = 0
    with open(path, "r") as fh:
        # skip leading '['
        buf = fh.read(1 << 20)
        i = buf.index("[") + 1
        buf = buf[i:]
        while True:
            buf = buf.lstrip()
            while buf.startswith(","):
                buf = buf[1:].lstrip()
            if buf.startswith("]") or buf == "":
                return
            while True:
                try:
                    obj, end = dec.raw_decode(buf)
                    break
                except ValueError:
                    chunk = fh.read(1 << 20)
                    if not chunk:
                        return
                    buf += chunk
            yield obj
            buf = buf[end:]


rows = []
kept = {}
n = 0
for fr in frames(PATH):
    n += 1
    it = fr["iteration"]
    L = cx(fr["L"])
    au = cx(fr["alpha_L_u"])
    al = cx(fr["alpha_L_l"])
    F = cx(fr["F"])
    wF = cx(fr["omega_F"])

    om_i = L[0].imag
    du = np.abs(np.diff(au))
    dl = np.abs(np.diff(al))

    # distance of each branch point to the F contour
    dFu = np.abs(au[:, None] - F[None, :]).min(axis=1)
    dFl = np.abs(al[:, None] - F[None, :]).min(axis=1)

    # pairwise branch distance (the pinch measure used in v4.6)
    dpair = np.abs(au[:, None] - al[None, :]).min()

    def seg(a, s, e):
        return a[s:e + 1]

    rows.append(dict(
        it=it, om_i=om_i,
        omF_max=wF.imag.max(),
        # step statistics per segment, upper branch
        du_hl=seg(du, HL0, HL1 - 1).max(), du_rl=seg(du, RL0, RL1 - 1).max(),
        du_ar=seg(du, AR0, AR1 - 1).max(), du_rr=seg(du, RR0, RR1 - 1).max(),
        du_hr=seg(du, HR0, HR1 - 1).max(),
        du_max=du.max(), du_argmax=int(du.argmax()),
        dl_max=dl.max(), dl_argmax=int(dl.argmax()),
        # junction jumps (segment boundaries)
        j_hl_rl_u=abs(au[RL0] - au[HL1]), j_rl_ar_u=abs(au[AR0] - au[RL1]),
        j_ar_rr_u=abs(au[RR0] - au[AR1]), j_rr_hr_u=abs(au[HR0] - au[RR1]),
        # key point values, upper
        au_left=au[0], au_botL=au[RL0], au_topL=au[RL1], au_arc=au[I_ARC_TOP],
        au_topR=au[RR0], au_botR=au[RR1], au_right=au[-1],
        al_left=al[0], al_botL=al[RL0], al_topL=al[RL1], al_arc=al[I_ARC_TOP],
        al_topR=al[RR0], al_botR=al[RR1], al_right=al[-1],
        # left/right horizontal asymmetry: same omega_i level, both horizontals
        med_dFu=np.median(dFu), med_dFl=np.median(dFl),
        med_dFu_hl=np.median(dFu[HL0:HL1 + 1]), med_dFu_hr=np.median(dFu[HR0:HR1 + 1]),
        med_dFl_hl=np.median(dFl[HL0:HL1 + 1]), med_dFl_hr=np.median(dFl[HR0:HR1 + 1]),
        max_dFu=dFu.max(), max_dFl=dFl.max(),
        absu_hl=np.median(np.abs(au[HL0:HL1 + 1])), absu_hr=np.median(np.abs(au[HR0:HR1 + 1])),
        absl_hl=np.median(np.abs(al[HL0:HL1 + 1])), absl_hr=np.median(np.abs(al[HR0:HR1 + 1])),
        dpair=dpair,
        # sign of imag part -- the physically decisive quantity for the ribbon
        au_arc_im=au[I_ARC_TOP].imag, al_arc_im=al[I_ARC_TOP].imag,
        n_au_imneg=int((au.imag < 0).sum()), n_al_impos=int((al.imag > 0).sum()),
        # diagnostics from the run
        sideclass_u=fr["sideclass_u"], sideclass_l=fr["sideclass_l"],
        nfallb_u=fr["nfallb_u"], nfallb_l=fr["nfallb_l"],
        loop_mono_u=fr["loop_mono_u"], loop_mono_l=fr["loop_mono_l"],
        F_span=np.abs(F).max(),
    ))
    if it in KEEP_FRAMES:
        kept[f"au_{it}"] = au
        kept[f"al_{it}"] = al
        kept[f"F_{it}"] = F
        kept[f"L_{it}"] = L
        kept[f"wF_{it}"] = wF
    if n % 100 == 0:
        print(n, "frames", flush=True)

print("total frames:", n)
keys = [k for k in rows[0]]
arr = {k: np.array([r[k] for r in rows]) for k in keys}
np.savez_compressed(OUT, **arr, **kept)
print("saved", OUT)
