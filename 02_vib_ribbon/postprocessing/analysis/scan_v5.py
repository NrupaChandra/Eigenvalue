import json, numpy as np
P="/sessions/magical-upbeat-hamilton/mnt/Eigenvalue/01_briggs/vib_ribbon/results/contour_iteration_v5_ribbon.json"
d=json.load(open(P))
print("frames:", len(d), " keys:", sorted(d[0].keys()))
cx=lambda L: np.array([c["re"] for c in L])+1j*np.array([c["im"] for c in L])
n_hl,n_ver,n_arc,n_hr = 287,160,60,287
print(f"n_L = {len(d[0]['L'])}  (expect 954)")
np.savez_compressed("v5.npz",
    **{f"au_{f['iteration']}":cx(f["alpha_L_u"]) for f in d},
    **{f"al_{f['iteration']}":cx(f["alpha_L_l"]) for f in d},
    **{f"F_{f['iteration']}":cx(f["F"]) for f in d},
    **{f"L_{f['iteration']}":cx(f["L"]) for f in d})
print(f"{'it':>4}{'om_i':>10}{'fixedw_u':>11}{'sideclass_u':>12}{'mono_u':>10}{'nfallb_u':>9}{'a_u(wf)':>22}")
for f in d:
    a=f["alpha_u_wf"]
    print(f"{f['iteration']:4d}{f['omega_i']:10.5f}{f['fixedw_u']:11.3e}{f['sideclass_u']:12.3e}"
          f"{f['loop_mono_u']:10.3e}{f['nfallb_u']:9d}   {a['re']:+.6f}{a['im']:+.6f}j")
