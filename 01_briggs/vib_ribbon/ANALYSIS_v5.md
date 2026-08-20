# v5 run — why refining the contour was not enough, and what fixes it

The first v5 run (N = 300, n_L = 954, Δω_r = 8.36e−4) still hopped, at the same
iteration and onto the same eigenvalue as v4.6. The fixed-ω check caught it and stopped
the run, which is what it is for — but my resolution diagnosis was incomplete. This is the
corrected analysis, from the 15-frame v5 JSON plus the independent Python solver
(`analysis/os_spatial.py`).

---

## 1. What the v5 run did

| frame | ω_i | fixedw_u | sideclass_u | α_u(ω_f) |
|---|---|---|---|---|
| 1–13 | 0 → −0.04579 | 0 | 0 | 0.413445 + 0.088846i |
| 14 | −0.04979 | 0 | **9.76e−2** | 0.413445 + 0.088846i |
| 15 | −0.05379 | **2.27e−1** | 2.58e−1 | **0.552821 + 0.255033i** → STOP |

Note the levels: the break is between **ω_i = −0.04979 and −0.05379**. My earlier tests
used −0.04730 and −0.05130, so the interval that actually fails was never tested. That is
the first reason the fix missed.

## 2. Neither frame jumps, and both are real eigenvalues

Max step along the left horizontal: **0.0037** (frame 14) and **0.0039** (frame 15),
against `DA_MAX = 0.012`. Zero side-fallbacks, zero flagged nodes. The two curves do not
jump apart — they **drift** apart, smoothly, from the seed onwards:

| contour idx | ω_r | α_u frame 14 | α_u frame 15 | \|diff\| |
|---|---|---|---|---|
| 148 (seed) | 0.12375 | +0.26367 −0.04141i | +0.26725 −0.05057i | 0.0098 |
| 174 | 0.14548 | +0.30621 −0.00401i | +0.32422 −0.01805i | 0.023 |
| 182 | 0.15217 | +0.30121 +0.00014i | +0.35233 −0.00913i | 0.052 |
| 203 | 0.16973 | +0.31548 +0.00084i | +0.41641 +0.00215i | 0.101 |
| 287 (ω_r = 0.24) | 0.24000 | +0.40660 +0.00568i | +0.65700 +0.07157i | **0.259** |

Spot-checking against the independent solver confirms both are genuine: at
ω = 0.20903 − 0.04979i the only eigenvalue near frame 14's value is 0.364293 + 0.003088i
(frame 14 has 0.36446 + 0.00303i); at ω = 0.23913 − 0.04979i there is exactly **one**
eigenvalue in the box \|Im α\| < 0.06, 0.2 < Re α < 0.8, namely 0.405307 + 0.005249i.

## 3. There is no branch point either

Separation of the two branches involved, over ω_r ∈ [0.12, 0.20] × ω_i ∈ [−0.046, −0.060]:

```
  om_i \ om_r   0.120   0.130   0.140   0.150   0.160   0.170   0.180   0.190   0.200
     -0.0480   0.0860  0.0701  0.0456  0.0385  0.0766  0.1001  0.0867  0.0899  0.0893
     -0.0500   0.0863  0.0704  0.0441  0.0361  0.0739  0.0954  0.0863  0.0887  0.0902
     -0.0540   0.0878  0.0717  0.0481  0.0405  0.0771  0.1015  0.0863  0.0881  0.0907
     -0.0580   0.0904  0.0768  0.0576  0.0551  0.0814  0.0810  0.0833  0.0875  0.0904
```

Minimum **0.036**, at ω ≈ 0.150 − 0.050i. They never collide, so there is no branch point
and the continuation is not ambiguous in principle. Consistent with the earlier
closed-rectangle test (monodromy 0.00000 once refined).

## 4. Which frame is right — settled by continuing in ω_i

Continue α vertically at fixed ω_r from frame 14's level down to frame 15's:

| ω_r | frame 14 value | continued to ω_i = −0.05379 | frame 15 value | verdict |
|---|---|---|---|---|
| 0.12375 (seed) | +0.26367 −0.04141i | **+0.26815 −0.05078i** | +0.26725 −0.05057i | agrees (9e−4) |
| 0.15050 | +0.30095 +0.00000i | **+0.30329 +0.00108i** | +0.34548 −0.01060i | differs 0.044 |
| 0.24000 | +0.40660 +0.00568i | **+0.40725 −0.00015i** | +0.65700 +0.07157i | differs 0.260 |

The vertical at ω_r = 0.24 has a max step of **0.0008** — utterly unambiguous.

**So the seed is right and the endpoint is wrong.** Frame 15 drifts onto the neighbouring
branch between the seed and ω_r ≈ 0.15, gradually, with every individual step comfortably
inside the guard. Refining Δω_r cannot stop it, because there is no large step to catch —
which is why N = 100 → 300 changed nothing.

Note also frame 14's curve nearly stalls in that window (steps of 0.0001–0.0005 around
ω_r = 0.150, α ≈ 0.301, Im α ≈ 0) while frame 15's moves normally at 0.0035. That window
is where the branch runs along the real α axis, on F, and where the discretisation noise
(0.005–0.015, from `analysis/conv.py`) is comparable to the 0.036 branch separation.

---

## 5. The fix: continue in ω_i, not only along the contour

The tracker had one piece of information and was not using it. Every frame re-seeded from
scratch (`classify_by_F` at the interior node) and tracked along the contour with **no
memory of the previous frame**. But between frames every node moves by at most
\|Δω_i\| ≈ 4e−3, and the branch moves by ≤ 1.1e−2 — against a branch separation of
≥ 3.6e−2. A factor-3 margin. So the previous frame's α **at the same node index** is a far
better predictor than the along-contour secant.

`continue_branch` now takes `prev` and uses `prev[j]` as the predictor when it is finite.
This is also the physically right notion: Briggs' method is about deforming the contour and
following the poles, so the branch should be defined by continuity in the deformation
parameter.

**Verified on the frame-15 data** (`analysis/` scripts, same nodes, same spectra):

| predictor | α at ω_r = 0.24, ω_i = −0.05379 | error vs truth |
|---|---|---|
| along-contour secant (as shipped) | +0.656941 +0.071618i | **0.2598** |
| previous frame (`CHANGE 5`) | **+0.406956 +0.000022i** | **0.00034** |

and the max deviation from the previous frame drops from 0.259 to 0.017.

Continuing on from that corrected endpoint, up the riser and round the arc:
riser top +0.401638 + 0.074318i, arc apex **+0.414044 + 0.089296i** against frame 1's
+0.413445 + 0.088846i — agreement to 7.5e−4, which is this solver's noise floor. In the
Julia run the arc node's ω is bit-identical between frames, so it will be exact.

### Keeping the check able to fail

`prev` is applied **only on the horizontals** (`mask_to_horizontals`). The risers and the
arc are deliberately left on the secant. Their ω is pinned, so if `prev` drove them too,
α there would reproduce frame 1 *by construction* and `fixedw_u/l` could never fail —
the exact trap `README_v4_series.md` warns about ("a good check is one that *can* fail").
As it stands they are reached by continuation from the prev-anchored end of the horizontal,
so `fixedw` still tests the thing that actually broke.

`framedev_u/l` (max \|α_new − α_prev\| on the horizontals) is reported with a WARN above
`DFRAME_MAX = 0.02`. That one is a **report, not a test** — it can be masked if the true
branch genuinely moves a long way in one ω_i step.

---

## 6. Residual risk, stated plainly

The branch is now *defined* by continuation from frame 1. If frame 1 is right and every
step is resolved, the answer is right. The checks verify "every step is resolved"
(`fixedw`, `gridconv`, `framedev`); nothing verifies "frame 1 is right" — that rests on the
spectrum being well separated at ω_i = 0, which it is.

Still not addressed (needs the rewrite in `ANALYSIS_v4.6.md` §5): the n⁸ conditioning of
the coefficient-space D3/D4 recursion, which is what makes the ω_r ≈ 0.15 window marginal
in the first place; F drifting bodily downwards; the radius-7 smoothing of F.
