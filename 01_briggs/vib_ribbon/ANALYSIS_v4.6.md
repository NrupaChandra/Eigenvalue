# v4.6 ribbon run — analysis of `contour_iteration_v4.6_ribbon.json` (1001 frames)

Everything below was checked against an **independently written Python spatial-OS solver**
(`os_spatial.py`, built from the Julia source's `:alpha_collocation` block only) plus a
streaming re-analysis of every one of the 1001 frames. Scripts are in the outputs folder:
`scan_v46.py`, `os_spatial.py`, `test_flip.py`, `rect_mono.py`, `resolution.py`, `conv.py`,
`make_figs.py`. Figures: `results/diag_v46_*.png`.

---

## 0. What the ribbon problem actually asks for

From Görtz's dissertation, §7.2.4 (*Harmonic source with finite extent leading to convective
instabilities*) — this is the theory section for our case:

- A ribbon forcing at real frequency ω₀ enters the Laplace transform as
  `s_t(ω) = 1 / (i(ω − ω₀))`: a **fixed, immovable simple pole on the real ω axis**.
- Provided there is **no absolute instability**, L can be lowered through the real axis,
  but causality requires the pole at ω₀ to stay **enclosed**. That is the only reason the
  Hankel/ribbon shape exists — horizontals below the axis, risers, semicircle over ω₀.
- Closing L in the lower half plane for t > 0, the response is the **residue at ω₀**:

  `P(x, t, y) ≈ P̂(x, ω₀, y) e^{iω₀t}` — harmonic in time, and spatial behaviour set by
  `α_u(ω₀)` for `x > d` and `α_l(ω₀)` for `x < −d` (his Eq. 7.46 & 7.49).
- The stability verdict is a **sign test**, not a collision test: as ω is lowered from L′
  onto the real axis, a branch that crosses the real α axis **from above** (ending with
  `Im α < 0`) is spatially amplifying downstream; one crossing **from below** (ending with
  `Im α > 0`) is amplifying upstream.

**Consequences for v4.6, and you are right about this:**

1. The pinch is not the target. It only appears as a *precondition* — you have to establish
   that no pinch with `Im ω > 0` occurs while lowering L. That is a one-off check, not the
   termination criterion of the loop.
2. **The answer lives at fixed real ω** — on the arc and on the real-axis ends of the risers.
   Those ω points never move. So lowering the horizontals further and further changes
   nothing about the ribbon answer. This is not a hypothesis; the run confirms it exactly
   (see §2).
3. `pinch_tol = 1e-4` on `branch_distance` is therefore the wrong stopping rule for this
   script, and it is also unreachable here (§4).

---

## 1. The tracking failure — located, and it is *not* the detour

Your reading of the video is right in its symptom and slightly off in its location.

**Where it happens.** The upper branch changes discontinuously between **iteration 13 and
iteration 14**, when the horizontal level passes `ω_i = −0.0473 → −0.0513`. Index-by-index
comparison of the two frames (contour index → ω_r on the left horizontal):

| idx | ω_r | α_u, it 13 | α_u, it 14 | \|Δ\| |
|---|---|---|---|---|
| 48 (seed) | 0.12121 | +0.2560 −0.0381i | +0.2593 −0.0470i | 0.0095 |
| 56 | 0.14141 | +0.2962 −0.0133i | +0.3051 −0.0206i | 0.0115 |
| **58** | **0.14646** | +0.2999 −0.0034i | +0.3230 −0.0038i | **0.0230** |
| 66 | 0.16667 | +0.3121 +0.0033i | +0.4043 +0.0089i | 0.0924 |
| 94 (end of horizontal) | 0.23737 | +0.4027 +0.0089i | +0.6440 +0.0771i | **0.2507** |
| 95 (riser bottom) | 0.24000 | +0.4064 +0.0091i | +0.6536 +0.0807i | 0.2574 |

The divergence starts at **index 58, ω_r ≈ 0.146 — on the left horizontal, well upstream of
the detour** (which lives at ω_r ∈ [0.24, 0.26]). By the time the branch reaches the riser it
is already on the wrong curve; the riser, arc and right horizontal merely carry the error.
So `loop_mono_u` switching on at it 14 is a **consequence**, not the cause — it compares
α at the two riser bottoms, and the left one had already moved.

**It is not the side-of-F filter.** Re-running the same continuation with the filter removed
gives bit-identical results (`test_flip.py`): pure continuation, side-filtered continuation
and the recorded run agree to 4 decimals at every node, with **zero** side fallbacks. So the
`(B)` fix from v4.5 is not implicated.

**It is not a real branch point either.** A closed rectangle in the ω-plane spanning the two
horizontal levels, ω_r ∈ [0.121, 0.24], ω_i ∈ [−0.0513, −0.0473]:

| nodes per horizontal leg | monodromy \|α_end − α_start\| |
|---|---|
| 48 (the run's own spacing, Δω_r = 0.00253) | **0.0844** |
| 170 (Δω_r = 0.00071) | **0.00000** |

The loop closes exactly once it is resolved. **There is no branch point of α(ω) between the
two levels.** The apparent monodromy is the under-resolution itself.

**It is a mode hop caused by under-resolving the horizontal.** Continuing α_u from the seed
to ω_r = 0.24 at each level, varying only the number of nodes (`resolution.py`):

| nodes | Δω_r | α(0.24) at ω_i = −0.04730 (it 13) | α(0.24) at ω_i = −0.05130 (it 14) |
|---|---|---|---|
| **48** (v4.6) | 0.00253 | +0.40636 +0.00911i | **+0.65360 +0.08051i** ← the run |
| 96 | 0.00125 | +0.40636 +0.00911i | **+0.40647 +0.00375i** |
| 192 | 0.00062 | +0.40636 +0.00911i | **+0.40647 +0.00375i** |

The it-13 value is converged; the it-14 value **is not**, and the converged answer
(0.40647 + 0.00375i) is the smooth continuation of the it-13 value — exactly what a
Δω_i = 0.004 step should give. **Iterations 14 → 1001 (98.7 % of the run) are on the wrong
α branch.**

Note the max accepted step at it 14 was **0.0151** against `da_max = 0.03`, and `sep_frac`
never tripped. Every guard in the tracker passed. Refining 4× does **not** reduce the max
step (0.0109 at 48 nodes vs 0.0109 at 192 on the it-13 level) — the step there is dominated
by hopping between near-degenerate eigenvalues, not by the branch's own variation.

### Why that region is hostile

Between ω_r ≈ 0.145 and 0.16 the upper branch runs almost exactly **along the real α axis at
α ≈ 0.30** (|Im α| ≲ 0.005) — i.e. straight through F, into a cluster of nearby eigenvalues.
Refining the sampling there produces *more* zero-crossings of Im α, not fewer (3 crossings at
48 nodes, 9 at 192): the curve is not smooth at that scale, it is noise.

And the noise has an identifiable source. The Chebyshev derivative matrices are built by the
coefficient-space recursion, so their entries grow like n⁸:

| num_modes | max\|D2\| | max\|D4\| |
|---|---|---|
| 100 | 1.3e8 | 1.4e15 |
| **150** (used) | 6.6e8 | **3.7e16** |
| 200 | 2.1e9 | 3.8e17 |

`eps · max|D4| ≈ 8`. Re-solving the spectrum near α ≈ 0.31 at ω = 0.1465 − 0.0513i:

| num_modes | the two eigenvalues near α = 0.31 |
|---|---|
| 100 | +0.32319 −0.00368i, +0.30341 +0.00218i |
| 130 | +0.31822 −0.00383i, +0.30897 +0.00259i |
| 150 | +0.32076 −0.01568i, +0.30594 +0.01590i |
| 180 | +0.32934 −0.01368i, +0.29502 +0.00658i |
| 220 | garbage (0.288, 0.357) |

**The uncertainty in these eigenvalues (~0.005–0.015) is the same size as their separation
(~0.015–0.023).** In that window the discretisation cannot resolve which branch is which —
no continuation scheme, however careful, can fix that. (Elsewhere it is fine: at the arc,
ω = 0.25 + 0.01i, the values are stable to ~3e-3 across num_modes 100–180.)

My solver agrees with the Julia one to ~4e-4 at the arc, consistent with this picture — not
the 1e-7 quoted for the earlier Python port in the README.

---

## 2. The ribbon answer was frozen from iteration 14 and never changed again

`α_u` at the arc top (ω = ω_f + i·r_arc) over the whole run takes **exactly two values**:

| iterations | α_u(ω_f) | α_l(ω_f) |
|---|---|---|
| 1 – 13 | +0.413445 +0.088846i | +0.040086 −3.303930i |
| 14 – 1001 | +0.552821 +0.255033i | +0.040086 −3.303930i |

Bit-identical within each block. This is correct behaviour — the arc sits at fixed ω, so the
spectrum there cannot change — but it means **987 of the 1001 iterations produced no new
information about the ribbon problem at all**, and the one change that did occur is the
tracking artefact of §1.

Combined with §0.2: the ribbon deliverable was already available at iteration ~2.

**Reading off the physics from the correct block (iterations ≤ 13):**
`α_u(ω_f) = 0.4134 + 0.0889i` → Im α > 0, has **not** crossed the real α axis → no downstream
spatial growth. `α_l(ω_f) = 0.0401 − 3.3039i` → Im α < 0 throughout, never crosses → no
upstream growth. Across all 1001 frames, `Im α_l > 0` at **0 of 570** contour points, always.

So: **no convective instability at ω_f = 0.25 for plane Couette at Re = 2000, β = 0** — which
is the expected answer (plane Couette is linearly stable at every Re), and it is reassuring
that the correct branch gives it.

---

## 3. The lower branch is not a pinch partner and never was

At ω = 0.12121 − 0.0473i the spectrum has 278 finite eigenvalues, **190 with Im α > 0 and 88
with Im α < 0**. Sorted by |α|, the ones with Im α < 0 are:

```
+0.257264 −0.038010i     <- this is the eigenvalue the tracker calls alpha_u
+0.025439 −3.325319i     <- the next one down
+0.024666 −6.565481i
+0.024635 −9.789097i          ... a near-vertical wall-mode ladder at Re(alpha) ~ 0.025
```

There is **nothing between them**. Once F has been pushed below the α_u curve, the "nearest
eigenvalue below F" is forced onto the −3.3i rung, and stays there for the entire run
(`α_l(ω_f)` is constant to all digits for 1001 frames).

That is why `branch_distance` starts at **3.31** and only creeps to **2.19** after 1000
iterations, against `pinch_tol = 1e-4` — five orders of magnitude away, decreasing at
~1e-3 per iteration. The loop cannot terminate on its own criterion; the 1000-iteration cap
is what stopped it.

This is a structural statement about the problem, not a bug: in this ω range there is no α
branch sitting just below the real axis, so the two branches that Briggs needs — one
approaching F from each side — **do not both exist on this contour**.

---

## 4. Two further things the frames show

**F is being dragged bodily downwards.** `Im F` goes 0 (it 1) → −0.084…0.000 (it 13) →
**−1.172…−1.126 (it 1000)**: by the end, the entire inverse-Fourier contour is a near-flat
line at Im α ≈ −1.15. In Briggs' construction F must remain a *deformation* of the real α
axis — homotopic to it, indenting around individual poles. A rigid translation of the whole
contour 1.15 below the axis has swept across a large part of the spectrum, changing which
poles count for x > 0 and which for x < 0. The deformation rule appears to keep F clear of
the two *tracked* branches only, not of the rest of the spectrum. See
`diag_v46_4_F_dragged.png`.

**The stall is confirmed and is exactly as the README describes.** From iteration ≈ 22 on,
`ω_i` of L and `max Im ω_F` are within ~1e-3 of each other and descend together
(`diag_v46_2_timeseries.png`, bottom right). L is riding on its own lower bound; the repair
step is doing the work, not the potential.

---

## 5. Suggestions for v5

Ordered by how much they change the answer.

1. **Fix the discretisation before anything else.** The n⁸ growth of the coefficient-space
   D3/D4 recursion is the root cause of §1. Options: (a) physical-space Chebyshev
   differentiation matrices (Weideman & Reddy `chebdif`), (b) a basis that satisfies the
   clamped BCs identically, e.g. `(1−η²)² T_k(η)`, which removes the boundary rows and the
   `−200i` trick as well, (c) drop `num_modes` to ~80–100, where D4 is 1e15 rather than
   4e16, and check convergence. Until the spectrum near α ≈ 0.3 is resolved to better than
   the branch separation, no tracker result in that window means anything.
2. **Add a resolution check to the tracker, not just a step check.** The existing guards
   (`da_max`, `sep_frac`) all passed while the branch hopped. A cheap decisive test: recompute
   each frame's branch on a 2× refined ω grid and compare the endpoint. If they disagree,
   the frame is not converged. That single check would have caught iteration 14.
3. **Change the target to the ribbon question.** Replace `pinch_tol` termination with:
   (i) a separate, one-off absolute-instability check; (ii) descend the horizontals only to
   just below the real axis; (iii) report `α_u(ω_f)`, `α_l(ω_f)` and the crossing history of
   `Im α = 0` along each branch. `branch_distance` becomes a diagnostic, not a stopping rule.
4. **Constrain F to a genuine deformation of the real axis.** Pin the endpoints, and forbid
   the mean of `Im F` from drifting; only local indentations around the tracked poles should
   be allowed.
5. **Re-check the v4.5 conclusion.** The "branch point at ω ≈ 0.25 − 0.05i inside the detour"
   was inferred from `loop_mono` switching on at ω_i = −0.051. In this run the same threshold
   is produced by an unresolved hop at ω_r ≈ 0.146 — nowhere near the detour — and the
   corresponding rectangle loop closes exactly once refined. The v4.5 finding may be the same
   artefact seen from a different angle; it is worth repeating that two-road test with the
   horizontal refined 4× before treating the monodromy as physical.
6. On the F smoothing already flagged in `README_v4_series.md`: with `radius=7` every
   iteration for 1000 iterations, and F reaching a nearly straight line at Im α ≈ −1.15, the
   smoothing looks like a plausible contributor to item 4. Worth testing by simply switching
   it off for a short run.
