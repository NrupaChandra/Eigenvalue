# vib_ribbon — vibrating-ribbon (signalling) problem for plane Couette flow

The ribbon is modelled as the NEAT paper's **interior time-periodic point forcing**

    f'(x, y, t) = δ(x) δ(y − y') e^(−i ω_f t),

with the channel walls keeping their homogeneous boundary conditions. Because the forcing
enters only the right-hand side, the **homogeneous Orr–Sommerfeld eigenvalue problems are
unchanged** — the ribbon just fixes ω to a real forcing frequency ω_f.

> **Note:** an earlier standalone prototype (`models.jl`, `ribbon_core.jl`, `run_classify.jl`,
> `run_sweep.jl`, plus `DESIGN_NOTES.md`/`SPEC.md` at the repository root) was removed because
> it did not work. The current scripts are instead ribbon variants of the main `briggsv4.jl`
> contour-relaxation solver in `../src/`.

## Folder structure

- `src/` — `briggsv4_ribbon.jl` … `briggsv4.6_ribbon.jl`, `briggsv5_ribbon.jl`
- `analysis/` — post-mortem scripts for the v4.6 run (`ANALYSIS_v4.6.md` at folder level)
- `postprocessing/` — `plot_contour_deformation_ribbon.m` (video renderer for the
  v4.2/v4.3 JSON; no exact-pinch markers, overlays the descent branches in gray)
- `results/` — output JSON (`contour_iteration_v4.1_ribbon.json`) and animation
  (`contour_video4.1ribbon.mp4`; videos are git-ignored, local only)

## Scripts

**`briggsv4_ribbon.jl`** — `briggsv4.jl` with one structural change: the temporal integration
contour L is the piecewise vibrating-ribbon contour of the NEAT paper (Eq. 24) — horizontal
parts at height ω_i, vertical risers, and a semicircular detour above the real forcing
frequency ω_f. The arc stays pinned to the real axis while the horizontal parts descend.
Everything else (matrices, tracking, guards, λ = 4.0, σ = 3e-5, Δt = 1e-3, ζ_common = 4e-4,
pinch_tol, loop length) is identical to `briggsv4.jl`, apart from three mechanical fixes needed
because L now has n_L = 306 points instead of N = 100 (piecewise `contour_L()`, endpoint
indexing, pairwise-minimum branch distances). Output: `contour_iteration_v4_ribbon.json`.

**`briggsv4.1_ribbon.jl`** — keeps the `briggsv4.jl` descent byte-identical (so the pinch forms
and dUL closes exactly as in the working Couette run / `contour_video13`), with one additive
change: the *saved* per-frame contour is the NEAT ribbon (Hankel) contour that detours around
ω_f, with its branches recomputed by robust per-point classification purely for the animation.
ω_f (default 0.25) only positions the drawn detour; it does not steer the descent — the pinch
it finds is the flow's intrinsic pinch (Im ω₀ < 0 ⇒ not absolutely unstable).
Output: `contour_iteration_v4.1_ribbon.json`.

> **Known v4.1 issue:** the saved ribbon branches were classified *independently* at each
> contour point (no continuity), so they mode-hopped — on the vertical risers essentially
> every point jumped to a different eigenvalue (jumps up to ~1.5 in |α|). That is the noise
> visible in `contour_video4.1ribbon.mp4`. The tracked descent branches, which do close
> smoothly to the pinch (dUL → 1e-4), were never saved. The old videos also carried an
> "exact pinch" marker hard-coded for the algebraic test case — meaningless here (no
> analytical pinch point exists for Couette flow).

**`briggsv4.2_ribbon.jl`** — same descent (still byte-identical to `briggsv4.jl`), animation
fixed: (1) display branches are continuity-tracked along the Hankel contour (spatial spectra
computed in parallel, then serial nearest-neighbour continuation with a damped linear
predictor, seeded from the descent's tracked branches at the left contour end); (2) riser
resolution 25 → 60 points each (n_L = 376); (3) the descent's straight contour and tracked
branches are saved per frame (`L_descent`, `alpha_L_u_descent`, `alpha_L_l_descent`) so the
real pinch closing is visible; (4) JSON frames kept in memory instead of re-parsing the file
every iteration. Render with `postprocessing/plot_contour_deformation_ribbon.m` (no
exact-pinch marker; descent branches overlaid in gray).
Output: `contour_iteration_v4.2_ribbon.json`.

> **Known v4.2 issue:** the display branches were continued from a *single* seed at the left
> end of the contour by plain nearest-neighbour matching on a uniform grid. On the vertical
> risers the **upper** branch legitimately travels ≈ 2.7 in |α| across only 60 points
> (≈ 0.06 per step — the same order as the local spacing of the spatial OS spectrum), so the
> continuation hopped modes going up the left riser and never recovered. Decisive test:
> ω = ω_f − r_arc + 0i is the *same* physical point in every frame, so α there must be
> frame-independent — frame 1 (descent on the real axis) gives α_u = 0.4012 + 0.0740i, while
> from frame 30 on v4.2 reported 0.556…0.572 + 0.215…0.227i, a different eigenvalue. Coming
> back down the right riser on the wrong mode then corrupted the entire right horizontal
> (median error vs. the descent branch grew to 2.43) — the "upper curve drifts away from F".
> The lower branch moves only ≈ 0.2 over the same riser (max step 0.044) and survived, which
> is why only the upper branch looked broken.

**`briggsv4.3_ribbon.jl`** — descent still byte-identical to `briggsv4.jl`; **no theory
changed**, only the sampling and continuation numerics:

1. the horizontal parts of the ribbon contour are placed *on the descent grid* `omega_r`, so
   they are the same ω points the descent already solved, and their α values are taken
   directly from the descent's tracked branches instead of being re-classified (same
   eigenvalue problem, same branch rule, no extra solves) — the horizontals now agree with
   the plain Briggs run exactly;
2. only the detour (risers + arc) is continued, anchored at **both** ends: the left riser is
   tracked upward from the descent value at ω_f − r_arc, the right riser upward from the
   descent value at ω_f + r_arc, so the two sides are independent and the out-and-back loop
   closes;
3. riser resolution 60 → 160 points, arc 40 → 60 (nominal step drops 0.045 → 0.017 in |α|),
   with a secant predictor, step-limited corrector, a separation guard (winner must be
   clearly nearer than the runner-up) and adaptive bisection of every rejected interval,
   re-solved in one parallel batch, up to 4 passes; output is resampled back onto the fixed
   display grid so the saved frame size stays constant (n_L = 476);
4. per-frame closure diagnostic `arc_mismatch_u` / `arc_mismatch_l` — the gap between the arc
   end and the independently tracked top of the right riser. Small ⇒ the loop closed on one
   mode; large ⇒ the detour is still hopping;
5. iteration cap 1500 → 300. In the v4.2 run ω_i was frozen at −0.57249 from iteration ≈ 276
   on, with `d_branch` flat at 1.16e−1 (`pinch_tol` = 1e−4 never reached), so iterations
   277–1102 produced no new information.

Output: `contour_iteration_v4.3_ribbon.json`.

> **Known v4.3 issue (methodological, not numerical):** to guarantee the horizontals agreed
> with the descent, v4.3 *copied* the descent's tracked α values onto them instead of
> recomputing. The values were correct — same ω, same spectrum — but 96 of 476 points
> (**36.5 % of the drawn curve length**) were then identical to the descent branch by
> construction, so the grey overlay could not disagree with the red curve there and stopped
> being a check at all.

**`briggsv4.4_ribbon.jl`** — descent still byte-identical to `briggsv4.jl`. **Nothing is
copied:** every point of the ribbon contour gets its own spectrum and its own branch
decision, so the grey descent overlay becomes a genuine independent cross-check.

- display branches are continuity-tracked over the **whole** contour (horizontals included);
- seeds are **not** the descent's stored values — they are the side-of-F classification (the
  descent's own rule) recomputed from the freshly solved spectrum at each end;
- horizontal resolution doubled (`n_sub = 2`), dropping the median α step there from 0.033 to
  ≈ 0.017, below `da_max = 0.03`. The guard stays strict rather than being loosened;
- descent grid points remain a **subset** of the ribbon horizontals, so red-vs-grey is
  compared at identical ω with no interpolation. n_L = 570.

Four diagnostics per frame (printed and stored in the JSON):

| field | meaning |
|---|---|
| `xcheck_u/l` | tracked ribbon branch vs descent branch (red vs grey) — the cross-check, not forced to zero |
| `sideclass_u/l` | tracked branch vs side-of-F classification evaluated independently from the same spectrum |
| `fwdbwd_u/l` | forward traversal (left seed) vs backward traversal (right seed) over the same nodes |
| `loop_mono_u/l` | **monodromy** of the closed detour: up the left riser + arc + down the right riser, vs straight along the bottom. Nonzero ⇒ the loop encloses a branch point of α(ω) |

`xcheck` is deliberately *not* forced to zero. The descent's rule is a selection ("nearest
eigenvalue on the correct side of F"), not an analytic continuation, so the two can
legitimately part company where the descent switches modes — that is a result, not a bug.

Output: `contour_iteration_v4.4_ribbon.json`.

> **Known v4.4 issue:** the *tracking* was clean — max step 0.017 against a 0.03 guard, zero
> guard trips, zero refinements, 100 % of points on the upper side of F at every frame — but
> the branch it followed **drifted away from F**: median |α − F| went 0.035 → 0.786 → **2.079**
> (it = 11 / 61 / 301), against 0.035 → 0.114 → 0.612 for the descent branch. Briggs needs the
> branches *adjacent* to F, the ones that collide at the pinch; a branch 2 units away is a
> valid eigenvalue branch and irrelevant.
>
> Two causes. **The seed:** v4.4 seeded at the contour's left *end*, ω_r = 0. At purely
> imaginary ω the spectrum is symmetric under α → −ᾱ, so the seed chooses between **mirror
> twins** (verified to 7 digits: +0.7836 − 1.9284i vs the descent's −0.7836 − 1.9284i at
> it = 301). The tiebreak is "closer to F" and the margin was **2.6 %** at it = 21 — a coin
> flip, which flipped. **The corrector:** it took the nearest eigenvalue to the predictor with
> no side-of-F constraint, so nothing kept it adjacent once it had left.

**`briggsv4.5_ribbon.jl`** — fixes **A + B**, both of which apply the *existing* descent rule.
No value is copied; every α still comes from this script's own freshly solved spectra.

- **A — seed at the interior node the descent seeds at** (`start_index = N/4`, ω_r ≈ 0.121)
  instead of the contour end, then march outward in both directions — exactly the protocol of
  `contour_alpha_L_conti`. At that node the branch sits next to F and there is no mirror twin.
- **B — the corrector filters candidates by side of F** before taking the nearest to the
  predictor, i.e. the rule of `couetteflow_spatial_sing_mode_comparison`. As there, an empty
  side falls back to nearest-overall — and every fallback is **counted and reported**
  (`nfallb_u/l`), so the filter cannot silently swallow a branch that genuinely leaves F.

Diagnostics per frame: `xcheck_u/l` (red vs grey at identical ω, not forced to zero),
`sideclass_u/l`, `nfallb_u/l`, `loop_mono_u/l` (monodromy of the closed detour).

Output: `contour_iteration_v4.5_ribbon.json`.

> **v4.5 findings (both branches verified with an independent Python OS solver, matching the
> Julia spectra to 1e-7):** the tracking is clean — red ≡ grey to 1e-6 on the left horizontal,
> zero side-fallbacks, and pure continuation over the detour reproduces v4.5's values to
> 1.5e-6. The red/grey disagreement after the detour (up to 2.85) is **real monodromy**: a
> branch point of α(ω) sits inside the detour strip at ω ≈ 0.25 − 0.05i (`loop_mono_u`
> switches on at exactly ω_i = −0.051). The ribbon contour passes *over* it, the straight
> contour *under* it — two different sheets, both correct. Three sheets checked; only one
> pair swaps. **Consequence:** the straight-line loop and the ribbon loop are genuinely
> different problems, which motivates v4.6.

**`briggsv4.6_ribbon.jl`** — **the closed loop.** The ribbon contour is *inside* the
spatial↔temporal iteration: F → (temporal OS) → ω_F bounds the descent of L = ribbon contour
→ (spatial OS, tracked) → α_u, α_l → potential Φ → deforms F → … The ribbon's branches —
detour excursions and sheet change included — now steer the pinch search; ω_f is no longer
cosmetic. Descent *theory* unchanged (same Φ, dynamics, acceptance rules, pinch criterion);
what changed is geometry + mechanism:

1. `contour_L()` **is** `hankel_L(omega_i)`; trial contours are ribbon contours;
2. in-loop branches from the v4.5 tracker via `ribbon_branches()` (interior seed by the
   side-of-F rule, side filter, secant predictor, separation guard, adaptive bisection) —
   `contour_alpha_L_conti` is dead code (it was the old per-point matcher that fed
   mode-hopping noise into Φ in the abandoned `briggsv4_ribbon.jl`);
3. ω-level update evaluated on the ribbon's horizontal nodes (the level belongs to them);
4. Φ_F endpoint indexing `N` → `end` (branches have n_L = 570 points);
5. `branch_distance` is the **pairwise** minimum (elementwise pairing is meaningless on the
   ribbon); `pinch_tol` applies to that;
6. diagnostics per frame: `sideclass_u/l`, `nfallb_u/l`, `loop_mono_u/l`, refinement counts.
   No descent overlay — there is no straight-line run to compare against.

Same parameters as v4.5 (ω_f = 0.25, r_arc = 0.01, 300 iterations). Output:
`contour_iteration_v4.6_ribbon.json`. Render with
`postprocessing/plot_contour_deformation_ribbon.m`.

**`briggsv5_ribbon.jl`** — v4.6 with the four fixes from the v4.6 post-mortem
(`../ANALYSIS_v4.6.md`). **No physics changed**: same potential functions, same descent
dynamics, same acceptance rules, same branch definition, same ribbon geometry. What changed
is resolution, verification, and what the run reports.

The v4.6 diagnosis in one line: the upper branch hopped at iteration 14 on the **left
horizontal** at ω_r ≈ 0.146 — upstream of the detour — and iterations 14–1001 (98.7 % of the
run) were on the wrong branch. Not the side-of-F filter (removing it reproduces the same
values to 4 decimals, zero fallbacks) and not a branch point (the closed ω-rectangle spanning
the two levels has monodromy 0.0844 at the run's spacing and exactly 0.00000 once refined) —
just under-resolution:

| nodes on the left horizontal | α at ω = 0.24 − 0.05130i |
|---|---|
| 48 (v4.6, Δω_r = 2.53e−3) | 0.65360 + 0.08051i ← what v4.6 got |
| 96 (Δω_r = 1.25e−3) | 0.40647 + 0.00375i |
| 192 (Δω_r = 6.2e−4) | 0.40647 + 0.00375i (converged) |

1. **Resolution.** `N` 100 → 300, so the ribbon horizontals sit at Δω_r = 8.4e−4, three times
   finer than the 1.25e−3 at which the continuation converged. `n_L` 570 → 954, F 100 → 300.
   **`num_modes` stays at 150** — raising it makes things worse: the coefficient-space
   Chebyshev recursion gives max|D4| ~ n⁸ (1.4e15 at 100, 3.7e16 at 150, 3.8e17 at 200) and at
   220 the spectrum near α ≈ 0.31 is already garbage. `da_max` 0.03 → `DA_MAX = 0.012`, since
   at 0.03 the guard sat ~12× above the typical step and never fired (the hop was accepted
   with a step of 0.0151).
2. **Fixed-ω consistency check** (`fixedw_u/l`). The riser tops and the whole arc are at ω
   values that never move, so α there must be identical in every frame. Measured against
   frame 1; above `fixed_tol = 1e-8` the run stops, because a hop has occurred upstream and
   every later frame is meaningless. This is the v4.2 "decisive test" reinstated — it would
   have caught iteration 14 at once, and it is worth more than all four v4.6 diagnostics,
   every one of which passed while the branch hopped.
3. **Short run.** 1000 → `n_iterations = 60`. In v4.6 α at the arc took exactly **two** values
   over 1001 frames, because the arc sits at fixed ω. A long descent adds nothing.
4. **Ribbon reporting instead of pinch termination.** Per Görtz §7.2.4, for a harmonic source
   at real ω₀ the response is the residue at that pole — harmonic in time, with the spatial
   behaviour set by α_u(ω₀) for x > d and α_l(ω₀) for x < −d — and the verdict is a **sign
   test** on Im α, not a collision test. The run now prints α_u(ω_f), α_l(ω_f), the
   Im α = 0 crossing counts and the resulting verdict. `pinch_tol` is kept only as the
   **absolute-instability alarm**: a pinch with Im ω > 0 would invalidate the ribbon argument
   entirely. (In v4.6 `branch_distance` went 3.31 → 2.19 against 1e−4 — unreachable, because
   below F there is nothing but the wall ladder at −3.33i, −6.57i, …)

Plus an optional **grid-convergence check** (`gridconv_u/l`, every `grid_check_every = 10`
iterations): the same contour retracked at `grid_check_mult` × the sampling and compared at
the coincident nodes — the direct measure of the v4.6 failure mode.

New per-frame JSON fields: `fixedw_u/l`, `gridconv_u/l`, `alpha_u_wf`, `alpha_l_wf`,
`ncross_u/l`, `omega_i`, `n_L`. Output: `contour_iteration_v5_ribbon.json`.

> **Not fixed in v5** (needs a bigger rewrite, see `ANALYSIS_v4.6.md` §5): the n⁸ conditioning
> of the D3/D4 recursion, which leaves the eigenvalues near α ≈ 0.31 uncertain by about their
> own separation; F being dragged bodily down to Im α ≈ −1.15 instead of staying a deformation
> of the real axis; the radius-7 rolling average applied to F every iteration.

All use Re = 2000, β = 0, num_modes = 150 and run in parallel via `Distributed`
(`addprocs(31)` — reduce for local runs).

## How to run

From `src/`:

```
julia briggsv4_ribbon.jl
julia briggsv4.1_ribbon.jl
```

The output JSON is written to the working directory — move it to `../results/` afterwards.
Videos are rendered from the JSON with the MATLAB scripts in `../postprocessing/`.
