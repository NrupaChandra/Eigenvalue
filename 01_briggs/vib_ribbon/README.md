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

- `src/` — `briggsv4_ribbon.jl`, `briggsv4.1_ribbon.jl`, `briggsv4.2_ribbon.jl`,
  `briggsv4.3_ribbon.jl`
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

Output: `contour_iteration_v4.3_ribbon.json`. Render with
`postprocessing/plot_contour_deformation_ribbon.m`.

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
