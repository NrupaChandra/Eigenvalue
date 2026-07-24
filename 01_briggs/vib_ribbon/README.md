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

- `src/` — `briggsv4_ribbon.jl`, `briggsv4.1_ribbon.jl`, `briggsv4.2_ribbon.jl`
- `postprocessing/` — `plot_contour_deformation_ribbon.m` (video renderer for the
  v4.2 JSON; no exact-pinch markers, overlays the descent branches in gray)
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
