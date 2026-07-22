# vib_ribbon — causal spatial-branch solver (prototype 1)

Vibrating-ribbon (signalling) problem for plane Couette flow, built on the bare-bones
`briggsv4.jl` Couette physics and the model-swap pattern of `briggs_test.jl`.

The ribbon is modelled as the NEAT paper's **interior time-periodic point forcing**

    f'(x, y, t) = δ(x) δ(y − y') e^(−i ω_f t),

with the channel walls keeping their homogeneous boundary conditions. Because the forcing
enters only the right-hand side, the **homogeneous Orr–Sommerfeld eigenvalue problems are
unchanged** — the ribbon just fixes ω to a real forcing frequency ω_f and asks for the
**causal spatial branches** α⁺(ω_f) (the waves that feed x > 0, downstream) and α⁻(ω_f)
(x < 0, upstream).

This prototype computes **only those branches**. Forced-response amplitudes/residues are not
computed yet, but the Couette spatial solve returns eigenvectors (top-N block of the companion
vector = the physical ṽ amplitude) so residues can be added later without restructuring.

Design rationale and the full comparison with the older contour-relaxation scripts are in
`../../vib_ribbon/DESIGN_NOTES.md` and `../../vib_ribbon/SPEC.md` at the repository root.

---

## What it does (the method, in one paragraph)

For a fixed real ω_f, start at ω = ω_f + i·σ₀ with σ₀ high enough that the two relevant
spatial roots sit on opposite sides of the real α-axis (σ₀ "above all relevant temporal
singularities"). At that height the spatial contour F_α is the real axis, so the causal side of
each root is simply the sign of Im α: the least-damped root with Im α > 0 is α⁺, the
least-damped with Im α < 0 is α⁻. Then lower Im ω toward 0 in fixed steps, at each step
re-solving the spatial spectrum and following the two leaders by continuity (nearest previous
value, kept distinct). If the upper leader collides with a root from the lower family while
Im ω is still > 0, a genuine Briggs pinch obstructs the descent — report `PINCH` and stop.
Otherwise the leaders reach ω = ω_f + i·0⁺ and are returned as α±(ω_f). No moving F-contour, no
potentials, no pseudo-time: for Couette the spatial roots never approach the real α-axis, so F
stays real and the whole calculation is a clean vertical continuation.

## Files

| file | responsibility |
|---|---|
| `models.jl` | Physics layer. `AnalyticModel` (quadratic test) + `CouetteModel` (OSE matrices copied from `briggsv4.jl`). Each provides `spatial_spectrum(model, ω)` and `temporal_spectrum(model, α)`; analytic models also `exact_roots`. This is the swappable interface. |
| `ribbon_core.jl` | Model-agnostic Briggs logic: `physical_spatial_roots` (spurious filter), `classify_leaders` (sign-of-Im), `descend` (the causal continuation + pinch test), `sweep` (ω_f loop with warm/cold checks). Knows nothing about Couette. |
| `run_classify.jl` | Driver for one ω_f. Runs the three validation checks below and prints α⁺/α⁻. |
| `run_sweep.jl` | Driver for the ω_f sweep; saves JSON to `results/`. |
| `results/` | Sweep output (one JSON array of per-ω_f records). |

## How to run

From this directory:

```
julia run_classify.jl      # single-frequency validation (all three cases)
julia run_sweep.jl         # frequency sweep for Case S (analytic) and Couette
```

`run_classify.jl` builds the N = 100 and N = 150 Couette operators (a few seconds) for the
two-resolution spurious filter. Everything is serial — no `addprocs`, no Distributed.

## Validation cases and expected output

**Case S — stable analytic** `C1 = −1 − 0.96i, C2 = 0.4 + 0.8i, C3 = 1 − 0.5i`.
Globally temporally stable (max Im ω over real α = −0.64), saddle Im ω₀ = −1.04 < 0.
Descent completes (`status = OK`) at every real ω_f and the returned α⁺/α⁻ match the closed-form
`exact_roots` to ≈ 1e-10, independent of σ₀ and the number of steps (path independence = the
causality check). Reference roots: ω_f = 0.5 → α⁺ = 1.04239 + 0.25618i, α⁻ = −1.04239 − 1.05618i.

**Case U — absolutely unstable analytic** `C1 = −1 + 0.96i, C2 = 0.4 + 0.8i, C3 = 1`
(the existing `briggs_test.jl` fixture). Saddle ω₀ = −0.88 + 0.8i (Im > 0). Descending the line
Re ω = Re ω₀ = −0.88 from σ₀ = 3.0 must return `status = PINCH` with ω_stop = −0.88 + 0.8i and
α⁺, α⁻ coinciding at α₀ = −0.2 − 0.4i. For an off-centre ω_f the descent instead ends
`NO_STRADDLE` (α⁺ has crossed the axis; both roots end below it) — the signature that the flow
is absolutely unstable and the signalling problem is ill-posed there.
(σ₀ must be high enough that the roots straddle at the top; for this fixture that means
σ₀ ≳ 1.3, so the driver uses 3.0.)

**Couette — Re = 2000, β = 0, U = y**, ω_f = 0.5. `status = OK`; least-damped
α⁺ = 0.7496 + 0.1016i (N = 100; a consistent N = 150 value with a small resolution delta),
least-damped α⁻ = 0.059 − 3.275i. Cross-check: the frozen arc value in
`../contour_iteration_v3.1.json` is α_u = 0.7367 + 0.1006i at ω = 0.49, consistent after the
dα_r/dω_r ≈ 1.3 offset. The sweep α⁺₁(ω_f) is smooth from 0.194 + 0.043i (ω_f = 0.1) to
1.254 + 0.137i (ω_f = 0.9); the min axis clearance stays well above 0 (closest ≈ 0.04 near
ω_f = 0.1), which documents that F may safely stay on the real axis.

The analytic descent logic and all of these numbers were checked independently (a Python mirror
of `descend`) before the Julia was written.

## Key parameters (in `descend`)

- `sigma0` — starting Im ω. Must be above all relevant temporal singularities AND high enough
  that the two branches straddle the real α-axis at the top.
- `nsteps` — number of descent steps (result must be independent of this; increase if a step
  produces a large `jump`).
- `pinch_tol`, `pinch_im_floor` — a pinch is declared if the α⁺/α⁻ gap drops below `pinch_tol`
  while Im ω > `pinch_im_floor`.
- `model_coarse` — second-resolution model for the spurious-mode filter (Couette only).
- `cutoff`, `rtol` (in `physical_spatial_roots`) — drop roots with |α| > cutoff or that fail
  the two-resolution agreement test.

## Status codes from `descend`

- `OK` — leaders tracked continuously to the real axis and still straddle it → causal α±(ω_f).
- `PINCH` — upper/lower branches collided at Im ω > 0 → genuine Briggs pinch obstructs the
  descent (absolute instability); `omega_stop` reports where.
- `NO_STRADDLE` — reached the axis but α⁺ crossed to the lower half-plane → absolutely unstable
  regime, signalling problem ill-posed at this ω_f.

## What is deliberately NOT here yet

Forced-response amplitudes and residues 1/(∂D/∂α); F_α deformation and repulsion potentials
(only needed if a spatial root ever approaches the real α-axis — monitored via `min_clear`);
oblique β ≠ 0 / complex-β extended-Squire runs; moving-frame saddles via `v_g`;
parallelisation over ω_f. Bring the older `briggsv4` / `contour_test_v3.1` contour machinery
back only if `min_clear` approaches 0.

## Next steps

1. Run `run_classify.jl`; confirm Case S exact match, Case U pinch, Couette reference numbers.
2. Run `run_sweep.jl`; confirm smooth Couette curves and warm/cold agreement.
3. Then: add residue/amplitude evaluation (phase 2) using the carried eigenvectors.
