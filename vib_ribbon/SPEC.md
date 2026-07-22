# vib_ribbon — implementation spec (prototype 1: causal spatial branches)

Planning document. No solver code committed yet. This spec is the hand-off for the build stage.

Scope fixed by the working assumptions:
- Ribbon = interior time-periodic point forcing f'(x,y,t)=δ(x)δ(y−y')e^{−iω_f t}. Walls keep
  homogeneous BCs. → The homogeneous OSE eigenvalue problems are UNCHANGED from `briggsv4.jl`;
  the ribbon only fixes ω to a real forcing frequency ω_f and asks for the causal spatial roots.
- Prototype computes only α⁺(ω_f) and α⁻(ω_f). No forced-response amplitudes yet. Eigenvectors
  are returned by the spatial solve and carried along so residues/receptivity can be added later.
- Two analytic validation cases (existing unstable + new verified stable) plus Couette.
- Serial only. No Distributed, no contour relaxation, no potentials in prototype 1.

The ribbon setup (ω_f, the Hankel L geometry) is taken from `Briggs/contour_test_v3.1.jl`.
Numerically, the only part of that Hankel contour that matters is the vertical descent at
Re ω = ω_f; the full Γ1–Γ5 object is kept only as an analysis/figure helper, not in the solver
core.

---

## 1. Verified constants (independently derived; re-derive in code comments)

Convention (from `briggs_test.jl`): D(α,ω) = ω − (C1 + C2·α + C3·α²) = 0.
- ω(α) = C1 + C2 α + C3 α²
- spatial roots: α±(ω) = [ −C2 ± √(C2² − 4 C3 (C1 − ω)) ] / (2 C3)
- saddle / branch point: α₀ = −C2/(2C3),  ω₀ = C1 − C2²/(4C3)

Case U (EXISTING, absolutely unstable — keep as-is):
  C1 = −1 + 0.96i, C2 = 0.4 + 0.8i, C3 = 1
  α₀ = −0.2 − 0.4i, ω₀ = −0.88 + 0.8i  (Im ω₀ = +0.80 > 0 ⇒ absolutely unstable)
  Descent along Re ω = ω_f = Re ω₀ = −0.88 MUST hit the pinch at Im ω = 0.8.

Case S (NEW, stable — add as second fixture; verified numerically this session):
  C1 = −1 − 0.96i, C2 = 0.4 + 0.8i, C3 = 1 − 0.5i
  α₀ = 0.0 − 0.4i, ω₀ = −0.84 − 1.04i  (Im ω₀ = −1.04 < 0 ⇒ NOT absolutely unstable)
  max Im ω(α real) = −0.64 < 0 (globally temporally stable)
  Exact roots straddle the real α-axis at every real ω_f, e.g.
    ω_f = 0.0 : α⁺ = 0.85198 + 0.28547i, α⁻ = −0.85198 − 1.08547i
    ω_f = 0.5 : α⁺ = 1.04239 + 0.25618i, α⁻ = −1.04239 − 1.05618i
  Descent to Im ω = 0 completes with no pinch; roots match the closed form above.

Couette reference (Re=2000, β=0, U=y; from independent solve, M=100), ω_f = 0.5:
  least-damped α⁺ = 0.7496 + 0.1016i ; least-damped α⁻ = 0.059 − 3.275i
  (consistent with the frozen arc value 0.7367 + 0.1006i at ω = 0.49 in contour_iteration_v3.1.json)
  sweep smoothness: α⁺₁ from 0.194 + 0.043i (ω_f=0.1) to 1.254 + 0.137i (ω_f=0.9)

---

## 2. Minimal file structure

    vib_ribbon/
    ├── DESIGN_NOTES.md          (rationale; already written)
    ├── SPEC.md                  (this file)
    ├── models.jl                model layer: the two-function interface + 3 concrete models
    ├── ribbon_core.jl           model-agnostic Briggs logic: filter, classify, descent, pinch test
    ├── run_classify.jl          driver: one ω_f, one model → prints α⁺/α⁻, saves JSON
    ├── run_sweep.jl             driver: ω_f sweep with warm-start + cold-start checks
    └── results/                 one JSON per run

`ribbon_core.jl` never mentions Couette or C1/C2/C3. `models.jl` never mentions descent logic.
That seam is the whole point (it is what `briggs_test.jl` already proved by swapping
`algebraic_alpha_roots` in for `couetteflow`).

---

## 3. Model interface (the contract `ribbon_core.jl` depends on)

A model is a struct + two functions:

    spatial_spectrum(model, ω) -> (alphas::Vector{ComplexF64}, vecs::Union{Nothing,Matrix})
        all spatial roots α solving D(α,ω)=0 at fixed ω; eigenvectors optional (nothing for analytic)
    temporal_spectrum(model, α) -> alphas... (omegas::Vector{ComplexF64})
        all temporal roots ω at fixed α; used ONLY to choose σ0 and for diagnostics

    exact_roots(model, ω)  (analytic models only) -> (α⁺, α⁻)   # for validation

Three concrete models in `models.jl`:
- `AnalyticModel(C1,C2,C3)` — spatial = quadratic formula (2 roots), temporal = evaluate ω(α),
  exact_roots = closed form. Instantiated as Case U and Case S.
- `CouetteModel(Re,β,num_modes)` — wraps the `couetteflow` matrix block copied verbatim from
  `briggsv4.jl`; spatial = 2N companion pencil eigen, temporal = N pencil eigen; no exact_roots.

---

## 4. Function responsibilities

`models.jl`
- build Couette D-matrices once (copy of briggsv4 lines 12–89), behind `CouetteModel`.
- `spatial_spectrum`, `temporal_spectrum`, `exact_roots` per model.
- finite-mask helper (drop Inf/NaN eigenvalues from the singular pencil).

`ribbon_core.jl`
- `physical_spatial_roots(model, ω; cutoff)` — solve at two resolutions (for Couette:
  num_modes 100 & 150; analytic: single solve), keep roots agreeing to `rtol` within |α|<cutoff.
  Returns filtered roots (+ vecs from the finer solve if present).
- `classify_at_top(roots)` — at ω = ω_f + iσ0 the F_α contour is the real axis, so split by
  sign(Im α): upper set = α⁺ candidates, lower set = α⁻ candidates. Pick the causal leaders:
  α⁺ = min Im among upper, α⁻ = max Im among lower.
- `match_sets(old, new)` — mutual nearest-neighbour matching of two root sets between adjacent
  descent steps; returns permutation + max match distance; flags any unmatched root.
- `descend(model, ω_f; σ0, nsteps)` — the core loop (pseudocode §5); returns
  (α⁺, α⁻, status, ω_stop, min_gap, min_axis_clearance, trace).
- `pinch_detected(...)` — gap between the tracked α⁺ leader and the α⁻ set falls below
  `pinch_tol` while Im ω > `pinch_im_floor`.
- diagnostics: per-step record of n⁺, n⁻, min |Im α|, leaders, match distance.

`run_classify.jl`
- pick model + ω_f + σ0; call `descend`; print result table; dump trace JSON to results/.

`run_sweep.jl`
- loop ω_f over a grid; warm-start each from previous leaders; every K-th ω_f cold-start a full
  `descend` and assert leader agreement; save per-ω_f record (leaders, status, clearances,
  resolution delta). Plot Re/Im α⁺₁ vs ω_f as PNG.

---

## 5. Pseudocode

Classification at the top of the descent:
```
roots0, vecs0 = physical_spatial_roots(model, ω_f + i·σ0)
upper = roots0[Im .> 0];  lower = roots0[Im .< 0]      # F_α = real axis at σ0
alpha_plus  = upper[argmin(Im(upper))]                 # least-damped downstream
alpha_minus = lower[argmax(Im(lower))]                 # least-damped upstream
```

Fixed-ω_f causal descent:
```
function descend(model, ω_f; σ0, nsteps, pinch_tol, pinch_im_floor):
    σ_grid = range(σ0, 0, length=nsteps+1)
    roots  = physical_spatial_roots(model, ω_f + i·σ0)
    (αp, αm) = classify_at_top(roots)
    min_gap = Inf;  min_clear = min(abs.(Im.(roots)))
    for σ in σ_grid[2:end]:
        ω = ω_f + i·σ
        roots_new = physical_spatial_roots(model, ω)
        # continuity tracking of the two leaders (nearest to previous leader),
        # keeping them distinct (mask αp's pick before choosing αm) — same rule as
        # contour_test_v3.1.jl couetteflow_spatial_pair_comparison
        αp = nearest(roots_new, αp)
        αm = nearest(roots_new \ {αp}, αm)
        gap = min(abs.(αp .- (roots_new[Im .< Im(αp)])))   # αp vs opposite family
        min_gap  = min(min_gap, gap)
        min_clear = min(min_clear, min(abs.(Im.(roots_new))))
        if gap < pinch_tol and σ > pinch_im_floor:
            return (αp, αm, status="PINCH", ω_stop=ω, min_gap, min_clear, trace)
        # sanity: per-side counts should be conserved; large match jump ⇒ warn/refine
    return (αp, αm, status="OK", ω_stop=ω_f + i·0⁺, min_gap, min_clear, trace)
```

Outer ω_f sweep:
```
prev = nothing
for (k, ω_f) in enumerate(ω_f_grid):
    guess = prev                                   # warm start from previous frequency
    (αp, αm, status, ...) = descend(model, ω_f; seed=guess)
    if k % K == 0:                                 # periodic independent check
        (αp_c, αm_c, ...) = descend(model, ω_f)    # cold start
        assert abs(αp - αp_c) < tol and abs(αm - αm_c) < tol
    store(ω_f => (αp, αm, status, min_gap, min_clear, resdelta))
    prev = (αp, αm)
```

---

## 6. Validation criteria (fixed in advance)

Case S (stable analytic), any σ0 ≥ 0.2 and any nsteps ≥ 10:
- status == OK for all real ω_f in the sweep;
- returned (α⁺, α⁻) equal to `exact_roots` to ≤ 1e-10;
- result independent of σ0 and nsteps (path independence = the causality check);
- sweep curves α±(ω_f) smooth and matching the closed form point-by-point.

Case U (unstable analytic), ω_f = Re ω₀ = −0.88, σ0 > 0.8:
- status == PINCH, with ω_stop ≈ −0.88 + 0.8i (Im to ≤ 1e-3), i.e. the descent halts at the
  branch point instead of reaching the axis;
- min_gap → 0 there; α⁺ and α⁻ coincide at α₀ = −0.2 − 0.4i to ≤ 1e-3;
- for ω_f well away from −0.88, descent may complete but the leaders must NOT straddle the axis
  at ω = ω_f (both Im < 0), flagged as "no causal upstream branch — absolutely unstable regime".

Couette (Re=2000, β=0), ω_f = 0.5:
- status == OK; least-damped α⁺ = 0.7496 + 0.1016i (M=100), consistent M=150 value with a
  reported resolution delta; least-damped α⁻ = 0.059 − 3.275i;
- min_clear stays well above 0 across the sweep (closest approach ≈ 0.04 near ω_f = 0.1) —
  documents that F may stay real and no deformation is needed;
- warm-start and cold-start leaders agree at every checked ω_f;
- sweep α⁺₁(ω_f) smooth, endpoints matching 0.194 + 0.043i (0.1) and 1.254 + 0.137i (0.9).

Spurious-mode / resolution guard (all models with matrices):
- every reported root survives the 100-vs-150 agreement filter;
- per-side counts n⁺, n⁻ conserved along each descent (a jump ⇒ tracking failure, not physics).

---

## 7. Deferred (not in prototype 1)

Forced-response amplitudes and residues 1/(∂D/∂α) and source projection at y' (needs the
carried eigenvectors); F_α deformation + potentials (only if min_clear ever approaches 0, i.e. a
regime with spatial roots reaching the real axis); oblique β≠0 / complex-β extended-Squire runs;
moving-frame saddles via v_g; parallelization over ω_f. The Hankel Γ1–Γ5 builder from
contour_test_v3.1.jl is retained only for figures.

---

## 8. First build task

`models.jl` (Case S + Case U) and `ribbon_core.jl` + `run_classify.jl`, validated on Case S
(machine-precision match to exact_roots, path independence) and Case U (pinch detected at
ω_f = −0.88). Only then add `CouetteModel` and hit the Couette numbers in §6. Serial throughout.
