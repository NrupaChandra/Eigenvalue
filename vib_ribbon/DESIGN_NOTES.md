# vib_ribbon — design notes (pre-code stage)

Status: planning only. No solver code exists in this folder yet.
Basis: analysis of `Briggs/briggs_test.jl` (analytic pinch validation — the run that worked),
`Briggs/briggsv4.jl` (bare-bones Couette), `Briggs/contour_test_v3.1.jl` (Hankel-contour
prototype), and the NEAT / TSFP14 / Huerre–Monkewitz literature in `Literature/`.

---

## 1. Verdict on the proposed starting point

**Instinct:** take the bare-bones analytic-test implementation (`briggs_test.jl`), which
demonstrably converged to the exact pinch point, and edit it into the vibrating-ribbon solver,
borrowing from `briggsv4.jl` / `contour_test_v3.1.jl`.

**Assessment: right skeleton, wrong main loop.**

Right, because `briggs_test.jl` is the cleanest code in the repo and it already demonstrates
the one architectural feature the ribbon solver must have: the physical model sits behind two
swappable functions (`algebraic_omega`, `algebraic_alpha_roots`) while the Briggs logic is
model-agnostic. That seam is exactly what lets us validate the ribbon machinery on an analytic
relation before touching Couette.

Wrong, because the loop `briggs_test.jl` runs — pseudo-time contour relaxation hunting a
*pinch point* — is not the computation the ribbon problem needs. The pinch hunt answers
"is the flow absolutely unstable?" (a one-time validity check). The ribbon/signalling problem
answers "what spatial waves does a harmonic source at real frequency ω_f produce?" — and for
plane Couette that computation needs **no moving contours at all**:

- Verified numerically (Re = 2000, β = 0): descending Im ω from 0.3 to 0 at ω_r = 0.5, no
  spatial root ever approaches the real α-axis (min |Im α| ≈ 0.10 at the axis). F can stay
  the real axis; branch identity is fixed once at the top of the descent by sign(Im α).
- The two relaxation pathologies already observed in the repo (L-update deadlock at k ≈ 51 and
  unbounded sinking of F in the `contour_test_v3.1.jl` run; fake-pinch branch collapses in the
  v2–v4 history) simply cannot occur in a fixed-contour design.

So: **keep the `briggs_test.jl` skeleton and its model seam; delete its relaxation loop;
replace it with a causal σ-descent classifier plus an ω_f sweep.** Borrow from
`contour_test_v3.1.jl` only the serial memory-based tracking style
(`couetteflow_spatial_pair_comparison` / `contour_alpha_L_conti(L, old_...)`) and the
polyline-intersection test (as a diagnostic); borrow from `briggsv4.jl` the Couette matrix
block unchanged. The Hankel Γ1–Γ5 contour builder from `contour_test_v3.1.jl` is kept as an
*analysis/figure* object — numerically the only part of it that matters is the vertical
descent at ω_r = ω_f.

## 2. Why the analytic test still comes first (your instinct, upgraded)

The existing analytic model (C1 = −1 + 0.96i, C2 = 0.4 + 0.8i, C3 = 1) is **absolutely
unstable** (Im ω₀ = +0.8), so the signalling problem is *meaningless* for it — transients
contaminate everything (Huerre & Monkewitz §3.1). It is therefore the right fixture for the
validity-check path, but not for the ribbon path.

Proposal: add a second analytic fixture with the same quadratic form but coefficients chosen
so the model is globally temporally stable, e.g.

    C1 = −1.0 − 0.96i,  C2 = 0.4 + 0.8i,  C3 = 1.0 − 0.5i

My derivation (re-verify at implementation time): Im ω(α real) = −0.96 + 0.8α − 0.5α² has
maximum −0.64 < 0 (stable for all real α); saddle at α₀ = −C2/(2C3) = −0.4i,
ω₀ = C1 − C2²/(4C3) = −0.84 − 1.04i, so Im ω₀ < 0 — no absolute instability, signalling
well-posed. The spatial roots α±(ω_f) are known in closed form (quadratic formula), and they
split cleanly across the real α-axis at real ω_f (spot-checked at ω_f = −0.9:
α⁺ ≈ 0.482 + 0.436i, α⁻ ≈ −0.482 − 1.236i).

This gives the ribbon solver what the pinch solver had: **an exact analytic answer to converge
to.** Every stage below is validated on this fixture before Couette is switched in.

## 3. Planned contents of this folder

    vib_ribbon/
    ├── DESIGN_NOTES.md            (this file)
    ├── model_analytic.jl          stable quadratic fixture: temporal_spectrum(α), spatial_spectrum(ω), exact α±(ω)
    ├── model_couette.jl           matrix block + couetteflow from briggsv4.jl, unchanged, behind the same two functions
    ├── ribbon_classify.jl         Stage A: σ-descent classifier at one ω_f (serial, no Distributed, no plots)
    ├── ribbon_sweep.jl            Stage B: ω_f sweep with warm-start continuation + cold-start checkpoints
    └── results/                   one JSON per run: ω_f grid, α⁺/α⁻ sets, guard flags, resolution deltas

Nothing from `Briggs/` is modified; the ten historical versions stay untouched as reference.

## 4. The two-stage algorithm (what replaces the relaxation loop)

Stage A — classify at one forcing frequency ω_f:
1. σ₀ chosen above all temporal growth (Couette: 0.3 is generous; analytic fixture: 0.1).
2. Full spatial spectrum at ω = ω_f + iσ₀; spurious filter = keep roots that agree between
   num_modes = 100 and 150 to tolerance, inside |α| < cutoff.
3. Split by sign(Im α) → α⁺ / α⁻ sets (causal identity, fixed once).
4. Descend σ → 0 in ~20 fixed steps; mutual nearest-neighbour matching of the whole set;
   guards: per-side count conservation, jump bound |Δα| ≤ c·|dα/dω|_est·Δσ, min |Im α| monitor.
5. Output: α⁺/α⁻ sets at ω = ω_f, the least-damped members, and the axis-clearance min |Im α|.

Stage B — sweep ω_f over the band:
- warm-start each ω_f from the previous root set (one solve + matching);
- every K-th frequency: independent cold-start descent, assert agreement;
- deliverable: curves Re α⁺₁(ω_f), Im α⁺₁(ω_f) (spatial decay rate of the ribbon signal), with
  resolution error bars from the 100/150-mode difference.

Validity certificate (once per Re, β): temporal scan max Im ω < 0 over the α-window, plus the
global min |Im α| over the sweep. If min |Im α| ever approaches 0, the fixed-contour
assumption fails and the F-deformation machinery (v4 / v3.1) gets reactivated — that is the
*only* trigger for bringing the relaxation code back.

## 5. Acceptance tests fixed in advance

Analytic fixture: classifier output must equal the closed-form α±(ω_f) to ~1e-12, for any σ₀
and step count (path independence = causality check).

Couette, ω_f = 0.5 (reference values computed independently in the design stage, M = 100):
- least-damped α⁺ = 0.7496 + 0.1016i
- least-damped α⁻ = 0.059 − 3.275i
- cross-check against the frozen arc values of the `contour_iteration_v3.1.json` run:
  α_u = 0.7367 + 0.1006i at ω = 0.49 (consistent after dα_r/dω_r ≈ 1.3 offset correction).
- sweep smoothness: α⁺₁ from 0.194 + 0.043i (ω_f = 0.1) to 1.254 + 0.137i (ω_f = 0.9).

## 6. Open decisions before the first line of code

1. Ribbon model: interior harmonic source δ(x)δ(y−y′)e^(−iω_f t) (NEAT §III.B) or
   oscillating-wall BC? If amplitudes are wanted, which y′?
2. Deliverable: branch curves α±(ω_f) only, or also receptivity amplitudes (residues, needs
   ∂D/∂α and eigenvectors — phase 2)?
3. ω_f band and Re list; β = 0 fixed for the first milestone?
4. Is `Eigenvalue_analysis/` the comparison dataset, and where is its generating solver?

## 7. First coding task

`model_analytic.jl` + `ribbon_classify.jl` run on the stable analytic fixture, matching the
closed-form roots to machine precision. Only after that is green: swap in `model_couette.jl`
and hit the Couette acceptance numbers above. Serial only; parallelization (pmap over ω_f) is
deferred until the serial sweep is validated.
