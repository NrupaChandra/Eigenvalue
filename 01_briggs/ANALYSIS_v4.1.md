# Briggs v4.1 — diagnosis of the near-pinch stagnation, and a design for v4.2

Run analysed: `results/json/contour_iteration_commonzeta_4e-4_overlapcheckv2.json`
(**1501 iterations** — the file accumulates a continued run, not the single `for k = 1:700` loop),
`results/videos/testv4.1.mp4`, source `src/briggsv4.1.jl`.

---

## Headline

**The run did not stagnate before the pinch. It converged to the pinch and then spent
~800 iterations unable to say so.**

The true branch point of this problem is

| quantity | value |
|---|---|
| pinch wavenumber | α_p = **0.572519492 − 3.038860757 i** |
| pinch frequency | ω_p = **0.298649115 − 0.572829695 i** |
| ∂²ω/∂α² at α_p | 0.406438 − 0.138791 i,  \|ω″\| = 0.42941 |
| absolute growth rate | σ_abs = Im ω_p = **−0.5728** → **absolutely stable** |

These are the *plateau* values: the fit was repeated for `num_modes` = 50, 60, 70, 80, 100,
120, 150, 180 and α_p, ω_p, ω″ are identical to 9 digits over 50–120. `num_modes = 150`
(what the run used) shifts α_p by 9 × 10⁻⁷ and ω_p by 1.3 × 10⁻⁸ — negligible here, but the
polynomial fit residual grows 5 × 10⁻¹¹ → 6.7 × 10⁻⁸ → 2.4 × 10⁻⁷ over nm = 50 → 150 → 180,
which is an independent, quantitative confirmation of the divergent-discretisation finding in
`vib_ribbon/ANALYSIS_v5.1.md` (it puts the edge of the plateau at nm ≈ 100–120).
**These numbers are a regression test**: any version of the solver, run on the straight
contour for Couette at Re = 2000, β = 0, must reproduce them.

and by the end of the run

* F passes within **3.2 × 10⁻³** of α_p — *one third of the F grid spacing* (h_F = 1.01 × 10⁻²);
* Im L is within **2.4 × 10⁻⁷** of Im ω_p;
* max_F Im ω_F agrees with Im ω_p to **6 × 10⁻⁷**.

Meanwhile `d_branch` was reporting 0.11–0.33 and the stopping test `d_branch < 1e-4`
was never going to fire. **`d_branch` is not a convergence measure — it is a measure of
the L-contour grid spacing.** That is cause #1. Cause #2 is a non-monotone L step that
puts the whole system into a ±4 × 10⁻³ limit cycle, which is what makes the video look
like the branches approach and then back off.

Everything below is verified against an independent re-implementation of your temporal
Orr–Sommerfeld operator (Chebyshev collocation, 150 modes, U(y)=y, Re=2000, β=0, the
−200i boundary-row trick). It reproduces the `omega_F` stored in your JSON to
**1.1 × 10⁻⁷ worst case, ~1 × 10⁻⁸ typical**, so agreement between the two is a real
cross-check, not a tautology.

---

## 1. How v4.1 works

### 1.1 The idea, in one paragraph

Briggs' question is: *can the ω-integration contour L be lowered from Im ω = +∞ toward
the real axis while the α-contour F is deformed so that the α⁺ and α⁻ spatial roots
stay on opposite sides of F?* The answer is blocked exactly when two roots from opposite
sides coalesce — a double root of D(α,ω)=0, i.e. ∂D/∂α = 0, i.e. **dω/dα = 0**
(zero group velocity: a saddle of ω(α)). The lowest height L can reach is therefore

> min over admissible contours F of max over α ∈ F of Im ω(α)  =  Im ω(α_p)

a **minimax problem**, whose solution is the steepest-descent contour through the saddle.
v4.1 solves this minimax by a coupled gradient flow: L is pushed down by a constant force
and pushed up by a repulsive potential from the ω_F curve; F is pushed away from the two
spatial branches by a second repulsive potential. The two pushes come into balance
exactly when F sits on the steepest-descent path through the saddle and L rests on
Im ω(α_p). **That is a correct formulation, and it worked.**

### 1.2 Code-level, in the order the loop executes

**ω_F.** `contour_omega_F(F)` calls `couetteflow_temporal_sing_mode` for each α ∈ F and
takes `argmax(imag(eigvals))` — the most unstable temporal mode. ω_F is therefore the
image of F under the (mode-selected) temporal dispersion relation. Its maximum imaginary
part is the ceiling that L may not cross.

**L motion.** L is always a *horizontal segment* `omega_r[j] + omega_i*im`, so only the
scalar `omega_i` evolves. The update is

```
omega_i_cache[j] = omega_i + omega_dt*( -lambda
                     + dω_i/dω_r · ∂_r Φ_L(L[j])     ← identically zero: L is horizontal
                     - ∂_i Φ_L(L[j]) )
```

so effectively `omega_i + dt·(−λ − ∂Φ_L/∂ω_i)`: gravity λ = 4 downward, exponential
repulsion `exp(ζ/|Δω|²) − 1` upward from every point of ω_F. This produces **100 candidate
heights** (one per L node — they differ because each L node is at a different distance
from the ω_F curve). The code then keeps only candidates above
`omega_lower_bound = max Im ω_F + 2e-7` (and, new in v4.1, further than
`min_useful_omega_jump = 1e-8` from the current value) and takes
`omega_candidate = minimum(greater_candidates)` — i.e. the candidate that hugs the ceiling
most closely. The candidate is accepted only if the jump is ≤ 40·omega_dt **and** the
retracked branches pass `branch_overlap_valid` (no upper/lower point pair closer than
1e-8); otherwise omega_dt is halved. A separate *repair* step lifts L back to the ceiling
if ω_F has risen above it.

**Branches.** For each ω ∈ L, `couetteflow(:alpha_collocation)` gives the spatial spectrum;
`dominant_eigvals` / `couetteflow_spatial_sing_mode_comparison` classify each root as
*upper* or *lower* by the **sign of its projection on the normal of the nearest F node** —
so "upper/lower" means "on this side of the current F", not "originating in the upper
half α-plane". `contour_alpha_L_conti` seeds at index N/4 and marches both directions with
`track_eigenvalue_simple(α_prev) = α_prev`, i.e. zeroth-order continuation ("nearest root
on the correct side of F to the previous one").

**F motion.** F is a *graph* α_i(α_r) on a fixed α_r grid. The update is the normal
gradient flow of Φ_F in graph form,

```
∂α_i/∂t = [ α_i′ ∂_r Φ_F − ∂_i Φ_F + σ α_i″ ] / (1 + α_i′²),   σ = 3e-5
```

with the RHS clamped to ±100, followed by an unconditional 15-point box filter
(`rolling_average_filter(·, 7)`) and Neumann ends.

**ζ.** `zeta_common = zeta_alpha = zeta_omega = 4e-4`, **fixed**. There is no adaptive ζ in
v4.1 (that was v2). So "ζ_alpha collapsing" cannot be the cause here.

**Adaptive/acceptance logic.** Three separate mechanisms:
`branch_slowdown_factor(d_branch) = clamp(d_branch/0.2, 0.1, 1)` scales Δt;
`acceptance_factor` shrinks from 10 to 0.01 as F becomes equally influenced by both
branches; the attempt loop halves `local_delta_t` up to 5 times when
`all(move_raw .≤ factor .* max.(|α_i|, 1e-12))` fails.

**Stopping.** `d_branch = min_i |α_u(ω_i) − α_l(ω_i)| < pinch_tol = 1e-4`.

### 1.3 What changes when the branches get close to F

**Nothing switches on.** This is worth stating plainly, because your description assumed a
second mechanism takes over. v4.1 has no near-pinch mode. The only things that change
continuously are: Δt is clamped at 0.1× by `branch_slowdown_factor`; `move_safety`
tightens from 2.0 → 0.8 as `d_contour` falls; and the ω_F ceiling becomes the binding
constraint on L instead of the potential. The same gradient flow runs from k=1 to k=1501.

---

## 2. What actually happened in the run

| phase | iterations | what happens |
|---|---|---|
| **descent** | 1 – ~160 | L falls from Im ω = 0 to −0.57 at the full rate λΔt = 4×10⁻³/step; F is dragged from the real axis down to Im α ≈ −3.0. Working exactly as intended. |
| **arrival** | ~160 – 212 | Im L comes within 10⁻³ of the true Im ω_p (k=212). The bulk deformation is finished. |
| **slow refinement** | ~212 – 1000 | F's shape converges onto the steepest-descent path: min\|F − α_p\| improves 0.10 (k=190) → 0.030 (k=518) → 0.010 (k=754) → 0.0040 (k=989). **Real progress, completely invisible in `d_branch`.** |
| **saturation** | ~1000 – 1501 | min\|F − α_p\| saturates at 3.2×10⁻³ ≈ h_F/3. Nothing further is available at N = 100. |

Throughout the last ~900 iterations L is in a hard limit cycle: **441 upward steps
(mean +3.4×10⁻³) and 460 downward steps (mean −3.3×10⁻³), net drift +6.4×10⁻⁴ — zero
progress.** `d_branch` flips between 0.112 (L on the ceiling) and 0.328 (L one bounce
above it). F is dragged along, oscillating by 2.4×10⁻² peak-to-peak at the pinch node;
its net motion over 900 iterations is only 19 % of its cumulative path length.

That bounce is what the video shows. It is not the branches failing to close — it is
L jumping away from the answer and back, ~450 times.

![diagnosis](fig1_diagnosis.png)

---

## 3. Why it stagnates

### 3.1 Primary cause — `d_branch` measures the L grid, not convergence

Near a simple branch point the two spatial roots are

> α_±(ω) = α_p ± √( 2(ω − ω_p) / ω″ ) + O(ω − ω_p)

so the gap is **square-root** in the ω-distance:

> d_branch = 2 √( 2 |ω − ω_p| / |ω″| )

Check against the run at its closest state (k=1034): the L node nearest ω_p is
ω = 0.297980 − 0.572816 i, so |ω − ω_p| = 6.694 × 10⁻⁴, giving a predicted gap of
**0.111674**. Observed `d_branch` = **0.111558**. Agreement to **0.1 %**. Over the whole
stalled phase the run's points lie on that curve (Fig. 1c).

The consequences are brutal:

* L is a fixed 100-point grid on ω_r ∈ [0, 0.5], spacing **5.05 × 10⁻³**. Even with
  Im L *exactly* on Im ω_p, the nearest node is generically ~10⁻³ away in Re ω, giving a
  **floor on d_branch of ~1.3–1.9 × 10⁻²** — and in this run the node happens to sit
  6.7 × 10⁻⁴ off, giving 0.11.
* To reach `pinch_tol = 1e-4` you need |ω − ω_p| < 5.4 × 10⁻¹⁰, i.e. an L grid of
  **≈ 5 × 10⁸ points**. The stopping test is unreachable by four to five orders of
  magnitude, for any run, at any resolution you would ever use.

Additionally, `d_branch` is *structurally* bounded away from zero: `dominant_eigvals`
assigns α_u and α_l to opposite sides of F, so they can never be reported as coincident
unless both lie exactly on F. That constraint is harmless at the present gap of 0.11 but
would bite immediately if you tried to drive d_branch to ~h_F.

### 3.2 Secondary cause — the L step is not monotone downward

I reproduced the ω-block of the loop exactly (same Φ_L, same ζ, same λ, same Δt, same
filters) and it predicts the accepted `omega_i` **exactly** in every iteration tested.
The mechanism:

* When the clearance `Im L − max Im ω_F` is larger than λΔt = 4 × 10⁻³, ~36–98 of the
  100 candidates lie above the ceiling and the minimum of them is a *descent*. Fine.
* When the clearance drops below λΔt, every descending candidate falls **below** the
  ceiling and is filtered out. The only survivors are the 3–5 nodes closest to the ω_F
  peak, which the exponential repulsion has pushed **upward**. `minimum(greater_candidates)`
  then returns a value *above* `omega_i`, and the code accepts it — there is no check
  that the step is downward.

Example, k=1442: clearance 8.0 × 10⁻⁷, 4 candidates above the bound, minimum
= −0.5715393 vs `omega_i` = −0.5728294 → **L jumps up 1.3 × 10⁻³**. Confirmed against
the stored run.

**This is precisely what v4.1 changed relative to v4.** The only functional differences
are `min_useful_omega_jump = 1e-8` and `minimum(greater_candidates)` in place of
`argmin(|x − bound|)` — which are the same thing for candidates above the bound. So the
sole behavioural change is the 1e-8 filter. In v4 the near-zero step was accepted (the
"technically accepts a step but the contour almost does not move" the README describes);
v4.1 rejects it and is forced onto an *upward* candidate instead. Both are symptoms of the
same missing rule: **the accepted L step must be downward and clamped to the ceiling.**

### 3.3 The other candidates, ruled out with evidence

| hypothesis | verdict | evidence |
|---|---|---|
| potential ineffective near the pinch | **not binding** | F reached h_F/3 from α_p; the flow was still moving F usefully to k≈1000 |
| adaptive ζ_alpha collapsing | **impossible in v4.1** | ζ is a fixed constant 4e-4; no adaptation code exists |
| acceptance criterion suppressing F motion | **inactive at the end**, active early | at the pinch \|α_i\| ≈ 3.03, so the RHS is 0.03–30 while moves are 5×10⁻⁴ — passes with 5 orders of margin. *Early* (k ≲ 20) α_i ≈ 0 so the RHS is ~10⁻⁹ and the test fails every attempt — but the code applies `alpha_i_smooth` regardless (see below), so the only effect is Δt/32. |
| smoothing/curvature cancelling the deformation | **secondary, real, but not the blocker** | the box filter still removes 3.5–4.7 × 10⁻⁴ per iteration from the converged F, comparable to the net drive (1.3–6.2 × 10⁻⁴) — it limits F's final accuracy, but F is already at h_F/3, so it cannot explain a 0.11 gap |
| F forced to stay between the branches | **structurally true, not binding here** | see §3.1; would bite only below d_branch ~ h_F |
| L-lowering rule hitting an artificial limit | **YES — but the limit is correct** | the ceiling `max Im ω_F + 2e-7` *is* the physical constraint; the bug is the bounce off it, §3.2 |
| branch tracking switching / losing branches | **NO, tracking is correct at the pinch** | at k=1034 both α_u and α_l are genuine roots of the temporal dispersion relation at ω = L[59] (residual 5×10⁻⁷); their midpoint is 2.4×10⁻³ from α_p and they match α_p ± √(2Δω/ω″) to 2.4×10⁻³ |
| discretization preventing resolution of the two roots | **YES, in ω — this is cause #1** | §3.1 |
| gradient flow only approaches asymptotically | **true in principle, never the binding constraint** | a repulsive-potential flow with F wedged between two branches does equilibrate at finite separation, and exp(ζ/d²) forbids contact — but F got to h_F/3 of α_p, so this was never what stopped you |

Two smaller code observations, neither causal:

* In the F block, `alpha_i_cache` and `accepted` are computed and never used —
  `alpha_i = alpha_i_smooth` runs unconditionally. The acceptance test therefore only
  controls Δt halving; it does not gate the update. Worth making explicit.
* The box filter is applied at **full strength once per iteration regardless of
  `local_delta_t`**, so as Δt shrinks the smoothing-to-drive ratio grows without bound.
  At the end of this run they are the same order.
* `dist_u = minimum(abs.(F .- alpha_L_u))` in the printout is *index-wise*, comparing
  F[j] with α_u[j] — two unrelated parameterisations. It is only logged, but it is
  misleading in the log.
* `argmax(Im ω_F)` landed on an endpoint of F in 23 of the last 1302 iterations. The
  ceiling is then set by an endpoint rather than by the interior ridge — worth a guard.

---

## 4. Against the actual Briggs criterion

A narrow gap is not a pinch. The three things that must hold, and what the run gives:

**(i) Double root.** Solved directly for dω/dα = 0 by a local polynomial fit of ω(α) on
circles of radius R = 0.06 / 0.03 / 0.015 about the run's own argmax node. The three radii
agree to **1.1 × 10⁻⁶** in α_p and the fit residual is 5–8 × 10⁻⁸:

> α_p = 0.572519351 − 3.038861399 i,  ω_p = 0.298649120 − 0.572829683 i,
> ω″ = 0.406438 − 0.138791 i (non-degenerate, so a simple square-root branch point).

*Practical warning:* a naive finite-difference Newton on dω/dα fails here. With h = 2×10⁻⁵
the eigen-solver's ~10⁻⁸ noise gives ω″ noise of order 4·10⁻⁸/h² ≈ 10², and I initially got
ω″ = 57.6 instead of 0.429 — a factor 134 error, with a plausible-looking convergence
history. Sample on a circle of radius ~10⁻² and fit; never difference at h ≪ 10⁻³.

**(ii) Square-root structure.** At ω = ω_p + iδ the measured root pair satisfies
|α₊ − α₋| = 2|h| with relative error 6.0 × 10⁻³ → 7.0 × 10⁻⁴ → 9.5 × 10⁻⁵ for
δ = 10⁻² → 10⁻³ → 10⁻⁴ (clean O(h²) = O(δ) convergence), and the midpoint offset divided
by |h|² is a constant 0.757 / 0.760 / 0.766. Textbook second-order branch point.

**(iii) Opposite origins — the test that actually makes it Briggs.** Continuing both roots
as Im ω is raised from Im ω_p to Im ω_p + 3.5:

* **root A** migrates from Im α = −2.98 up through the real axis at **Im ω ≈ −0.083** and
  on to Im α = +11.3 → genuine **α⁺**;
* **root B** stays at Im α ≈ −3.09 → −3.21 for the whole path, never crossing → genuine **α⁻**.

They end on opposite sides of Im α = 0. **This is a genuine α⁺/α⁻ Briggs pinch**, and the
crossing at Im ω ≈ −0.083 is exactly why F had to be deformed below the real axis at all.

![physics](fig2_physics.png)

**Physical answer.** σ_abs = Im ω_p = −0.5728 < 0: plane Couette at Re = 2000, β = 0 is
**absolutely stable**, and strongly so — max over real α of Im ω is −0.0191 (also stable),
so the absolute mode is far more damped than the least-damped convective one. This is the
expected result and it is a good sign that the machinery reproduces it.

---

## 5. What must not be changed

These are working and should be left alone:

1. **The whole Stage-1 formulation** — the coupled L/F gradient flow, the exponential
   repulsive potentials Φ_L and Φ_F, ζ = 4 × 10⁻⁴, λ = 4, σ = 3 × 10⁻⁵, the graph-form
   normal flow for F. It delivered the correct minimax contour to a third of a grid cell.
2. **`omega_lower_bound = max Im ω_F + clearance`** — this is the physical constraint, not
   an artificial one. Keep it; only fix how the step approaches it.
3. **Branch classification by signed normal projection onto F.** This is the correct
   Briggs bookkeeping once F is deformed; it was verified correct at the pinch. (Verify
   the α⁺/α⁻ identity once at the end by the origin test — not every iteration.)
4. **`branch_overlap_valid`** as a guard against the two branches selecting the same
   eigenvalue.
5. **The JSON history format.** It contains everything needed to recover the pinch
   *post hoc*, which is what let this diagnosis be done without re-running anything.

---

## 6. Fixes, least invasive first

| # | change | effort | risk | what it buys |
|---|---|---|---|---|
| **1** | **Post-process the existing JSON** with a local-fit pinch solver. No solver change at all. | 0 | none | You already have the answer — α_p, ω_p, ω″, σ_abs — from runs on disk. |
| **2** | **Replace the stopping test.** Drop `d_branch < 1e-4`; stop on clearance + a converged parabolic cap in ω_F (§7). | ~10 lines | none | The loop terminates when it is actually done, at k≈700 instead of never. |
| **3** | **Make the L step monotone.** Require `omega_candidate < omega_i`; if no descending candidate is above the bound, halve `omega_dt` instead of accepting an ascent; clamp to the bound. Make `min_useful_omega_jump` relative to the clearance. | ~10 lines | low | Removes the ±4×10⁻³ limit cycle; L converges geometrically onto the ceiling and F settles instead of being shaken. |
| **4** | **Add Stage 2** as a routine called *after* the loop (or when the gate fires). Purely additive. | ~60 lines | low | Turns "close to the pinch" into α_p, ω_p, ω″ at machine-ish accuracy plus a verification report. |
| **5** | **Fix the F-update bookkeeping**: use `alpha_i_cache` on acceptance, and scale the smoothing with `local_delta_t` (or use only the σ α_i″ term, which is already dt-scaled). | ~15 lines | medium — changes Stage-1 trajectories, needs re-validation | Improves F's converged accuracy below the current ~4×10⁻⁴ smoothing floor. |
| **6** | **Adaptive ω-resolution**: refine L locally around ω_p, or drop L entirely once Stage 2 starts (Stage 2 needs only ω(α)). | larger | medium | Only needed if you want d_branch itself to be small, which — see §3.1 — you should not. |

Fixes 1–4 are sufficient. I would do 1 immediately, then 2 + 3 + 4 together as v4.2, and
treat 5 as a separate experiment.

---

## 7. Recommended design — Briggs v4.2

Two stages, with a quantitative gate between them. Stage 1 is v4.1 with two small
corrections; Stage 2 is new and additive.

### Stage 1 — global contour deformation (unchanged physics)

Keep everything in §5. Two corrections:

```julia
# --- L step: monotone descent, clamped to the ceiling ---
descending = filter(x -> isfinite(x) && x > omega_lower_bound && x < omega_i, omega_i_cache)

if length(descending) >= min_valid_omega_candidates
    omega_candidate = minimum(descending)        # hug the ceiling from above
elseif (omega_i - omega_lower_bound) > omega_settle_tol
    omega_dt *= 0.5                              # step too large: shrink, never ascend
    continue
else
    omega_settled   = true                       # L is resting on the ceiling
    omega_candidate = omega_i
end
```

and replace `min_useful_omega_jump = 1e-8` by something relative, e.g.
`1e-3 * (omega_i - omega_lower_bound)`, so that when the clearance is 10⁻⁶ a step of
10⁻⁹ still counts as useful. Note this changes no physics — it only removes steps that
move *away* from the target.

### The gate — a quantitative pinch-neighbourhood test

Enter Stage 2 when **all** of the following hold (all quantities already computed):

1. `omega_settled` — L is resting on the ceiling with clearance < ε_c (ε_c ≈ 10⁻⁶,
   or 10⁻⁴ × the total descent);
2. `i = argmax(imag.(omega_F))` is interior (3 ≤ i ≤ N−2) and is a genuine local maximum,
   `Im ω_F[i−1] < Im ω_F[i] > Im ω_F[i+1]` — this rules out the endpoint case seen in 23
   iterations of this run;
3. the branches straddle F near node i (already available from the signed projections);
4. `max Im ω_F` has changed by less than ε_c over the last M iterations.

Do **not** gate on `d_branch`.

### Stage 2 — local pinch solver

Of the options you listed, the one most compatible with the existing code and most robust
numerically is **direct solution of dω/dα = 0 using only the temporal operator you already
have**, by local polynomial fit rather than by differencing. Reasons: (a) it reuses
`couetteflow(:omega_collocation)` and needs no new operator; (b) it never differentiates an
eigenvalue at small h, which is where the naive Newton fails (§4); (c) it returns ω″ for
free, which is what you need for the verification and for the error bars; (d) the spatial
solver's own accuracy at the coalescence is only ~5 × 10⁻⁷, so anything built on the
spatial roots (`∂D/∂α = 0` directly, or minimum-distance-between-branches) inherits a worse
conditioning.

**Algorithm.**

1. **Seed.** α_c = ½(F[i] + ½(α_u(ω*) + α_l(ω*))) where i is the gate's argmax node and
   ω* is the L node nearest ω_F[i]. In this run the two seeds are 8 × 10⁻³ apart and either
   works.
2. **Sample.** M = 24 points on |α − α_c| = R, R ≈ 5 h_F ≈ 0.05. At each point take the
   eigenvalue **nearest the previous one** (continuation), *not* `argmax(imag)` — the
   argmax selection is discontinuous and will corrupt the fit.
3. **Fit.** Least squares, degree 6, in s = (α − α_c)/R. Report the residual; it should be
   ~10⁻⁸ (it is the built-in health check).
4. **Solve.** p′(s) = 0, take the root nearest s = 0. Then
   α_p = α_c + R s₀, ω_p = p(s₀), ω″ = p″(s₀)/R².
5. **Repeat** with R → R/2, three levels; require the α_p values to agree to < 10⁻⁵.
   (In this run: 1.1 × 10⁻⁶.)

Cost: 3 × 24 = 72 temporal eigen-solves — about one Stage-1 iteration.

Use the branch data you already have as an **independent seed and cross-check**, not as the
answer: from any two L nodes near the turning point, ½(α_u + α_l) → α_p and
(α_u − α_l)² = 8(ω − ω_p)/ω″ give α_p, ω_p, ω″ from two linear fits at zero extra cost. In
this run that route gives α_p to 2.4 × 10⁻³ and Im ω_p to ~10⁻⁶ immediately.

### Verification and reporting (replaces the pinch tolerance)

Report, and require:

1. |dω/dα|(α_p) below the fit noise floor, and **ω″ ≠ 0** (non-degenerate).
2. The square-root table: |α₊ − α₋| vs 2√(2δ/|ω″|) for δ = 10⁻², 10⁻³, 10⁻⁴, with the
   relative error falling like δ.
3. **The α⁺/α⁻ origin test** — continue both roots to Im ω = Im ω_p + O(1) with a ladder
   fine enough that no single step exceeds ~0.5 in α, and confirm they end on opposite
   sides of Im α = 0. This is the step that distinguishes a Briggs pinch from two same-side
   roots touching, and it is the one thing v4.1 never checked.
4. σ_abs = Im ω_p, ω_r,abs = Re ω_p, α_p, ω″.
5. Stage-1 quality, none of which involves d_branch: min|F − α_p| relative to h_F;
   |Im L − Im ω_p|; |max Im ω_F − Im ω_p|.

Optionally still log a *derived* d_branch, but label it for what it is:
2√(2|ω_nearest − ω_p|/|ω″|) — a measure of the ω grid.

### Two items on the project's open list that this closes

`vib_ribbon/src/README_v4_series.md` records, under "Things to know before writing a v5":

> *"The calculation stalls before it pinches. This is the biggest open problem and it is older
> than all of these versions. ω_i freezes around −0.5725 and the gap between the curves
> flattens at about 0.12, never reaching the 1e-4 target."*

−0.5725 **is** Im ω_p = −0.5728297, and 0.12 **is** the L-grid floor 2√(2|Δω|/|ω″|). There is
nothing wrong with the descent. Note also that raising N (v5's 100 → 300, Δω_r = 8.4 × 10⁻⁴)
only moves that floor from ~0.19 to ~0.089 — still ~900× above `pinch_tol`. **No achievable
resolution will make the pinch test fire.**

The same list records *"F drifting bodily downwards rather than staying a deformation of the
real α axis"*. For the straight-contour Couette case that drift is **correct**: the α⁺ root
crosses Im α = 0 at Im ω ≈ −0.083 (§4), so F *must* go below the real axis to stay below it,
and Im α ≈ −3.04 is exactly where the saddle sits. The test of whether the drift is right is
not how far F moved, but whether it ends on the steepest-descent path through α_p — which
`pinch_refine.py` now measures directly as min|F − α_p| / h_F.

### One caveat

This diagnosis characterises the pinch that *this* homotopy found. In principle another
saddle with larger Im ω_p could be encountered by a different deformation path. Nothing in
the run suggests one, but a coarse scan of stationary points of ω(α) over the region F
sweeps would be a cheap addition to the v4.2 verification.

---

## Files

* `pinch_refine.py` — read-only Stage-2 solver/verifier. Post-processes any
  `contour_iteration*.json`; prints α_p, ω_p, ω″, σ_abs, the Stage-1 quality metrics, and
  (with `--verify`) the square-root table and the α⁺/α⁻ origin test. Contains an
  independent re-implementation of the temporal operator and self-checks it against the
  `omega_F` in the file. This is fix #1 and the reference implementation of Stage 2.

  ```
  python3 pinch_refine.py results/json/contour_iteration_commonzeta_4e-4_overlapcheckv2.json --verify
  ```

* `fig1_diagnosis.png`, `fig2_physics.png` — the figures above.
