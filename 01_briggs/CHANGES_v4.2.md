# briggsv4.2.jl — change specification

Scope: `src/briggsv4.1.jl` → `src/briggsv4.2.jl`. Roughly **90 lines changed or added**,
in five places. No physics is modified. The diff is deliberately small so that the
regression test in §7 can attribute any change in the answer to a specific edit.

---

## 0. What is byte-identical to v4.1

Do not touch any of this — it produced a correct minimax contour to a third of a grid cell:

* `couetteflow` (both modes), `couetteflow_temporal_sing_mode`,
  `couetteflow_spatial_sing_mode_comparison`, `dominant_eigvals`
* `contour_omega_F`, `contour_normals`, `contour_alpha_L_conti`, `track_branch_pmap`
* `phi_L`, `Phi_L`, `phi_F`, `Phi_F` and all four gradient routines
* `zeta_common = 4e-4`, `lambda = 4.0`, `sigma = 3e-5`, `delta_t = 1e-3`
* the F update itself (the `rhs_j` expression), `branch_slowdown_factor`,
  `directional_dt_check`, `acceptance_factor`, `influence_ratio`
* `branch_overlap_valid`, the repair step, the JSON history format
* `num_modes = 150` — see §6, this is a *separate* commit

---

## 1. L step: monotone descent, clamped to the ceiling  (~15 lines)

**Where:** inside `for omega_attempt in 1:50`, replacing the `greater_candidates` /
`omega_candidate` block.

**Why:** the `j`-loop runs `for j in 2:(length(omega_i_cache)-1)`, so entries `[1]` and `[N]`
are never updated and stay exactly equal to `omega_i`. In v4 those two produced the
zero-length step the README complains about; v4.1's `min_useful_omega_jump = 1e-8` removes
them but leaves only *upward* candidates once the clearance drops below `lambda*omega_dt`,
and `minimum(greater_candidates)` then returns a value **above** `omega_i` which the code
accepts. That is the ±4e-3 limit cycle.

```julia
# --- v4.2 -------------------------------------------------------------
# the endpoints of omega_i_cache are never written by the j-loop; drop them
# explicitly instead of relying on a magnitude filter
interior_cache = omega_i_cache[2:end-1]

descending = filter(x -> isfinite(x) && x > omega_lower_bound && x < omega_i,
                    interior_cache)

if length(descending) >= min_valid_omega_candidates
    omega_candidate = minimum(descending)          # hug the ceiling from above
    # ... unchanged from here: jump check, L_trial, contour_alpha_L_conti,
    #     branch_overlap_valid, accept / halve omega_dt ...

elseif (omega_i - omega_lower_bound) > omega_settle_tol
    omega_dt *= 0.5                                # step overshoots: shrink, never ascend
    continue

else
    omega_settled   = true                         # L is resting on the ceiling
    omega_accepted  = true                         # do NOT break the outer k-loop
    omega_candidate = omega_i
    omega_status    = "settled"
    break
end
```

Delete `min_useful_omega_jump` (superseded). Add near the other tolerances:

```julia
omega_settle_tol = 1e-8      # clearance below which L is considered resting
omega_settled    = false     # reset each k
```

**Behaviour this produces:** while the clearance exceeds `lambda*omega_dt` the descent is
unchanged. Below that, `omega_dt` halves until the step fits, so L converges geometrically
onto the ceiling instead of bouncing. When F later lowers the ceiling, the clearance grows,
descending candidates reappear, and L resumes on its own — no extra logic needed. The repair
step is untouched and still lifts L if ω_F rises above it.

---

## 2. The gate: a quantitative pinch-neighbourhood test  (~20 lines, new function)

Replaces `d_branch < pinch_tol` as the thing that ends the loop.

```julia
function pinch_gate(omega_F, omega_i, ceiling_hist;
                    clear_tol = 1e-6, hold = 25)
    N   = length(omega_F)
    i   = argmax(imag.(omega_F))
    interior = 3 <= i <= N - 2
    cap = interior && imag(omega_F[i-1]) < imag(omega_F[i]) > imag(omega_F[i+1])
    clearance = omega_i - imag(omega_F[i])
    settled   = 0.0 <= clearance < clear_tol
    steady    = length(ceiling_hist) >= hold &&
                (maximum(ceiling_hist[end-hold+1:end]) -
                 minimum(ceiling_hist[end-hold+1:end])) < clear_tol
    return (interior && cap && settled && steady), i, clearance
end
```

Four conditions, all cheap and all from quantities already computed:

1. the ω_F maximum is **interior** — this rules out the endpoint case, which occurred in 23
   of the last 1302 iterations of the v4.1 run and would otherwise set a bogus ceiling;
2. it is a genuine **local maximum** with a resolved cap (the saddle crossing);
3. L is **resting** on it (`omega_settled` from §1);
4. the ceiling has been **steady** for `hold` iterations.

`d_branch` appears nowhere in it.

---

## 3. Stopping and reporting  (~10 lines)

**Where:** the `if d_branch < pinch_tol` block.

```julia
push!(ceiling_hist, maximum(imag.(omega_F)))
gate_ok, i_peak, clearance = pinch_gate(omega_F, omega_i, ceiling_hist)

if gate_ok
    stop_after_save = true
    stop_reason     = "pinch neighbourhood reached"
end
```

Keep `pinch_tol` **only** as an absolute-instability alarm, exactly as v5 already proposed:
if `imag(omega_p) > 0` after §4 runs, flag absolute instability. Delete it as a convergence
test.

Log line: replace `dUL/dUF/dLF` with the quantities that actually converge —
`clearance`, `i_peak`, `max Im ω_F`, and `omega_status`. Keep `d_branch` if you like, but
label it for what it is. (Also fix `dist_u = minimum(abs.(F .- alpha_L_u))`: that compares
`F[j]` with `alpha_L_u[j]`, two unrelated parameterisations, and is meaningless as printed.)

---

## 4. Stage 2: the local pinch solver  (~45 lines, new, called after the loop)

Uses only `couetteflow(:omega_collocation)`, which you already have. No new operator, no new
physics, no dependency beyond `LinearAlgebra` (already imported).

```julia
function poly_roots(c::Vector{ComplexF64})          # roots of sum c[k+1] s^k
    n = length(c)
    while n > 1 && abs(c[n]) < 1e-30; n -= 1; end
    n <= 1 && return ComplexF64[]
    a = c[1:n] ./ c[n]
    C = zeros(ComplexF64, n-1, n-1)
    for i in 1:n-2; C[i+1, i] = 1.0; end
    C[:, end] .= -a[1:n-1]
    return eigvals(C)
end

function omega_branch(alpha, ref)                   # continuation, NOT argmax(imag)
    ev, _ = couetteflow(alpha, nothing, Val(:omega_collocation), Val(:eigen))
    ev = ev[[isfinite(real(e)) && isfinite(imag(e)) for e in ev]]
    return ev[argmin(abs.(ev .- ref))]
end

function refine_pinch(alpha_c, omega_ref; R = 0.05, M = 24, deg = 6, levels = 3)
    ap, wp, w2 = ComplexF64(alpha_c), ComplexF64(omega_ref), 0.0 + 0im
    history = Tuple{Float64,ComplexF64,ComplexF64,Float64}[]
    for lev in 1:levels
        pts  = [alpha_c + R * cis(2pi * (m - 1) / M) for m in 1:M]
        vals = ComplexF64[]
        ref  = omega_ref
        for z in pts
            ref = omega_branch(z, ref); push!(vals, ref)
        end
        closure = abs(omega_branch(pts[1], vals[end]) - vals[1])   # check that CAN fail
        closure > 1e-9 && @warn "continuation did not close ($closure); raise M or lower R"

        S = [(z - alpha_c) / R for z in pts]
        V = [S[m]^k for m in 1:M, k in 0:deg]
        c = V \ vals
        resid = maximum(abs.(V * c .- vals))

        dc  = ComplexF64[k * c[k+1]           for k in 1:deg]        # p'
        ddc = ComplexF64[k * (k-1) * c[k+1]   for k in 2:deg]        # p''
        rts = filter(z -> abs(z) < 1.2, poly_roots(dc))
        isempty(rts) && error("no stationary point of omega(alpha) inside R=$R")
        s0 = rts[argmin(abs.(rts))]

        ap = alpha_c + R * s0
        wp = evalpoly(s0, c)
        w2 = evalpoly(s0, ddc) / R^2
        push!(history, (R, ap, wp, resid))
        alpha_c, omega_ref, R = ap, wp, R / 2
    end
    spread = maximum(abs(history[i][2] - history[j][2]) for i in 1:levels, j in 1:levels)
    return ap, wp, w2, spread, history
end
```

**Seed** (from the gate's `i_peak`, and the branch midpoint as an independent second opinion):

```julia
j_star   = argmin(abs.(L .- omega_F[i_peak]))
seed_mid = 0.5 * (alpha_L_u[j_star] + alpha_L_l[j_star])
alpha_c  = 0.5 * (F[i_peak] + seed_mid)
alpha_p, omega_p, w2, spread, hist = refine_pinch(alpha_c, omega_F[i_peak])
```

**Why this variant and not the alternatives you listed:** it reuses the temporal solver;
it never differentiates an eigenvalue at small `h` (with `h = 2e-5` the ~1e-8 eigen-noise
gives `ω″` noise of order 1e2 — I got 57.6 instead of 0.429, with a convincing-looking
convergence history); the fit residual and the closure check are built-in failure modes; and
`ω″` comes out for free, which §5 needs. Anything built on the *spatial* roots — solving
`∂D/∂α = 0` directly, or minimising the branch-to-branch distance — inherits the spatial
solver's ~5e-7 accuracy at the coalescence, which is 50× worse.

Cost: `levels*M = 72` temporal eigen-solves ≈ one Stage-1 iteration.

---

## 5. Verification block  (~30 lines, new, runs once after §4)

Report, and require:

| check | requirement | v4.1 run gives |
|---|---|---|
| radius consistency | `spread < 1e-5` over the three levels | 1.1e-6 |
| fit residual | ~1e-8 (grows with `num_modes`, see §6) | 6.7e-8 at nm=150 |
| non-degenerate | `abs(w2)` well away from 0 | 0.42941 |
| square-root law | `abs(a₊-a₋)` vs `2*sqrt(2δ/abs(w2))` at δ = 1e-2, 1e-3, 1e-4; relative error must fall like δ | 6.0e-3 → 7.0e-4 → 9.5e-5 |
| midpoint offset | `abs(mid-α_p)/abs(h)^2` constant | 0.757 / 0.760 / 0.766 |
| **α⁺/α⁻ origin test** | continue both roots to Im ω = Im ω_p + O(1), ladder fine enough that no single step exceeds ~0.5 in α; they must end on **opposite** sides of Im α = 0 | A → Im α = +11.3, B → −3.21 ✔ |
| Stage-1 quality | `min abs(F .- α_p) / h_F`, `abs(Im L - Im ω_p)`, `abs(max Im ω_F - Im ω_p)` | 0.32, 2.4e-7, 6e-7 |

The origin test is the one that makes it Briggs rather than "two curves got close", and v4.1
never did it. `postprocessing/pinch_refine.py --verify` is the reference implementation of
this whole block; port it or just call it on the JSON.

**Output:** write these to a **separate** file `*_pinch.json`, not appended to the history
array. Adding fields to only the last entry breaks MATLAB's `jsondecode`, which needs every
element of the array to have the same fields — `plot_contour_deformation.m` would return a
cell array instead of a struct array and fall over.

---

## 6. Separate commit — `num_modes`

`vib_ribbon/ANALYSIS_v5.1.md` established that `num_modes = 150` is past the divergent
threshold. My independent check on the pinch point: α_p, ω_p and ω″ are identical to 9
digits for `num_modes` = 50…120, `num_modes = 150` costs 9e-7 in α_p and 1.3e-8 in ω_p, and
the fit residual grows 5.4e-11 → 6.7e-8 → 2.4e-7 over nm = 50 → 150 → 180.

Recommendation: **`num_modes = 100`** (safely inside the plateau, no loss of resolved
eigenvalues, ~3× cheaper). But this changes ω_F and therefore the whole Stage-1 trajectory,
so commit it **separately** from §1–§5, run the regression either side, and record the
difference. Do not bundle it.

---

## 7. Acceptance criteria

v4.2 is done when, on the straight contour, Couette, Re = 2000, β = 0:

1. **Regression.** Stage 2 returns
   `α_p = 0.572519492 − 3.038860757i`, `ω_p = 0.298649115 − 0.572829695i`,
   `ω″ = 0.406438 − 0.138791i` to within 1e-6 / 1e-8. This is the ground truth this project
   has never had — `README_v4_series.md` says "there is no analytical pinch point for Couette
   flow"; there is now a numerical one.
2. **No upward L steps** except via the repair step. Count them; the count must be 0.
   (v4.1: 441 upward steps out of 901 in the stalled phase.)
3. **The clearance decreases monotonically** once L is within `lambda*delta_t` of the ceiling.
4. **The loop terminates** at k ≈ 700–900 with `stop_reason = "pinch neighbourhood reached"`.
   (v4.1: never; 1501 iterations and counting, ~51 % of them after the answer was reached.)
5. **Every check in §5 passes**, and each of them is capable of failing.
6. **Time to gate is not worse than v4.1's k ≈ 700** — §1 must not slow the descent down.

---

## 8. Explicitly *not* in v4.2

* **The F-update bookkeeping.** `alpha_i_cache` and `accepted` are computed and discarded
  (`alpha_i = alpha_i_smooth` runs unconditionally), and the radius-7 box filter is applied
  at full strength once per iteration regardless of `local_delta_t`, so the smoothing-to-drive
  ratio grows without bound as Δt shrinks. Both are real; the filter still removes 3.5–4.7e-4
  per iteration from the converged F, comparable to the net drive, and caps F's accuracy near
  4e-4. But F still reached h_F/3 of the saddle, so neither is why v4.1 stalled — and both
  change Stage-1 trajectories. Do them **after** §7 is passing, as a measured experiment
  against the regression number, one at a time.
* **Adaptive ω-resolution / refining L near ω_p.** Only worth it if you want `d_branch` itself
  to be small, which — see `ANALYSIS_v4.1.md` §3.1 — you should not.
* **Anything in the ribbon line.** Port §1–§5 to `briggsv4.6/v5.x_ribbon.jl` only once v4.2
  passes on the straight contour.
