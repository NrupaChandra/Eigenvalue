# briggs v4.3 — diagnosis of the `LAPACKException(150)` crash

Run: `julia briggsv4.3.jl`, crashed after 121 printed iterations.
Data used: `results/json/contour_iteration_v4.3.json` (122 entries; entry 1 is the
pre-loop initial state written at src line 623–640, so entry `it` = loop step `k = it-1`).

## Short version

`α` on the F contour became **NaN** during the alpha update of loop step k = 122.
`eigen(A,B)` was then called on matrices full of NaN and LAPACK's QZ iteration
returned `INFO = 150`, which Julia wraps as `LAPACKException(150)`.

The NaN came from **overflow of the repulsion potential** `exp(ζ/(|α_F−α_L|² + ε))`
in `d_d_alpha_r_Phi_F` / `d_d_alpha_i_Phi_F`. It is not a bad eigenvalue problem,
not the ill-conditioning of the D4 collocation matrices, and not a Distributed/pmap
issue. Worker 6 simply happened to draw one of the poisoned α values.

## 1. What `LAPACKException(150)` means

The stack is `eigen → ggev3! → chklapackerror(150)`. `ZGGEV3` documents

- `1 ≤ INFO ≤ N` — the QZ iteration in `ZHGEQZ` failed; eigenvalues `INFO+1 … N` are
  usable, the rest are not,
- `INFO = N+1` — other failure in `ZHGEQZ`,
- `INFO = N+2` — failure in `ZTGEVC`.

Here the temporal problem is `Val(:omega_collocation)`, so `N = num_modes = 150`.
`INFO = 150 = N` means the QZ iteration failed on the very first deflation — nothing
was computed at all.

## 2. It is NaN input, not a genuine QZ convergence failure

Reproduced outside Julia by rebuilding the identical `A`, `B` (same Chebyshev
recursion, same `-200i` boundary rows, `Re = 2000`, `β = 0`, `num_modes = 150`) and
calling LAPACK directly:

| α passed in | `ZGGEV3` INFO |
|---|---|
| `0.7 - 2.09im` (a normal contour point) | 0 |
| `0.7 + NaN·im` | **150** |
| `NaN + NaN·im` | **150** |

Exactly the observed code. To exclude the alternative (well-formed but nasty pencil),
`ZGGEV3` — the same routine Julia calls, eigenvectors included — was run on **2800**
α values: the 600 α actually used in the last six saved iterations plus 2200 random
points in `Re α ∈ [0,1]`, `Im α ∈ [−2.30, −1.75]`, the strip the contour occupies.
**Zero failures.** (`ZGGEV` via SciPy on a separate 1189-point grid: also zero.)
The pencil is badly scaled — `‖A‖ ≈ 7·10¹²`, `cond(A) ≈ 7·10¹⁰`, `cond(B) ≈ 2·10⁸` —
but it solves reliably.

## 3. Where the NaN is produced

`s_alpha = 2`, `ζ = 4e-4`, `ε = 1e-10`, so with `d = |α_F − α_L|`

```
term  =  −ζ·s·d^(s−2)/(d^s + ε)²  ·  exp(ζ/(d² + ε))  ·  Δ  ·  w
```

`exp` overflows at `ζ/(d²+ε) > 709.78`, i.e. **d < 7.506e-4**; including the
prefactor and a typical quadrature weight the whole term overflows at
**d ≈ 7.55e-4**.

State at the last successful step (k = 121, JSON entry 122):

| quantity | value |
|---|---|
| closest interior F point to `alpha_L_u` | `d = 1.1848e-3` (F #76 ↔ `alpha_L_u` #198, 1-based) |
| exponent `ζ/(d²+ε)` | **284.9** (overflow at ~702) |
| `max|rhs_j|` before `clamp` | **3.28e126** |
| margin to overflow | factor **1.57** in distance |

And the sign structure at that point is already the fatal one:

```
alpha_i_r = −3.086e-1
gr = d_d_alpha_r_Phi_F = −1.954e126      alpha_i_r·gr = +6.030e125
gi = d_d_alpha_i_Phi_F = +4.194e126             −gi  = −4.194e126
numerator = alpha_i_r·gr − gi
```

Both halves share the same `exp` factor, so they overflow together, and they have
**opposite signs**: once `d` drops below 7.5e-4 the numerator is `(+Inf) + (−Inf) = NaN`.

The guard on line 1038 does not catch it — `clamp(NaN, −100, 100) === NaN` in Julia.
`clamp` only rescues the `±Inf` case, which is why the run survived the earlier close
approaches.

Then:
`rhs_j = NaN` → `alpha_i_trial[j] = NaN` → `rolling_average_filter(·, 7)` smears it
over ~15 points → line 1082 assigns `alpha_i = alpha_i_smooth` **unconditionally**
(the accept/reject loop cannot block it) → `F = contour_F()` contains NaN →
`pmap` → `eigen` → `LAPACKException(150)`.

## 4. Why v4.3 and not v4.2

The eigenvalue code, the potential, the clamp and the smoothing are byte-identical
between the two files. The only substantive change is `omega_r`: 100 → **478** points
(200 each over the base cells 47–57 and 57–67).

`alpha_L_u`, `alpha_L_l` are the point sources of `Φ_F`, and they inherit `length(L)`.
Sampling the same branch curve 5× more finely makes the *discrete* nearest distance
from F to the branch point cloud systematically smaller, and `exp(ζ/d²)` is
hypersensitive to that. Direct test at k = 121, subsampling the same 478-point branch
back to ~100 points:

| branch sampling | min distance | `|∇Φ_F|` at the closest F point |
|---|---|---|
| 478 pts (v4.3) | 1.185e-3 | 4.63e126 |
| ~100 pts (v4.2-like) | 1.409e-3 | 8.75e90 |

A 19 % change in distance, **36 orders of magnitude** in the gradient. The near field
of `Φ_F` is not a convergent quadrature — refining the branch makes it diverge.

Run comparison:

| | v4.2 | v4.3 |
|---|---|---|
| saved iterations | 690 | 121 (crash) |
| `length(L)` | 100 | 478 |
| closest interior F→branch over the whole run | 1.585e-3 | **1.185e-3** |
| max exponent `ζ/(d²+ε)` reached | 159 | **285** |
| iterations with d < 2e-3 | 3 | 6 |
| final `ω_i` | −0.5728 (pinch) | −0.4819 |

v4.3 reached a closer approach in 121 steps than v4.2 did in 690, at 84 % of the
pseudo-time.

## 5. Secondary findings

- **Branch-tracking discontinuity.** `alpha_L_u` at JSON entry 122 has one jump of
  **1.227** between 0-based indices 199→200 (0.7604−2.0878i → 0.5957−0.8715i), inside
  the first refined block where neighbouring steps are 1.7e-3, i.e. ~700× the local
  step and 20× the median `|dα/dω|`. v4.2's final `alpha_L_u`/`alpha_L_l` have no such
  jump. The crash-adjacent closest approach is at index 197, right beside it.
  `couetteflow_spatial_sing_mode_comparison` accepts only eigenvalues strictly on one
  side of F (`proj > 0`); when the branch runs within ~1e-3 of F the side test is
  ill-conditioned and the tracker can be forced onto a different eigenvalue.
- **The accept/reject loop at lines 1016–1080 has no effect on the result.**
  `alpha_i_cache` is filled on acceptance and never read; line 1082 always applies
  `alpha_i_smooth`. The loop only shrinks `local_delta_t`.
- **`local_delta_t` is local.** The global `delta_t` stays 1e-3 for the whole run
  (only the `omega_repaired` branch can change it), so a rejection never carries a
  reduced step into the next iteration.
- **`pinch_tol = 1e-4` remains unreachable.** With `d_branch ≈ 2√(2|ω−ω_p|/|ω''|)` and
  `|ω''| ≈ 0.4295`, the floor set by the L grid is ≈ 0.217 at Δω_r = 5.05e-3 and
  ≈ 0.0486 at the refined Δω_r = 2.54e-4. Reaching 1e-4 needs `|ω−ω_p| ≤ 5.4e-10`,
  i.e. Δω_r ~ 1e-9. The refinement buys a factor 4.5, not the factor 500 needed. This
  is the same objection recorded in the v4.2 header.
- **Re-running the script wipes the results.** Lines 623–640 execute before the loop
  and truncate `contour_iteration_v4.3.json` to a single entry. Back up the 122-entry
  file before the next run.

## 6. Fix applied (v4.3a, in `src/briggsv4.3.jl`)

Physics, discretization, potential parameters and time stepping are untouched. The
200/200 refinement of L over `omega_r_base[47..57..67]` stays exactly as intended.

1. **Overflow guard on the alpha-side potential** — the root cause.

   ```julia
   const EXP_ARG_MAX = 400.0
   expc(x) = exp(x < EXP_ARG_MAX ? x : EXP_ARG_MAX)
   ```

   used in `phi_F`, `d_d_alpha_r_phi_F`, `d_d_alpha_i_phi_F`. The cap starts to
   engage at d = 1.0e-3, below the 1.185e-3 the run ever reached, and `rhs_cap = 100`
   already saturates from d ≈ 9.5e-3, so the cap is invisible to the dynamics.
   The **omega-side functions are deliberately left uncapped**: their exponent runs to
   ~1e4 by design (531/691 of v4.2's iterations), and the resulting NaN is already
   sanitised by the `isfinite(x) && abs(x) < 10.0` filter on the omega candidates.
   Capping them would change which candidates survive.

2. **NaN-safe clamp**, line ~1098: `rhs_j = isnan(rhs_j) ? 0.0 : clamp(rhs_j, -rhs_cap, rhs_cap)`.
3. **Non-finite alpha updates are discarded**, not applied: `global alpha_i = copy(alpha_i_smooth)`
   is now guarded by `all(isfinite, alpha_i_smooth)`, with a warning line if it fires.
   Nothing non-finite can reach `contour_F()` and hence `eigen(A,B)`.
4. **Iteration budget 700 → 2000** (v4.2 needed 690).
5. **`const RESUME`** — `false` reproduces the original behaviour exactly (the JSON is
   overwritten and the run restarts); `true` reloads the last entry of `filename` and
   appends from there, so a long run survives a restart. Back up
   `contour_iteration_v4.3.json` before any run with `RESUME = false`.

### Verification

- Replaying all 121 historical iterations with and without the cap: the *applied*
  `rhs_j` (after `clamp`) is **bit-identical**, `max|Δ| = 0.000e+00`. The fix cannot
  have changed the trajectory that already ran.
- Failure scenario, shrinking the closest interior gap at k = 121 below the overflow
  radius:

  | d | exponent | old `rhs_j` | applied | new `rhs_j` | applied |
  |---|---|---|---|---|---|
  | 1.185e-3 | 285 | −3.28e126 | −100 | −3.28e126 | −100 |
  | 8.89e-4 | 507 | −1.32e223 | −100 | −7.37e176 | −100 |
  | 7.11e-4 | 791 | **NaN** | crash | −1.44e177 | −100 |
  | 1.19e-4 | 28293 | **NaN** | crash | −3.07e179 | −100 |

  Note the third row: the NaN-safe clamp alone would map NaN → 0, i.e. it would switch
  the repulsion **off** at exactly the moment the contour touches the branch. The
  exponent cap is what keeps the intended −100 push. Both are needed, in that order.
- Block/bracket structure of the patched file checked against the unmodified v4.2 and
  v4.3 sources.

### Still open (not changed here)

- `pinch_tol = 1e-4` cannot fire (§5), so the run will use its whole iteration budget
  rather than reporting convergence.
- The branch-tracking discontinuity in `alpha_L_u` (§5) is untouched.
- The near field of `Φ_F` is still `length(L)`-dependent. The structural fix, if the
  contour is meant to keep approaching the branch, is to measure the distance from
  `F[j]` to the nearest *segment* of the `alpha_L_u`/`alpha_L_l` polyline rather than
  to the sampled points, or to floor `d²` at a fixed value. Otherwise every further
  refinement of L moves the blow-up threshold closer — the guard only makes that
  survivable, it does not remove it.
