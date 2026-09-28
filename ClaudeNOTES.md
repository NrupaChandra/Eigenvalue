# Project notes: what happened, in short

One-page summary of all the analysis notes written with Claude during this project
(Aug–Sep 2026). It replaces the three `claude_notes/` folders. The full, detailed notes
(numbers, derivations, figures) are kept in the Claude project **"Briggs"**; this file is the
overview.

---

## 1. What the project does

`01_briggs_couette/src/briggs*.jl` finds the **absolute-instability pinch point** of plane
Couette flow (Re = 2000, β = 0) with Briggs' method: two contours are deformed together —
L in the ω-plane is pushed down, F in the α-plane is pushed away from the spatial branches —
until the upper (α⁺) and lower (α⁻) branches pinch.

**Answer (verified independently, analytic dispersion relation, 22 digits):**

```
alpha_p = 0.572519491739732725 - 3.038860757129136972 i
omega_p = 0.298649114756735772 - 0.572829694808979836 i
omega'' = 0.406311314 - 0.138945509 i
```

Im ω_p = −0.5728 < 0 → **absolutely stable**. It is a genuine α⁺/α⁻ pinch (origin test done).

**The physics has never changed.** Every version from Kilian's original `briggs.jl` to
v7.1 uses the same Orr–Sommerfeld operator (150 Chebyshev modes, same boundary rows, same
temporal and spatial pencils). What changed is numerics, the contour-moving algorithm,
bookkeeping and infrastructure.

---

## 2. Version history (Couette)

| version | what changed | what we learned |
|---|---|---|
| `briggs.jl` | Kilian's original: coupled L/F gradient flow, exponential repulsion | the formulation is correct and was never replaced |
| v2 – v3.3 | adaptive vs fixed ζ, memory tracking, overlap guard | branches could select the same eigenvalue → guard added |
| v4 / v4.1 | fixed ζ = 4e-4, λ = 4, branches classified by side of F; v4.1 no zero-length steps | converged to the pinch but could not say so: `d_branch` measured the L-grid spacing; L step oscillated ±4e-3 |
| v4.2 | monotone L step | first run from which the pinch was extracted |
| v4.3 | L refined 100 → 478 points | crashed (`LAPACKException(150)`): `exp(ζ/d²)` overflowed to NaN → overflow cap added |
| v4.4 / v4.5 | bisection endgame for ω_i; adaptive L grid; eigvals-only, parallel spectra | refinement switched itself off too early (89 % of run wasted) |
| v4.6 | refinement reaches level 6 | grid solved; now ζ too large for the gap, F squeezed, ω_i drifted up |
| v4.7 | ζ ∝ d², checkpoint/resume, branch-tracking continuity guard | all three kept; best "reported" d ≈ 2.3e-3 |
| v4.8 | F refined (no separate .jl) | force clamp fired on 25 % of nodes → wiggles |
| v4.9 | tanh limiter, ε tied to ζ, line search | smoothest F ever; F horizontal in a tilted saddle costs ~3× in d |
| v5 | descent term on Im ω, trust region | best ω_i until v6.1; filter switched off, F went jagged |
| v5.1 | spatial pencil equilibration + fixes | **key finding:** `d_branch` was never measurable below ~3e-3 in any earlier version (badly scaled pencil); v5.1 run itself regressed |
| v6 – v6.3 | graded F, self-extending α_r window, smoothing gates | v6.1 got closest (|ω_i − Im ω_p| ≈ 1e-9); v6.2/6.3 unstable |
| **v7** | "reset": v4.7 + equilibration of **both** pencils, nothing else | temporal error 6.7e-9 → 1.9e-14; run can now see d below 3e-3 |
| **v7.1** | v7 + every iteration logged to JSONL | **best working version**; branches reach ~1e-3 … 4e-4 |

---

## 3. Lessons worth remembering

- **`d_branch` is not a convergence measure.** Use `|omega_i − Im omega_p|`
  (true d = 2·sqrt(2·that/|ω''|)).
- **Pencil scaling matters more than precision.** Equilibration (R·P·S, same eigenvalues)
  fixed errors that looked like "double-precision limits". 150 modes unequilibrated was the
  worst choice; 60 equilibrated is as accurate and ~12× faster.
- **Anything tuned in grid-index units breaks when the grid is refined** (filters, clamps,
  stencils).
- **The exponential barrier needs ζ scaled with the gap** (ζ ∝ d²) and an overflow cap.
- **Keep logs append-only and checkpoint** — a crash once destroyed a whole run.

---

## 4. Current status and open issues (v7.1)

- **Upper branch crosses F near the end (from iteration ~276, clearly from ~490).** The only
  guard is `omega_i ≥ peak(Im ω_F) + 1e-9`, but the peak is estimated from 100 F nodes
  (spacing 0.01). Near the pinch the true peak lies between nodes and is underestimated by
  ~1e-7, so F really rises above ω_i there. The last ~15 frames' small d is partly bought by
  this violation; ω_i still stays above Im ω_p, so the bound on ω_p holds.
  Fixes, smallest first: refine the peak with a few extra eigen-solves between nodes; add a
  side test (upper-branch nodes must stay above F) to the ω step; refine F near the pinch.
- The α-step acceptance check is not binding since v4 (the smoothed step is applied even if
  all 5 attempts fail).
- The side-of-F branch classification is only as fine as F's node spacing; the continuity
  guard (v4.7) covers it.
- Everything is Couette, Re = 2000, β = 0; saddle constants must be recomputed elsewhere.

---

## 5. Vibrating ribbon (`02_vib_ribbon`)

- **v4.6:** set up from Görtz §7.2.4 — the ribbon frequency ω₀ is a fixed pole on the real
  axis that L must keep enclosed; stability is a **sign test** on how branches cross the real
  α axis, not a collision test.
- **v5:** refining L was not enough; the branch drifted onto another eigenvalue between
  ω_i = −0.0498 and −0.0538.
- **v5.1:** the "tear" was the discretisation, not the flow — `num_modes = 150` was past the
  divergent threshold and carried 5.5e-3 error; `num_modes = 80` removes it.

---

## 6. Repository

Reorganised on 2026-09-28 into `01_briggs_couette`, `02_vib_ribbon`,
`03_eigenvalue_analysis`, `04_literature`, `05_briggs_ASBL`. Figures, videos and large JSON
logs are git-ignored and stay local.
