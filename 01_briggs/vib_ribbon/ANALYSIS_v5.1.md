# v5.1 — the tear was the discretisation, not the flow

The first v5.1 run stopped at iteration 14 on the tear detector, at ω_i = −0.05422, with a
max along-contour step of 3.587e−2 against tear_tol = 0.024. The detector was right to fire.
What it was detecting was the eigenvalue solver.

**`num_modes = 150` was past the divergent threshold of the discretisation and carried
5.5e−3 of error. That error is the entire "unresolvable branch crossing".** Setting
`num_modes = 80` removes it. Everything below is measured with the independent Python
solver in `analysis/os_spatial.py`.

---

## 1. What the run showed, and the one thing that was new

Frames 1–12 clean (maxstep ≈ 1.75e−3). Frame 13: `sideclass_u` 1.03e−1, maxstep 5.1e−3.
Frame 14: maxstep 3.587e−2 → stop.

New in v5.1, and worth recording: **the symmetric contour produces two mirror-image trouble
spots**, at ω_r = +0.1444 and ω_r = −0.1452. That is exactly the Orr–Sommerfeld symmetry
for a real base flow — if (α, ω) solves it then (−ᾱ, −ω̄) does — so a feature at
ω = 0.148 − 0.050i must have a twin at ω = −0.148 − 0.050i. Checked numerically: at
ω_r = −0.14524 the run gives α = −0.297896 + 0.005341i, and −conj(0.297896 + 0.005341i)
= −0.297896 + 0.005341i. Exact. The symmetric contour is doing its job and the solver is
self-consistent.

(The tear happened to open on the negative side first, simply because the tracker marches
outward from the seed at ω_r = +0.12 and that direction reaches its crossing first.)

## 2. The decisive test: refine and watch

Track α across the crossing at frame 15's level, ω_i = −0.05422, ω_r ∈ [0.135, 0.160],
refining Δω_r. If the jump were an under-resolved sharp turn, the max step must fall like
1/factor.

**At `num_modes = 150`:**

| Δω_r | factor | max step | median step | max/median |
|---|---|---|---|---|
| 8.35e−4 | 1× | 0.00637 | 0.00326 | 2.0 |
| 2.09e−4 | 4× | 0.00654 | 0.00083 | 7.8 |
| 1.04e−4 | 8× | 0.00863 | 0.00048 | 17.9 |

The median converges perfectly — the smooth part of the curve is fine. The **max step does
not fall; it grows.** On the finest grid the values jitter ±0.004 with no trend
(0.3479, 0.3455, 0.3445, 0.3501, 0.3509, 0.3439, 0.3524, …). That is not a curve being
sampled, it is noise being sampled.

## 3. Where the noise comes from

`max|D4|` grows like n⁸ (1.4e15 at n = 100, 3.7e16 at 150), so past a threshold the
discretisation diverges — extra modes buy round-off, not accuracy. One eigenvalue of the
crossing pair at ω = 0.15075 − 0.05422i:

| num_modes | α | error vs plateau |
|---|---|---|
| 50 | 0.34692327 − 0.01213617i | |
| 60 | 0.34692303 − 0.01213601i | plateau, spread **4e−7** |
| 70 | 0.34692327 − 0.01213612i | |
| 80 | 0.34692285 − 0.01213514i | |
| 100 | 0.34692320 − 0.01212789i | |
| 110 | 0.34719170 − 0.01184683i | 3.9e−4 |
| 130 | 0.34728424 − 0.01193706i | 4.1e−4 |
| 140 | 0.35001044 − 0.01209022i | 3.1e−3 |
| **150** | 0.35221286 − 0.01057364i | **5.5e−3** ← what the code used |

The boundary-row scaling (`−200i`) is not implicated: varying it 20 → 2000 moves the
eigenvalues by ~1e−4.

## 4. With `num_modes = 80`

Same refinement test:

| Δω_r | factor | max step | median | max/median |
|---|---|---|---|---|
| 8.35e−4 | 1× | **0.00375** | 0.00321 | 1.2 |
| 2.09e−4 | 4× | 0.00099 | 0.00081 | 1.2 |
| 1.04e−4 | 8× | 0.00053 | 0.00040 | 1.3 |
| 5.22e−5 | 16× | 0.00030 | 0.00020 | 1.5 |

Falls like 1/factor, max/median steady at 1.2–1.5 — the curve is uniformly smooth, with no
sharp feature anywhere. At v5.1's own spacing the worst step is 0.00375, against
`DA_MAX = 0.012` and `tear_tol = 0.024`. The tear detector will stay quiet.

**The avoided crossing is real but ordinary.** Minimum branch separation over the same
ω-window is **0.0149** at num_modes = 80 (it read 0.0027 at 150). Against ~1e−5 accuracy
that is a margin of about 1500. There is no branch point on the contour, no sheet exchange,
and no need to indent L.

**80 loses nothing.** Identical count of eigenvalues with |α| < 5 (23) at num_modes
80/100/120/150. α at the arc: 0.41381153 (nm 80) vs 0.41381657 (nm 150) — 5e−6, so the
ribbon answer is unchanged. At the contour ends ω = ±0.5 − 0.05i, nm 80 and 120 agree to
1.2e−5. And the 160×160 pencil is roughly 6× cheaper to solve than the 300×300 one, so the
run gets faster as well as more accurate.

---

## 5. What this invalidates

I have to withdraw several earlier conclusions, all of which rested on num_modes = 150:

- **"The branches close to 0.0027, below the noise floor, so no tracker can resolve it."**
  Wrong. The gap is 0.0149 and the accuracy is 1e−5.
- **"There is a branch point / near-collision at ω\* ≈ 0.148 − 0.050i."** There is an
  avoided crossing there, but nothing pathological, and the "0.0027" was noise on top of it.
- **"Keep num_modes at 150; raising it makes things worse."** (`ANALYSIS_v4.6.md` §5,
  and the v5 header.) Raising it is indeed bad — but so is 150. I tested 100/130/150/180/220
  back then, saw the values scatter, and read the scatter as an irreducible floor instead of
  testing *downward* to find the plateau. That was the error, and it cost v5 and the first
  v5.1 run.

What survives: the fixed-ω check, the grid-convergence check, the frame filter, the tear
detector, and the symmetric contours. The tear detector in particular did its job — it
refused to produce 985 more frames of a torn curve, which is how this got found at all.

## 6. Still open

- The `alpha_i_trial` computed, acceptance-tested and then discarded in favour of
  `alpha_i_smooth` — a radius-7 rolling average applied to F every iteration. Nobody chose
  this; `README_v4_series.md` has flagged it since v4.
- F drifting bodily downwards rather than staying a deformation of the real α axis.
- Whether the ribbon should have a second detour at −ω_f (a real forcing cos ω_f t has poles
  at both ±ω_f; Görtz's e^{iω_f t} has only one). Now visible as an asymmetry on the
  symmetric contour.
