# The v4 ribbon series, in plain English

A reference for what each version does, what broke, and why. Written so it can be used as
the starting point for a v5.

---

## What the whole thing is trying to do

We want to know whether a disturbance in the flow **grows in place** (absolute instability)
or **gets washed downstream** (convective instability). The test for this is the Briggs
method, and it works like this:

- There are two planes. The **ω-plane** (frequency) and the **α-plane** (wavenumber).
- We draw a line in the ω-plane, called **L**, and a line in the α-plane, called **F**.
- For every point on L we solve an eigenvalue problem. That gives us many α values. Two of
  them matter: the one just **above** F (`alpha_L_u`, drawn red) and the one just **below**
  F (`alpha_L_l`, drawn green).
- We slowly push L downwards. F bends out of the way so it never gets crossed.
- Eventually the upper and lower α curves **collide**. That collision is the **pinch**. Where
  it happens tells us the answer.

The vibrating ribbon adds one thing: the flow is being forced at a fixed real frequency ω_f.
The correct L for that problem is not a straight line — it is a line with a **detour**: it
runs flat, goes up to the real axis, hops over ω_f in a little semicircle, and comes back
down. That shape is the "ribbon contour" or "Hankel contour".

---

## The single most important thing to understand

There are **two completely different designs** in this folder, and mixing them up causes
endless confusion.

**Design 1 — ribbon is decoration (v4.1 to v4.5).**
The calculation runs on a plain **straight** L, exactly like `briggsv4.jl`. The ribbon
contour is drawn separately once per frame, purely for the video. It has **no effect on
anything**. You could change ω_f to any value and every number in the result would be
identical. Only the picture would move.

**Design 2 — ribbon is real (`briggsv4_ribbon.jl`, and v4.6).**
The ribbon contour *is* L. Its α curves feed back and push F around. ω_f actually matters.

Versions v4.1–v4.5 are Design 1. Everything we fixed in them was fixing **the drawing**, not
the physics. That was deliberate (the descent was to be left alone), but it is easy to forget.

---

## Version by version

### `briggsv4_ribbon.jl` — the first attempt (abandoned)

Put the ribbon inside the loop straight away. Design 2.

It was dropped because the α curves were computed by asking, at each point separately,
"which eigenvalue is nearest?" — with no memory of the previous point. On the steep vertical
parts of the ribbon that jumps between different curves constantly. That noise then fed into
the calculation every single iteration, so the whole thing misbehaved.

Worth reading before writing a v5: it already solved several bookkeeping problems that come
back as soon as L stops being a straight line.

### v4.1 — ribbon becomes decoration

The calculation was reset to the plain straight-line version so it would behave, and the
ribbon was demoted to a drawing.

**What broke:** the drawing used the same no-memory method as above. The red curve jumped
around wildly — on the vertical parts, essentially every point landed on a different curve.
Jumps of up to 1.5. That is the noise visible in the v4.1 video.

Also, the old videos carried a black "exact pinch" dot. That belonged to a completely
different test case and means nothing here. There is **no analytical pinch point** for
Couette flow, and none has been drawn since.

### v4.2 — give the drawing a memory

Now each point starts from where the previous point was, and picks the nearest eigenvalue to
that. Vertical parts got more sample points (25 → 60).

**Much better, but still wrong.** The upper curve has to travel a long way up the vertical
part — about 2.7 — and 60 steps meant each step was about 0.06. The eigenvalues themselves
sit about 0.06 apart. So "nearest to where I was" stopped being reliable, and it stepped onto
a neighbouring curve without noticing.

Proof it was wrong: the point ω = 0.24 on the real axis is the *same point* in every frame,
so the answer there can never change. Frame 1 said 0.4012 + 0.0740i. Frame 30 onwards said
0.556 + 0.215i. Different eigenvalue, so it had hopped.

The lower curve only moves about 0.2 over the same stretch, so its steps stayed tiny and it
was fine. That is why only the upper curve looked broken.

### v4.3 — smaller steps, and a safety check

Three changes:

- Vertical parts 60 → 160 points, so each step is ~0.017 instead of 0.06.
- Start from **both ends** instead of one, so neither side has far to go.
- A guard at each step: if the step is too big, or if the nearest eigenvalue is not *clearly*
  nearer than the runner-up, cut the step in half and try again.

That worked. But v4.3 also did something it should not have: on the flat parts, it **copied**
the numbers from the main calculation instead of working them out. The numbers were correct
(same input, same answer), but it meant the grey reference line and the red line were
identical *by construction* over 36% of the drawn curve. So they could not disagree, and the
comparison stopped meaning anything.

**Separate problem, worth remembering:** the v4.3 run crashed at the end. The file was being
completely rewritten from scratch every iteration — about 4 GB of writing to save 26 MB — and
one write came up short, cutting the file off at exactly 16 MB. The fix was to write each
frame onto the end of the file instead, keeping a valid `]` on disk at all times. Never go
back to rewriting the whole file.

### v4.4 — nothing copied

Every point works out its own answer. Flat parts got twice as many points so the steps stay
small. The starting point is chosen using the same rule the main calculation uses, but
recalculated fresh.

**What broke:** the red curve slowly **wandered away from F**. Its typical distance from F
went 0.035 → 0.79 → 2.08 while the grey line stayed around 0.6.

That matters because the curves we need are the ones sitting *next to* F — those are the ones
that eventually collide. A curve two units away is a real eigenvalue curve and completely
irrelevant.

Two causes:

1. **Bad starting point.** It started at the far-left edge, ω_r = 0. At that exact spot the
   problem is mirror-symmetric, so the eigenvalues come in identical twin pairs — one at
   +0.78 and one at −0.78, same height. The tiebreak was "which is slightly closer to F", and
   the margin was **2.6%**. Basically a coin flip. It flipped, and every later frame faithfully
   followed the wrong twin.
2. **No constraint.** Nothing stopped it drifting once it had left.

### v4.5 — start in a sensible place, and keep checking

Two fixes, both using rules that already existed in the main calculation:

- **Start in the middle** of the contour (the same place the main calculation starts) instead
  of at the edge. There the curve sits right next to F and there is no twin to confuse it.
  Then walk outwards in both directions.
- **At each step, check the side of F** before accepting a point. If nothing is on the right
  side, fall back — and **count** every fallback, so the check can never quietly hide
  something.

**Result: the left half is now perfect** — red matches grey to 1e-6, zero fallbacks all run.

**But the right half (after the detour) still differs, by up to 2.85.** This turned out not
to be a bug.

---

## The finding: it is a real branch point

To check this properly, the eigenvalue solver was **rebuilt from scratch in Python**,
independently of the Julia code, and agreed with it to 1e-7.

Then the same curve was followed from the same starting point along two different roads to
the same destination:

| road | answer |
|---|---|
| up the riser, over the arc, down the other side | 0.53691 + 0.15612i |
| straight along the bottom | 0.76600 + 0.00051i |
| what v4.5 produced | 0.53691 + 0.15612i (matches road 1) |
| what the main calculation has | 0.77649 + 0.00372i (matches road 2) |

Both are correct. They are answers to different questions.

Think of a spiral staircase. Walk around the central pole and come back to the same spot on
the map — you are one floor up. There is a **branch point** (the "pole") sitting inside the
gap between the two vertical parts, at roughly **ω ≈ 0.25 − 0.05i**. The ribbon contour goes
*over* it; the straight contour goes *under* it. Different side, different floor.

Confirmation: the `loop_mono` number switches on at exactly ω_i = −0.051, which is the depth
at which the detour first wraps around that point.

So the red/grey disagreement is not an error. It is the most interesting number in the run —
and it is exactly the kind of structure the ribbon detour exists to reveal.

### v4.6 — close the loop

Because of that finding, the straight-line problem and the ribbon problem are **genuinely
different problems**. So v4.6 puts the ribbon contour back inside the calculation, the way
`briggsv4_ribbon.jl` originally intended — but now with the tracking machinery that actually
works.

The physics functions are untouched: the potential, the descent rules, the acceptance tests,
the pinch test are all byte-identical to `briggsv4.jl`. Six functions changed, and all six
are bookkeeping:

- the contour is now the ribbon shape;
- branches come from the v4.5 tracker instead of the old no-memory one;
- the pinch distance is measured **every point against every point**, because on a ribbon
  contour matching point 5 to point 5 is meaningless;
- three places had `100` hardcoded as the number of points; the ribbon has 570.

---

### v5 — refine the contour, and add a check that can fail

The v4.6 post-mortem (`../ANALYSIS_v4.6.md`) found the upper branch hopping at iteration 14
on the **left horizontal** at ω_r ≈ 0.146 — before the detour, not in it — so iterations
14–1001 were on the wrong curve. Two candidate explanations were tested and ruled out:

- **Not the side-of-F filter.** Remove it and the continuation gives the same numbers to four
  decimals, with zero fallbacks.
- **Not a branch point.** Walk a closed rectangle in the ω-plane between the two horizontal
  levels: at the run's own spacing you come back 0.084 away from where you started, but at
  3.5× finer spacing you come back to **exactly** the same value. There is nothing enclosed.

It was under-resolution. 48 nodes on the horizontal gives 0.65360 + 0.08051i at ω = 0.24;
96 and 192 both give 0.40647 + 0.00375i, which is the smooth continuation of the previous
frame. And the tracker never noticed: the biggest accepted step was 0.0151 against a guard of
0.03.

So v5 changes four things and no physics:

1. **N 100 → 300.** Horizontal spacing 8.4e−4, three times finer than the 1.25e−3 where the
   answer converged. `da_max` 0.03 → 0.012 to match. **Do not touch `num_modes`** — that is a
   different knob and it goes the other way: max|D4| grows like n⁸, and at 220 the spectrum
   near α ≈ 0.31 is already meaningless.
2. **A check that can fail.** The riser tops and the arc are at fixed ω, so α there must be
   the same in every frame. `fixedw_u/l` measures it; the run stops if it moves. This is the
   v4.2 test, which existed and was dropped. It catches iteration 14 instantly, where all four
   of v4.6's diagnostics were clean.
3. **60 iterations, not 1000.** α at the arc took two values in 1001 frames. Everything after
   the hop was noise.
4. **Report the ribbon answer, not the pinch.** α_u(ω_f), α_l(ω_f), and whether either crossed
   Im α = 0. `pinch_tol` stays, but only as the absolute-instability alarm.

Also added: an optional grid-convergence check that retracks the same contour at double
sampling and compares. That is the thing that would have measured the v4.6 error directly.

**Still open after v5:** the eigenvalue conditioning (α ≈ 0.31 is uncertain by about its own
branch separation), F drifting bodily down to Im α ≈ −1.15, and the radius-7 smoothing.

---

## Things to know before writing a v5

**The calculation stalls before it pinches.** This is the biggest open problem and it is
older than all of these versions. The contour stops descending — ω_i freezes around −0.5725
and the gap between the curves flattens at about 0.12, never reaching the 1e-4 target. Look
for `jump=0.000e+00` in the log: it appears from about iteration 137 and gets more frequent.
The cause is that L cannot go below ω_F, and ω_F keeps rising to meet it. Nothing in v4.1–v4.6
addresses this, because it lives in the descent, which was off-limits.

**F is heavily smoothed, and this was never asked for.** In the original `briggsv4.jl`, F gets
a rolling average of radius 7 applied **every iteration**, and the accepted trial value is
computed and then thrown away — the smoothed version is always used regardless. Over 300
iterations that is a lot of smoothing. Worth examining; it may be related to the stall.

**Small `v_g` inconsistency.** With `v_g = 0` the spatial solver matches Orr–Sommerfeld
exactly (verified symbolically). With `v_g ≠ 0` there is a leftover term. Harmless now, would
matter in a moving frame.

**What the diagnostics mean.**

| name | meaning | healthy value |
|---|---|---|
| `xcheck` | red vs grey at the same ω | small, but not forced to zero |
| `sideclass` | tracked curve vs the plain "nearest on the correct side" rule | small |
| `nfallb` | times the side check found nothing and fell back | 0 |
| `loop_mono` | going over the detour vs along the bottom | non-zero means a branch point is enclosed |
| `refined` | steps that were cut in half and retried | small |

**Lessons that cost real time:**

- Never copy numbers from one part of the calculation into another for display. It destroys
  the ability to check anything. (v4.3)
- The starting point matters more than the tracking. v4.4's tracking was flawless and the
  answer was still wrong. (v4.4)
- When two things disagree, find out which is right before making them agree. Twice the
  "broken" curve turned out to be correct. (v4.5)
- Write results incrementally, never by rewriting the whole file. (v4.3 crash)
- A good check is one that *can* fail. If two numbers cannot disagree, comparing them tells
  you nothing.
