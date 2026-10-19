# Point-cloud spatial branches in `pointcloud_test.jl`

Session of 2026-10-13. This replaces seed-and-walk tracking of `alpha_L_u` / `alpha_L_l` with a frontier built from pooled spatial spectra.

## Motivation

The previous method (`contour_alpha_L_conti`, v4.5–v4.7):

1. took the eigenvalue closest to F at one seed omega (4/6 of L, or `omega_r_base[25]` in v4.5),
2. walked forwards and backwards along L, at each omega taking the eigenvalue nearest to the previous one on the requested side of F.

Two problems:

- **Jumps.** A single wrong pick at one omega makes every later omega follow the wrong mode. CHANGE 9 (`TRACK_SIDE_SLACK`) reduced this but did not remove it.
- **Cost and structure.** The selection is an inherently serial chain. The v4.4 reference path additionally ran one `pmap` per step.

## New approach

1. **Parallel spectra.** One `pmap` over all of L. Each worker solves one spatial eigenproblem and classifies the eigenvalues as above/below F with the existing normal-projection test.
2. **Two clouds.** All classified eigenvalues of all omegas are pooled into an upper and a lower cloud. Implausible eigenvalues are removed.
3. **Frontier.** In each cloud, the outer line facing F is taken as the branch trajectory.

No seed, no chain: the result at one omega cannot affect any other.

### Coordinates

Each eigenvalue `e` is described relative to F by

- `s`: arclength position of its orthogonal foot point on F,
- `d`: distance from `e` to the F polyline.

In the (s, d) plane F is the axis d = 0, and the outer line facing F is the **lower envelope of d over s**.

## Functions

### Worker side

| Function | Purpose |
|---|---|
| `arclength(F)` | Cumulative arclength of the polyline F, `s_F[1] = 0`. |
| `foot_on_F(e, F, s_F, k)` | Projects `e` onto the two segments adjacent to its nearest node `k`. Returns `(d, s)`. Interior segments clamp the projection parameter to [0, 1]; the two end segments are not clamped outwards, so eigenvalues beyond F's ends get s < 0 or s > len(F) instead of collapsing onto the end nodes. |
| `spatial_payload(omega, F, normals_F)` | Unchanged eigensolve and side test (sign of `Re(conj(n_k)·(e − F[k]))`). Now also returns `up_s` / `lo_s`, and `up_d` / `lo_d` are the distance to the polyline instead of to the nearest node. |

Projection parameter on segment a→b: `t = Re(conj(b−a)·(e−a)) / |b−a|²`, foot `a + t(b−a)`, `s = s_F[i] + t|b−a|`.

### Main process

| Function | Purpose |
|---|---|
| `contour_alpha_L_cloud(L)` | Driver: `pmap` of `spatial_payload`, then `gather_cloud` and `cloud_frontier` for each side. Stores the clouds in `cloud_u` / `cloud_l`, sets `track_jump_u/l` (largest step along each frontier), resets `track_overrides` to 0. |
| `gather_cloud(payloads, side)` | Pools one side over all omegas into flat arrays `pts, s, d`, applying the filters below. The omega of each point is discarded. |
| `cloud_frontier(pts, s, d; nbins)` | Builds the trajectory (algorithm below). |
| `smooth_polyline(z, w)` | Optional moving average with half-width `w`, shrinking near the ends so the endpoints stay fixed. `w = 0` returns the input unchanged. |
| `contour_alpha_L_init(L)` | Calls the cloud method, then clears the worker copies of `alpha_L_u/l` as before. |
| `contour_alpha_L_conti(L)` | Alias for `contour_alpha_L_cloud(L)`, so the iteration loop needs no change. |
| `contour_alpha_L_walk(L)` | The previous seed-and-walk method, kept for comparison. |
| `plot_cloud()` | Scatter of both clouds with F, `alpha_L_u`, `alpha_L_l`. |

`contour_alpha_L_conti_ref` (v4.4 reference path) is untouched.

## Filters (`gather_cloud`)

A point is kept only if all hold:

- `|α| ≤ CLOUD_ABS_MAX`: removes spurious modes of unreasonable magnitude,
- `d ≤ CLOUD_DIST_MAX`: points this far from F cannot be on the frontier,
- `−CLOUD_S_PAD·len(F) ≤ s ≤ (1 + CLOUD_S_PAD)·len(F)`: foot point on F. Beyond F's ends the near branch runs out and far modes would form the frontier.

## Frontier algorithm (`cloud_frontier`)

**A. Minimum per bin.** The s-range is cut into `nbins` equal bins. Each bin keeps its eigenvalue with the smallest d; empty bins are skipped. The result is an s-ordered sequence of actual eigenvalues.

**B. Spike removal (two-sided valley test).** A bin with no sample of the near branch picks up a farther mode, which appears as a spike or zig-zag away from F. For each frontier point, take the smallest d among frontier points within `CLOUD_WIN` bin widths to the left (`mL`) and to the right (`mR`, `Inf` if none). Then

```
keep  ⇔  d ≤ (1 + CLOUD_SPIKE_REL) · max(mL, mR)
```

i.e. drop a point only if **both** sides hold a clearly closer point. Repeat until nothing is removed.

| Situation | Verdict |
|---|---|
| far mode in a gap between two branch samples (lower on both sides) | dropped |
| flank of a dip towards F (lower on one side only) | kept |
| bottom of a dip | kept |
| gentle hump on a smooth branch (neighbours < 10 % lower) | kept |
| end point (no neighbour on one side) | kept |

A one-sided test (compare against the overall local minimum) cannot tell a spike from the flank of a steep dip; it was tried first and either cut real flanks or let the zig-zag through.

**C. Optional smoothing** with `smooth_polyline(·, CLOUD_SMOOTH)`.

## Parameters

| Constant | Default | Meaning |
|---|---|---|
| `CLOUD_ABS_MAX` | 10.0 | drop \|α\| above this |
| `CLOUD_DIST_MAX` | 2.0 | drop eigenvalues farther than this from F |
| `CLOUD_NBINS` | 0 | number of bins along F; 0 → `length(L)` |
| `CLOUD_WIN` | 3 | spike-filter half-window, in bins |
| `CLOUD_SPIKE_REL` | 0.1 | relative margin of the valley test |
| `CLOUD_SMOOTH` | 0 | moving-average half-width; 0 = raw eigenvalues |
| `CLOUD_S_PAD` | 0.0 | allowed foot-point overshoot beyond F's ends, as a fraction of len(F) |

## Test results

Setup: Re = 5000, β = 0, 150 Chebyshev modes, straight F on [−1.5, 1.5], 150 omegas on [−1.5, 1.5] + 1e-3 i, 5 workers. Total runtime of the file ≈ 30 s.

| | cloud size | frontier points | max step |
|---|---|---|---|
| upper | 5210 | 116 | 0.054 |
| lower | 150 | 48 | 0.063 |
| old walk, upper / lower | – | 150 / 150 | 0.044 / 0.063 |

Observations:

- **Lower branch** is identical to the old walk.
- **Upper branch** matches the old walk except that it dips to the **Tollmien–Schlichting mode**, α ≈ 1.06 + 0.0065i at ω ≈ 0.29 (and its mirror at α_r ≈ −1.06). The walk never followed this mode and stayed on the wall-mode branch, which is exactly the kind of mode jump the new method avoids.
- Near α_r ≈ 0 the old walk had a small loop; the frontier is clean there.
- Iterations during development:
  1. no s-filter → frontier climbed onto far modes beyond F's ends (max step 0.70); fixed with the s-range filter,
  2. one-sided spike filter → zig-zag at |α_r| ≈ 0.8 and 1.25, where TS samples are sparser than the bins; fixed with the two-sided valley test.

## Caveats and porting notes (briggsv4.4)

- **No omega correspondence.** The frontier is ordered by s along F. `alpha_L_u[j]` no longer belongs to `L[j]`, and the branches can differ in length.
- **`branch_distance`** uses `minimum(abs.(alpha_L_u .- alpha_L_l))`, which pairs by index. It must become a set distance (minimum over all pairs). The line integrals over `alpha_L_u` / `alpha_L_l` (`Phi_F` and derivatives, `influence_ratio`) work unchanged.
- **Envelope, not a single mode.** Frontier vertices are actual eigenvalues, but where the closest mode changes (wall mode → TS at |α_r| ≈ 0.75) the connecting segment belongs to neither mode.
- **Smoothing near a pinch.** Averaging rounds off a dip towards F, so the line moves away from F there. This overestimates the distance to the nearest eigenvalue and could let F cross a real one undetected. Hence the default `CLOUD_SMOOTH = 0`.
- **Bin width vs. omega sampling.** Bins narrower than the spacing of the branch samples leave gaps (repaired by step B); much wider bins blur a sharp dip, since one point per bin survives. `nbins = length(L)` keeps both in proportion under omega refinement.
- **Not yet tested on a deformed F** or near the pinch.
