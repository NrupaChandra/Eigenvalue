# Kilian Vinzenz Wilhelm

# ---------------------------------------------------------------------------
# v5  --  F stops being a straight line, and stops being crossed.
#
# v4.9 fixed the STEP (the barrier can no longer diverge, the displacement is
# limited smoothly, and every step is a line search on max Im(omega_F)).  Over
# 1501 iterations it ran clean: f_ripple 1.25e-5 against v4.7's 3.7e-4 and
# v4.8's 1.3e-3, 3 of 98 nodes at the limiter, d_branch 3.61e-3.  All of that
# is carried forward unchanged.  What it exposed is two faults in the SHAPE of
# F, both measured on contour_iteration_v4.9.json (1501 entries).
#
# FAULT 1 -- F IS CROSSED BY THE UPPER BRANCH, AND IT IS GETTING WORSE.
#
#   iterations    1- 400 : crossing on   0.0 % of frames
#   iterations  400- 800 : crossing on  28.2 %      worst penetration -1.95e-4
#   iterations  800-1100 : crossing on  42.7 %      worst -3.56e-4
#   iterations 1100-1300 : crossing on  59.0 %      median clearance -4.5e-5
#   iterations 1300-1502 : crossing on  83.2 %      median clearance -1.38e-4
#
# First crossing at iteration 507.  The median LATE frame is crossed, so the
# last third of the v4.9 run is not a strictly admissible Briggs state: a pole
# has passed through the inversion contour.  The cause is that
# phi_F = exp(zeta/r^2) - 1 is a function of DISTANCE ONLY -- it has no notion
# of side, so once a branch point is on the wrong side of F the barrier holds
# it there just as firmly.  branch_overlap_valid() does not catch this: it
# compares alpha_u against alpha_l, not the branches against F.
#
# FAULT 2 -- F IS HORIZONTAL IN A SADDLE THAT IS TILTED, AND THAT ALONE COSTS
# A FACTOR 3 IN d.
#
# Write alpha - alpha_p = a + ib.  To second order
#     Im omega - Im omega_p = A a b + (B/2)(a^2 - b^2),   A = Re w'', B = Im w''
# and for a straight F of slope m passing a height b0 above the saddle,
#     max (Im omega - Im omega_p) = k(m) b0^2,  minimised at m = A/B, k = -B/2.
# For this flow k(0) = 0.6636 and k(m_opt) = 0.0695: a factor 9.45 in the peak,
# and since d ~ sqrt(peak), a factor 3.07 in d -- at the SAME b0.  The model is
# not a sketch; using F's own logged b0 and slope it reproduces the measured
# peak to 0-3 % at every checkpoint from iteration 400 to 1501.
#
# F's slope at the pinch over the v4.9 run: -0.347 (k=400), -0.107 (600),
# -0.034 (800), -0.012 (1000), -0.0019 (1501).  It is not merely flat, it is
# being FLATTENED by a factor 180, and the culprit is measurable:
#   rolling_average_filter(., 7) is a 15-point box.  Applied once it has
#   variance (M^2-1)/12 h^2 = 18.7 h^2; applied on 1500 iterations that is a
#   diffusion length of 167 h = 1.7 in alpha_r -- longer than the whole
#   contour.  Over the same 1500 steps sigma = 3e-5 gives sqrt(sigma t) =
#   2.1e-3.  The filter out-smooths the equation's own regulariser by ~800x.
#
# ---------------------------------------------------------------------------
# WHAT v5 CHANGES
#
# CHANGE 13.  alpha_i_rr on a non-uniform grid.  The old stencil
# (y+ - 2y + y-)/(h+ h-) is inconsistent when h+ != h-: it leaves a spurious
# (h+ - h-)/(h+ h-) * y' that does not vanish for a straight line and grows as
# 1/h.  Replaced by 2((y+ - y)/h+ - (y - y-)/h-)/(h+ + h-), and alpha_i_r by
# the second-order unequal-spacing weights.  Both reduce EXACTLY to the v4.9
# formulas on the uniform grid this version still runs, so this is a no-op
# today and a prerequisite for any later refinement of F.
#
# CHANGE 14.  The smoothing gets a width in alpha_r instead of a width in
# nodes.  SMOOTH_HALF_W = min(SMOOTH_FRAC * d_branch, SMOOTH_W_MAX).  Early,
# with d ~ 3, that is 0.15 -- the v4.9 15-node window, so nothing changes
# during the descent.  Late, with d = 3.6e-3, it is 1.8e-4, narrower than one
# cell, so the filter becomes the identity exactly when it was doing the
# damage.  SMOOTH_MODE = :index restores v4.9 behaviour, :off disables it.
#
# CHANGE 15.  sigma * y_rr is treated IMPLICITLY -- one tridiagonal solve on
# 100 nodes, free.  This removes the dt <= h^2/(2 sigma) limit entirely, which
# is what makes the filter dispensable and what will make refinement of F
# survivable.  The operator is the CHANGE 13 non-uniform one, so it is
# consistent on any grid.  DIFFUSION_IMPLICIT = false restores v4.9.
#
# CHANGE 16.  THE DESCENT TERM.  This is the substantive one.
#
# What Briggs' method actually asks for is
#     Im omega_p = min over admissible F  of  max over alpha in F  of Im omega(alpha)
# and the code has always computed the inner max (omegaF_peak_imag) while
# minimising it only by accident -- "repel from the poles" lowers the peak, it
# does not minimise it.  v5 adds the missing half.
#
# Because omega(alpha) is analytic, with omega = u + iv and alpha = x + iy the
# Cauchy-Riemann relations give dv/dy = Re(domega/dalpha) exactly.  F is a
# graph and alpha_i is its only degree of freedom, so the steepest-descent
# direction for the peak, in the one direction F can move, is simply
#
#       dy  =  -lambda * Re(domega/dalpha)
#
# and domega/dalpha costs NOTHING: omega is analytic, so its derivative is the
# same in every direction, and the chord derivative along F,
# (omega_F[j+1] - omega_F[j-1])/(F[j+1] - F[j-1]), is the whole complex
# derivative.  Both arrays are already computed every iteration.  Checked
# against the exact w''(alpha - alpha_p) on the final v4.9 frame: median
# agreement 9.3 % over the 17 nodes around the peak, 2-4 % on the nodes
# nearest it -- far better than a descent direction needs.
#
# Applied with a softmax weight exp((v_j - v_max)/T) so that only the
# neighbourhood of the peak is driven, T taken from the spread over the top
# DESC_NW nodes of the run's own omega_F.  Replayed on the final v4.9 frame it
# pushes the nodes left of the peak UP and those right of it DOWN -- i.e. it
# rotates F clockwise, from -0.0019 toward the -2.92 the geometry wants.
#
# NO PINCH INFORMATION ENTERS.  The slope is FOUND, never imposed: the code
# contains no omega'', no alpha_p, no omega_p, and no constant derived from
# them.  The -2.9242 above is a diagnostic used to size the prize, and it is
# deliberately absent from the source.  The same code transfers unchanged to a
# different flow, where the optimal slope is different and unknown -- which is
# the test that it cannot be encoding this pinch.
#
# CHANGE 17.  A HARD NON-CROSSING CONSTRAINT, which is what FAULT 1 needs.
# Every trial contour is checked: every upper-branch point must lie above F and
# every lower-branch point below it, F being extended flat beyond its ends as
# the boundary condition already implies.  A trial that crosses is rejected and
# the step is halved.  Because the v4.9 state this may resume from is ALREADY
# crossed, the rule accepts a trial that is merely LESS crossed than the
# current state, so a violated configuration can climb back out instead of
# deadlocking.
#
# CHANGE 18.  A PER-NODE trust region.  v4.9 limited every node by the same
# theta = MOVE_FRAC * d_branch.  At the throat that allowed a node with 2.36e-4
# of clearance to move 1.8e-4 -- 76 % of the way onto the branch -- while the
# wings, which have ~0.1 of room and are where the rotation has to happen, were
# held to the same tiny step.  v5 sets theta_j = MOVE_FRAC * (distance from
# node j to the nearest branch point), capped by global_max_move.  This is
# simultaneously far safer at the throat and ~500x faster in the wings, and it
# is what lets CHANGE 16 rotate F in a reasonable number of iterations.
#
# ---------------------------------------------------------------------------
# NOT IN v5, deliberately
#   * F is still 100 uniform nodes.  The v4.9 note argued refinement was
#     needed for a "corner" 3.2e-3 wide; that corner is largely an artefact of
#     F being flat.  A correctly tilted F crosses the branch structure
#     transversally, and at m = m_opt the linear term vanishes so Im omega
#     along F is a pure parabola about the crossing -- which is exactly what
#     omegaF_peak_imag's three-point fit reproduces, independent of spacing.
#     Fix the shape first; refine only if the diagnostics still ask for it.
#   * The eigensolver floor (K = 2.56e-6, d_floor 1.6e-3) is untouched.  The
#     pinch-free backward-error test on P(alpha) = alpha^2 B11 - alpha A11 - A12
#     decides whether that is a pencil-scaling fix or a Newton-on-D fix.
#
# EXPECTED, and how to tell it worked
#   1. n_cross -> 0 and min_clear > 0 on every frame.  This is pass/fail, and
#      it is the one that makes the result admissible.
#   2. f_slope (F's slope at the peak) moves from -0.002 toward -1 .. -3.
#      Do NOT tune anything to make it hit a particular number.
#   3. d_branch falls below the v4.9 plateau of ~3.4e-3.  The saddle model
#      says a factor up to 3 is available from orientation alone; treat
#      anything in that direction as confirmation, not as a target.
#   4. f_ripple stays at or below 1e-5 and n_limited stays small.  If either
#      degrades, CHANGE 14/15 went too far -- set SMOOTH_MODE = :index.
#
# ---------------------------------------------------------------------------
# CORRECTION after the first 100 frames (contour_iteration_v5.0.json)
#
# CHANGE 18 as first written froze the run.  theta_j = MOVE_FRAC * r_j with no
# floor looks reasonable late, but at the START F *is* the real axis and the
# two root families sit on it, so r_j ~ 1e-5 at every node by construction.
# Measured: theta at the peak node was 8.355e-7 against v4.9's 1.670e-1 at the
# same iteration -- 199,839x smaller -- and d moved 3.32 -> 3.29 in 100
# iterations where v4.9 reached 0.997.  CHANGE 18b floors theta_j at the v4.9
# global value, so the per-node rule can only LOOSEN.  Non-crossing safety is
# CHANGE 17's job, which is where it belongs.
#
# CHANGE 17b: branch_clearance's 0.5 gap bound meant no pair qualified while d
# was large, so the crossing test was dormant for the entire descent.  The
# bound is now max(0.5, 20 d_branch).  Verified on the v5.0 frames: 0 points on
# the wrong side at iterations 4, 20, 50, 80, 100, so the test engages without
# rejecting anything valid.
#
# Everything else in the v5.0 run looked right: f_slope moved +0.0008 ->
# -0.0088 (CHANGE 16 rotating F the correct way), f_ripple ~2e-5, and the line
# search accepted on attempt 1 on every frame.
#
# Writes contour_iteration_v5.json.  The 100 frames taken with the frozen theta
# are in contour_iteration_v5.0.json and the 912 deadlocked ones in _v5.1.json;
# neither is this run.
# ---------------------------------------------------------------------------
#
# ---------------------------------------------------------------------------
# SECOND CORRECTION, after contour_iteration_v5.1.json (912 entries)
#
# The v5.1 run descended cleanly to d = 1.756 by iteration 68, slightly AHEAD of
# v4.9 at the same point, with f_slope reaching -0.28 (v4.9 ended at -0.0019),
# so CHANGE 16 works.  Then at iteration 70 it froze completely: identical
# numbers for 843 iterations, max|dF| = 0, omega_i pinned, and the L grid burned
# twelve refinement levels on a geometry that was not moving (191 points,
# h = 2.50e-11) until the 39 MB log failed to write.
#
# The cause was CHANGE 17b, which I added in the previous correction.  Widening
# the clearance test's gap bound to max(0.5, 20 d_branch) pulled ~100 far-field
# branch pairs into the test.  At iteration 70 all 28 reported "crossings" sat at
# alpha_r = 1.22 .. 2.19 -- OUTSIDE F's span of [0,1] -- with branch gaps of 1.8
# to 2.5.  They were being compared against interp_flat's level extension of F,
# which is a fiction: in the construction F is the real axis out there, not a
# horizontal line at F's end height.  The test became unsatisfiable, every trial
# was vetoed, and nothing could ever move again.
#
# The way I checked 17b was also wrong, and that is the more useful lesson: I
# validated it on the v5.0 frames, but that run was frozen by the CHANGE 18 theta
# bug, so F never left the real axis and the fictitious extension happened to lie
# above the far-field points.  A constraint validated on a degenerate run.
#
#   FIX 1.  branch_clearance only considers pairs BOTH of whose points lie
#   inside F's own alpha_r span.  This is what makes the test mean something.
#   FIX 1b.  gbound reverts to the plain 0.5.  Replayed over the v5.1 frames, the
#   v5 rule reports 0 crossings on every frame including the frozen one, so the
#   run would have continued instead of deadlocking.
#
#   FIX 2.  A veto that cannot be overridden is worse than no veto.  After
#   CROSS_ESCAPE_AFTER = 8 consecutive fully-discarded iterations the crossing
#   veto is suspended for one step and says so on stdout.
#
#   FIX 3.  ADAPT_MIN_H 1e-11 -> 5e-7 (below K|omega''|/2 the L grid resolves
#   round-off, not the branch), and a discarded step is no longer fed to the
#   stall detector, which was manufacturing a fake stall every 20 iterations.
#
# Writes contour_iteration_v5.json.
# ---------------------------------------------------------------------------
#
# NOT PARSED BY JULIA.  No Julia was available where this was written.  Run
#     julia -e 'ex = Meta.parseall(read("src/briggsv5.jl", String));
#               bad = filter(e -> e isa Expr && (e.head === :error ||
#                                                e.head === :incomplete), ex.args);
#               isempty(bad) ? println("parse OK") : (foreach(println, bad); exit(1))'
# and then one iteration, before starting a segment.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# v4.9  --  the F step.  Three changes, one block, nothing else touched.
#
# WHY.  Measured on contour_iteration_v4.7.json and _v4.8.json (501 entries
# each), reconstructing the repulsion field from the logged F, alpha_L_u/l and
# zeta_alpha.  The reconstruction reproduces the run's own exp_arg_max on 60 of
# 60 frames in both runs, so the force model below is the one that ran.
#
#   * rhs_cap = 100 fires on 27 of 146 interior F nodes on v4.8 iteration 501.
#     Every node it fires on moves by exactly local_delta_t * 100 = 1.0e-2 per
#     iteration -- 2.6 * d_branch -- and only the SIGN of the true force
#     survives.  Nine sign flips across 41 nodes.  That is the ripple.
#   * The clamped fraction tracks n_F step for step: 0-1 % at 100 nodes, 6 % at
#     112, 14 % at 124, 25 % at 136, 27 % at 148.  Refining F is what put the
#     nodes inside the saturation zone; it did not cause the wiggle directly.
#   * The cause of the saturation is the barrier's dynamic range.  With
#     epsilon = 1e-10 fixed, phi_F = exp(zeta/(r^2+eps)) - 1 reaches an exponent
#     of 1574 at r = 8.75e-5 and |grad Phi| spans 1e-4 .. 3.05e177 across 40
#     neighbouring nodes.  No step-size rule can integrate that.
#   * Consequence.  F's ripple is 1.3e-3 (v4.8) against 3.7e-4 (v4.7).  It lifts
#     max Im(omega_F) by ~0.5*|omega''|*ripple^2 ~ 3.6e-7, and the measured
#     frame-to-frame jitter of that peak is 7.5e-7 while the descent step is
#     9e-8 -- a ratio of 8.  omega_i has been random-walking since ~iteration
#     380, and d_branch is now set by the omega_i offset ALONE (median predicted
#     vs measured over the last 150 frames: 7.86e-3 vs 7.74e-3, log-corr 0.99).
#
# CHANGE 10.  epsilon_alpha = zeta_alpha / EXP_ARG_TARGET, i.e. the softening
# in phi_F is tied to zeta rather than fixed at 1e-10.  The exponent is then
# capped at EXP_ARG_TARGET = 10 BY CONSTRUCTION at every stage of the run, early
# and late, whatever the relation between zeta and d happens to be.  (Tying it
# to d_branch instead does not work: zeta is pinned at ZETA_REF until iteration
# ~246 while d is O(1), so the cap would be ~2e-3 at iteration 50 and the wall
# would vanish.)  The omega-side epsilon is untouched -- phi_L keeps 1e-10.
#
#   early (zeta = 4.0e-4): hard core sqrt(eps) = 6.3e-3 vs wall range
#                          sqrt(zeta) = 2.0e-2  -- the wall still works
#   late  (zeta = 1.2e-5): hard core sqrt(eps) = 1.1e-3 vs d = 3.8e-3
#
# The barrier no longer has to be infinite because branch_overlap_valid() and
# directional_dt_check() are already hard non-crossing guards.  A finite wall
# plus a hard constraint is the sound combination; an infinite wall as the only
# guard is what forced rhs_cap in the first place.
#
# CHANGE 11.  The F step limits the DISPLACEMENT, smoothly, instead of clamping
# the force:
#
#     theta = min(MOVE_FRAC * d_branch, global_max_move)
#     step  = theta * tanh(local_delta_t * rhs_j / theta)
#
# tanh is odd, monotone and C-infinity.  Small forces pass through unchanged
# (tanh x ~ x), large ones level off gradually, and two neighbours can never be
# handed identical magnitudes with opposite signs -- which is what clamp() did.
# theta is tied to d_branch, so the allowed step tightens automatically as the
# gap closes.  At d = 2.3 (iteration 50) theta = 0.115 against the old clamp's
# local_delta_t*100 = 0.1, so this is a near-no-op early and only bites late.
#
# CHANGE 12.  The alpha attempt loop becomes a real line search on the quantity
# the run is actually minimising, max Im(omega_F).  Previously: the loop tested
# whether the MOVE was small, alpha_i_cache and `accepted` were assigned and
# never read anywhere in the file, and alpha_i_smooth from the last attempt was
# applied whether the test passed or failed -- so the loop halved delta_t and
# never rejected anything.  Now each attempt evaluates the trial peak and the
# best trial seen is the one applied, so the peak can only ever go down.
#
# COST is zero on an accepted first attempt: the loop already called
# contour_omega_F(F) once after the update, and attempt 1's evaluation replaces
# it.  Only REJECTED attempts cost extra, bounded by ALPHA_MAX_ATTEMPTS.
#
# The tolerance carries no pinch information: it is DESCENT_TOL_FRAC times the
# run's own recent average descent rate of that peak (half-window medians, the
# same shape as adapt_stalled()).  Early, when the peak is falling fast, it is
# loose enough never to reject; late, when the run has stalled, it goes to zero
# and the rule becomes strictly monotone.
#
# NOT IN THIS VERSION, deliberately -- one variable at a time:
#   * F refinement stays OFF.  N = 100, uniform, as in v4.7.
#   * the non-uniform alpha_i_rr stencil, implicit sigma diffusion, and the
#     15-point index-space box filter are unchanged.  They only matter on a
#     graded F grid, which this version does not build.
#   * ADAPT_MIN_H, pinch_fit and omegaF_peak_imag are unchanged.  The
#     clean-window sqrt fit is LOGGED alongside pinch_fit (omega_pr_clean,
#     omega2_clean, clean_spread, clean_n) so the replacement can be validated
#     against the existing one without changing any behaviour.
#
# WHAT TO WATCH, in this order, over 30 iterations:
#   1. n_limited -> 0 or near it          (27 of 146 clamped in v4.8)
#   2. f_ripple  -> below 3e-4            (1.3e-3 in v4.8, 3.7e-4 in v4.7)
#   3. dpk (max Im(omega_F) step) <= 0 on nearly every iteration.  It is NOT
#      guaranteed: when no attempt lowers the peak the least-bad trial is still
#      applied, so the rule is "take the best available", not "never rise".
#      A run of positive dpk with alpha_accepted = false means the line search
#      has nothing to work with and the geometry, not the step, is the problem.
#   4. d_branch p90/p10 below 1.3         (2.0 in v4.7, 3.9 in v4.8)
#
# SEEDING.  This file writes contour_iteration_v4.9.json and checkpoint_v4.9.json.
# To continue from the v4.7 geometry rather than from scratch, copy
# contour_iteration_v4.7.json to contour_iteration_v4.9.json before the first
# launch and leave RESUME = true: with no checkpoint_v4.9.json present the
# legacy geometry-only path picks up F (100 nodes, so it is compatible),
# omega_i, the omega_r grid, zeta_alpha and the zeta median window from the last
# entry.  Do NOT seed from contour_iteration_v4.8.json -- its F has 148 nodes
# and N is 100 here.  A constants diff against a v4.7 checkpoint is expected and
# correct; v4.9 is a different experiment and gets its own log.
#
# NOT PARSED BY JULIA.  Same caveat as CHANGE 8 and CHANGE 9: no Julia was
# available where this was written.  Run
#     julia -e 'Meta.parse(read("src/briggsv4.9.jl", String))'
# and then one iteration as a smoke test before starting a segment.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# v4.5  --  (A) make one iteration cheap, (B) put the L points where they help.
#
# WHAT v4.4 ESTABLISHED (measured on contour_iteration_v4.4.json, 567 entries)
#   The bisection endgame works: omega_gap sits at the 1e-9 clearance from
#   iteration ~275 on, and omega_i ends at -0.572830165888 against
#   Im(omega_p) = -0.572829694809, i.e. 4.7e-7 past the true pinch height.
#   The residual branch gap is then PURELY the L discretisation:
#
#       d = 2 |sqrt( 2 (L[i] - omega_p) / omega'' )|      (matches to 0.1-0.3 %)
#
#   with Re(omega_p) = 0.2986491 sitting 49.4 % into the cell
#   [0.2985982, 0.2987013] of width 1.031e-4.  Half-cell miss 5.089e-5 gives
#   d = 3.079e-2; observed 3.083e-2.  Same check on v4.3 (h = 2.538e-4):
#   predicted 4.86e-2, observed best 4.73e-2.  Both runs plateaued at exactly
#   the half-cell value, so d_min has been a grid diagnostic, not a physical
#   one, since iteration ~275.  Those last ~290 iterations produced nothing.
#
# CHANGES
#   1. eigvals INSTEAD OF eigen.  No eigenvector is used anywhere in the loop:
#      dominant_eigvals returned them and every caller discarded them
#      (`_, _`), couetteflow_temporal_sing_mode likewise.  Val(:eigvals) skips
#      the Schur-vector back-substitution in ggev.
#   2. ONE EIGENSOLVE PER omega, NOT TWO.  v4.4 swept the upper branch and the
#      lower branch separately, so every omega on L was factorised twice for
#      the same spectrum.  Both branches are now selected from one spectrum.
#   3. PARALLEL SPECTRA, SERIAL SELECTION.  The spectrum at omega does not
#      depend on alpha_prev -- only the final argmin does.  So all N_L
#      eigenproblems are solved by one balanced pmap and the tracking sweep
#      afterwards is pure arithmetic on the main process.  v4.4 ran a
#      460-long SERIAL chain (pmap over a ONE-element array) on 2 of 5
#      workers.  Results are unchanged by construction; verify_fast_path()
#      checks that against the v4.4 code path, which is kept below as
#      contour_alpha_L_conti_ref.
#   4. ADAPTIVE L.  The base 100-point grid is kept and a small window around
#      the current argmin is subdivided by ADAPT_FACTOR whenever the run
#      stalls (same winning omega_r and d_min flat to ADAPT_STALL_TOL over
#      ADAPT_STALL_WINDOW iterations).  

#   5. PINCH_TOL 1e-4 -> 3e-3.  The sqrt-law residual on the v4.4 run is white
#      noise (lag-1 autocorrelation -0.02, uncorrelated between neighbouring
#      grid points) of size delta_d = K/d with K = 2.56e-6, so signal equals
#      noise at d = sqrt(K) = 1.6e-3.  1e-4 is below the double-precision
#      floor of this 300x300 pencil and could never fire.

# ---------------------------------------------------------------------------

begin
    using Distributed, Plots, BenchmarkTools, FFTW, JSON, Statistics, Printf
    addprocs(5)
    w = workers()
end
begin
    @everywhere using LinearAlgebra, Statistics
    ###############
    # EIGENVALUES #
    ###############
    @everywhere begin
        Re = 2000.0
        beta = 0.0 + 0.0 * im
        num_modes = 150
        start = 0
        terminate = 1
        v_g = 0.0 + 0.0 * im
    end
    @everywhere begin
        y_colloc_points = [cos((j - 1) * pi / (num_modes - 1)) for j = 1:num_modes]
        y_colloc_points_new = ((start + terminate) / 2) .- y_colloc_points * ((terminate - start) / 2)
        D0_static = zeros(Float64, num_modes, num_modes)
        for j = 1:num_modes
            D0_static[:, j] .= cos.((j - 1) * acos.(y_colloc_points))
        end
        D1_static = [zeros(num_modes, 1)    D0_static[:, 1]         4 * D0_static[:, 2]]
        D2_static = [zeros(num_modes, 1)    zeros(num_modes, 1)     4 * D0_static[:, 1]]
        D3_static = [zeros(num_modes, 1)    zeros(num_modes, 1)     zeros(num_modes, 1)]
        D4_static = [zeros(num_modes, 1)    zeros(num_modes, 1)     zeros(num_modes, 1)]
        D1_static_V2 = zeros(Float64, num_modes, num_modes)
        D2_static_V2 = zeros(Float64, num_modes, num_modes)
        D3_static_V2 = zeros(Float64, num_modes, num_modes)
        D4_static_V2 = zeros(Float64, num_modes, num_modes)
        D0_static_V2 = D0_static
        D1_static_V2[:, 1:3] .= D1_static
        D2_static_V2[:, 1:3] .= D2_static
        D3_static_V2[:, 1:3] .= D3_static
        D4_static_V2[:, 1:3] .= D4_static
        for j = 4:num_modes
            D1_static_V2[:, j] .= 2 * (j - 1) * D0_static_V2[:, j - 1] + (j - 1) * D1_static_V2[:, j - 2] / (j - 3)
            D2_static_V2[:, j] .= 2 * (j - 1) * D1_static_V2[:, j - 1] + (j - 1) * D2_static_V2[:, j - 2] / (j - 3)
            D3_static_V2[:, j] .= 2 * (j - 1) * D2_static_V2[:, j - 1] + (j - 1) * D3_static_V2[:, j - 2] / (j - 3)
            D4_static_V2[:, j] .= 2 * (j - 1) * D3_static_V2[:, j - 1] + (j - 1) * D4_static_V2[:, j - 2] / (j - 3)
        end
        D1_static_V3 = D1_static_V2 / (-(terminate - start) / 2)^1
        D2_static_V3 = D2_static_V2 / (-(terminate - start) / 2)^2
        D3_static_V3 = D3_static_V2 / (-(terminate - start) / 2)^3
        D4_static_V3 = D4_static_V2 / (-(terminate - start) / 2)^4
        D0 = D0_static
        D1 = D1_static_V3
        D2 = D2_static_V3
        D3 = D3_static_V3
        D4 = D4_static_V3
        u = y_colloc_points_new
        d2u = 0.0
    end
    # -----------------------------------------------------------------------
    # v4.5 CHANGE 1: mode2 === Val(:eigvals) returns the spectrum only.
    # A, B are built exactly as before; only the LAPACK call differs
    # (ggev with jobvl = jobvr = 'N').  Val(:eigen) is untouched so the
    # reference path below still runs.
    # -----------------------------------------------------------------------
    @everywhere function couetteflow(alpha, omega, mode::Val, mode2::Val)
        setprecision(53) do
            if mode === Val(:omega_collocation)
                A11 = - im * alpha * (u * ones(Complex{Float64}, 1, length(u))) .* D2 + im * alpha * (u * ones(Complex{Float64}, 1, length(u))) * (alpha^2 + beta^2) .* D0 + im * alpha * (d2u * ones(1, length(u))) .* D0 + 1 / Re .* D4 - 2 / Re * (alpha^2 + beta^2) .* D2 + 1 / Re * (alpha^2 + beta^2)^2 .* D0 + alpha * v_g .* D0
                A11 = [-200 * im * [D0[1:1, :]; D1[1:1, :]];   A11[3:num_modes - 2, :];    -200 * im * [D1[num_modes:num_modes, :]; D0[num_modes:num_modes, :]]]
                A = A11
                B11 = - im .* D2 + im * (alpha^2 + beta^2) .* D0
                B11 = [[D0[1:1, :]; D1[1:1, :]];   B11[3:num_modes - 2, :];    [D1[num_modes:num_modes, :]; D0[num_modes:num_modes, :]]]
                B = B11
            elseif mode === Val(:alpha_collocation)
                A11 = -2 * im * omega * D1 - 4 / Re * D3 + 4 / Re * beta^2 * D1 - im * (u * ones(Complex{Float64}, 1, length(u))) .* D2 + im * beta^2 * (u * ones(1, length(u))) .* D0 + im * (d2u * ones(1, length(u))) .* D0 - im * v_g .* D2 + im * v_g * beta^2 .* D0
                A12 = im * omega * D2 - im * omega * beta^2 * D0 + 1 / Re * D4 - 2 / Re * beta^2 * D2 + 1 / Re * beta^4 * D0
                A11 .= [zeros(Complex{Float64}, 2, num_modes);                    A11[3:num_modes - 2, :];    zeros(Complex{Float64}, 2, num_modes)]
                A12 .= [-200 * im * [D0[1:1, :]; D1[1:1, :]];   A12[3:num_modes - 2, :];    -200 * im * [D1[num_modes:num_modes, :]; D0[num_modes:num_modes, :]]]
                A21 = 1 * Matrix{Complex{Float64}}(I, num_modes, num_modes)
                A22 = zeros(Complex{Float64}, num_modes, num_modes)
                A = [A11 A12; A21 A22]
                B11 = - 4 / Re * D2 - 2 * im * (u * ones(Complex{Float64}, 1, length(u))) .* D1 + 2 * im * v_g .* D1
                B11 = [zeros(Complex{Float64}, 2, num_modes);  B11[3:num_modes - 2,:];    zeros(Complex{Float64}, 2, num_modes);]
                B12 = zeros(Complex{Float64}, num_modes, num_modes)
                B21 = zeros(Complex{Float64}, num_modes, num_modes)
                B22 = 1 * Matrix{Complex{Float64}}(I, num_modes, num_modes)
                B = [B11 B12; B21 B22]
            end
            if mode2 === Val(:matrix)
                return A, B
            elseif mode2 === Val(:eigen)
                eigvals_, eigvecs_ = eigen(A, B)
                return eigvals_, eigvecs_
            elseif mode2 === Val(:eigvals)
                return eigvals(A, B)
            end
        end
    end
end
############
# CONTOURS #
############
begin
    global L = ComplexF64[]
    global F = ComplexF64[]
    global omega_F = ComplexF64[]
    global alpha_L_u = ComplexF64[]
    global alpha_L_l = ComplexF64[]

    function load_on_workers()
        @sync begin
            for (name, data) in [(:L, L), (:F, F), (:omega_F, omega_F), (:alpha_L_u, alpha_L_u), (:alpha_L_l, alpha_L_l)]
                for pid in workers()
                    @async remotecall_wait(Core.eval, pid, Main, :($name = $(deepcopy(data))))
                end
            end
        end
    end

    # v4.5: omega_r is now mutated at run time by the adaptive refinement, so
    # the workers' stale copy from the @everywhere block below has to be
    # replaced whenever it changes.  Nothing on a worker currently reads
    # omega_r (L is always passed as an argument), but keeping the copies
    # consistent costs nothing and removes a trap.
    function push_omega_r()
        @sync for pid in workers()
            @async remotecall_wait(Core.eval, pid, Main, :(omega_r = $(deepcopy(omega_r))))
        end
    end

    # ALPHA contour: F
    @everywhere begin
        alpha_r_start = 0.0
        alpha_r_end = 1.0
        N = 100
        alpha_r = range(alpha_r_start, alpha_r_end, length=N)
        alpha_i = fill(0.0, N)
    end
    function contour_F()
        F = [alpha_r[j] + alpha_i[j] * im for j in 1:N]
        return F
    end
    @everywhere F = Vector{Complex{Float64}}[]
    F = contour_F()

    # -----------------------------------------------------------------------
    # v4.5 CHANGE 1: no eigenvector.  The caller only ever took the first
    # return value.
    # -----------------------------------------------------------------------
    @everywhere function couetteflow_temporal_sing_mode(alpha)
        ev = couetteflow(alpha, nothing, Val(:omega_collocation), Val(:eigvals))
        ev = ev[[isfinite(real(e)) && isfinite(imag(e)) for e in ev]]
        return ev[argmax(imag.(ev))]
    end

    # -----------------------------------------------------------------------
    # v4.4 REFERENCE PATH.  Unused by the loop, kept so verify_fast_path()
    # can check the new selection against it.  Do not delete.
    # -----------------------------------------------------------------------
    @everywhere function couetteflow_spatial_sing_mode_comparison(
        omega, alpha_approximation, branch_side::Symbol, F, normals_F; side_tol = 0.0)
        eigvals_, _ = couetteflow(nothing, omega, Val(:alpha_collocation), Val(:eigen))
        mask = [isfinite(real(e)) && isfinite(imag(e)) for e in eigvals_]
        eigvals_ = eigvals_[mask]
        candidates = ComplexF64[]
        for eigval in eigvals_
            distances = [abs(eigval - f) for f in F]
            idx_min = argmin(distances)
            f_near = F[idx_min]
            normal = normals_F[idx_min]
            proj = real(conj(normal) * (eigval - f_near))
            if branch_side == :upper && proj > side_tol
                push!(candidates, eigval)
            elseif branch_side == :lower && proj < -side_tol
                push!(candidates, eigval)
            end
        end
        if isempty(candidates)
            diffs = abs.(ComplexF64(alpha_approximation) .- ComplexF64.(eigvals_))
            return eigvals_[argmin(diffs)]
        end
        diffs = abs.(ComplexF64(alpha_approximation) .- candidates)
        return candidates[argmin(diffs)]
    end

    function contour_omega_F(F)
        @everywhere omega_F = Complex{Float64}[]
        omega_F = pmap(alpha -> couetteflow_temporal_sing_mode(alpha), F)
        return ComplexF64.(omega_F)
    end
    omega_F = contour_omega_F(F)

    # OMEGA contour: L
    @everywhere begin
        omega_r_start = 0.0
        omega_r_end = 0.5
        omega_i = 0.0

        # v4.5: this is the ONLY grid built up front.  v4.4's REFINE_LO /
        # REFINE_HI / PTS_PER_CELL are gone -- the refinement below is placed
        # by the run, not guessed before it.
        omega_r_base = collect(range(omega_r_start, omega_r_end, length=100))
        omega_r = copy(omega_r_base)
        N_L = length(omega_r)
    end
    function contour_L()
        L = Complex{Float64}[
            omega_r[j] + omega_i * im
            for j in eachindex(omega_r)
        ]
        return L
    end
    @everywhere L = Vector{Complex{Float64}}[]
    L = contour_L()
    #
    @everywhere begin
        alpha_L_u = Vector{Complex{Float64}}[]
        alpha_L_l = Vector{Complex{Float64}}[]
    end
    #
    load_on_workers()
    @everywhere function contour_normals(F)
        normals = Complex{Float64}[]
        for j in 2:(length(F) - 1)
            tangent = F[j+1] - F[j-1]
            normal = im * tangent / abs(tangent)
            push!(normals, normal)
        end
        insert!(normals, 1, normals[1])
        push!(normals, normals[end])
        return normals
    end
    function plot_normals()
        x = real.(F)
        y = imag.(F)
        u_vec = real(contour_normals(F))
        v_vec = imag(contour_normals(F))
        quiver(x, y, quiver=(u_vec, v_vec), aspect_ratio=1; xlims=(-1.0, 1.0), ylims=(-1.0, 1.0))
    end
    @everywhere begin
        normals_F = contour_normals(F)
    end

    # -----------------------------------------------------------------------
    # v4.5 CHANGES 2 + 3.
    #
    # spatial_payload does the expensive part for ONE omega: one eigensolve,
    # then the side classification, which depends only on (omega, F).  It
    # returns both sides at once, so the upper and lower branches share a
    # single factorisation instead of forcing two.
    #
    # Everything that depends on alpha_prev -- i.e. the tracking -- is a
    # nearest-neighbour argmin over at most a few hundred numbers, so it is
    # done afterwards on the main process.  This is what breaks v4.4's serial
    # chain without changing any result: the spectrum at omega_j is the same
    # number whether it was computed first, last, or in parallel.
    #
    # Payload fields
    #   all   every finite eigenvalue at this omega        (fallback pool)
    #   up    those on the +normal side of F,  up_d their distance to F
    #   lo    those on the -normal side of F,  lo_d likewise
    # -----------------------------------------------------------------------
    @everywhere function spatial_payload(omega, F, normals_F)
        ev = couetteflow(nothing, omega, Val(:alpha_collocation), Val(:eigvals))
        ev = ev[[isfinite(real(e)) && isfinite(imag(e)) for e in ev]]
        up = ComplexF64[];  up_d = Float64[]
        lo = ComplexF64[];  lo_d = Float64[]
        for e in ev
            dmin = Inf
            k = 0
            @inbounds for t in eachindex(F)
                dd = abs(e - F[t])
                if dd < dmin
                    dmin = dd
                    k = t
                end
            end
            proj = real(conj(normals_F[k]) * (e - F[k]))
            if proj > 0.0
                push!(up, e);  push!(up_d, dmin)
            elseif proj < 0.0
                push!(lo, e);  push!(lo_d, dmin)
            end
        end
        return (all = ev, up = up, up_d = up_d, lo = lo, lo_d = lo_d)
    end

    # Same rule as v4.4 dominant_eigvals: on each side, the eigenvalue whose
    # nearest F point is nearest.  Used for the initial contour and for the
    # single seed point of the tracking sweep.
    function dominant_from(pl)
        eu = isempty(pl.up) ? nothing : pl.up[argmin(pl.up_d)]
        el = isempty(pl.lo) ? nothing : pl.lo[argmin(pl.lo_d)]
        return eu, el
    end

    # -----------------------------------------------------------------------
    # v4.7 CHANGE 9: continuity wins when the side test contradicts it.
    #
    # THE BUG.  spatial_payload decides "upper" or "lower" from the sign of
    # the projection of (eigenvalue - nearest F node) onto that node's normal.
    # The resolution of that test is set by F's node spacing.  F is 100 points
    # over alpha_r in [0,1], so the spacing is ~1.0e-2, while the two branch
    # points it is being asked to separate are d_branch apart -- 3.2e-3 on the
    # last frame of the 501-iteration run, and shrinking every time the run
    # succeeds.  Measured there, at omega index 82 of 141:
    #
    #     both near-pinch roots land on the SAME F node (k = 57)
    #     upper root: |e - F[k]| = 2.07e-3, of which the normal part is
    #                 -1.64e-4 and the tangential part +2.06e-3
    #     so the sign that decides the side is set by 1.6e-4 against a node
    #     spacing of 1.0e-2 -- a margin 62x below the resolution of the
    #     polyline whose normal is being used.
    #
    # Both roots then came out classified `lo`, pl.up held nothing within 0.02
    # of the pinch, and this function returned the nearest OTHER upper-classed
    # eigenvalue: 0.685013 - 1.142555i, a jump of 1.899 from alpha_prev.  The
    # sweep is a chain, so alpha_prev was then wrong for every remaining omega
    # and the branch stayed latched on that mode to the end of L.  187 of the
    # 501 logged frames have an upper-branch step > 0.5; it first appears at
    # iteration 32 and is present on nearly every frame from 101 on.  Whether
    # it fires depends on where F's nodes happen to fall, which is why it looks
    # intermittent.
    #
    # THE FIX.  Keep the side test -- it is what separates the branches when
    # they are far apart -- but do not let it move the tracked root by an
    # amount that pure continuity says is absurd.  If honouring the side
    # classification costs TRACK_SIDE_SLACK times more movement than the
    # nearest eigenvalue overall, the classification is what is wrong, not the
    # continuation, so take the continuation.
    #
    # This cannot silently merge the two branches: a swap would need the other
    # branch to be the nearest root, and near the tip the branches are
    # d_branch = 3.2e-3 apart while consecutive omega samples move a root by
    # ~1.7e-4, a ratio of 19.  branch_overlap_valid() is still the backstop.
    #
    # Replayed on the recorded iteration-501 geometry, recomputing all 141
    # spectra: the guard fires 7 times out of 282 selections, the largest
    # upper-branch step falls from 1.8983 to 0.0857, no step exceeds 0.5, the
    # branch becomes monotone in Re past the pinch, and d_branch is unchanged
    # to within the resolution of the check.  Only indices 80-140 -- the tail
    # after the pinch -- differ; nothing before the pinch moves.
    #
    # TRACK_SIDE_SLACK = Inf restores the v4.6/v4.7 behaviour exactly.
    # -----------------------------------------------------------------------
    const TRACK_SIDE_SLACK = 10.0   # Inf -> exactly the old behaviour
    const TRACK_MIN_STEP   = 1e-9   # floor, so a root that barely moved cannot
                                    # make the ratio blow up on its own

    global track_overrides = 0      # reset by contour_alpha_L_conti
    global track_jump_u    = 0.0    # largest |d alpha| along the upper branch
    global track_jump_l    = 0.0    # ... and the lower, set by the same sweep

    # Same rule as v4.4 couetteflow_spatial_sing_mode_comparison with
    # side_tol = 0: nearest to alpha_prev among the requested side, falling
    # back to nearest among all finite eigenvalues if that side is empty --
    # plus the CHANGE 9 continuity guard above.
    function select_tracked(pl, alpha_prev, side::Symbol)
        isempty(pl.all) && error("spatial_payload: no finite eigenvalues")
        ap = ComplexF64(alpha_prev)
        all_best = pl.all[argmin(abs.(ap .- pl.all))]
        cand = side === :upper ? pl.up : pl.lo
        isempty(cand) && return all_best
        side_best = cand[argmin(abs.(ap .- cand))]
        if isfinite(TRACK_SIDE_SLACK)
            d_side = abs(ap - side_best)
            d_all  = abs(ap - all_best)
            if d_side > TRACK_SIDE_SLACK * max(d_all, TRACK_MIN_STEP)
                global track_overrides += 1
                return all_best
            end
        end
        return side_best
    end

    # Kept for callers that want one omega on the main process.
    function dominant_eigvals(omega, F, normals_F)
        return dominant_from(spatial_payload(omega, F, normals_F))
    end

    function contour_alpha_L_init(L)
        payloads = pmap(w -> spatial_payload(w, F, normals_F), L)
        au = Vector{ComplexF64}(undef, length(L))
        al = Vector{ComplexF64}(undef, length(L))
        for j in eachindex(L)
            eu, el = dominant_from(payloads[j])
            (eu === nothing || el === nothing) &&
                error("contour_alpha_L_init: empty branch side at omega = $(L[j])")
            au[j] = eu
            al[j] = el
        end
        @everywhere alpha_L_u = Complex{Float64}[]
        @everywhere alpha_L_l = Complex{Float64}[]
        return au, al
    end

    # -----------------------------------------------------------------------
    # v4.5 contour_alpha_L_conti.
    # One balanced pmap over the whole of L (this is >99 % of the cost), then
    # the same seed-and-walk selection v4.4 performed, done locally.
    # -----------------------------------------------------------------------
    function contour_alpha_L_conti(L)
        payloads = pmap(w -> spatial_payload(w, F, normals_F), L)

        # Same physical starting omega as the old 100-point grid
        omega_start_target = omega_r_base[25]
        s = argmin(abs.(real.(L) .- omega_start_target))

        au = Vector{ComplexF64}(undef, length(L))
        al = Vector{ComplexF64}(undef, length(L))

        eu, el = dominant_from(payloads[s])
        (eu === nothing || el === nothing) &&
            error("contour_alpha_L_conti: empty branch side at the seed omega = $(L[s])")
        au[s] = eu
        al[s] = el

        # v4.7 CHANGE 9 diagnostics.  Counted per sweep; this function is also
        # called for rejected omega trials, so what ends up in the log is the
        # last sweep before the log write, i.e. the accepted geometry.
        global track_overrides = 0

        for j in (s + 1):length(L)
            au[j] = select_tracked(payloads[j], au[j - 1], :upper)
            al[j] = select_tracked(payloads[j], al[j - 1], :lower)
        end
        for j in (s - 1):-1:1
            au[j] = select_tracked(payloads[j], au[j + 1], :upper)
            al[j] = select_tracked(payloads[j], al[j + 1], :lower)
        end

        # Largest step along each branch between neighbouring omega samples.
        # This is the direct measure of the artefact CHANGE 9 removes: on the
        # broken frames it was ~1.9, on healthy ones it is ~0.09.
        global track_jump_u = length(au) < 2 ? 0.0 :
                              maximum(abs.(au[2:end] .- au[1:end-1]))
        global track_jump_l = length(al) < 2 ? 0.0 :
                              maximum(abs.(al[2:end] .- al[1:end-1]))
        return au, al
    end

    #######################
    # v4.4 REFERENCE PATH #
    #######################
    # Kept verbatim so verify_fast_path() has something to compare against.
    @everywhere function track_eigenvalue_simple(alpha_0)
        return alpha_0
    end
    @everywhere function track_branch_pmap(L, start_index, alpha_0, direction, branch_side)
        N = length(L)
        alpha_current = alpha_0
        current_index = start_index
        results = [(start_index, alpha_current)]
        while 1 <= current_index + direction <= N && current_index + direction >= 1
            indices_segment = current_index + direction:direction:current_index + direction
            L_segment = L[indices_segment]
            args = [(idx, omega, alpha_current) for (idx, omega) in zip(indices_segment, L_segment)]
            tracked = pmap(arg -> begin
                index, omega, alpha_prev = arg
                alpha_tracked = track_eigenvalue_simple(alpha_prev)
                alpha_corrected = couetteflow_spatial_sing_mode_comparison(
                    omega, alpha_tracked, branch_side, F, normals_F)
                (index, alpha_corrected)
            end, args)
            result = tracked[1]
            push!(results, result)
            alpha_current = result[2]
            current_index = result[1]
        end
        return results
    end
    function contour_alpha_L_conti_ref(L)
        omega_start_target = omega_r_base[25]
        start_index = argmin(abs.(real.(L) .- omega_start_target))
        alpha_u_start, alpha_l_start = dominant_eigvals(L[start_index], F, normals_F)
        future_u_fwd = @spawn track_branch_pmap(L, start_index, alpha_u_start, +1, :upper)
        future_u_bwd = @spawn track_branch_pmap(L, start_index, alpha_u_start, -1, :upper)
        future_l_fwd = @spawn track_branch_pmap(L, start_index, alpha_l_start, +1, :lower)
        future_l_bwd = @spawn track_branch_pmap(L, start_index, alpha_l_start, -1, :lower)
        results_u = vcat(fetch(future_u_bwd), fetch(future_u_fwd)[2:end])
        results_l = vcat(fetch(future_l_bwd), fetch(future_l_fwd)[2:end])
        sort!(results_u, by = x -> x[1])
        sort!(results_l, by = x -> x[1])
        return [x[2] for x in results_u], [x[2] for x in results_l]
    end

    alpha_L_u, alpha_L_l = contour_alpha_L_init(L)
    load_on_workers()

    #########
    # PLOTS #
    #########
    function plot_omega()
        plot(L)
        plot!(omega_F)
    end
    function plot_alpha()
        plot(F)
        plot!(alpha_L_u)
        plot!(alpha_L_l, color=3)
    end
    function reset()
        global alpha_i = fill(0.0, N)
        global F = contour_F()
        global omega_F = contour_omega_F(F)
        global omega_i = 0.0
        global omega_r = copy(omega_r_base)
        global adapt_level  = 0
        global adapt_nogain = 0
        global adapt_done   = false
        adapt_reset_history!()
        push_omega_r()
        global L = contour_L()
        load_on_workers()
        @everywhere begin
            global normals_F = contour_normals(F)
        end
        global alpha_L_u, alpha_L_l = contour_alpha_L_init(L)
        load_on_workers()
        iteration_step = 1
    end
end

# ---------------------------------------------------------------------------
# EQUIVALENCE CHECK.  Run once after loading, before starting a long run:
#
#     verify_fast_path()
#
# It compares the v4.5 selection against the v4.4 code path on the current L
# and prints the largest absolute difference on each branch.  Expect exactly
# 0.0 -- the eigenvalues come from the same LAPACK routine on the same
# matrices and only the selection was reorganised.  A difference at the 1e-16
# level would still be acceptable (ggev with and without vectors may order a
# degenerate pair differently); anything larger means the refactor is wrong.
# ---------------------------------------------------------------------------
function verify_fast_path(L_test = L)
    @printf("verify_fast_path: %d points, this runs BOTH paths and is slow.\n", length(L_test))
    flush(stdout)
    t_fast = @elapsed (au_f, al_f) = contour_alpha_L_conti(L_test)
    t_ref  = @elapsed (au_r, al_r) = contour_alpha_L_conti_ref(L_test)
    du = maximum(abs.(au_f .- au_r))
    dl = maximum(abs.(al_f .- al_r))
    @printf("  max |upper_fast - upper_ref| = %.3e\n", du)
    @printf("  max |lower_fast - lower_ref| = %.3e\n", dl)
    @printf("  fast %.2f s | ref %.2f s | speedup %.1fx\n", t_fast, t_ref, t_ref / t_fast)
    flush(stdout)
    return du, dl
end

######################
# POTENTIAL FUNCTION #
######################
begin
    s_omega = 2.0
    s_alpha = 2.0
    epsilon = 1e-10

    global zeta_common = 4e-4
    global zeta_omega = zeta_common
    global zeta_alpha = zeta_common

    # -----------------------------------------------------------------------
    # v4.3a OVERFLOW GUARD  (fixes the LAPACKException(150) crash at k = 122)
    # exp(zeta/(d^s + epsilon)) overflows once d = |alpha_F - alpha_L| < 7.5e-4.
    # Both d_d_alpha_r_Phi_F and d_d_alpha_i_Phi_F carry the same exp factor and
    # enter rhs_j with opposite signs, so the numerator becomes (+Inf) + (-Inf)
    # = NaN, and clamp(NaN, -100, 100) === NaN passes it straight through into
    # alpha_i -> F -> eigen(A, B) -> LAPACK INFO = N = 150.
    #
    # v4.5 note: with zeta fixed at 4e-4 the repulsion length is sqrt(zeta) =
    # 0.02, and F has to thread the gap between the two branches, so it must
    # sit within ~d/2 of each.  On the v4.4 run F sat 1.3e-2 from the upper
    # branch at d = 3.08e-2, giving zeta/r^2 ~ 2 -- harmless.  At d = 1e-3 the
    # same geometry gives zeta/r^2 = 1600, i.e. past EXP_ARG_MAX with rhs_j
    # pinned at rhs_cap over the whole neighbourhood.  The scaling that
    # reproduces the current hand-tuned value is zeta ~ d^2/2 (d = 3.08e-2
    # gives 4.7e-4 against the 4e-4 set here).  NOT implemented in v4.5: it is
    # the next change, and it should not be mixed into the same run as the
    # grid change.  Watch d_contour and the exp arguments once d_branch drops
    # below ~1e-2.
    # -----------------------------------------------------------------------
    const EXP_ARG_MAX = 400.0
    expc(x) = exp(x < EXP_ARG_MAX ? x : EXP_ARG_MAX)

    # -----------------------------------------------------------------------
    # v4.9 CHANGE 10.  The alpha-side softening is tied to zeta, so the
    # exponent phi_F can ever see is bounded by EXP_ARG_TARGET at every stage
    # of the run.  epsilon (1e-10) stays as it is for the omega side: phi_L
    # runs against zeta_omega, which is frozen, and nothing there saturates.
    # epsilon_alpha is refreshed in the main loop wherever zeta_alpha is.
    # -----------------------------------------------------------------------
    const EXP_ARG_TARGET = 10.0
    global epsilon_alpha = zeta_alpha / EXP_ARG_TARGET

    function phi_L(omega_L, omega)
        phi_L = 0.0
        phi_L = exp(zeta_omega / (abs(omega_L - omega)^s_omega + epsilon)) - 1.0
        return phi_L
    end
    function Phi_L(omega_L)
        Phi_L = 0.0
        d_omega_1 = omega_F[2] - omega_F[1]
        Phi_L += phi_L(omega_L, omega_F[1]) * abs(d_omega_1)
        for j in 2:(length(omega_F) - 1)
            d_omega_j = 0.5 * (omega_F[j+1] - omega_F[j-1])
            Phi_L += phi_L(omega_L, omega_F[j]) * abs(d_omega_j)
        end
        d_omega_N = omega_F[N] - omega_F[N-1]
        Phi_L += phi_L(omega_L, omega_F[N]) * abs(d_omega_N)
        return Phi_L
    end
    function phi_F(alpha_F, alpha)
        phi_F = 0.0
        phi_F = expc(zeta_alpha / (abs(alpha_F - alpha)^s_alpha + epsilon_alpha)) - 1.0
        return phi_F
    end
    # NOTE on the adaptive grid: Phi_F and its gradients weight every alpha_L
    # point by the local arc length |0.5 (alpha[j+1] - alpha[j-1])|, so
    # clustering L points near the pinch does not change the potential -- it
    # only resolves the same charge distribution better.  Near the cusp
    # |dalpha/domega| ~ |omega - omega_p|^(-1/2) diverges, so the arc length
    # per omega cell GROWS towards the pinch and v4.4's uniform spacing was
    # under-resolving exactly the part that carries the most weight.
    function Phi_F(alpha_F)
        Phi_F = 0.0
        d_alpha_u_1 = alpha_L_u[2] - alpha_L_u[1]
        Phi_F += phi_F(alpha_F, alpha_L_u[1]) * abs(d_alpha_u_1)
        for j in 2:(length(alpha_L_u) - 1)
            d_alpha_u_j = 0.5 * (alpha_L_u[j+1] - alpha_L_u[j-1])
            Phi_F += phi_F(alpha_F, alpha_L_u[j]) * abs(d_alpha_u_j)
        end
        d_alpha_u_N = alpha_L_u[end] - alpha_L_u[end-1]
        Phi_F += phi_F(alpha_F, alpha_L_u[end]) * abs(d_alpha_u_N)
        ######
        d_alpha_l_1 = alpha_L_l[2] - alpha_L_l[1]
        Phi_F += phi_F(alpha_F, alpha_L_l[1]) * abs(d_alpha_l_1)
        for j in 2:(length(alpha_L_l) - 1)
            d_alpha_l_j = 0.5 * (alpha_L_l[j+1] - alpha_L_l[j-1])
            Phi_F += phi_F(alpha_F, alpha_L_l[j]) * abs(d_alpha_l_j)
        end
        d_alpha_l_N = alpha_L_l[end] - alpha_L_l[end-1]
        Phi_F += phi_F(alpha_F, alpha_L_l[end]) * abs(d_alpha_l_N)
        return Phi_F
    end
    ###############################
    # POTENTIAL FUNCTION GRADIENT #
    ###############################
    function d_d_omega_r_phi_L(omega_L, omega)
        d_d_omega_r_phi_L = 0.0
        d_d_omega_r_phi_L = -zeta_omega * (real(omega_L) - real(omega)) * s_omega * abs(omega_L - omega)^(s_omega - 2) / (abs(omega_L - omega)^s_omega + epsilon)^2.0 * exp(zeta_omega / (abs(omega_L - omega)^s_omega + epsilon))
        return d_d_omega_r_phi_L
    end
    function d_d_omega_r_Phi_L(omega_L)
        d_d_omega_r_Phi_L = 0.0
        d_omega_1 = omega_F[2] - omega_F[1]
        d_d_omega_r_Phi_L += d_d_omega_r_phi_L(omega_L, omega_F[1]) * abs(d_omega_1)
        for j in 2:(length(omega_F) - 1)
            d_omega_j = 0.5 * (omega_F[j+1] - omega_F[j-1])
            d_d_omega_r_Phi_L += d_d_omega_r_phi_L(omega_L, omega_F[j]) * abs(d_omega_j)
        end
        d_omega_N = omega_F[N] - omega_F[N-1]
        d_d_omega_r_Phi_L += d_d_omega_r_phi_L(omega_L, omega_F[N]) * abs(d_omega_N)
        return d_d_omega_r_Phi_L
    end
    #
    function d_d_omega_i_phi_L(omega_L, omega)
        d_d_omega_i_phi_L = 0.0
        d_d_omega_i_phi_L = -zeta_omega * (imag(omega_L) - imag(omega)) * s_omega * abs(omega_L - omega)^(s_omega - 2) / (abs(omega_L - omega)^s_omega + epsilon)^2.0 * exp(zeta_omega / (abs(omega_L - omega)^s_omega + epsilon))
        return d_d_omega_i_phi_L
    end
    function d_d_omega_i_Phi_L(omega_L)
        d_d_omega_i_Phi_L = 0.0
        d_omega_1 = omega_F[2] - omega_F[1]
        d_d_omega_i_Phi_L += d_d_omega_i_phi_L(omega_L, omega_F[1]) * abs(d_omega_1)
        for j in 2:(length(omega_F) - 1)
            d_omega_j = 0.5 * (omega_F[j+1] - omega_F[j-1])
            d_d_omega_i_Phi_L += d_d_omega_i_phi_L(omega_L, omega_F[j]) * abs(d_omega_j)
        end
        d_omega_N = omega_F[N] - omega_F[N-1]
        d_d_omega_i_Phi_L += d_d_omega_i_phi_L(omega_L, omega_F[N]) * abs(d_omega_N)
        return d_d_omega_i_Phi_L
    end
    ###
    function d_d_alpha_r_phi_F(alpha_F, alpha)
        d_d_alpha_r_phi_F = 0.0
        d_d_alpha_r_phi_F = -zeta_alpha * (real(alpha_F) - real(alpha)) * s_alpha * abs(alpha_F - alpha)^(s_alpha - 2) / (abs(alpha_F - alpha)^s_alpha + epsilon_alpha)^2.0 * expc(zeta_alpha / (abs(alpha_F - alpha)^s_alpha + epsilon_alpha))
        return d_d_alpha_r_phi_F
    end
    function d_d_alpha_r_Phi_F(alpha_F)
        d_d_alpha_r_Phi_F = 0.0
        d_alpha_u_1 = alpha_L_u[2] - alpha_L_u[1]
        d_d_alpha_r_Phi_F += d_d_alpha_r_phi_F(alpha_F, alpha_L_u[1]) * abs(d_alpha_u_1)
        for j in 2:(length(alpha_L_u) - 1)
            d_alpha_u_j = 0.5 * (alpha_L_u[j+1] - alpha_L_u[j-1])
            d_d_alpha_r_Phi_F += d_d_alpha_r_phi_F(alpha_F, alpha_L_u[j]) * abs(d_alpha_u_j)
        end
        d_alpha_u_N = alpha_L_u[end] - alpha_L_u[end-1]
        d_d_alpha_r_Phi_F += d_d_alpha_r_phi_F(alpha_F, alpha_L_u[end]) * abs(d_alpha_u_N)
        ######
        d_alpha_l_1 = alpha_L_l[2] - alpha_L_l[1]
        d_d_alpha_r_Phi_F += d_d_alpha_r_phi_F(alpha_F, alpha_L_l[1]) * abs(d_alpha_l_1)
        for j in 2:(length(alpha_L_l) - 1)
            d_alpha_l_j = 0.5 * (alpha_L_l[j+1] - alpha_L_l[j-1])
            d_d_alpha_r_Phi_F += d_d_alpha_r_phi_F(alpha_F, alpha_L_l[j]) * abs(d_alpha_l_j)
        end
        d_alpha_l_N = alpha_L_l[end] - alpha_L_l[end-1]
        d_d_alpha_r_Phi_F += d_d_alpha_r_phi_F(alpha_F, alpha_L_l[end]) * abs(d_alpha_l_N)
        return d_d_alpha_r_Phi_F
    end
    #
    function d_d_alpha_i_phi_F(alpha_F, alpha)
        d_d_alpha_i_phi_F = 0.0
        d_d_alpha_i_phi_F = -zeta_alpha * (imag(alpha_F) - imag(alpha)) * s_alpha * abs(alpha_F - alpha)^(s_alpha - 2) / (abs(alpha_F - alpha)^s_alpha + epsilon_alpha)^2.0 * expc(zeta_alpha / (abs(alpha_F - alpha)^s_alpha + epsilon_alpha))
        return d_d_alpha_i_phi_F
    end
    function d_d_alpha_i_Phi_F(alpha_F)
        d_d_alpha_i_Phi_F = 0.0
        d_alpha_u_1 = alpha_L_u[2] - alpha_L_u[1]
        d_d_alpha_i_Phi_F += d_d_alpha_i_phi_F(alpha_F, alpha_L_u[1]) * abs(d_alpha_u_1)
        for j in 2:(length(alpha_L_u) - 1)
            d_alpha_u_j = 0.5 * (alpha_L_u[j+1] - alpha_L_u[j-1])
            d_d_alpha_i_Phi_F += d_d_alpha_i_phi_F(alpha_F, alpha_L_u[j]) * abs(d_alpha_u_j)
        end
        d_alpha_u_N = alpha_L_u[end] - alpha_L_u[end-1]
        d_d_alpha_i_Phi_F += d_d_alpha_i_phi_F(alpha_F, alpha_L_u[end]) * abs(d_alpha_u_N)
        ######
        d_alpha_l_1 = alpha_L_l[2] - alpha_L_l[1]
        d_d_alpha_i_Phi_F += d_d_alpha_i_phi_F(alpha_F, alpha_L_l[1]) * abs(d_alpha_l_1)
        for j in 2:(length(alpha_L_l) - 1)
            d_alpha_l_j = 0.5 * (alpha_L_l[j+1] - alpha_L_l[j-1])
            d_d_alpha_i_Phi_F += d_d_alpha_i_phi_F(alpha_F, alpha_L_l[j]) * abs(d_alpha_l_j)
        end
        d_alpha_l_N = alpha_L_l[end] - alpha_L_l[end-1]
        d_d_alpha_i_Phi_F += d_d_alpha_i_phi_F(alpha_F, alpha_L_l[end]) * abs(d_alpha_l_N)
        return d_d_alpha_i_Phi_F
    end
    function plot_omega_potential()
        x = real.(L)
        y = imag.(L)
        u = [d_d_omega_r_Phi_L(omega) for omega in L]
        v = [d_d_omega_i_Phi_L(omega) for omega in L]
        quiver(x, y, quiver=(u, v); xlims=(omega_r_start, omega_r_end))
        plot!(omega_F)
    end
    function plot_alpha_potential()
        x = real.(F)
        y = imag.(F)
        u = [d_d_alpha_r_Phi_F(alpha) for alpha in F]
        v = [d_d_alpha_i_Phi_F(alpha) for alpha in F]
        quiver(x, y, quiver=(u, v); xlims=(alpha_r_start, alpha_r_end))
        plot!(alpha_L_u)
        plot!(alpha_L_l)
    end
    function influence_ratio(alpha_F, alpha_L_u, alpha_L_l)
        infl_u = 0.0
        infl_l = 0.0
        d_alpha_u_1 = alpha_L_u[2] - alpha_L_u[1]
        infl_u += phi_F(alpha_F, alpha_L_u[1]) * abs(d_alpha_u_1)
        for j in 2:(length(alpha_L_u)-1)
            d_alpha_u_j = 0.5 * (alpha_L_u[j+1] - alpha_L_u[j-1])
            infl_u += phi_F(alpha_F, alpha_L_u[j]) * abs(d_alpha_u_j)
        end
        infl_u += phi_F(alpha_F, alpha_L_u[end]) * abs(alpha_L_u[end] - alpha_L_u[end-1])
        d_alpha_l_1 = alpha_L_l[2] - alpha_L_l[1]
        infl_l += phi_F(alpha_F, alpha_L_l[1]) * abs(d_alpha_l_1)
        for j in 2:(length(alpha_L_l)-1)
            d_alpha_l_j = 0.5 * (alpha_L_l[j+1] - alpha_L_l[j-1])
            infl_l += phi_F(alpha_F, alpha_L_l[j]) * abs(d_alpha_l_j)
        end
        infl_l += phi_F(alpha_F, alpha_L_l[end]) * abs(alpha_L_l[end] - alpha_L_l[end-1])
        ratio = infl_l / (infl_u + infl_l + eps())
        return ratio
    end
    function acceptance_factor(alpha_F, alpha_L_u, alpha_L_l)
        r = influence_ratio(alpha_F, alpha_L_u, alpha_L_l)
        closeness = 4 * r * (1 - r)
        factor = 10.0 - 9.99 * closeness
        return factor
    end
end
#################
# TIME-STEPPING #
#################
begin
    lambda = 4.0       # much stronger downward push
    sigma = 3e-5       # weaker smoothing/diffusion on F
    delta_t = 1e-3     # aggressive initial pseudo-time step
    function spectral_filter(alpha_i, cutoff_fraction)
        x = alpha_i
        X = fft(x)
        N = length(X)
        cutoff = floor(Int, cutoff_fraction * N/2)
        X[cutoff+1:end-cutoff] .= 0
        x_filtered = real(ifft(X))
        return x_filtered
    end
    function rolling_average_filter(x, window_radius)
        N = length(x)
        x_smooth = similar(x)
        for i in 1:N
            # Determine window bounds (handle edges safely)
            left = max(1, i - window_radius)
            right = min(N, i + window_radius)
            x_smooth[i] = mean(x[left:right])
        end
        return x_smooth
    end
    function smooth_alpha_curve(x; radius=3, passes=2)
        y = copy(x)
        for _ in 1:passes
            y = rolling_average_filter(y, radius)
            y[1] = y[2]
            y[end] = y[end-1]
        end
        return y
    end
    function complexvec_to_json(vec)
        return [Dict("re" => real(x), "im" => imag(x)) for x in vec]
    end
    # v4.6: JSON.jl does not round-trip NaN/Inf, and this file is re-parsed on
    # every iteration (read - push - write), so one NaN would break the run at
    # the next save.  Non-finite diagnostics are stored as null.
    jsonnum(x) = (x isa Real && isfinite(x)) ? x : nothing
    filename = "contour_iteration_v5.json"
end

# ---------------------------------------------------------------------------
# v4.7 CHANGE 5:  zeta_alpha scaled with the branch gap.
#
# WHY.  phi_F = exp(zeta_alpha / r^2) - 1 is what keeps F off the branch
# points, and sqrt(zeta_alpha) is the RANGE of that repulsion.  But F has to
# thread BETWEEN the two branch points, so it must sit at r ~ d_branch/2 from
# each, where the exponent is
#
#     zeta_alpha / (d_branch/2)^2  =  4 zeta_alpha / d_branch^2
#
# With zeta_alpha frozen at 4e-4 that exponent grows as 1/d^2.  Harmless at the
# d = 3.1e-2 the value was hand-tuned at (exponent ~1.7); past EXP_ARG_MAX =
# 400 by d ~ 2e-3, where expc clips and the force is wrong over the whole
# neighbourhood.  sqrt(4e-4) = 0.02 is the repulsion range, and by then the gap
# F must fit through is 5e-3 -- four times narrower.
#
# THE v4.6 RUN HIT EXACTLY THIS, and section 5 of the v4.5 design note called
# it in advance ("watch d_contour and the exp arguments once d_branch drops
# below ~1e-2").  Over its last 156 iterations, at level 6:
#
#     zeta/r^2 past EXP_ARG_MAX on 26 of them, past 100 far more often
#     min |F - branch| down to 7.6e-4
#     threading ratio (|F-au| + |F-al|)/d degraded from 1.00 to 1.1 - 2.5
#     max Im(omega_F) therefore drifted UP, pushing omega_i from 1.2e-6 to
#       4.9e-6 ABOVE the pinch -- and that, not the grid, is what now sets
#       d_branch: the omega_i offset alone forces d = 9.6e-3 against the
#       9.3e-3 observed, while the grid (h = 3.8e-7, miss 2.3e-7) would
#       allow 2.1e-3.
#     frame-to-frame jitter 2.7e-3 against the K/d round-off prediction of
#       3.6e-4 -- 7.5x too large to be arithmetic.  That is F being kicked,
#       not the double-precision floor, so there is real headroom left.
#
# THE FIX.  Hold the exponent at the geometry F actually has to occupy fixed,
# rather than the length sqrt(zeta):
#
#     zeta_target(d) = ZETA_REF * (d / ZETA_D_REF)^2
#
# calibrated on the one setting known to work -- zeta = 4e-4 at d = 3.08e-2,
# the v4.4/v4.5 plateau, where the threading ratio was 1.00 and the exponent
# ~1.7.  Capped at ZETA_REF so it can only ever shrink from the hand-tuned
# value as the gap closes; cap and formula meet exactly at the calibration
# point, so there is no discontinuity.
#
# STABILITY.  This closes a loop: zeta -> F -> omega_F -> omega_i -> d_branch
# -> zeta.  Two guards.  The input is a running MEDIAN of d_branch, because the
# raw value swings +-40 % frame to frame at this level and zeta ~ d^2 would
# square that into a +-80 % kick every iteration.  And zeta is rate-limited to
# ZETA_RATE per iteration, so the potential never jumps under F: 4e-4 down to
# the 4e-7 that d = 1e-3 asks for takes ~230 iterations at 3 %.
#
# NOTHING HERE KNOWS WHERE THE PINCH IS.  d_branch is a measured separation
# between two tracked branches; ZETA_REF and ZETA_D_REF are a working setting
# taken from this project's own run history; and the d^2 law comes from the
# requirement that F fit between two branch points, not from where those
# branch points are.  Set ZETA_ADAPTIVE = false and zeta_alpha stays at 4e-4,
# i.e. exactly v4.6.
# ---------------------------------------------------------------------------
const ZETA_ADAPTIVE = true
const ZETA_REF      = 4.0e-4     # the hand-tuned value ...
const ZETA_D_REF    = 3.08e-2    # ... and the d_branch it was measured healthy at
const ZETA_MIN      = 1.0e-12    # numerical floor; not expected to bind
const ZETA_RATE     = 0.03       # maximum fractional change per iteration
const ZETA_WINDOW   = 25         # iterations in the running median of d_branch

global zeta_hist_d = Float64[]

function zeta_push_d!(d)
    isfinite(d) || return nothing
    push!(zeta_hist_d, d)
    while length(zeta_hist_d) > ZETA_WINDOW
        popfirst!(zeta_hist_d)
    end
    return nothing
end

# Target zeta for a given (smoothed) gap, capped at the hand-tuned value.
zeta_target(d) = clamp(ZETA_REF * (d / ZETA_D_REF)^2, ZETA_MIN, ZETA_REF)

# One rate-limited step of zeta_alpha towards that target.  Returns zeta
# unchanged while the adaptation is off or the median window is not yet full,
# so the first ZETA_WINDOW iterations run at exactly the v4.6 value.
function zeta_step(zeta_now)
    (!ZETA_ADAPTIVE || length(zeta_hist_d) < ZETA_WINDOW) && return zeta_now
    dm = median(zeta_hist_d)
    (isfinite(dm) && dm > 0) || return zeta_now
    zt = zeta_target(dm)
    return clamp(zt, zeta_now * (1.0 - ZETA_RATE), zeta_now * (1.0 + ZETA_RATE))
end

# ---------------------------------------------------------------------------
# v4.9 CHANGE 11/12 CONTROLS
#
# MOVE_FRAC.  The largest displacement any F node may take in one iteration,
# as a fraction of d_branch.  0.05 means "no node moves more than 5 % of the
# gap per step".  The old behaviour was local_delta_t * rhs_cap, which at
# d = 3.8e-3 was 1.0e-2, i.e. 2.6 gaps.
#
# LIMITER_SAT_AT.  Diagnostic only: |local_delta_t * rhs / theta| above this
# counts as saturated for n_limited.  tanh(2) = 0.964, i.e. within 4 % of the
# cap.  This is the direct successor of "how many nodes hit rhs_cap".
#
# DESCENT_*.  The line search of CHANGE 12.  The tolerance is a fraction of the
# run's own recent descent RATE of max Im(omega_F), measured as the difference
# of two half-window medians -- the same shape as adapt_stalled(), and for the
# same reason: the quantity is noisy and a fixed threshold is either always or
# never satisfied.  No pinch quantity enters it.  While the window is filling
# the tolerance is Inf, i.e. the check is inactive, which is what you want for
# the first iterations of a fresh run; a resume seeds the window from the log.
# ---------------------------------------------------------------------------
const MOVE_FRAC          = 0.05     # max node displacement, as a fraction of d_branch
const LIMITER_ON         = true     # false -> clamp(rhs, +-rhs_cap), i.e. v4.7
const LIMITER_SAT_AT     = 2.0      # |u| above this counts towards n_limited
const RHS_CAP            = 100.0    # only used when LIMITER_ON = false

const DESCENT_ON         = true     # false -> apply attempt 1 unconditionally, i.e. v4.7
const DESCENT_WINDOW     = 24       # iterations in the peak-descent-rate window (even)
const DESCENT_TOL_FRAC   = 0.10     # allowed rise, as a fraction of the recent rate
const ALPHA_MAX_ATTEMPTS = 5        # halvings of local_delta_t before giving up

global peak_hist = Float64[]

function peak_push!(p)
    isfinite(p) || return nothing
    push!(peak_hist, p)
    while length(peak_hist) > DESCENT_WINDOW
        popfirst!(peak_hist)
    end
    return nothing
end

# Allowed rise in max Im(omega_F) for one accepted step.  Inf while the window
# fills; -> 0 as the run stalls, which makes the rule strictly monotone exactly
# when the random walk is the thing being fought.
function descent_tolerance()
    length(peak_hist) < DESCENT_WINDOW && return Inf
    h = DESCENT_WINDOW >> 1
    older = median(peak_hist[1:h])
    newer = median(peak_hist[h+1:end])
    rate  = (older - newer) / h          # > 0 while descending
    return isfinite(rate) ? DESCENT_TOL_FRAC * max(rate, 0.0) : Inf
end

# ---------------------------------------------------------------------------
# v4.9 DIAGNOSTICS (logged, nothing acts on them)
#
# f_chord_ripple.  Deviation of y_j from the chord through its two neighbours.
# For a resolved curve this is 0.5*h+*h-*y'', so it must fall as h^2 under
# refinement.  On v4.8 it ROSE 3.5x while h fell 256x -- the direct evidence
# that the refined block was resolving nothing.  3.7e-4 on v4.7, 1.3e-3 on v4.8.
#
# clean_pinch_fit.  Least squares of gap^2 = (8/|omega''|)*(omega_r - omega_p)
# over the part of the d profile where the sqrt law actually holds, one side at
# a time.  On the last frame of either run this recovers omega_p to 1.4e-7 and
# |omega''| to 0.06 %, where pinch_fit's +-2-node d^4 parabola returns
# omega2_fit = 8.8e-7 against the true 0.4294 and grid_share = 468.  The spread
# between the two sides is a free, pinch-independent error bar.  LOGGED ONLY in
# v4.9 -- pinch_fit still drives the adaptive controller, so that this can be
# validated against it on the same frames before anything switches over.
# ---------------------------------------------------------------------------
const CLEAN_GAP_LO = 1.2e-2
const CLEAN_GAP_HI = 9.0e-2

function f_chord_ripple(Fc)
    n = length(Fc)
    n < 3 && return NaN
    x = real.(Fc); y = imag.(Fc)
    m = 0.0
    for j in 2:(n - 1)
        den = x[j+1] - x[j-1]
        den == 0 && continue
        t = (x[j] - x[j-1]) / den
        r = y[j] - (y[j-1] + (y[j+1] - y[j-1]) * t)
        isfinite(r) && (m = max(m, abs(r)))
    end
    return m
end

function clean_pinch_fit(wr, d; lo = CLEAN_GAP_LO, hi = CLEAN_GAP_HI)
    n = length(d)
    (n < 8 || n != length(wr)) && return (NaN, NaN, NaN, 0)
    i  = argmin(d)
    wp = Float64[]; w2 = Float64[]; ntot = 0
    for side in (-1, 1)
        idx = [j for j in 1:n if sign(wr[j] - wr[i]) == side &&
                                 isfinite(d[j]) && lo <= d[j] <= hi]
        length(idx) < 3 && continue
        x = Float64[wr[j] for j in idx]
        y = Float64[d[j]^2 for j in idx]
        xm = sum(x) / length(x); ym = sum(y) / length(y)
        sxx = sum((x .- xm) .^ 2)
        sxx > 0 || continue
        m = sum((x .- xm) .* (y .- ym)) / sxx
        (isfinite(m) && m != 0) || continue
        c = ym - m * xm
        push!(wp, -c / m)
        push!(w2, abs(8.0 / m))
        ntot += length(idx)
    end
    isempty(wp) && return (NaN, NaN, NaN, 0)
    spread = length(wp) == 2 ? abs(wp[1] - wp[2]) : NaN
    return (sum(wp) / length(wp), sum(w2) / length(w2), spread, ntot)
end

# ---------------------------------------------------------------------------
# v5 CONTROLS
#
# CHANGE 14.  Smoothing width, in alpha_r rather than in nodes.  SMOOTH_FRAC is
# set so that at the START of a descent (d ~ 3) the window matches v4.9's
# 15-node box, and shrinks with the gap thereafter.  SMOOTH_W_MAX caps it at
# that same width.
#
# CHANGE 16.  DESC_FRAC is the share of a node's trust region the descent term
# may use on the most strongly driven node; the rest is left to the barrier.
# DESC_NW sets the softmax temperature from the spread of Im(omega_F) over the
# top DESC_NW nodes -- a run-measured quantity, so the weight adapts to the
# sharpness of the peak instead of carrying a fixed scale.
#
# CHANGE 17.  CLEAR_FRAC is the clearance a trial must leave between F and the
# nearest branch point, as a fraction of d_branch.  0.0 means "simply do not
# cross".  Raise it if the run keeps grazing.
#
# CHANGE 18.  MOVE_FRAC is reinterpreted: it is now a fraction of each node's
# OWN distance to the nearest branch, not of the global gap.
# ---------------------------------------------------------------------------
const SMOOTH_MODE    = :width     # :width (v5) | :index (v4.9) | :off
const SMOOTH_FRAC    = 0.05       # half-width in alpha_r, as a fraction of d_branch
const SMOOTH_W_MAX   = 0.1515     # cap: the v4.9 15-node window at h = 1.0101e-2
const SMOOTH_RADIUS  = 7          # nodes, used only by SMOOTH_MODE = :index

const DIFFUSION_IMPLICIT = true   # CHANGE 15; false -> explicit, i.e. v4.9

const DESC_ON        = true       # CHANGE 16; false -> v4.9 (barrier only)
const DESC_FRAC      = 0.30       # share of theta_j the descent term may use
const DESC_NW        = 9          # nodes defining the softmax temperature

const CROSS_CHECK    = true       # CHANGE 17; false -> v4.9 (no side test)
const CLEAR_FRAC     = 0.0        # required clearance, as a fraction of d_branch
# v5.  A veto that can never be overridden is worse than no veto: the v5.1 run
# hit an unsatisfiable crossing test at iteration 70 and then sat on the SAME
# numbers for 843 iterations, burning twelve L-refinement levels on a frozen
# geometry until the log reached 39 MB and the write failed.  After this many
# consecutive fully-discarded iterations the veto is suspended for one step so
# the run can move, and it says so loudly.
const CROSS_ESCAPE_AFTER = 8

global discard_run = 0            # consecutive iterations with no admissible trial

const THETA_PER_NODE = true       # CHANGE 18; false -> the single v4.9 theta

# Linear interpolation of the polyline (xs, ys) at xq, held flat outside the
# ends -- the same boundary condition alpha_i[1] = alpha_i[2] already imposes.
function interp_flat(xs, ys, xq)
    n = length(xs)
    xq <= xs[1]   && return ys[1]
    xq >= xs[n]   && return ys[n]
    j = searchsortedlast(xs, xq)
    j >= n && return ys[n]
    h = xs[j+1] - xs[j]
    h <= 0 && return ys[j]
    t = (xq - xs[j]) / h
    return ys[j] + t * (ys[j+1] - ys[j])
end

# CHANGE 17.  Signed clearance between F and the two branches: positive means
# every upper point is above F and every lower point is below it.  Returns the
# worst (smallest) clearance and how many points are on the wrong side.
function branch_clearance(Fv, alpha_u, alpha_l; gap_bound = 0.5)
    xs = real.(Fv); ys = imag.(Fv)
    lo = xs[1]; hi = xs[end]
    worst = Inf; ncross = 0; any_pair = false
    for k in eachindex(alpha_u)
        (k <= length(alpha_l)) || break
        abs(alpha_u[k] - alpha_l[k]) < gap_bound || continue
        # v5.  ONLY points that lie inside F's own alpha_r span count.  F is a
        # finite polyline on [0,1]; interp_flat holds it level outside that, and
        # that extension is a fiction -- in the construction F is the real axis
        # out there, not a horizontal line at F's end height.  Testing against
        # the fiction is what deadlocked the v5.1 run: at iteration 70 all 28
        # "crossings" sat at alpha_r = 1.22 .. 2.19, outside F entirely, with
        # branch gaps of 1.8 to 2.5 -- far-field points nowhere near a pinch.
        (lo <= real(alpha_u[k]) <= hi) || continue
        (lo <= real(alpha_l[k]) <= hi) || continue
        any_pair = true
        cu = imag(alpha_u[k]) - interp_flat(xs, ys, real(alpha_u[k]))
        cl = interp_flat(xs, ys, real(alpha_l[k])) - imag(alpha_l[k])
        cu < worst && (worst = cu)
        cl < worst && (worst = cl)
        cu <= 0 && (ncross += 1)
        cl <= 0 && (ncross += 1)
    end
    return any_pair ? (worst, ncross) : (Inf, 0)
end

# CHANGE 16.  domega/dalpha at every F node, from the two arrays the loop
# already has.  omega is analytic, so the chord derivative along F is the full
# complex derivative -- no extra eigenproblem is solved anywhere for this.
function domega_dalpha(Fv, wF)
    n = length(Fv)
    g = zeros(ComplexF64, n)
    n < 2 && return g
    for j in 2:(n-1)
        dz = Fv[j+1] - Fv[j-1]
        g[j] = dz == 0 ? 0.0 + 0.0im : (wF[j+1] - wF[j-1]) / dz
    end
    dz1 = Fv[2] - Fv[1];       g[1] = dz1 == 0 ? 0.0+0.0im : (wF[2] - wF[1]) / dz1
    dzn = Fv[n] - Fv[n-1];     g[n] = dzn == 0 ? 0.0+0.0im : (wF[n] - wF[n-1]) / dzn
    for j in eachindex(g)
        isfinite(real(g[j])) && isfinite(imag(g[j])) || (g[j] = 0.0 + 0.0im)
    end
    return g
end

# CHANGE 16.  Softmax weights on Im(omega_F), peaked at the maximum.  The
# temperature is the spread over the top DESC_NW nodes of this same array, so
# it carries no external scale.
function descent_weights(wF)
    v = imag.(wF)
    n = length(v)
    vmax = maximum(v)
    k = min(DESC_NW, n)
    vs = sort(v; rev = true)
    T = vs[1] - vs[k]
    (isfinite(T) && T > 0) || (T = max(abs(vmax), 1.0) * 1e-12)
    w = similar(v)
    @inbounds for j in eachindex(v)
        w[j] = exp(clamp((v[j] - vmax) / T, -50.0, 0.0))
    end
    return w
end

# CHANGE 14.  A filter whose window is a width in alpha_r.  Returns the input
# unchanged when the window is narrower than one cell, which is what happens
# once d_branch is small -- the point of the change.
function width_average(xs, ys, half_w)
    n = length(ys)
    (half_w <= 0 || n < 3) && return copy(ys)
    out = similar(ys)
    @inbounds for j in 1:n
        lo = j; while lo > 1 && xs[j] - xs[lo-1] <= half_w; lo -= 1; end
        hi = j; while hi < n && xs[hi+1] - xs[j] <= half_w; hi += 1; end
        if hi == lo
            out[j] = ys[j]
        else
            acc = 0.0; wsum = 0.0
            for m in lo:hi
                wm = 1.0 - abs(xs[m] - xs[j]) / half_w      # triangular, C0 at the edge
                wm <= 0 && continue
                acc += wm * ys[m]; wsum += wm
            end
            out[j] = wsum > 0 ? acc / wsum : ys[j]
        end
    end
    return out
end

# CHANGE 13 + 15.  The consistent non-uniform second-derivative operator, and
# one implicit step of  y_t = sigma * m_j * L y  with the zero-slope end
# conditions the explicit code already applies.
function diffusion_implicit(xs, ys, mscale, dt, sig)
    n = length(ys)
    (n < 3 || dt <= 0 || sig <= 0) && return copy(ys)
    dl = zeros(Float64, n-1); dg = ones(Float64, n); du = zeros(Float64, n-1)
    rhs = copy(ys)
    @inbounds for j in 2:(n-1)
        hp = xs[j+1] - xs[j]; hm = xs[j] - xs[j-1]
        (hp > 0 && hm > 0) || continue
        a = 2.0 / (hm * (hp + hm))
        c = 2.0 / (hp * (hp + hm))
        f = dt * sig * mscale[j]
        dl[j-1] = -f * a
        du[j]   = -f * c
        dg[j]   = 1.0 + f * (a + c)
    end
    dg[1] = 1.0; du[1] = -1.0; rhs[1] = 0.0          # y_1 = y_2
    dg[n] = 1.0; dl[n-1] = -1.0; rhs[n] = 0.0        # y_n = y_{n-1}
    out = try
        Tridiagonal(dl, dg, du) \ rhs
    catch
        copy(ys)
    end
    return all(isfinite, out) ? out : copy(ys)
end

# ---------------------------------------------------------------------------
# v4.4 omega-descent controls (unchanged)
# ---------------------------------------------------------------------------
const OMEGA_BISECT        = true    # false -> exactly the v4.3 omega update
const OMEGA_GAP_TARGET    = 1e-8    # wanted omega_i - max(imag(omega_F))
# v4.6 CHANGE 3.  Was 1e-3.  On the v4.5 run omega_i overshoots UPWARD on 563
# of the 1752 iterations after k = 250 (32 %), landing 2e-3 to 4.3e-3 above the
# bound, and 561 of those 563 had a gap ABOVE 1e-3 -- so the bisection, which
# exists to pull omega_i back down, refused to engage on exactly the iterations
# that needed it.  Max gap after k = 100 is 4.9e-3, so 1e-2 covers every
# overshoot.  Cost is near zero: probe 1 tries omega_lower_bound itself, which
# on a recovery frame is admissible, so the loop breaks after one probe (v4.5
# averaged 0.35 probes/iteration with a maximum of 1).
const OMEGA_BISECT_ENGAGE = 1e-2    # only bisect once the gap is below this
const OMEGA_BISECT_MAX    = 25      # probes per iteration (each = 1 branch track)
const OMEGA_BISECT_TOL    = 1e-12   # stop when the bracket is this narrow
const OMEGA_BISECT_MIN_DFC = 0.0

# ---------------------------------------------------------------------------
# v4.6 CHANGE 2:  parabolic peak of Im(omega_F) for omega_lower_bound.
#
# omega_lower_bound was maximum(imag.(omega_F)) -- a maximum over 100 discrete
# F nodes.  On the v4.5 final frames the true peak sits 0.28 of a node spacing
# away from the nearest node, so the discrete maximum under-reports it by
# 9.7e-7 and omega_i is allowed that far BELOW Im(omega_p), past the optimum.
#
# The contour itself is not the problem: a 3-point parabola through the peak
# gives -0.572829718733 against the analytic Im(omega_p) = -0.572829694809,
# i.e. 2.4e-8.  Measured on v4.5 frames 700..2000:
#
#     discrete max  - Im(omega_p) :  -8.0e-7 ... -1.2e-6
#     parabolic max - Im(omega_p) :  -3.1e-9 ... -4.2e-8
#
# That moves the ceiling this error puts on d_min from 4.4e-3 to 6.7e-4.
# Fitting in index rather than arclength is fine here: F node spacing is
# 1.0100e-2 to 1.0330e-2, and the two give the same answer to 4e-9.
#
# The estimate is always >= the discrete maximum, so the bound only ever moves
# UP -- the conservative direction.  Guarded: a non-peak (den >= 0), a vertex
# more than one node away, or a non-finite result falls back to the discrete
# maximum, i.e. to exact v4.5 behaviour.
# ---------------------------------------------------------------------------
const OMEGA_LB_PARABOLIC = true

function omegaF_peak_imag(omega_F)
    y = imag.(omega_F)
    j = argmax(y)
    (!OMEGA_LB_PARABOLIC || j == firstindex(y) || j == lastindex(y)) && return y[j]
    y0, y1, y2 = y[j-1], y[j], y[j+1]
    den = y0 - 2 * y1 + y2
    den >= 0 && return y1                       # not a proper interior peak
    t = 0.5 * (y0 - y2) / den                   # vertex offset, in node units
    (!isfinite(t) || abs(t) > 1.0) && return y1
    ypk = y1 - (y2 - y0)^2 / (8 * den)          # >= y1 because den < 0
    return (isfinite(ypk) && ypk >= y1) ? ypk : y1
end

# ---------------------------------------------------------------------------
# v4.5 ADAPTIVE L CONTROLS
#
# TRIGGER.  Not a fixed d_min threshold: the level at which d_min stalls is
# itself set by the current spacing (base grid -> 0.22, v4.4's grid -> 0.031),
# so a fixed threshold either fires immediately or never.  Instead: refine
# when the argmin has stopped moving and d_min has stopped drifting.  That is
# the signal that the contour has parked and only the grid is holding d_min
# up.  On the v4.4 run it first becomes true around iteration 295 and stays
# true for the remaining ~270 iterations.
#
# Both halves of that test have to tolerate noise, and a first version that
# did not was silently useless.  The eigensolver noise on d is delta_d = K/d,
# so its RELATIVE size K/d^2 GROWS as the pinch closes: at d = 3e-2 it is
# 0.3 %, at d = 5e-3 it is 10 %.  Hence
#   - the winner is allowed to jitter by one cell (it flips between the two
#     points bracketing the tip), not required to be identical; and
#   - "d_min stopped changing" is a drift-versus-noise test, not a spread
#     test: compare the two half-window means against the standard error of
#     their difference.  A fixed 1 % spread test stops firing near d = 5e-3
#     and switches the adaptivity off exactly when one more round is still
#     wanted.  Simulated against the measured noise, the fixed test died at
#     level 4 (d = 5.9e-3); this one reaches level 7 (d = 2.2e-3).
#
# ACTION.  Subdivide the ADAPT_HALF_CELLS cells either side of the winner by
# ADAPT_FACTOR.  The tip of the V is bracketed by the winner's two immediate
# neighbours (d is monotone on each side of it), so +-2 cells contains it with
# margin against noise.  Cost: (2*ADAPT_HALF_CELLS)*(ADAPT_FACTOR-1) = 12 new
# points per round, and one extra contour_alpha_L_conti to rebuild the
# branches on the new grid.
#
# STOPPING.  Four independent guards, whichever bites first:
#   ADAPT_D_FLOOR      d_branch is already at the noise floor
#   ADAPT_MIN_GAIN     a round failed to buy anything, ADAPT_NOGAIN_MAX times
#                      in a row -- the grid is no longer the limiter and the
#                      omega_i error (i.e. F's resolution) has taken over
#   ADAPT_MAX_POINTS   point budget
#   ADAPT_MIN_H        spacing floor
# The no-gain rule needs TWO consecutive failures, not one: a single round can
# land the new points badly and gain nothing purely by luck.  In simulation
# level 3 gained 0.4 % and level 4 then gained a factor 6, so a one-strike
# rule would have stopped at d = 2.7e-2, i.e. no better than v4.4.
# ---------------------------------------------------------------------------
const ADAPT_ON           = true
const ADAPT_FACTOR       = 4        # spacing divided by this per round
const ADAPT_HALF_CELLS   = 2        # subdivide omega_r[i-2 .. i+2]
const ADAPT_STALL_WINDOW = 20       # iterations the argmin must hold still
const ADAPT_STALL_TOL    = 1e-2     # relative drift of d_min allowed on top of noise
# v4.6: was 10.  With the score gate a round taken during the descent is no
# longer counted as a failure, but it IS still a level, so the descent can now
# spend levels that v4.5 never reached.  20 levels x 5 nodes = 200 points, well
# inside ADAPT_MAX_POINTS; ADAPT_MIN_H and ADAPT_D_FLOOR remain the real stops.
const ADAPT_MAX_LEVEL    = 20
const ADAPT_MAX_POINTS   = 400
# v5.  Was 1e-11.  Two reasons.  (a) Consecutive alpha values stop being
# distinguishable from eigensolver noise below h = K|omega''|/2 = 5.5e-7, so
# anything finer resolves round-off, not the branch.  (b) With the v5.1 geometry
# frozen, the stall detector fired every 20 iterations and burned levels down to
# h = 2.50e-11 on 191 points, which is what inflated the log to 39 MB.
const ADAPT_MIN_H        = 5e-7
const ADAPT_MIN_GAIN     = 0.90     # d_after/d_before above this counts as no gain
const ADAPT_NOGAIN_MAX   = 2        # consecutive no-gain rounds before giving up
# Below this d_branch, refining resolves eigensolver noise rather than the
# pinch: the v4.4 residual is white with delta_d = K/d, K = 2.56e-6, so
# signal = noise at d = 1.6e-3.
const ADAPT_D_FLOOR      = 2.0e-3

# ---------------------------------------------------------------------------
# v4.6 CHANGE 1 (the one that unblocks the run) and CHANGE 4.
#
# WHAT WENT WRONG IN v4.5.  The adaptive grid fired exactly twice, at k = 193
# and k = 227, and then switched itself off for the remaining 1774 iterations.
# The cause is the no-gain rule being applied at a moment when its answer was
# predetermined.  d_min has two independent causes,
#
#     d = 2 sqrt( 2 |Delta omega_r + i Delta omega_i| / |omega''| )
#
# and refining omega_r can only ever remove the first.  At both refinements
# omega_i was still 1.2e-3 ABOVE the pinch while the omega_r miss was ~1e-4,
# so d was pinned by the vertical offset.  Since
#
#     d(omega_r) >= 2 sqrt( 2 |omega_i - Im omega_p| / |omega''| ) = d_floor
#
# the achievable gain was capped at d_floor/d_branch = 0.9296 (k = 193) and
# 0.9517 (k = 227) BEFORE any refinement ran.  Both exceed ADAPT_MIN_GAIN =
# 0.90, so both rounds were scored as failures no matter what the grid did,
# adapt_nogain hit ADAPT_NOGAIN_MAX, and adapt_done was set -- permanently,
# there is no re-enable path.  The stall detector itself was fine: replayed
# against the recorded history it would have fired at k = 248 and 1607 times
# afterwards, with the winner index steady at 70 in 100 % of 20-iteration
# windows after k = 250.  It simply was never called again.
#
# CHANGE 1 -- SCORE GATE.  Before the no-gain counter is touched, ask how much
# of the current d_branch the horizontal miss actually explains.  Near the
# pinch
#
#     d(omega_r)^4 = (64/|omega''|^2) [ (omega_r - omega_pr)^2
#                                     + (omega_i - omega_pi)^2 ]
#
# so d^4 is an exact parabola in omega_r at leading order, and its vertex gives
# omega_pr without any knowledge of the vertical offset (which only shifts the
# parabola up).  With
#
#     miss    = |omega_r[i_pinch] - omega_pr_fit|
#     d_horiz = 2 sqrt(2 miss / |omega''|_fit)
#
# the ratio d_horiz/d_branch is the share of d the grid is responsible for.
# On the v4.5 frames: 0.47 (k = 192), 0.66 (k = 226), 0.81 (k = 300), then
# 0.97-0.99 from k = 700 on, and 0.04 on an omega_i overshoot frame.  A round
# below ADAPT_GRID_SHARE carries no information about the grid and is not
# scored -- neither killing round would have counted, adapt_nogain stays 0, and
# the run keeps refining.
#
# Failure direction is safe: a failed fit gives grid_limited = false, the round
# is not scored, and adaptivity can only stay alive LONGER.  The level cap, the
# point budget and ADAPT_D_FLOOR still bound the run.
#
# CHANGE 4 -- VERTEX PLACEMENT.  The same fit locates omega_pr far better than
# the grid does, so put the new nodes there instead of subdividing blindly.
# On the v4.5 frames (n = 2, h = 3.16e-4, 1006 quiet frames):
#
#     grid argmin  - Re(omega_p) :  -3.80e-5      (fixed, it is the grid)
#     fitted vertex - Re(omega_p) :  -7.18e-7 mean, 1.55e-7 s.d.
#
# a factor 53 in the miss, i.e. a factor 7 in d, in a SINGLE round against the
# factor 2 per round that subdivision buys.  The residual is a systematic
# stencil bias, not noise (mean/s.d. ~ 5, and it grows monotonically with the
# stencil: -0.72, -1.48, -2.50, -3.77 e-6 for n = 2..5), so it shrinks as the
# cluster tightens.
#
# The noise floor on the vertex is width-independent: an antisymmetric
# perturbation of size 4 d^2 K (from sigma_d = K/d) on a parabola of curvature
# 64/|omega''|^2 shifts the vertex by ~ K |omega''| / 4 = 2.7e-7, matching the
# observed 1.6e-7 s.d.  So the fit cannot locate omega_pr better than ~3e-7,
# which corresponds to d ~ 2.4e-3 -- the same place sqrt(K) = 1.6e-3 puts the
# wall.  Expect this run to converge to d ~ 2-4e-3 in two or three rounds and
# then genuinely stop; that is the float64 limit, not a tuning failure.
#
# NOT DONE HERE, deliberately: zeta ~ d^2/2, the non-uniform alpha_i_rr, and
# the 15-point rolling filter.  They start to matter below d ~ 3e-3.  A
# non-uniform L grid is safe for the omega update: omega_i_vectorization is
# constant, so the d_d_omega_r_Phi_L central difference -- the only place
# omega_r spacing enters -- is identically zero.
# ---------------------------------------------------------------------------
const ADAPT_SCORE_GATE   = true     # false -> v4.5 no-gain bookkeeping
const ADAPT_VERTEX_PLACE = true     # false -> v4.5 refine_omega_r subdivision
const ADAPT_FIT_HALF     = 2        # points either side of the argmin in the d^4 fit
const ADAPT_GRID_SHARE   = 0.90     # horizontal miss must explain this much of d
const ADAPT_CLUSTER_HALF = 2        # 2*this+1 nodes inserted around the vertex
const ADAPT_VERTEX_MAXMOVE = 2.0    # reject a vertex further than this * h_local

# v4.7 CHANGE 6 (fixing a v4.6 bug of mine).  grid_share had a lower bound only.
# Once the vertical offset dominates, the varying part of d^4 across the stencil
# is a couple of per cent of the constant part, the fitted curvature is then
# dominated by noise and collapses, and grid_share comes out at 200-400 --
# which sails past a ">= 0.90" test and reads as "definitely grid-limited".
# On the v4.6 run |omega''|_fit had a median of 0.19 against the true 0.4294
# after k = 412, and grid_share exceeded 1.5 on 48 of 156 iterations.  A share
# far above 1 is a broken fit, not a strong signal: bound it on both sides.
# The degeneracy is itself informative -- the fit falling apart is the run
# telling you the grid has stopped being the limiter.
const ADAPT_GRID_SHARE_MAX = 1.50

# v4.7 CHANGE 7: measure "the winner has parked" in omega_r against the local
# cell width, instead of in grid-index units.  See adapt_stalled().
const ADAPT_STALL_BY_OMEGA = true
const ADAPT_WINDOW_CELLS   = 2.0    # how far the winning omega_r may wander

# v4.4 used 1e-4, which is under the double-precision floor above and so could
# never fire.  v4.5 used 3e-3, just clear of it.
#
# v4.6: 5e-4, i.e. deliberately unreachable.  The point of this run is to see
# WHERE it plateaus, and 3e-3 would stop it exactly in the interesting band --
# and could stop it spuriously, since at d = 2e-3 the jitter K/d is 1.3e-3, so
# a single frame can dip through a 3e-3 threshold by luck alone.  Raise this
# back to ~2e-3 once the plateau is known and you want the run to terminate on
# success.
const PINCH_TOL          = 5.0e-4

global adapt_level  = 0
global adapt_nogain = 0
global adapt_done   = false
global adapt_hist_i = Int[]
global adapt_hist_d = Float64[]
global adapt_hist_w = Float64[]   # v4.7: the WINNING omega_r, not just its index
global adapt_hist_h = Float64[]   # v4.7: h_local at that moment

function adapt_push!(i, d, w, h)
    push!(adapt_hist_i, i)
    push!(adapt_hist_d, d)
    push!(adapt_hist_w, w)
    push!(adapt_hist_h, h)
    while length(adapt_hist_i) > ADAPT_STALL_WINDOW
        popfirst!(adapt_hist_i)
        popfirst!(adapt_hist_d)
        popfirst!(adapt_hist_w)
        popfirst!(adapt_hist_h)
    end
    return nothing
end

function adapt_reset_history!()
    empty!(adapt_hist_i)
    empty!(adapt_hist_d)
    empty!(adapt_hist_w)
    empty!(adapt_hist_h)
    return nothing
end

# True when the argmin has parked AND d_min's drift is inside the noise.
# The history is cleared on every refinement, so the entries in a window are
# always comparable (omega_r cannot change inside one).
#
# v4.7 CHANGE 7: the "has it parked" test is now in omega_r, not in index.
# v4.6 required the winning INDEX to hold within +-1.  That is the same thing
# as +-1 cell on a uniform grid, but once vertex placement has packed a dozen
# nodes inside one old cell the winner hops several indices while barely
# moving in omega_r -- and the test then can never pass.  That is what froze
# the v4.6 run at level 6 (winner spread 6 to 9 indices over every 20-iteration
# window after k = 420).  Measuring the same thing in omega_r, against the
# local cell width, is index-count-independent and reduces to the old test on
# a uniform grid.  Set ADAPT_STALL_BY_OMEGA = false for v4.6 behaviour.
function adapt_stalled()
    n = length(adapt_hist_d)
    n < ADAPT_STALL_WINDOW && return false
    if ADAPT_STALL_BY_OMEGA
        href = minimum(adapt_hist_h)
        (isfinite(href) && href > 0) || return false
        (maximum(adapt_hist_w) - minimum(adapt_hist_w)) > ADAPT_WINDOW_CELLS * href && return false
    else
        # the winner may flip between the two points bracketing the tip
        (maximum(adapt_hist_i) - minimum(adapt_hist_i)) > 1 && return false
    end
    m = abs(mean(adapt_hist_d))
    m <= 0.0 && return false
    h = n ÷ 2
    drift = abs(mean(adapt_hist_d[(h + 1):end]) - mean(adapt_hist_d[1:h]))
    noise = 2.0 * std(adapt_hist_d) / sqrt(n)   # s.e. of the half-mean difference
    return drift <= ADAPT_STALL_TOL * m + noise
end

# Local cell width at index i (the smaller of the two neighbouring cells).
function local_h(wr, i)
    n = length(wr)
    n < 2 && return Inf
    if i == 1
        return wr[2] - wr[1]
    elseif i == n
        return wr[n] - wr[n-1]
    end
    return min(wr[i+1] - wr[i], wr[i] - wr[i-1])
end

# Subdivide the cells [i-half, i+half] by `factor`.  Everything outside is
# untouched, so the grid stays sorted and strictly increasing.
function refine_omega_r(wr, i; factor = ADAPT_FACTOR, half = ADAPT_HALF_CELLS)
    n = length(wr)
    lo = max(1, i - half)
    hi = min(n, i + half)
    hi <= lo && return copy(wr)
    out = Float64[]
    append!(out, wr[1:(lo - 1)])
    for c in lo:(hi - 1)
        seg = collect(range(wr[c], wr[c + 1], length = factor + 1))
        append!(out, seg[1:(end - 1)])   # right endpoint comes from the next cell
    end
    push!(out, wr[hi])
    append!(out, wr[(hi + 1):n])
    @assert issorted(out)
    @assert minimum(diff(out)) > 0
    return out
end

# Least-squares slope of y on x.
function _slope(x, y)
    n = length(x)
    n < 2 && return NaN
    mx = mean(x); my = mean(y)
    sxx = sum((x .- mx) .^ 2)
    sxx <= 0 && return NaN
    return sum((x .- mx) .* (y .- my)) / sxx
end

# Diagnostic only, never used to place points.  Near a pinch
# d^2 = (8/|omega''|) |omega_r - omega_r_p|, so the slope of d^2 against
# omega_r on either arm gives |omega''|.  On the v4.4 data this returns 0.477
# against the true 0.429 (11 % high) -- fine as a sanity number, and the
# reason the vertex is NOT extrapolated: the same fit locates the tip only
# 1.8x better than the grid already does, because the cubic term bends the
# arms.  Bracket and subdivide instead.
function omega2_estimate(wr, d, i; n = 4)
    lo = max(firstindex(d), i - n)
    hi = min(lastindex(d), i + n)
    (i - lo < 2 || hi - i < 2) && return NaN
    sl = _slope(wr[lo:i], d[lo:i] .^ 2)
    sr = _slope(wr[i:hi], d[i:hi] .^ 2)
    (isnan(sl) || isnan(sr)) && return NaN
    s = 0.5 * (abs(sl) + abs(sr))
    return s > 0 ? 8.0 / s : NaN
end

# Predicted reachable d for a half-cell miss on a grid of spacing h.
adapt_predicted_d(h, w2) = (isnan(w2) || w2 <= 0) ? NaN : 2 * sqrt(2 * (h / 2) / w2)

# ---------------------------------------------------------------------------
# v4.6: locate the pinch in omega_r from the shape of the V.
#
# Near the pinch  d^4 = (64/|omega''|^2) [ (omega_r - omega_pr)^2 + v^2 ]  with
# v = omega_i - Im(omega_p).  So d^4 is a parabola in omega_r whose vertex is
# omega_pr and whose curvature is 64/|omega''|^2; v only shifts it upward and
# drops out of both.  Returns (omega_pr, |omega''|), or (NaN, NaN) if the fit
# is not usable.
#
# This is NOT omega2_estimate.  That fits d^2 linearly along each arm, where
# the |.| kink and the cubic term cost it 11 % on |omega''| -- the reason v4.5
# kept it diagnostic and refused to extrapolate the vertex.  The d^4 form has
# neither problem: on the v4.5 data it returns |omega''| = 0.4295 against the
# true 0.42941 (0.02 %) and the vertex to 7.2e-7 mean / 1.6e-7 s.d. over 1006
# frames, against the grid argmin's fixed 3.8e-5.
# ---------------------------------------------------------------------------
function pinch_fit(wr, d, i; n = ADAPT_FIT_HALF)
    lo = max(firstindex(d), i - n)
    hi = min(lastindex(d), i + n)
    (i - lo < 2 || hi - i < 2) && return (NaN, NaN)
    x = wr[lo:hi] .- wr[i]
    y = d[lo:hi] .^ 4
    (all(isfinite, x) && all(isfinite, y)) || return (NaN, NaN)
    xs = maximum(abs, x)
    ys = maximum(abs, y)
    (isfinite(xs) && isfinite(ys) && xs > 0 && ys > 0) || return (NaN, NaN)
    u = x ./ xs                       # both O(1), so the Vandermonde stays tame
    z = y ./ ys
    c = hcat(u .^ 2, u, ones(length(u))) \ z
    (all(isfinite, c) && c[1] > 0) || return (NaN, NaN)
    w_pr = wr[i] - xs * c[2] / (2 * c[1])
    curv = ys * c[1] / xs^2           # = 64 / |omega''|^2
    return (w_pr, curv > 0 ? 8.0 / sqrt(curv) : NaN)
end

# Insert 2m+1 nodes of spacing h centred on w_c, keeping every existing node.
# New nodes outside the grid, or within h/8 of a node already present, are
# dropped.  Returns a sorted, strictly increasing grid.
function insert_cluster(wr, w_c, h; m = ADAPT_CLUSTER_HALF)
    (isfinite(w_c) && isfinite(h) && h > 0) || return copy(wr)
    (w_c <= wr[1] || w_c >= wr[end]) && return copy(wr)
    out = collect(wr)
    for j in -m:m
        w = w_c + j * h
        (w <= wr[1] || w >= wr[end]) && continue
        minimum(abs.(out .- w)) < h / 8 && continue
        push!(out, w)
    end
    sort!(out)
    @assert issorted(out)
    @assert minimum(diff(out)) > 0
    return out
end

# ---------------------------------------------------------------------------
# v4.7 CHANGE 8:  full-state checkpoint / resume.
#
# WHY.  The v4.5 resume block restored the GEOMETRY -- L, F, the two branches,
# the adapted grid -- and nothing else.  That was complete for v4.5, where
# every control constant was frozen for the whole run.  It is NOT complete for
# v4.7, because CHANGE 5 made zeta_alpha a STATE VARIABLE.  A resume that
# leaves it out restarts the continuation at ZETA_REF = 4e-4.  From a segment
# that ended at zeta_alpha = 2.3e-5 that is a 17.7x jump, held for
# ZETA_WINDOW = 25 iterations while the median window refills (zeta_step
# returns its argument unchanged until then), and only then walked back down at
# ZETA_RATE = 3 %/iteration -- ln(17.7)/ln(1.03) ~ 97 iterations more.  So
# about 120 iterations of every resumed segment would run at the wrong
# repulsion strength with exp_arg 17.7x too large: straight back into the
# EXP_ARG_MAX clipping that CHANGE 5 exists to remove.  adapt_done,
# adapt_nogain and the stall history were dropped as well, so a run that had
# already retired the refinement would silently switch it back on.
#
# WHAT IS SAVED.  Everything the iteration carries across k, stored as the
# PRIMITIVES rather than the derived arrays so the two cannot disagree:
#
#     omega_i, omega_r          ->  L = contour_L()
#     alpha_i                   ->  F = contour_F()
#     alpha_L_u, alpha_L_l          (a full re-track if lost)
#     omega_F                       (a pmap if lost)
#     zeta_alpha, zeta_hist_d       CHANGE 5 state
#     adapt_level, adapt_nogain, adapt_done, adapt_hist_{i,d,w,h}
#     delta_t
#
# WHERE.  A separate small file, rewritten in full every iteration.  NOT the
# JSON log: that log is a diagnostic record, it was 19 MB at k = 567 on v4.6,
# and its entries are written BEFORE the adaptive refinement block, so the last
# entry always describes the PRE-refinement grid.  The checkpoint is written at
# the very END of the iteration, so it describes the state the next iteration
# actually starts from, refinement included.
#
# WHEN THE TWO DISAGREE.  The log entry for iteration j is written before the
# checkpoint that says "next is j+1", so a kill in between leaves the log one
# entry ahead of the state.  On resume the log is truncated back to agree with
# the checkpoint and the original is kept once as <log>.pretruncate.
#
# NOTHING HERE KNOWS WHERE THE PINCH IS.  It stores and restores values the run
# produced, and compares constants against this same source file.
# ---------------------------------------------------------------------------
const CKPT_FILE     = "checkpoint_v5.json"
const CKPT_VERSION  = "v5-ckpt-1"
const CKPT_EVERY    = 1        # write a checkpoint every N steps; 0 disables

# RESUME = true   ->  continue an existing run if a checkpoint (or, failing
#                     that, a log) is there.  This is the default: a run that
#                     stopped at 500 picks up at 501 on the next launch.
# RESUME = false  ->  start from scratch.  The existing log and checkpoint are
#                     MOVED ASIDE to .bak1/.bak2/... rather than overwritten.
const RESUME        = true
# Refuse to resume across a change to any constant in ckpt_constants().  A
# continuation under different numerics is a different experiment and must not
# be written into one log as if it were not.
const RESUME_STRICT = true

# Absolute stopping point, counted in the JSON "iteration" field, NOT in k.
# Segment a run by raising this between launches: 500, then 1000, then 1500.
# MAX_STEPS is only a safety bound on how many steps ONE process will take.
const ITER_TARGET = 1000
const MAX_STEPS   = 2000

# ---------------------------------------------------------------------------
# Constants that define WHAT IS BEING COMPUTED.  Run control (ITER_TARGET,
# MAX_STEPS, file names, CKPT_EVERY) is deliberately NOT in here -- those are
# meant to change from segment to segment.
# ---------------------------------------------------------------------------
function ckpt_constants()
    return Dict{String,Any}(
        "ZETA_ADAPTIVE"        => ZETA_ADAPTIVE,
        "ZETA_REF"             => ZETA_REF,
        "ZETA_D_REF"           => ZETA_D_REF,
        "ZETA_MIN"             => ZETA_MIN,
        "ZETA_RATE"            => ZETA_RATE,
        "ZETA_WINDOW"          => ZETA_WINDOW,
        "OMEGA_BISECT"         => OMEGA_BISECT,
        "OMEGA_GAP_TARGET"     => OMEGA_GAP_TARGET,
        "OMEGA_BISECT_ENGAGE"  => OMEGA_BISECT_ENGAGE,
        "OMEGA_BISECT_MAX"     => OMEGA_BISECT_MAX,
        "OMEGA_BISECT_TOL"     => OMEGA_BISECT_TOL,
        "OMEGA_BISECT_MIN_DFC" => OMEGA_BISECT_MIN_DFC,
        "OMEGA_LB_PARABOLIC"   => OMEGA_LB_PARABOLIC,
        "ADAPT_ON"             => ADAPT_ON,
        "ADAPT_FACTOR"         => ADAPT_FACTOR,
        "ADAPT_HALF_CELLS"     => ADAPT_HALF_CELLS,
        "ADAPT_STALL_WINDOW"   => ADAPT_STALL_WINDOW,
        "ADAPT_STALL_TOL"      => ADAPT_STALL_TOL,
        "ADAPT_MAX_LEVEL"      => ADAPT_MAX_LEVEL,
        "ADAPT_MAX_POINTS"     => ADAPT_MAX_POINTS,
        "ADAPT_MIN_H"          => ADAPT_MIN_H,
        "ADAPT_MIN_GAIN"       => ADAPT_MIN_GAIN,
        "ADAPT_NOGAIN_MAX"     => ADAPT_NOGAIN_MAX,
        "ADAPT_D_FLOOR"        => ADAPT_D_FLOOR,
        "ADAPT_SCORE_GATE"     => ADAPT_SCORE_GATE,
        "ADAPT_VERTEX_PLACE"   => ADAPT_VERTEX_PLACE,
        "ADAPT_FIT_HALF"       => ADAPT_FIT_HALF,
        "ADAPT_GRID_SHARE"     => ADAPT_GRID_SHARE,
        "ADAPT_GRID_SHARE_MAX" => ADAPT_GRID_SHARE_MAX,
        "ADAPT_CLUSTER_HALF"   => ADAPT_CLUSTER_HALF,
        "ADAPT_VERTEX_MAXMOVE" => ADAPT_VERTEX_MAXMOVE,
        "ADAPT_STALL_BY_OMEGA" => ADAPT_STALL_BY_OMEGA,
        "ADAPT_WINDOW_CELLS"   => ADAPT_WINDOW_CELLS,
        "PINCH_TOL"            => PINCH_TOL,
        "EXP_ARG_MAX"          => EXP_ARG_MAX,
        "TRACK_SIDE_SLACK"     => TRACK_SIDE_SLACK,
        "TRACK_MIN_STEP"       => TRACK_MIN_STEP,
        "EXP_ARG_TARGET"       => EXP_ARG_TARGET,
        "MOVE_FRAC"            => MOVE_FRAC,
        "LIMITER_ON"           => LIMITER_ON,
        "RHS_CAP"              => RHS_CAP,
        "DESCENT_ON"           => DESCENT_ON,
        "DESCENT_WINDOW"       => DESCENT_WINDOW,
        "DESCENT_TOL_FRAC"     => DESCENT_TOL_FRAC,
        "ALPHA_MAX_ATTEMPTS"   => ALPHA_MAX_ATTEMPTS,
        "SMOOTH_MODE"          => String(SMOOTH_MODE),
        "SMOOTH_FRAC"          => SMOOTH_FRAC,
        "SMOOTH_W_MAX"         => SMOOTH_W_MAX,
        "SMOOTH_RADIUS"        => SMOOTH_RADIUS,
        "DIFFUSION_IMPLICIT"   => DIFFUSION_IMPLICIT,
        "DESC_ON"              => DESC_ON,
        "DESC_FRAC"            => DESC_FRAC,
        "DESC_NW"              => DESC_NW,
        "CROSS_CHECK"          => CROSS_CHECK,
        "CLEAR_FRAC"           => CLEAR_FRAC,
        "THETA_PER_NODE"       => THETA_PER_NODE,
        "num_modes"            => num_modes,
        "N"                    => N,
        "s_omega"              => s_omega,
        "s_alpha"              => s_alpha,
        "epsilon"              => epsilon,
        "lambda"               => lambda,
        "sigma"                => sigma,
    )
end

function ckpt_constants_diff(stored)
    now = ckpt_constants()
    diffs = String[]
    for key in sort(collect(keys(now)))
        if !haskey(stored, key)
            push!(diffs, string(key, ": absent from the checkpoint, now ", now[key]))
            continue
        end
        a = stored[key]
        b = now[key]
        same = (a isa Number && b isa Number) ? (Float64(a) == Float64(b)) : (a == b)
        same || push!(diffs, string(key, ": checkpoint ", a, "  ->  now ", b))
    end
    return diffs
end

# ---------------------------------------------------------------------------
# open(path, "w") TRUNCATES FIRST and only then spends seconds writing tens of
# MB.  A kill inside that window destroys the whole run history -- which is
# exactly what a resume mechanism must not depend on.  Write to a temporary
# file and rename instead, so the failure window is the rename, not the write.
# ---------------------------------------------------------------------------
function atomic_write(path, str; tries = 5)
    tmp = path * ".tmp"
    open(tmp, "w") do io
        write(io, str)
    end
    # On Windows the rename has to delete the target first, and that fails if
    # anything else holds the file open -- an editor, a plotting script, a
    # virus scanner.  Losing a 10-hour run to a transient lock is not
    # acceptable, so retry, and fall back to a direct write rather than throw.
    for attempt in 1:tries
        try
            mv(tmp, path; force = true)
            return nothing
        catch err
            if attempt == tries
                @printf("[io] rename onto %s failed %d times (%s); writing in place instead\n",
                        path, tries, sprint(showerror, err))
                flush(stdout)
                open(path, "w") do io
                    write(io, str)
                end
                try
                    rm(tmp; force = true)
                catch
                end
                return nothing
            end
            sleep(0.2 * attempt)
        end
    end
    return nothing
end

function ckpt_state_ok()
    isfinite(omega_i)    || return false, "omega_i"
    isfinite(zeta_alpha) || return false, "zeta_alpha"
    isfinite(delta_t)    || return false, "delta_t"
    all(isfinite, omega_r)   || return false, "omega_r"
    all(isfinite, alpha_i)   || return false, "alpha_i"
    all(isfinite, alpha_L_u) || return false, "alpha_L_u"
    all(isfinite, alpha_L_l) || return false, "alpha_L_l"
    all(isfinite, omega_F)   || return false, "omega_F"
    return true, ""
end

function write_checkpoint!()
    ok, bad = ckpt_state_ok()
    if !ok
        @printf("[ckpt] NOT written before iteration %d: %s holds non-finite values; the previous checkpoint stands\n",
                iteration_step, bad)
        flush(stdout)
        return false
    end
    # JSON.jl does not round-trip NaN/Inf, and the four adapt_hist arrays have
    # to stay index-aligned, so if any one of them is polluted the whole window
    # is dropped.  That costs ADAPT_STALL_WINDOW iterations of stall detection
    # -- exactly what every refinement already costs via adapt_reset_history!.
    hist_ok = all(isfinite, adapt_hist_d) && all(isfinite, adapt_hist_w) &&
              all(isfinite, adapt_hist_h)
    ckpt = Dict{String,Any}(
        "ckpt_version"   => CKPT_VERSION,
        "written_at"     => Base.Libc.strftime("%Y-%m-%d %H:%M:%S", time()),
        "next_iteration" => iteration_step,
        "constants"      => ckpt_constants(),
        # geometry, as primitives
        "omega_i"        => omega_i,
        "omega_r"        => collect(Float64, omega_r),
        "alpha_i"        => collect(Float64, alpha_i),
        "alpha_L_u"      => complexvec_to_json(alpha_L_u),
        "alpha_L_l"      => complexvec_to_json(alpha_L_l),
        "omega_F"        => complexvec_to_json(omega_F),
        # v4.7 CHANGE 5 state -- the part the v4.5 resume used to drop
        "zeta_alpha"     => zeta_alpha,
        "zeta_hist_d"    => collect(Float64, zeta_hist_d),
        # v4.9 CHANGE 12 state.  Without this the descent tolerance is Inf for
        # DESCENT_WINDOW iterations after every resume, i.e. the line search is
        # off exactly across the joins.  Same reasoning as CHANGE 8b for zeta.
        "peak_hist"      => all(isfinite, peak_hist) ? collect(Float64, peak_hist) : Float64[],
        "delta_t"        => delta_t,
        # adaptive-refinement state
        "adapt_level"    => adapt_level,
        "adapt_nogain"   => adapt_nogain,
        "adapt_done"     => adapt_done,
        "adapt_hist_ok"  => hist_ok,
        "adapt_hist_i"   => hist_ok ? collect(Int, adapt_hist_i)     : Int[],
        "adapt_hist_d"   => hist_ok ? collect(Float64, adapt_hist_d) : Float64[],
        "adapt_hist_w"   => hist_ok ? collect(Float64, adapt_hist_w) : Float64[],
        "adapt_hist_h"   => hist_ok ? collect(Float64, adapt_hist_h) : Float64[],
    )
    atomic_write(CKPT_FILE, JSON.json(ckpt))
    return true
end

# Drop log entries at or past the checkpoint's next_iteration so the log and
# the restored state describe the same run.  The original is kept once.
function truncate_log_to!(next_iter)
    if !isfile(filename)
        # A checkpoint with no log: the loop's read-push-write would throw on
        # its first save.  Start an empty array so the run can continue; the
        # history before this point is gone and the operator should know.
        @printf("!! %s is missing but %s is not.  Starting an EMPTY log -- every entry before iteration %d is lost.\n",
                filename, CKPT_FILE, next_iter)
        atomic_write(filename, JSON.json(Any[]))
        return 0
    end
    data = JSON.parse(read(filename, String))
    keep = [e for e in data if e["iteration"] < next_iter]
    dropped = length(data) - length(keep)
    if dropped > 0
        backup = filename * ".pretruncate"
        isfile(backup) || cp(filename, backup)
        atomic_write(filename, JSON.json(keep))
    end
    return dropped
end

# RESUME = false must not silently delete a finished run.
function archive_existing!()
    moved = String[]
    for p in (filename, CKPT_FILE)
        isfile(p) || continue
        n = 1
        while isfile(string(p, ".bak", n))
            n += 1
        end
        dest = string(p, ".bak", n)
        mv(p, dest)
        push!(moved, dest)
    end
    return moved
end

# Iteration stamped "resumed" in the log, so every join is visible in the data
# rather than having to be remembered.  0 = this process started a fresh run.
global resume_join_at = 0

begin
    ckpt = nothing
    if RESUME && isfile(CKPT_FILE)
        ckpt = JSON.parse(read(CKPT_FILE, String))
        diffs = ckpt_constants_diff(get(ckpt, "constants", Dict{String,Any}()))
        if !isempty(diffs)
            println()
            println("!! CHECKPOINT CONSTANTS DIFFER FROM THIS SOURCE FILE:")
            for d in diffs
                println("!!     ", d)
            end
            if RESUME_STRICT
                error("Refusing to resume across a constant change -- the continuation " *
                      "would not be the same run as the segment already in $(filename). " *
                      "Either restore the constants, set RESUME_STRICT = false if the " *
                      "hybrid is what you want, or set RESUME = false to start fresh.")
            else
                println("!! RESUME_STRICT = false: continuing anyway.  THE SEGMENTS IN")
                println("!! $(filename) ARE NOT THE SAME EXPERIMENT.")
                println()
            end
        end
    end

    if ckpt !== nothing
        # ---------------- full-state resume ----------------
        global iteration_step = Int(ckpt["next_iteration"])
        global omega_i   = Float64(ckpt["omega_i"])
        global omega_r   = Float64[x for x in ckpt["omega_r"]]
        global alpha_i   = Float64[x for x in ckpt["alpha_i"]]
        global alpha_L_u = ComplexF64[complex(x["re"], x["im"]) for x in ckpt["alpha_L_u"]]
        global alpha_L_l = ComplexF64[complex(x["re"], x["im"]) for x in ckpt["alpha_L_l"]]
        global omega_F   = ComplexF64[complex(x["re"], x["im"]) for x in ckpt["omega_F"]]
        # v4.7 state
        global zeta_alpha   = Float64(ckpt["zeta_alpha"])
        global delta_t      = Float64(ckpt["delta_t"])
        global adapt_level  = Int(ckpt["adapt_level"])
        global adapt_nogain = Int(ckpt["adapt_nogain"])
        global adapt_done   = ckpt["adapt_done"] === true
        empty!(zeta_hist_d)
        append!(zeta_hist_d, Float64[x for x in ckpt["zeta_hist_d"]])
        empty!(peak_hist)
        if haskey(ckpt, "peak_hist")
            append!(peak_hist, Float64[x for x in ckpt["peak_hist"]])
        end
        adapt_reset_history!()
        if get(ckpt, "adapt_hist_ok", false) === true
            append!(adapt_hist_i, Int[x     for x in ckpt["adapt_hist_i"]])
            append!(adapt_hist_d, Float64[x for x in ckpt["adapt_hist_d"]])
            append!(adapt_hist_w, Float64[x for x in ckpt["adapt_hist_w"]])
            append!(adapt_hist_h, Float64[x for x in ckpt["adapt_hist_h"]])
        end
        # derived from the primitives above, never stored, so they cannot
        # disagree with them
        global L = contour_L()
        global F = contour_F()
        push_omega_r()
        load_on_workers()
        @everywhere begin
            normals_F = contour_normals(F)
        end
        dropped = truncate_log_to!(iteration_step)
        global resume_join_at = iteration_step
        println()
        println("====================================================================")
        @printf("RESUMED from %s (written %s)\n", CKPT_FILE, get(ckpt, "written_at", "?"))
        @printf("  next iteration   : %d   (target %d)\n", iteration_step, ITER_TARGET)
        @printf("  omega_i          : %.15e\n", omega_i)
        @printf("  N_L / adapt_level: %d / %d   (adapt_done = %s, nogain = %d)\n",
                length(omega_r), adapt_level, adapt_done, adapt_nogain)
        @printf("  zeta_alpha       : %.6e   (median window %d/%d full)\n",
                zeta_alpha, length(zeta_hist_d), ZETA_WINDOW)
        @printf("  stall history    : %d/%d entries\n",
                length(adapt_hist_d), ADAPT_STALL_WINDOW)
        @printf("  delta_t          : %.6e\n", delta_t)
        if dropped > 0
            @printf("  log truncated    : %d entr%s past the checkpoint dropped; original kept as %s.pretruncate\n",
                    dropped, dropped == 1 ? "y" : "ies", filename)
        end
        println("====================================================================")
        println()
        flush(stdout)

    elseif RESUME && isfile(filename)
        # ---------------- legacy geometry-only resume (v4.5 path) ----------------
        println()
        println("!! NO $(CKPT_FILE) FOUND -- falling back to the v4.5 log-based resume.")
        println("!! v4.7 CHANGE 8b: zeta_alpha and the zeta median window ARE now")
        println("!! recovered from the log below, so the repulsion strength is")
        println("!! continuous across the join.  What is still lost: the stall")
        println("!! history (20 iterations of no refinement, the same cost every")
        println("!! refinement round already pays), adapt_nogain, and adapt_done.")
        println("!! Use this only to pick up a run started before checkpointing")
        println("!! existed; every launch after this one will find the checkpoint.")
        println()
        resume_array = JSON.parse(read(filename, String))
        resume_entry = resume_array[end]
        global iteration_step = resume_entry["iteration"] + 1
        global L         = ComplexF64[complex(x["re"], x["im"]) for x in resume_entry["L"]]
        global alpha_L_u = ComplexF64[complex(x["re"], x["im"]) for x in resume_entry["alpha_L_u"]]
        global alpha_L_l = ComplexF64[complex(x["re"], x["im"]) for x in resume_entry["alpha_L_l"]]
        global F         = ComplexF64[complex(x["re"], x["im"]) for x in resume_entry["F"]]
        global omega_F   = ComplexF64[complex(x["re"], x["im"]) for x in resume_entry["omega_F"]]
        global omega_i   = imag(L[1])
        global alpha_i   = imag.(F)
        # v4.5: the grid IS the stored L, not the constant built at load time.
        global omega_r      = real.(L)
        global adapt_level  = get(resume_entry, "adapt_level", 0)
        global adapt_nogain = 0
        global adapt_done   = false

        # v4.7 CHANGE 8b: the legacy path was written as if the log carried no
        # zeta information.  The v4.7 log does -- every entry has zeta_alpha,
        # and d_branch is the only thing zeta_hist_d is ever fed -- so the
        # CHANGE 5 state does NOT have to be thrown away here.  Dropping it
        # cost, resuming the 501-iteration run:
        #     zeta_alpha 1.0029e-05 -> ZETA_REF 4.0e-04           = 39.9x
        #     exp_arg    2.75       -> 121 on the first step
        #     recovery   25 frozen + ln(39.9)/ln(1.03) ~ 125      ~ 150 iters
        # i.e. 150 iterations back inside the EXP_ARG_MAX clipping CHANGE 5
        # exists to remove.  With this block, zeta_step has a full window on
        # the first step and moves 1.0029e-05 -> 9.73e-06 instead.
        z_log = get(resume_entry, "zeta_alpha", nothing)
        if z_log isa Real && isfinite(z_log) && z_log > 0
            global zeta_alpha = float(z_log)
            @printf("   restored zeta_alpha = %.6e from the log entry\n", zeta_alpha)
        else
            @printf("   log entry has no usable zeta_alpha; keeping ZETA_REF = %.3e\n",
                    ZETA_REF)
        end
        empty!(zeta_hist_d)
        for e_prev in resume_array[max(1, length(resume_array) - ZETA_WINDOW + 1):end]
            db_prev = get(e_prev, "d_branch", nothing)
            db_prev isa Real && zeta_push_d!(float(db_prev))
        end
        @printf("   restored zeta median window: %d of %d entries%s\n",
                length(zeta_hist_d), ZETA_WINDOW,
                length(zeta_hist_d) < ZETA_WINDOW ?
                    "  (short -- zeta stays frozen until it fills)" : "")

        # v4.9: rebuild the peak-descent window from the stored omega_F arrays,
        # so the CHANGE 12 line search is live on the first step after a join
        # instead of running with an infinite tolerance for DESCENT_WINDOW
        # iterations.
        empty!(peak_hist)
        for e_prev in resume_array[max(1, length(resume_array) - DESCENT_WINDOW + 1):end]
            wF_prev = get(e_prev, "omega_F", nothing)
            if wF_prev isa AbstractVector && !isempty(wF_prev)
                peak_push!(omegaF_peak_imag(
                    ComplexF64[complex(x["re"], x["im"]) for x in wF_prev]))
            end
        end
        @printf("   restored peak window: %d of %d entries%s\n",
                length(peak_hist), DESCENT_WINDOW,
                length(peak_hist) < DESCENT_WINDOW ?
                    "  (short -- descent tolerance stays Inf until it fills)" : "")

        adapt_reset_history!()
        push_omega_r()
        load_on_workers()
        @everywhere begin
            normals_F = contour_normals(F)
        end
        global resume_join_at = iteration_step
        @printf("RESUMED (geometry only) from iteration %d: omega_i = %.9e, N_L = %d, adapt_level = %d, %d entries in %s\n",
                resume_entry["iteration"], omega_i, length(omega_r), adapt_level,
                length(resume_array), filename)
        flush(stdout)

    else
        # ---------------- fresh run ----------------
        if isfile(filename) || isfile(CKPT_FILE)
            moved = archive_existing!()
            println("RESUME = false: existing run moved aside -> ", join(moved, ", "))
        end
        iteration_step = 1
        dict_to_JSON = Dict(
            "iteration" => iteration_step,
            "L" => complexvec_to_json(L),
            "alpha_L_u" => complexvec_to_json(alpha_L_u),
            "alpha_L_l" => complexvec_to_json(alpha_L_l),
            "F" => complexvec_to_json(F),
            "omega_F" => complexvec_to_json(omega_F),
            "omega_F_at_L" => complexvec_to_json(omega_F),
            "omega_gap" => omega_i - maximum(imag.(omega_F)),
            "omega_bisect_status" => "initial",
            "omega_bisect_tries" => 0,
            "d_branch" => nothing,
            "d_contour" => nothing,
            # v4.5: the grid varies between entries now, so every entry says
            # which grid it was computed on.
            "n_L" => length(L),
            "adapt_level" => adapt_level,
            "h_local" => nothing,
            # v4.7 CHANGE 9b: the seed entry used to carry a SHORTER key list
            # than every later entry.  MATLAB's jsondecode builds a struct
            # array only when every object in the array has the same field
            # NAMES -- one odd entry and it returns a cell array of structs
            # instead, so every reader doing jsonData(k).field dies on the
            # first pass.  That is why the v4.5 log animates and the v4.6 and
            # v4.7 logs do not.  Differing value sizes and nulls are fine; it
            # is only the names that matter.  Keep these in step with the
            # per-iteration Dict below.
            "omega_pr_fit" => nothing,
            "omega2_fit"   => nothing,
            "pinch_miss"   => nothing,
            "grid_share"   => nothing,
            "zeta_alpha"   => nothing,
            "zeta_d_med"   => nothing,
            "exp_arg_max"  => nothing,
            "dist_u"       => nothing,
            "dist_l"       => nothing,
            "track_overrides" => nothing,
            "track_jump_u"    => nothing,
            "track_jump_l"    => nothing,
            # v4.9: keep this list in step with the per-iteration Dict, for the
            # jsondecode reason spelled out above.
            "n_limited"      => nothing,
            "alpha_accepted" => nothing,
            "alpha_attempt"  => nothing,
            "theta_move"     => nothing,
            "f_ripple"       => nothing,
            "peak_move"      => nothing,
            "descent_tol"    => nothing,
            "exp_arg_post"   => nothing,
            "exp_arg_bare"   => nothing,
            "dist_u_post"    => nothing,
            "dist_l_post"    => nothing,
            "epsilon_alpha"  => nothing,
            "omega_pr_clean" => nothing,
            "omega2_clean"   => nothing,
            "clean_spread"   => nothing,
            "clean_n"        => nothing,
            # v5: keep in step with the per-iteration Dict, for the jsondecode
            # reason spelled out above.
            "f_slope"     => nothing,
            "min_clear"   => nothing,
            "n_cross"     => nothing,
            "theta_peak"  => nothing,
            "desc_scale"  => nothing,
            "smooth_hw"   => nothing,
        )
        current_array = Any[]
        push!(current_array, dict_to_JSON)
        atomic_write(filename, JSON.json(current_array))
        iteration_step += 1
        write_checkpoint!()
        @printf("FRESH RUN: log %s, checkpoint %s, target iteration %d\n",
                filename, CKPT_FILE, ITER_TARGET)
        flush(stdout)
    end
end

function branch_distance(alpha_L_u, alpha_L_l)
    return minimum(abs.(alpha_L_u .- alpha_L_l))
end

function contour_L_at(omega_i_trial)
    return ComplexF64[
        omega_r[j] + omega_i_trial * im
        for j in eachindex(omega_r)
    ]
end

function branch_overlap_valid(alpha_u, alpha_l; overlap_tol = 1e-8)
    min_dist = Inf
    min_i = 0
    min_j = 0
    for i in eachindex(alpha_u)
        for j in eachindex(alpha_l)
            d = abs(alpha_u[i] - alpha_l[j])
            if d < min_dist
                min_dist = d
                min_i = i
                min_j = j
            end
        end
    end
    if min_dist < overlap_tol
        return false, @sprintf(
            "upper/lower overlap: d=%.3e at upper[%d], lower[%d]", min_dist, min_i, min_j)
    end
    return true, @sprintf("no overlap: min upper/lower distance=%.3e", min_dist)
end

function omega_trial_ok(omega_i_trial)
    L_try = contour_L_at(omega_i_trial)
    au_try, al_try = contour_alpha_L_conti(L_try)
    if !all(isfinite, au_try) || !all(isfinite, al_try)
        return false, "non-finite branch", L_try, au_try, al_try
    end
    ok, reason = branch_overlap_valid(au_try, al_try)
    return ok, reason, L_try, au_try, al_try
end

function local_contour_distances(F, alpha_L_u, alpha_L_l)
    all_branches = vcat(alpha_L_u, alpha_L_l)
    return [minimum(abs.(f .- all_branches)) for f in F]
end

function contour_distance(F, alpha_L_u, alpha_L_l)
    return minimum(local_contour_distances(F, alpha_L_u, alpha_L_l))
end

function nearest_branch_info(f, alpha_L_u, alpha_L_l)
    all_branches = vcat(alpha_L_u, alpha_L_l)
    distances = abs.(all_branches .- f)
    idx = argmin(distances)
    return all_branches[idx], distances[idx]
end

function directional_dt_check(F_old, F_trial, alpha_L_u, alpha_L_l;
                          move_safety = 0.5, global_max_move = 0.02)
    N = length(F_old)
    toward_amounts = zeros(Float64, N)
    local_allowed_toward = zeros(Float64, N)
    full_moves = abs.(F_trial .- F_old)
    dt_ok = true
    for j in 1:N
        f_old = F_old[j]
        f_new = F_trial[j]
        move_vec = f_new - f_old
        nearest_branch, dist = nearest_branch_info(f_old, alpha_L_u, alpha_L_l)
        if dist < 1e-12
            toward_amounts[j] = abs(move_vec)
            local_allowed_toward[j] = 0.0
            dt_ok = false
            continue
        end
        direction_to_branch = nearest_branch - f_old
        unit_to_branch = direction_to_branch / dist
        toward_amount = real(conj(unit_to_branch) * move_vec)
        allowed_toward = move_safety * dist
        toward_amounts[j] = toward_amount
        local_allowed_toward[j] = allowed_toward
        if toward_amount > allowed_toward
            dt_ok = false
        end
        if abs(move_vec) > global_max_move
            dt_ok = false
        end
    end
    return dt_ok, toward_amounts, local_allowed_toward, full_moves
end

function branch_slowdown_factor(d_branch)
    d_safe = 0.2      # start slowing below this
    min_factor = 0.1  # never go below  10% speed
    return clamp(d_branch / d_safe, min_factor, 1.0)
end
function print_iteration_header(k)
    println()
    println("====================================================================")
    println("                      ITERATION k = $k STARTED")
    println("====================================================================")
    flush(stdout)
end
function print_iteration_footer(k)
    println("--------------------------------------------------------------------")
    println("                      ITERATION k = $k FINISHED")
    println("--------------------------------------------------------------------")
    println()
    flush(stdout)
end
function print_block(title)
    println()
    println("---- $title ----")
    flush(stdout)
end
#
for k = 1:MAX_STEPS
    global omega_i, L, alpha_L_u, alpha_L_l, alpha_i, F, omega_F, iteration_step
    global omega_r, adapt_level, adapt_nogain, adapt_done
    local dict_to_JSON, current_array

    # v4.7 CHANGE 8: stop on the ABSOLUTE iteration count, not on how many
    # steps this process happens to have taken.  k still counts steps in this
    # segment; iteration_step is the number that goes into the log and that
    # every previous segment also counted in.
    if iteration_step > ITER_TARGET
        @printf("TARGET REACHED: next iteration would be %d, ITER_TARGET = %d.  Stopping after %d step(s) this run.\n",
                iteration_step, ITER_TARGET, k - 1)
        flush(stdout)
        break
    end

    omega_i_old = omega_i
    omega_status = "none"
    omega_attempt_used = 0
    omega_valid_count = 0
    omega_jump = NaN
    omega_dt_used = NaN

    alpha_status = "not-run"
    alpha_attempt_used = 0
    max_raw_move = NaN
    max_smooth_move = NaN

    begin
        # omega_F as it stands NOW, i.e. the one L is about to be placed
        # against.  Stored in the JSON as "omega_F_at_L".
        omega_F_at_L = copy(omega_F)

        # v4.6 CHANGE 2: parabolic peak instead of the discrete node maximum.
        # Falls back to maximum(imag.(omega_F)) whenever the fit is not a clean
        # interior peak, so OMEGA_LB_PARABOLIC = false reproduces v4.5 exactly.
        omegaF_max_i = omegaF_peak_imag(omega_F)
        omega_clearance = 1e-9
        omega_lower_bound = omegaF_max_i + omega_clearance

        # ------------------------------------------------------------
        # REPAIR STEP
        # ------------------------------------------------------------
        omega_repaired = false
        omega_accepted = false

        if omega_i <= omega_lower_bound
            global omega_i = omega_lower_bound
            global L = contour_L()
            load_on_workers()
            omega_repaired = true
            omega_accepted = true
            global delta_t = min(delta_t, 1.5e-3)
            omega_status = "repaired"
            omega_jump = abs(omega_i - omega_i_old)
            omega_dt_used = 0.0
            omega_valid_count = 0
        end

        omega_dt = delta_t
        omega_dt_min = 1e-12
        min_valid_omega_candidates = 2
        max_omega_jump_factor = 40.0
        min_useful_omega_jump = 1e-12

        if !omega_repaired
            for omega_attempt in 1:50
                # NOTE (v4.5, unchanged on purpose): omega_i_vectorization is
                # a constant vector and is never written inside this loop, so
                # the central difference below is identically zero and the
                # d_d_omega_r_Phi_L coupling contributes nothing.  That term
                # is the only place omega_r spacing enters the omega update,
                # so a non-uniform grid cannot perturb the descent.
                omega_i_vectorization = fill(omega_i, length(omega_r))
                omega_i_cache = copy(omega_i_vectorization)

                for j in 2:(length(omega_i_cache) - 1)
                    omega_i_cache[j] =
                        omega_i_vectorization[j] +
                        omega_dt * (
                            -lambda
                            + (omega_i_vectorization[j+1] - omega_i_vectorization[j-1]) /
                            (omega_r[j+1] - omega_r[j-1]) *
                            d_d_omega_r_Phi_L(L[j])
                            - d_d_omega_i_Phi_L(L[j])
                        )
                end

                omega_i_cache = [isfinite(x) && abs(x) < 10.0 ? x : -Inf for x in omega_i_cache]
                greater_candidates = filter(x -> isfinite(x) && x > omega_lower_bound &&  abs(x - omega_i) > min_useful_omega_jump, omega_i_cache)

                if length(greater_candidates) >= min_valid_omega_candidates
                    omega_candidate = minimum(greater_candidates)
                    omega_jump = abs(omega_candidate - omega_i)
                    max_allowed_omega_jump = max_omega_jump_factor * omega_dt
                    if omega_jump <= max_allowed_omega_jump
                        L_trial = contour_L_at(omega_candidate)
                        alpha_L_u_trial, alpha_L_l_trial = contour_alpha_L_conti(L_trial)
                        overlap_ok, overlap_reason =
                            branch_overlap_valid(alpha_L_u_trial, alpha_L_l_trial)
                        if overlap_ok
                            global omega_i = omega_candidate
                            global L = L_trial
                            global alpha_L_u = alpha_L_u_trial
                            global alpha_L_l = alpha_L_l_trial
                            omega_accepted = true
                            omega_status = "accepted"
                            omega_attempt_used = omega_attempt
                            omega_valid_count = length(greater_candidates)
                            omega_dt_used = omega_dt
                            break
                        else
                            @printf(
                                "[k=%d] omega rejected: branch tracking failed (%s) | reducing dt_omega %.2e -> %.2e\n",
                                k, overlap_reason, omega_dt, 0.5 * omega_dt)
                            omega_dt *= 0.5
                        end
                    else
                        @printf(
                            "[k=%d] omega rejected: jump %.3e > allowed %.3e | reducing dt_omega %.2e -> %.2e\n",
                            k, omega_jump, max_allowed_omega_jump, omega_dt, 0.5 * omega_dt)
                        omega_dt *= 0.5
                    end
                else
                    @printf(
                        "[k=%d] omega rejected: valid/useful=%d < required=%d | lower=%.6e | omega_i=%.6e | reducing dt_omega %.2e -> %.2e\n",
                        k, length(greater_candidates), min_valid_omega_candidates,
                        omega_lower_bound, omega_i, omega_dt, 0.5 * omega_dt)
                    omega_dt *= 0.5
                end

                if omega_dt < omega_dt_min
                    @printf("[k=%d] STOP: omega_dt below minimum. No safe omega step found.\n", k)
                    omega_accepted = false
                    break
                end
            end
        end

        if !omega_accepted
            println("No safe omega accepted.")
            break
        end
        if omega_repaired
            global L = contour_L()
            load_on_workers()
            @everywhere begin
                normals_F = contour_normals(F)
            end
            global alpha_L_u, alpha_L_l = contour_alpha_L_conti(L)
            load_on_workers()
        else
            load_on_workers()
            @everywhere begin
                normals_F = contour_normals(F)
            end
        end

        # -------------------------------------------------------------------
        # v4.4 BISECTION ENDGAME (unchanged)
        # -------------------------------------------------------------------
        omega_bisect_status = OMEGA_BISECT ? "idle" : "off"
        omega_bisect_tries  = 0
        omega_gap_before    = omega_i - omega_lower_bound
        omega_gap_after     = omega_gap_before
        omega_wall_reason   = ""

        if OMEGA_BISECT && omega_accepted &&
           omega_gap_before > OMEGA_GAP_TARGET &&
           omega_gap_before < OMEGA_BISECT_ENGAGE

            a_lo = omega_lower_bound     # deepest legal height, validity unknown
            b_hi = omega_i               # current height, known admissible

            best_w = omega_i
            best_L = L
            best_u = alpha_L_u
            best_l = alpha_L_l

            for probe in 1:OMEGA_BISECT_MAX
                omega_bisect_tries += 1
                w_try = (probe == 1) ? a_lo : 0.5 * (a_lo + b_hi)
                ok, reason, L_try, au_try, al_try = omega_trial_ok(w_try)
                if ok && OMEGA_BISECT_MIN_DFC > 0.0
                    dfc_try = contour_distance(F, au_try, al_try)
                    if !isfinite(dfc_try) || dfc_try < OMEGA_BISECT_MIN_DFC
                        ok = false
                        reason = @sprintf("F-branch distance %.3e < %.3e",
                                          dfc_try, OMEGA_BISECT_MIN_DFC)
                    end
                end
                if ok
                    best_w, best_L, best_u, best_l = w_try, L_try, au_try, al_try
                    b_hi = w_try
                else
                    a_lo = w_try
                    omega_wall_reason = reason
                end
                if (best_w - omega_lower_bound) <= OMEGA_GAP_TARGET
                    omega_bisect_status = "target"
                    break
                end
                if (b_hi - a_lo) < OMEGA_BISECT_TOL
                    omega_bisect_status = "wall"
                    break
                end
                omega_bisect_status = "budget"
            end

            if best_w < omega_i
                global omega_i   = best_w
                global L         = best_L
                global alpha_L_u = best_u
                global alpha_L_l = best_l
                load_on_workers()
                omega_jump = abs(omega_i - omega_i_old)
                omega_status = omega_status * "+bisect"
            end

            omega_gap_after = omega_i - omega_lower_bound

            if omega_bisect_status != "target"
                @printf(
                    "[k=%d] bisect %s after %d probes: gap %.3e -> %.3e | dFC=%.3e | wall: %s\n",
                    k, omega_bisect_status, omega_bisect_tries,
                    omega_gap_before, omega_gap_after,
                    contour_distance(F, alpha_L_u, alpha_L_l),
                    isempty(omega_wall_reason) ? "-" : omega_wall_reason)
                flush(stdout)
            end
        end
        # v4.9: the stray `alpha_i_cache = copy(alpha_i)` that stood here is
        # gone with the variable -- it was assigned in two places and read in
        # none, which is what let the v4.7 attempt loop apply a rejected step.

        d_vec     = abs.(alpha_L_u .- alpha_L_l)
        i_pinch   = argmin(d_vec)
        d_branch  = d_vec[i_pinch]
        h_local   = local_h(omega_r, i_pinch)
        d_contour = contour_distance(F, alpha_L_u, alpha_L_l)

        # v4.6: pinch location from the d^4 parabola, and the share of d_branch
        # that the omega_r miss accounts for.  Computed here so it lands in the
        # JSON for every iteration; the adaptive block at the end of the
        # iteration reuses these, it does not refit.
        #   grid_share -> 1   : the grid is what limits d_branch
        #   grid_share << 1   : omega_i is still above the pinch
        w_pr_fit, w2_fit = pinch_fit(omega_r, d_vec, i_pinch)
        pinch_miss = abs(omega_r[i_pinch] - w_pr_fit)
        d_horiz    = (isfinite(pinch_miss) && isfinite(w2_fit) && w2_fit > 0) ?
                     2 * sqrt(2 * pinch_miss / w2_fit) : NaN
        grid_share = (isfinite(d_horiz) && d_branch > 0) ? d_horiz / d_branch : NaN
        dist_u = minimum(minimum(abs.(f .- alpha_L_u)) for f in F)
        dist_l = minimum(minimum(abs.(f .- alpha_L_l)) for f in F)

        # ---------------------------------------------------------------
        # v4.7 CHANGE 5: retune the repulsion range to the gap F now has to
        # thread.  Placed here, after d_branch and dist_u/dist_l exist and
        # before the alpha update that is the only consumer of zeta_alpha.
        # exp_arg is the largest exponent phi_F will see this iteration --
        # the direct saturation monitor.  Watch it stay well under
        # EXP_ARG_MAX = 400; in v4.6 it went past 400 on 26 of the last 156
        # iterations.
        # ---------------------------------------------------------------
        zeta_push_d!(d_branch)
        global zeta_alpha = zeta_step(zeta_alpha)
        # v4.9 CHANGE 10: epsilon_alpha is a function of zeta_alpha and must be
        # refreshed with it, here, before the alpha update that is the only
        # consumer of either.
        global epsilon_alpha = max(zeta_alpha / EXP_ARG_TARGET, 1e-300)
        r_min   = min(dist_u, dist_l)
        # v4.9: deliberately the BARE 1e-10 softening, not epsilon_alpha, so the
        # console `xarg` and the logged exp_arg_max stay on the same scale as
        # the v4.6/v4.7/v4.8 history.  It is still the PRE-update value -- the
        # post-update pair (exp_arg_post, exp_arg_bare) is logged after the F
        # step, and those are the ones to read.
        exp_arg = zeta_alpha / (r_min^2 + epsilon)

        branch_factor = branch_slowdown_factor(d_branch)
        local_delta_t = delta_t * branch_factor

        dt_min = 1e-10
        dt_max = 2e-2
        dt_start = min(delta_t, dt_max)

        if d_contour > 0.01
            move_safety = 2.0
        elseif d_contour > 0.003
            move_safety = 1.2
        else
            move_safety = 0.8
        end

        global_max_move = 0.20
        stop_after_save = false
        stop_reason = ""

        if !isfinite(d_branch) || !isfinite(d_contour)
            @printf("[k=%d] STOP: non-finite distance | d_branch=%.3e | d_contour=%.3e\n",
                    k, d_branch, d_contour)
            break
        end

        if d_branch < PINCH_TOL
            stop_after_save = true
            stop_reason = "pinch tolerance reached"
        end

    # -------------------------------------------------------------------
    # v4.9 CHANGE 11 + 12.  The F step.
    #
    # CHANGE 11 replaces clamp(rhs, +-rhs_cap) with a smooth limit on the
    # DISPLACEMENT.  theta is the largest move any node may take this
    # iteration; tanh maps the raw step into (-theta, theta) monotonically, so
    # a node pushed at 1e177 and its neighbour pushed at 1e5 no longer receive
    # the same magnitude with opposite signs.  Set LIMITER_ON = false for the
    # v4.7 clamp.
    #
    # CHANGE 12 makes the attempt loop a line search on max Im(omega_F), the
    # quantity the descent is actually minimising.  The best trial seen is the
    # one applied, so the peak can only go down.  The v4.7 loop tested the size
    # of the move, wrote alpha_i_cache and `accepted` and read neither, and
    # applied the last attempt pass or fail.  Set DESCENT_ON = false to apply
    # attempt 1 unconditionally, i.e. the v4.7 behaviour with the new limiter.
    #
    # omega_F for the accepted trial is computed here, so the
    # contour_omega_F(F) that used to follow this block is gone -- an accepted
    # first attempt costs exactly what v4.7 cost.  Only rejections cost extra.
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    # v5 THE F STEP
    #   barrier (v4.7 CHANGE 10 softening)  +  descent on Im(omega)  [CHANGE 16]
    #   per-node trust region               [CHANGE 18]
    #   implicit sigma                      [CHANGE 15]
    #   width-based smoothing               [CHANGE 14]
    #   hard non-crossing                   [CHANGE 17]
    #   line search on max Im(omega_F)      [v4.9 CHANGE 12, unchanged]
    # -------------------------------------------------------------------
    peak_old    = omegaF_peak_imag(omega_F)
    descent_tol = descent_tolerance()

    # CHANGE 18: each node is limited by ITS OWN clearance, not by the global
    # gap.  At the throat that is far tighter than v4.9 (which let a node with
    # 2.4e-4 of room move 1.8e-4); in the wings, where the rotation has to
    # happen, it is a few hundred times looser.
    dist_node_u = [minimum(abs.(f .- alpha_L_u)) for f in F]
    dist_node_l = [minimum(abs.(f .- alpha_L_l)) for f in F]
    r_node = min.(dist_node_u, dist_node_l)
    theta  = (isfinite(d_branch) && d_branch > 0) ?
             min(MOVE_FRAC * d_branch, global_max_move) : global_max_move
    # CHANGE 18b (correction, measured on contour_iteration_v5.0.json, 100 frames).
    # The first cut of CHANGE 18 used theta_j = MOVE_FRAC * r_j with NO floor, and
    # that is wrong at the start of a run for a structural reason: F begins AS the
    # real axis and the alpha+ / alpha- roots begin sitting right on it, so r_j is
    # ~1e-5 for every node by construction.  Measured at iteration 4: v4.9 allowed
    # theta = 1.670e-1, v5.0 allowed 8.355e-7 -- a factor 199,839 -- and the run
    # moved d from 3.32 to 3.29 in 100 iterations where v4.9 reached 0.997.
    # The floor below makes theta_j >= the v4.9 global value everywhere, so the
    # per-node rule can only ever give a node MORE room, never less.  Safety at
    # the throat is CHANGE 17's job (reject any trial that crosses), not theta's;
    # that is what theta was being asked to do here, and it could not.
    theta_j = THETA_PER_NODE ?
              [clamp(MOVE_FRAC * rj, theta, global_max_move) for rj in r_node] :
              fill(theta, length(F))

    # CHANGE 16: the descent direction.  dv/dy = Re(domega/dalpha) exactly, by
    # Cauchy-Riemann, and domega/dalpha comes from F and omega_F alone.
    dwda       = DESC_ON ? domega_dalpha(F, omega_F) : ComplexF64[]
    wdesc      = DESC_ON ? descent_weights(omega_F)  : Float64[]
    desc_raw   = DESC_ON ? [wdesc[j] * real(dwda[j]) for j in eachindex(F)] : Float64[]
    desc_scale = (DESC_ON && !isempty(desc_raw)) ? maximum(abs, desc_raw) : 0.0
    dt0        = local_delta_t

    # CHANGE 17: where the present geometry stands.  A v4.9 state may already
    # be crossed, so the acceptance rule below also lets a trial through on
    # the grounds that it is LESS crossed than what we have.
    # v5.  CHANGE 17b (gbound = max(0.5, 20 d_branch)) is REVERTED.  It was
    # wrong, and the way it was checked was wrong too: I validated it on the v5.0
    # frames, but that run was frozen by the CHANGE 18 theta bug, so F never left
    # the real axis and its flat extension happened to sit above the far-field
    # branch points.  The moment F actually descended (v5.1) the extension dropped
    # below them and the widened bound dragged ~100 far-field pairs into the test.
    # The narrow bound is correct: it keeps branch-switch garbage (gap ~ 1.9) out,
    # and the pinch-region pairs it is meant to police all have gap << 0.5 by the
    # time there is anything to police.  The span restriction added in
    # branch_clearance is the change that actually makes the test meaningful.
    gbound = 0.5
    clear_now, ncross_now = CROSS_CHECK ?
                            branch_clearance(F, alpha_L_u, alpha_L_l; gap_bound = gbound) :
                            (Inf, 0)
    cross_veto = CROSS_CHECK && (discard_run < CROSS_ESCAPE_AFTER)
    if CROSS_CHECK && !cross_veto
        @printf("[k=%d] CROSS VETO SUSPENDED after %d discarded iterations (clear=%.3e, n_cross=%d)\n",
                k, discard_run, clear_now, ncross_now)
        flush(stdout)
    end
    clear_req = CLEAR_FRAC * (isfinite(d_branch) ? d_branch : 0.0)

    accepted     = false
    best_peak    = Inf
    best_alpha   = Float64[]
    best_omega   = ComplexF64[]
    best_attempt = 0
    n_limited    = 0
    min_clear    = clear_now
    n_cross      = ncross_now

    for attempt in 1:ALPHA_MAX_ATTEMPTS
        alpha_i_trial = copy(alpha_i)
        n_lim_try = 0
        for j in 2:(length(alpha_i_trial) - 1)
            hp = alpha_r[j+1] - alpha_r[j]
            hm = alpha_r[j]   - alpha_r[j-1]
            # CHANGE 13.  Both reduce exactly to the v4.9 expressions when
            # hp == hm, which is every node of the uniform grid this version
            # still runs -- so this edit changes nothing today and is what
            # makes a graded F survivable later.
            alpha_i_r =
                (-hp / (hm * (hp + hm))) * alpha_i[j-1] +
                ((hp - hm) / (hp * hm))  * alpha_i[j]   +
                ( hm / (hp * (hp + hm))) * alpha_i[j+1]
            alpha_i_rr =
                2.0 * ((alpha_i[j+1] - alpha_i[j]) / hp -
                       (alpha_i[j] - alpha_i[j-1]) / hm) / (hp + hm)
            sig_expl = DIFFUSION_IMPLICIT ? 0.0 : sigma
            rhs_j =
                (
                    alpha_i_r * d_d_alpha_r_Phi_F(F[j])
                    - d_d_alpha_i_Phi_F(F[j])
                    + sig_expl * alpha_i_rr
                ) / (1.0 + alpha_i_r^2)
            # CHANGE 16.  Sign: moving y by dy changes Im(omega) by
            # Re(domega/dalpha) * dy, so the descent direction is minus that.
            # Scaled so the most strongly driven node spends DESC_FRAC of its
            # own trust region on descent and leaves the rest to the barrier.
            if DESC_ON && desc_scale > 0
                rhs_j -= DESC_FRAC * (theta_j[j] / dt0) * desc_raw[j] / desc_scale
            end
            rhs_j = isnan(rhs_j) ? 0.0 : rhs_j
            th = theta_j[j]
            if LIMITER_ON
                u_lim = local_delta_t * rhs_j / th
                abs(u_lim) > LIMITER_SAT_AT && (n_lim_try += 1)
                alpha_i_trial[j] = alpha_i[j] + th * tanh(u_lim)
            else
                rc = clamp(rhs_j, -RHS_CAP, RHS_CAP)
                abs(rhs_j) >= RHS_CAP && (n_lim_try += 1)
                alpha_i_trial[j] = alpha_i[j] + local_delta_t * rc
            end
        end
        alpha_i_trial[1] = alpha_i_trial[2]
        alpha_i_trial[end] = alpha_i_trial[end-1]

        # CHANGE 15.  The implicit half of the split.  Unconditionally stable,
        # so sigma no longer constrains dt and the filter is free to stand down.
        if DIFFUSION_IMPLICIT
            metric = similar(alpha_i_trial)
            @inbounds for j in eachindex(alpha_i_trial)
                jm = max(j - 1, 1); jp = min(j + 1, length(alpha_i_trial))
                dx = alpha_r[jp] - alpha_r[jm]
                sl = dx > 0 ? (alpha_i_trial[jp] - alpha_i_trial[jm]) / dx : 0.0
                metric[j] = 1.0 / (1.0 + sl^2)
            end
            alpha_i_trial = diffusion_implicit(collect(Float64, alpha_r),
                                               alpha_i_trial, metric,
                                               local_delta_t, sigma)
        end

        # CHANGE 14.  Width in alpha_r, not width in nodes.
        alpha_i_try =
            if SMOOTH_MODE === :index
                rolling_average_filter(alpha_i_trial, SMOOTH_RADIUS)
            elseif SMOOTH_MODE === :width
                hw = isfinite(d_branch) ?
                     min(SMOOTH_FRAC * d_branch, SMOOTH_W_MAX) : SMOOTH_W_MAX
                width_average(collect(Float64, alpha_r), alpha_i_trial, hw)
            else
                copy(alpha_i_trial)
            end
        alpha_i_try[1] = alpha_i_try[2]
        alpha_i_try[end] = alpha_i_try[end-1]

        max_raw_move    = maximum(abs.(alpha_i_trial .- alpha_i))
        max_smooth_move = maximum(abs.(alpha_i_try   .- alpha_i))
        alpha_attempt_used = attempt

        if !all(isfinite, alpha_i_try)
            local_delta_t *= 0.5
            alpha_status = "nonfinite/reduced"
            continue
        end

        F_try = ComplexF64[alpha_r[j] + alpha_i_try[j] * im for j in 1:N]

        # CHANGE 17.  No pole may sit on the wrong side of F.
        c_try = Inf; nx_try = 0
        if CROSS_CHECK
            c_try, nx_try = branch_clearance(F_try, alpha_L_u, alpha_L_l; gap_bound = gbound)
            if cross_veto && !(c_try >= clear_req || c_try > clear_now)
                local_delta_t *= 0.5
                alpha_status = "crossing/reduced"
                continue
            end
        end

        omega_try = contour_omega_F(F_try)
        peak_try  = omegaF_peak_imag(omega_try)

        if isfinite(peak_try) && peak_try < best_peak
            best_peak    = peak_try
            best_alpha   = copy(alpha_i_try)
            best_omega   = copy(omega_try)
            best_attempt = attempt
            n_limited    = n_lim_try
            min_clear    = c_try
            n_cross      = nx_try
        end

        if !DESCENT_ON || !isfinite(peak_old) || !isfinite(descent_tol) ||
           (isfinite(peak_try) && peak_try <= peak_old + descent_tol)
            accepted = true
            alpha_status = "accepted"
            break
        else
            local_delta_t *= 0.5
            alpha_status = "descent/reduced"
        end
    end

        # Always take the BEST trial seen, never a rejected one.  best_attempt
        # is 0 only if every attempt was non-finite or crossing, in which case
        # the geometry is left alone and the next iteration retries with the
        # halved local_delta_t.
        if best_attempt == 0
            global discard_run += 1
            @printf("[k=%d] WARNING: no admissible alpha trial in %d attempts; step discarded (run of %d)\n",
                    k, ALPHA_MAX_ATTEMPTS, discard_run)
            flush(stdout)
            alpha_status = "discarded"
        else
            global discard_run = 0
            global alpha_i = copy(best_alpha)
            global omega_F = copy(best_omega)
            peak_push!(best_peak)
        end
        peak_move = (best_attempt == 0 || !isfinite(peak_old)) ? NaN : best_peak - peak_old
        global F = contour_F()
        load_on_workers()

        # ---- diagnostics, all on the POST-update geometry -----------------
        # alpha_L_u/l here are the pre-update branches, the convention
        # dist_u/dist_l already use.
        dist_u_post = minimum(minimum(abs.(f .- alpha_L_u)) for f in F)
        dist_l_post = minimum(minimum(abs.(f .- alpha_L_l)) for f in F)
        r_post       = min(dist_u_post, dist_l_post)
        exp_arg_post = zeta_alpha / (r_post^2 + epsilon_alpha)
        exp_arg_bare = zeta_alpha / (r_post^2 + 1e-10)
        f_ripple     = f_chord_ripple(F)
        if CROSS_CHECK && best_attempt != 0
            min_clear, n_cross = branch_clearance(F, alpha_L_u, alpha_L_l; gap_bound = gbound)
        end
        # v5: F's slope where its omega-image peaks.  This is the number
        # CHANGE 16 exists to move; watch it leave -0.00x.  Nothing reads it.
        jpk = argmax(imag.(omega_F))
        jlo = max(jpk - 2, firstindex(F)); jhi = min(jpk + 2, lastindex(F))
        f_slope = (jhi > jlo && real(F[jhi]) != real(F[jlo])) ?
                  (imag(F[jhi]) - imag(F[jlo])) / (real(F[jhi]) - real(F[jlo])) : NaN
        smooth_hw = SMOOTH_MODE === :width ?
                    (isfinite(d_branch) ? min(SMOOTH_FRAC * d_branch, SMOOTH_W_MAX) :
                                          SMOOTH_W_MAX) :
                    (SMOOTH_MODE === :index ?
                     SMOOTH_RADIUS * (alpha_r[2] - alpha_r[1]) : 0.0)
        theta_pk = (jpk >= 1 && jpk <= length(theta_j)) ? theta_j[jpk] : theta
        w_pr_clean, w2_clean, clean_spread, clean_n = clean_pinch_fit(omega_r, d_vec)
        @printf(
            "[%04d] jump=%9.3e | dUL=%9.3e | dUF=%9.3e | dLF=%9.3e | dt=%9.3e | z=%9.3e xarg=%8.1f | gap=%9.3e (%s,%d) | NL=%3d L%d h=%8.2e gs=%5.2f miss=%8.2e | jmpU=%6.3f ovr=%2d | lim=%3d rip=%8.2e dpk=%+9.2e | slp=%+8.4f clr=%+9.2e nx=%2d | %s/%s\n",
            k, omega_jump, d_branch, dist_u, dist_l, local_delta_t, zeta_alpha, exp_arg,
            omega_gap_after, omega_bisect_status, omega_bisect_tries,
            length(omega_r), adapt_level, h_local, grid_share, pinch_miss,
            track_jump_u, track_overrides,
            n_limited, f_ripple, peak_move,
            f_slope, min_clear, n_cross,
            omega_status, alpha_status)
        flush(stdout)
        load_on_workers()
        local dict_to_JSON = Dict(
            "iteration" => iteration_step,
            "L" => complexvec_to_json(L),
            "alpha_L_u" => complexvec_to_json(alpha_L_u),
            "alpha_L_l" => complexvec_to_json(alpha_L_l),
            "F" => complexvec_to_json(F),
            "omega_F" => complexvec_to_json(omega_F),
            "omega_F_at_L" => complexvec_to_json(omega_F_at_L),
            "omega_gap" => omega_gap_after,
            "omega_bisect_status" => omega_bisect_status,
            "omega_bisect_tries" => omega_bisect_tries,
            "d_branch" => d_branch,
            "d_contour" => d_contour,
            "n_L" => length(L),
            "adapt_level" => adapt_level,
            "h_local" => h_local,
            # v4.6 diagnostics
            "omega_pr_fit" => jsonnum(w_pr_fit),
            "omega2_fit"   => jsonnum(w2_fit),
            "pinch_miss"   => jsonnum(pinch_miss),
            "grid_share"   => jsonnum(grid_share),
            # v4.7 diagnostics
            "zeta_alpha"   => jsonnum(zeta_alpha),
            "zeta_d_med"   => jsonnum(length(zeta_hist_d) < ZETA_WINDOW ?
                                      NaN : median(zeta_hist_d)),
            "exp_arg_max"  => jsonnum(exp_arg),
            "dist_u"       => jsonnum(dist_u),
            "dist_l"       => jsonnum(dist_l),
            # v4.7 CHANGE 9.  track_overrides: how many times the continuity
            # guard overrode the F-normal side test this sweep (7 of 282 on
            # the replayed iteration 501).  track_jump_u/l: the largest step
            # along each branch between neighbouring omega samples -- the
            # direct measure of the branch-switching artefact.  Watch
            # track_jump_u stay in the 0.0x range; ~1.9 means it switched.
            "track_overrides" => track_overrides,
            "track_jump_u"    => jsonnum(track_jump_u),
            "track_jump_l"    => jsonnum(track_jump_l),
            # v4.9 diagnostics.  n_limited is the direct successor of "how many
            # nodes hit rhs_cap" (27 of 146 on v4.8 iteration 501) and is the
            # first thing to check.  f_ripple is the max chord deviation of F,
            # which must fall as h^2 under any honest refinement.  peak_move is
            # the step taken by max Im(omega_F) this iteration: under CHANGE 12
            # it is never positive.  The clean_* fields are the sqrt-law fit on
            # the resolved part of the d profile, logged for comparison against
            # omega_pr_fit / omega2_fit -- nothing acts on them in v4.9.
            "n_limited"      => n_limited,
            "alpha_accepted" => accepted,
            "alpha_attempt"  => alpha_attempt_used,
            "theta_move"    => jsonnum(theta),
            "f_ripple"      => jsonnum(f_ripple),
            "peak_move"     => jsonnum(peak_move),
            "descent_tol"   => jsonnum(descent_tol),
            "exp_arg_post"  => jsonnum(exp_arg_post),
            "exp_arg_bare"  => jsonnum(exp_arg_bare),
            "dist_u_post"   => jsonnum(dist_u_post),
            "dist_l_post"   => jsonnum(dist_l_post),
            "epsilon_alpha" => jsonnum(epsilon_alpha),
            "omega_pr_clean" => jsonnum(w_pr_clean),
            "omega2_clean"   => jsonnum(w2_clean),
            "clean_spread"   => jsonnum(clean_spread),
            "clean_n"        => clean_n,
            # v5.  f_slope is the number CHANGE 16 exists to move -- F's slope
            # where its omega-image peaks; v4.9 ended at -0.0019 and the
            # geometry wants O(1).  min_clear and n_cross are CHANGE 17: the
            # worst signed clearance between F and the branches, and how many
            # points are on the wrong side.  min_clear MUST stay positive; it
            # was negative on 83 % of the last 200 v4.9 frames.
            "f_slope"     => jsonnum(f_slope),
            "min_clear"   => jsonnum(min_clear),
            "n_cross"     => n_cross,
            "theta_peak"  => jsonnum(theta_pk),
            "desc_scale"  => jsonnum(desc_scale),
            "smooth_hw"   => jsonnum(smooth_hw),
        )
        # v4.7 CHANGE 8: stamp the first entry after a resume, so a join is
        # visible in the data instead of having to be remembered.  Any
        # discontinuity at a join is then attributable rather than mysterious.
        if iteration_step == resume_join_at
            dict_to_JSON["resumed"] = true
        end
        json_str = open(filename, "r") do file
            read(file, String)
        end
        local current_array = JSON.parse(json_str)
        push!(current_array, dict_to_JSON)
        atomic_write(filename, JSON.json(current_array))
        global iteration_step += 1

        if stop_after_save
            write_checkpoint!()
            @printf(
                "[k=%d] STOP AFTER SAVE: %s | d_branch=%.3e | d_contour=%.3e | omega_i=%.6e\n",
                k, stop_reason, d_branch, d_contour, omega_i)
            flush(stdout)
            break
        end

        # -------------------------------------------------------------------
        # v4.5 ADAPTIVE L REFINEMENT
        #
        # Placed at the very end of the iteration, AFTER the JSON save, so the
        # stored entry always describes the grid its numbers were computed on
        # and the next iteration simply starts on the finer grid.
        # -------------------------------------------------------------------
        if ADAPT_ON && !adapt_done
            # v5: a discarded step means the geometry did not move, so feeding
            # it to the stall detector manufactures a fake stall and spends a
            # refinement level on nothing.
            if alpha_status != "discarded"
                adapt_push!(i_pinch, d_branch, omega_r[i_pinch], h_local)
            end

            if adapt_stalled()
                if adapt_level >= ADAPT_MAX_LEVEL
                    @printf("[k=%d] ADAPT OFF: level cap %d reached | d_branch=%.4e\n",
                            k, ADAPT_MAX_LEVEL, d_branch)
                    flush(stdout)
                    global adapt_done = true
                elseif d_branch <= ADAPT_D_FLOOR
                    @printf("[k=%d] ADAPT OFF: d_branch=%.4e is at or below the noise floor %.1e\n",
                            k, d_branch, ADAPT_D_FLOOR)
                    flush(stdout)
                    global adapt_done = true
                elseif i_pinch <= 1 || i_pinch >= length(omega_r)
                    @printf("[k=%d] adapt: argmin is at a grid endpoint (i=%d) -- not refining\n",
                            k, i_pinch)
                    flush(stdout)
                    adapt_reset_history!()
                else
                    h_new = h_local / ADAPT_FACTOR
                    n_add = 2 * ADAPT_HALF_CELLS * (ADAPT_FACTOR - 1)
                    if h_new < ADAPT_MIN_H
                        @printf("[k=%d] ADAPT OFF: spacing floor reached (h=%.2e)\n", k, h_local)
                        flush(stdout)
                        global adapt_done = true
                    elseif length(omega_r) + n_add > ADAPT_MAX_POINTS
                        @printf("[k=%d] ADAPT OFF: point budget reached (N_L=%d)\n", k, length(omega_r))
                        flush(stdout)
                        global adapt_done = true
                    else
                        w2   = omega2_estimate(omega_r, d_vec, i_pinch)
                        dpre = adapt_predicted_d(h_new, w2)
                        n_before = length(omega_r)
                        w_lo = omega_r[max(1, i_pinch - ADAPT_HALF_CELLS)]
                        w_hi = omega_r[min(n_before, i_pinch + ADAPT_HALF_CELLS)]

                        # ----------------------------------------------------
                        # v4.6 CHANGES 1 and 4.  w_pr_fit / pinch_miss /
                        # grid_share were computed with d_vec earlier in this
                        # iteration and stored in the JSON entry above.
                        #   (a) grid_share says whether this round is even
                        #       about the grid  -> CHANGE 1, the score gate.
                        #   (b) w_pr_fit says where the new nodes belong
                        #                       -> CHANGE 4, vertex placement.
                        # ----------------------------------------------------
                        grid_limited = !ADAPT_SCORE_GATE ||
                                       (isfinite(grid_share) &&
                                        grid_share >= ADAPT_GRID_SHARE &&
                                        grid_share <= ADAPT_GRID_SHARE_MAX)
                        h_cluster = max(h_new, 2 * pinch_miss)
                        vertex_ok = ADAPT_VERTEX_PLACE && isfinite(w_pr_fit) &&
                                    isfinite(pinch_miss) &&
                                    pinch_miss <= ADAPT_VERTEX_MAXMOVE * h_local

                        placement = "subdivide"
                        h_placed  = h_new
                        if vertex_ok
                            omega_r_try = insert_cluster(omega_r, w_pr_fit, h_cluster)
                            if length(omega_r_try) > n_before
                                global omega_r = omega_r_try
                                placement = "cluster"
                                h_placed  = h_cluster
                            else
                                # every new node collided with an existing one
                                global omega_r = refine_omega_r(omega_r, i_pinch)
                            end
                        else
                            global omega_r = refine_omega_r(omega_r, i_pinch)
                        end
                        global adapt_level += 1
                        push_omega_r()
                        global L = contour_L()
                        load_on_workers()
                        # F was replaced a few lines above (contour_F after the
                        # alpha update) but normals_F still belongs to the OLD
                        # F -- in v4.4 nothing re-tracked after that point so it
                        # never mattered.  spatial_payload takes both from the
                        # main process, so they have to agree before the
                        # re-track, or the side classification uses one F's
                        # points with another F's normals.
                        @everywhere begin
                            normals_F = contour_normals(F)
                        end
                        global alpha_L_u, alpha_L_l = contour_alpha_L_conti(L)
                        load_on_workers()
                        adapt_reset_history!()

                        d_after = branch_distance(alpha_L_u, alpha_L_l)
                        gain    = d_branch > 0 ? d_after / d_branch : NaN
                        @printf(
                            "[k=%d] ADAPT -> level %d (%s) | window [%.9f, %.9f] | h %.3e -> %.3e | N_L %d -> %d | |w''| arms~%.4f fit~%.4f | vertex %.9f (miss %.3e, %.0f %% of d) | d %.4e -> %.4e (gain %.3f, predicted %.4e)\n",
                            k, adapt_level, placement, w_lo, w_hi, h_local, h_placed,
                            n_before, length(omega_r), w2, w2_fit, w_pr_fit, pinch_miss,
                            100 * grid_share, d_branch, d_after, gain, dpre)
                        flush(stdout)

                        # A round that buys nothing means the grid is no longer
                        # what limits d_min -- from here it is the omega_i error
                        # (i.e. F's resolution near alpha_p) or the noise floor.
                        # Two in a row, because one can be bad luck in where the
                        # new points land.
                        #
                        # v4.6 CHANGE 1: only score the round when the grid is
                        # what is being judged.  While omega_i is still above
                        # the pinch, d_branch is set by the vertical offset and
                        # the gain is capped near 1 whatever the grid does --
                        # that is what killed v4.5 at k = 227.
                        if !grid_limited
                            gate_reason = !isfinite(grid_share) ? "the d^4 fit failed" :
                                          grid_share > ADAPT_GRID_SHARE_MAX ?
                                             "the d^4 fit has degenerated (|w''| unusable)" :
                                             "omega_i is still the limiter"
                            @printf("[k=%d] adapt: round NOT scored -- horizontal miss explains %.0f %% of d_branch, outside [%.0f, %.0f] %%: %s | adapt_nogain stays %d\n",
                                    k, 100 * grid_share, 100 * ADAPT_GRID_SHARE,
                                    100 * ADAPT_GRID_SHARE_MAX, gate_reason, adapt_nogain)
                            flush(stdout)
                        elseif isfinite(gain) && gain > ADAPT_MIN_GAIN
                            global adapt_nogain += 1
                        else
                            global adapt_nogain = 0
                        end
                        if adapt_nogain >= ADAPT_NOGAIN_MAX
                            @printf("[k=%d] ADAPT OFF: %d refinements in a row gained less than %.0f%% -- the L grid is no longer the limiter (look at omega_i / F resolution / zeta next) | d_branch=%.4e\n",
                                    k, adapt_nogain, 100 * (1 - ADAPT_MIN_GAIN), d_after)
                            flush(stdout)
                            global adapt_done = true
                        end
                    end
                end
            end
        end

        # -------------------------------------------------------------------
        # v4.7 CHANGE 8: CHECKPOINT.  The last statement of the iteration, on
        # purpose.  The JSON entry above was written BEFORE the refinement
        # block, so it still describes the pre-refinement grid; this describes
        # the state iteration `iteration_step` will actually start from.
        # Resuming from the log instead would silently throw away the last
        # refinement round and then be unable to re-earn it for
        # ADAPT_STALL_WINDOW iterations.
        # -------------------------------------------------------------------
        if CKPT_EVERY > 0 && k % CKPT_EVERY == 0
            write_checkpoint!()
        end
    end
end

##############
# READING IN #
##############
begin
    function json_to_complexvec(arr)
        return ComplexF64[complex(x["re"], x["im"]) for x in arr]
    end
    function load_step(filename; offset=0)
        json_str = open(filename, "r") do file
            read(file, String)
        end
        data = JSON.parse(json_str)
        n = length(data)
        idx = n + offset
        entry = data[idx]
        iteration_step = entry["iteration"]
        L = json_to_complexvec(entry["L"])
        alpha_L_u = json_to_complexvec(entry["alpha_L_u"])
        alpha_L_l = json_to_complexvec(entry["alpha_L_l"])
        F = json_to_complexvec(entry["F"])
        omega_F = json_to_complexvec(entry["omega_F"])
        return iteration_step, L, alpha_L_u, alpha_L_l, F, omega_F
    end
end
# v4.7 CHANGE 8: this block runs AFTER the loop and OVERWRITES the live
# globals with the last logged entry.  That is harmless at the end of a batch
# run, but it would make a manual write_checkpoint!() from the REPL store a
# state one refinement behind the one the run actually ended on.  Off by
# default; set POSTLOAD_AFTER_RUN = true to get the old interactive behaviour.
const POSTLOAD_AFTER_RUN = false
if POSTLOAD_AFTER_RUN && isfile(filename)
    iteration_step, L, alpha_L_u, alpha_L_l, F, omega_F = load_step(filename; offset=0)
    println("Loaded iteration: ", iteration_step)
    global omega_i = imag(L[1])
    global alpha_i = imag.(F)
    # v4.5: the grid comes back from the stored L, not from the constant.
    global omega_r = real.(L)
    push_omega_r()
    load_on_workers()
end
#plot_omega()
#plot_alpha()

# ---------------------------------------------------------------------------
# v4.7 CHANGE 8: EQUIVALENCE TEST for the checkpoint.
#
# A resumed run is only the same object as a continuous one if the checkpoint
# carries ALL the state.  The way to know rather than hope:
#
#   1.  RESUME = false, ITER_TARGET = 60.  Run.  Keep the log as run_A.json.
#   2.  RESUME = false, ITER_TARGET = 30.  Run.
#   3.  RESUME = true,  ITER_TARGET = 60.  Run again -- it resumes at 31.
#       Keep the log as run_B.json.
#   4.  compare_runs("run_A.json", "run_B.json")
#
# A complete checkpoint gives max|A-B| = 0.0 on every field, because nothing in
# this code is stochastic and the arithmetic is identical.  Anything non-zero
# names the first iteration that diverges and the field that did it, which is
# the piece of state the checkpoint is still missing.  Do this once after any
# change to what the iteration carries across k.
# ---------------------------------------------------------------------------
# JSON.jl stores non-finite numbers as null; read defensively so a diagnostic
# never throws on a damaged run.  A NaN in the report is itself the signal.
ckpt_num(x) = x isa Number ? Float64(x) : NaN

function compare_runs(file_a, file_b)
    a = JSON.parse(read(file_a, String))
    b = JSON.parse(read(file_b, String))
    n = min(length(a), length(b))
    @printf("compare_runs: %d vs %d entries, comparing the first %d\n",
            length(a), length(b), n)
    scalar_fields = ["d_branch", "d_contour", "omega_gap", "zeta_alpha",
                     "zeta_d_med", "exp_arg_max", "dist_u", "dist_l",
                     "h_local", "pinch_miss", "grid_share", "omega_pr_fit",
                     "omega2_fit", "n_L", "adapt_level"]
    vec_fields = ["L", "alpha_L_u", "alpha_L_l", "F", "omega_F", "omega_F_at_L"]
    worst = Dict{String,Float64}()
    first_bad = 0
    first_bad_field = ""
    for j in 1:n
        ea = a[j]
        eb = b[j]
        it = get(ea, "iteration", j)
        for f in scalar_fields
            xa = get(ea, f, nothing)
            xb = get(eb, f, nothing)
            (xa isa Number && xb isa Number) || continue
            d = abs(Float64(xa) - Float64(xb))
            worst[f] = max(get(worst, f, 0.0), d)
            if d > 0 && first_bad == 0
                first_bad = it
                first_bad_field = f
            end
        end
        for f in vec_fields
            va = get(ea, f, nothing)
            vb = get(eb, f, nothing)
            (va isa AbstractVector && vb isa AbstractVector) || continue
            if length(va) != length(vb)
                worst[f] = Inf
                if first_bad == 0
                    first_bad = it
                    first_bad_field = f * " (length)"
                end
                continue
            end
            d = 0.0
            for i in eachindex(va)
                d = max(d, abs(complex(ckpt_num(va[i]["re"]), ckpt_num(va[i]["im"])) -
                               complex(ckpt_num(vb[i]["re"]), ckpt_num(vb[i]["im"]))))
            end
            worst[f] = max(get(worst, f, 0.0), d)
            if d > 0 && first_bad == 0
                first_bad = it
                first_bad_field = f
            end
        end
    end
    for f in sort(collect(keys(worst)))
        @printf("  %-14s max|A-B| = %.6e\n", f, worst[f])
    end
    if first_bad == 0
        println("  IDENTICAL over the compared range -- the checkpoint is complete.")
    else
        @printf("  FIRST DIVERGENCE at iteration %d, field %s -- state is missing from the checkpoint.\n",
                first_bad, first_bad_field)
    end
    flush(stdout)
    return worst
end
function truncate_json!(filename; offset=0)
    json_str = open(filename, "r") do file
        read(file, String)
    end
    data = JSON.parse(json_str)
    n = length(data)
    idx = n + offset
    truncated = data[1:idx]
    # v4.7 CHANGE 8: atomic, like every other write to this file.
    atomic_write(filename, JSON.json(truncated))
    println("Truncated JSON to step with iteration=", truncated[end]["iteration"], " (kept $idx entries).")
end

#truncate_json!("contour_iteration.json"; offset=-1)
load_on_workers()
#plot_alpha()
