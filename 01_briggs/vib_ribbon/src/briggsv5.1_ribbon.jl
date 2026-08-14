# briggsv5.1_ribbon.jl -- symmetric contours + a tracker that respects BOTH
# continuities. Built on briggsv5_ribbon.jl; see ../ANALYSIS_v5.md for the data
# behind every number quoted here.
#
# WHAT THE v5 RUN SHOWED (1001 frames)
#   v5's frame-to-frame predictor did fix the global branch swap -- fixedw_u = 0
#   for all 1001 frames, arc value constant at 0.413445+0.088846i, loop_mono = 0,
#   zero fallbacks. But it introduced a TEAR: a discontinuity in alpha_u at
#   contour index 172|173 (omega_r = 0.1438|0.1447) that opens at iteration 15
#   and never heals, growing 0.0035 -> 0.0347 -> 0.0508 -> ... -> 0.1824 by
#   iteration 1000, always at the same node.
#
#   Cause: two alpha-branches nearly collide at omega* ~ 0.148 - 0.050i, and the
#   descending horizontal sweeps straight through it between iterations 14 and
#   15. Their separation:
#        omega_i = -0.04579 (frame 13)   min gap over omega_r  0.040
#        omega_i = -0.04979 (frame 14)                         0.0059
#        omega_i = -0.05055                                    0.0027
#        omega_i = -0.05379 (frame 15)                         0.033
#   At frame 14 the tracked curve does a fake "turn" at idx 173-176 with the
#   steps collapsing 0.0037 -> 0.0003 -- that is the tracker sliding across the
#   crossing. From frame 15 on, prev copies that shape forward node by node:
#   nodes <= 172 stay anchored to the pre-crossing branch, nodes >= 173 to the
#   post-crossing one, and since the prev-predictor has NO along-contour
#   coupling, nothing ever joins them again.
#
#   The v5 mistake was replacing along-contour continuity with frame-to-frame
#   continuity instead of imposing both. The old secant tracker had the first
#   and not the second (it extrapolates past the turn and flips the whole
#   downstream curve -- v4.6); v5 had the second and not the first (it tears the
#   curve at the crossing). Both fail in the same neighbourhood.
#
#   HONEST LIMIT: at closest approach the branches are 0.0027 apart while the
#   eigenvalue accuracy of this discretisation there is 5e-3..1.5e-2 (see
#   ../analysis/conv.py). The separation is BELOW the noise floor, so no
#   continuation scheme can decide which branch is which at omega*. v5.1 makes
#   the curve continuous and makes the conflict VISIBLE and FATAL rather than
#   silent; it does not pretend to resolve it. Resolving it needs the
#   discretisation rewrite (physical-space chebdif, or a clamped
#   (1-eta^2)^2 T_k basis) or an L contour indented around omega*.
#
# THE v5.1 CHANGES
#   A. SYMMETRIC CONTOURS. Goertz Sec. 7.2.2: the inverse Fourier contour F is
#      the WHOLE real alpha axis and the inverse Laplace contour L the whole
#      line Im omega = sigma. v5 truncated both to one side (alpha_r in [0,1],
#      omega_r in [0,0.5]), which is a half-line truncation that breaks the
#      symmetry of the problem. v5.1 uses
#           alpha_r in [-1, 1],   omega_r in [-0.5, 0.5].
#      N is doubled 300 -> 600 so the node SPACING is unchanged: the horizontals
#      stay at d omega_r = 1.0/599/n_sub = 8.35e-4 and F at d alpha_r = 3.34e-3,
#      exactly as in v5. n_L goes 954 -> 1554.
#      Keep N EVEN: with N even neither alpha_r = 0 nor omega_r = 0 is a grid
#      point, which is what you want -- omega_r = 0 is the imaginary omega axis,
#      where the spectrum is symmetric under alpha -> -conj(alpha) (the mirror
#      twins that flipped v4.4's seed), and alpha = 0 is the degenerate
#      no-streamwise-variation limit of the temporal problem.
#      NOTE, not changed here: the ribbon still has ONE detour, over +omega_f.
#      That is correct for the complex source s_t(t) = exp(i omega_f t) of
#      Goertz Eq. (7.44), which has a single pole. A physically real forcing
#      cos(omega_f t) would have poles at BOTH +omega_f and -omega_f and would
#      need a second detour -- now visible as an asymmetry, since the contour
#      covers omega_r < 0 with no detour there. Worth deciding deliberately.
#   B. BOTH CONTINUITIES IN THE CORRECTOR. The prev-frame value is no longer the
#      predictor. It is now a FILTER, in the same two-stage pattern the code
#      already uses for the side of F:
#          1. keep eigenvalues on the required side of F   (as before)
#          2. keep those within dframe_max of prev[j]      (NEW, frame identity)
#          3. of what survives, take the nearest to the along-contour secant
#      Step 3 restores along-contour smoothness AWAY from the crossing; step 2
#      prevents the drift onto a neighbouring branch across frames. Each
#      fallback (empty filter) is counted in nframefb_u/l and flags the node.
#      MEASURED, on the v5 frame-15 nodes and spectra (../analysis/):
#          v4.6 secant only        alpha(0.24) = 0.653742+0.070422i  (WRONG
#                                  branch, err 0.256), max step 0.0093
#          v5   prev as predictor  alpha(0.24) = 0.405795-0.000051i  (right
#                                  branch), max step 0.0352 at idx 173  <- TEAR
#          v5.1 filter + secant    alpha(0.24) = 0.405795-0.000051i  (right
#                                  branch), max step 0.0344 at idx 174
#      So B keeps the correct branch AND the correct endpoint, but it does NOT
#      remove the tear -- it moves it by one node. That is not a shortcoming of
#      the rule, it is the physics: at omega* the two branches are 0.0027 apart
#      and the discretisation resolves 5e-3..1.5e-2, so the frame constraint and
#      the along-contour constraint genuinely disagree and NO local rule can
#      arbitrate. That is precisely why C is fatal rather than advisory: the run
#      stops and says it cannot resolve the crossing, instead of silently
#      producing a torn curve for 1000 iterations as v5 did.
#      The real cure is either the discretisation rewrite, or indenting L around
#      omega* the way Briggs indents F around poles. Neither belongs in v5.1.
#   C. TEAR DETECTOR. maxstep_u/l = max along-contour |d alpha| per frame, AFTER
#      refinement. Above tear_tol = 2*DA_MAX the run stops. v5 had no such check:
#      track_path bisects flagged intervals up to max_refine times and then
#      silently returns whatever it has, discarding the flags -- which is why a
#      discontinuity survived 1000 iterations unreported. On the v5 data this
#      fires at frame 15 (0.0347 > 0.024) with frame 14 clean (0.0037).
#   D. SEED BY PHYSICAL omega_r. v4.6/v5 used start_index = N/4, which on the
#      one-sided grid [0, 0.5] happened to mean omega_r ~ 0.12. On the symmetric
#      grid the same fraction lands at omega_r ~ -0.25. The seed is now given as
#      seed_omega_r = 0.12 and the nearest horizontal node is used, so it stays
#      put whatever the range is.
#   E. side_projections and classify_by_F rewritten allocation-free. They are
#      O(n_evs * N) per node and N has doubled; the old versions allocated a
#      length-N temporary per eigenvalue per node.
#
# ---------------------------------------------------------------------------
# Inherited description (v5), unchanged in substance:
#
# v4.6 with the four fixes the v4.6 post-mortem called
# for. The PHYSICS IS UNCHANGED: same potential functions, same descent
# dynamics, same acceptance rules, same branch definition (side of F, nearest),
# same ribbon geometry. Everything changed below is resolution, verification,
# and what the run is asked to report.
#
# WHY  (see ../ANALYSIS_v4.6.md and ../analysis/ for scripts and data)
#   In the v4.6 run the upper branch hopped onto a neighbouring eigenvalue at
#   iteration 14, on the LEFT HORIZONTAL at omega_r ~ 0.146 -- upstream of the
#   detour, not in it. Iterations 14..1001 (98.7 % of the run) were on the wrong
#   branch. It was NOT the side-of-F filter (removing it reproduces the same
#   values to 4 decimals, zero fallbacks) and NOT a branch point: the closed
#   omega-rectangle spanning the two horizontal levels has monodromy 0.0844 at
#   the run's own spacing and exactly 0.00000 once refined. It was plain
#   under-resolution of the horizontal --
#
#     nodes on the left horizontal   alpha at omega_r=0.24, omega_i=-0.05130
#       48  (v4.6,  d_omega_r=2.53e-3)     0.65360 + 0.08051i   <- what v4.6 got
#       96  (       d_omega_r=1.25e-3)     0.40647 + 0.00375i
#      192  (       d_omega_r=6.2e-4 )     0.40647 + 0.00375i   <- converged
#
#   -- and every guard passed while it happened: max accepted step 0.0151
#   against da_max = 0.03, zero side-fallbacks, sideclass and loop_mono flat.
#
# THE FOUR CHANGES
#   1. RESOLUTION.  N = 100 -> 300, so the ribbon horizontals are sampled at
#      d_omega_r = 0.5/299/n_sub = 8.4e-4, three times finer than the 1.25e-3 at
#      which the continuation converged.  n_L: 570 -> 954, F: 100 -> 300 points.
#      num_modes STAYS at 150 -- raising that makes things worse, because the
#      coefficient-space Chebyshev recursion gives max|D4| ~ n^8 (1.4e15 at 100,
#      3.7e16 at 150, 3.8e17 at 200) and at 220 the spectrum near alpha ~ 0.31
#      is already garbage.  da_max is retightened 0.03 -> DA_MAX = 0.012 to
#      match the new spacing; at 0.03 the guard was ~12x the typical step.
#   2. FIXED-OMEGA CONSISTENCY CHECK.  The riser tops and the whole arc sit at
#      omega values that never move, so alpha there MUST be identical in every
#      frame.  fixedw_u/l measures exactly that against frame 1; above
#      fixed_tol the run stops, because a hop has occurred upstream and every
#      later frame is meaningless.  This is the v4.2 "decisive test" reinstated:
#      it would have caught iteration 14 at once, and it is worth more than all
#      four v4.6 diagnostics, every one of which passed while the branch hopped.
#   3. SHORT RUN.  1000 -> n_iterations = 60.  The arc values froze at iteration
#      14 and never moved again; for the ribbon there is nothing to gain from a
#      long descent, because the answer lives at fixed real omega.
#   4. RIBBON REPORTING instead of pinch termination.  Goertz, Sec. 7.2.4: for a
#      harmonic source at real omega_0 the response is the residue at that pole,
#      P ~ P_hat(x, omega_0, y) e^{i omega_0 t}, with the spatial behaviour set
#      by alpha_u(omega_0) for x > d and alpha_l(omega_0) for x < -d.  The
#      verdict is a SIGN test -- which branches cross Im alpha = 0 as L is
#      lowered -- not a collision test.  So the run reports alpha_u(omega_f),
#      alpha_l(omega_f) and the crossing counts, and branch_distance/pinch_tol
#      is kept only as the ABSOLUTE-INSTABILITY ALARM (a pinch with Im omega > 0
#      would invalidate the whole ribbon argument).  In v4.6 branch_distance
#      went 3.31 -> 2.19 against pinch_tol = 1e-4: unreachable here, because
#      below F there is nothing but the wall ladder at -3.33i, -6.57i, ...
#   5. FRAME-TO-FRAME CONTINUITY (added after the first v5 run, which still
#      hopped at iteration 14 even at N=300 -- refining the contour was not
#      enough).  Diagnosis, from the 15-frame v5 JSON:
#        - frames 14 and 15 are BOTH smooth (max step 0.0037 / 0.0039 against
#          DA_MAX = 0.012) and both consist of genuine eigenvalues.  Neither
#          jumps: the two curves drift apart gradually from the seed onwards.
#        - the two branches involved NEVER collide.  Minimum separation 0.036
#          over omega_r in [0.12,0.20] x omega_i in [-0.046,-0.060].  So there
#          is no branch point and the continuation is not ambiguous in
#          principle -- the tracker simply loses it.
#        - continuing alpha VERTICALLY from frame 14's value down to frame 15's
#          level settles which is right:
#             omega_r = 0.12375  0.26367-0.04141i -> 0.26815-0.05078i
#                                frame 15: 0.26725-0.05057i   agrees to 9e-4
#             omega_r = 0.15050  0.30095+0.00000i -> 0.30329+0.00108i
#                                frame 15: 0.34548-0.01060i   differs by 0.044
#             omega_r = 0.24000  0.40660+0.00568i -> 0.40725-0.00015i
#                                frame 15: 0.65700+0.07157i   differs by 0.260
#          with a max step of 0.0008 on that vertical.  The seed is right, the
#          endpoint is wrong.
#      The fix uses the one piece of information the tracker had and ignored:
#      each frame re-seeded from scratch and tracked along the contour with no
#      memory of the previous frame.  But between frames every node moves by at
#      most |d omega_i| ~ 4e-3 (horizontals shift by exactly that, riser nodes by
#      less, the arc not at all), and the branch moves by <= 1.1e-2 -- against a
#      branch separation of >= 3.6e-2.  So the previous frame's alpha at the SAME
#      node index is a much better predictor than the along-contour secant, with
#      a factor-3 margin.  continue_branch() uses it where available and falls
#      back to the secant otherwise.  This also removes the v4.4-style seed-flip
#      risk, since the seed node is predicted the same way.
#      IMPORTANT: prev is applied ONLY on the horizontals (mask_to_horizontals).
#      The risers and the arc are deliberately left on the secant, because their
#      omega is pinned -- if prev drove them too, alpha there would reproduce
#      frame 1 BY CONSTRUCTION and fixedw_u/l could never fail.  As it stands
#      they are reached by continuation from the prev-anchored end of the
#      horizontal, so fixedw still tests the thing that actually broke.
#      Verified against the v5 frame-15 data: with prev on the horizontal the
#      endpoint at omega_r = 0.24 comes out 0.406956+0.000022i (truth
#      0.40725-0.00015i, error 3.4e-4) instead of 0.656941+0.071618i, and the
#      riser+arc traverse from there lands on 0.414044+0.089296i against
#      frame 1's 0.413445+0.088846i.
#      Also reported: framedev_u/l = max |alpha_new - alpha_prev| on the
#      horizontals, with a WARN above DFRAME_MAX.  That is a report, not a test
#      -- it can be masked if the true branch genuinely moves a long way.
#
#   Plus, cheap and optional: a GRID-CONVERGENCE check (grid_check_every) that
#   retracks the same contour with every segment refined by grid_check_mult and
#   compares at the coincident nodes.  gridconv_u/l measures the v4.6 failure
#   mode directly.
#
# NOT fixed here -- these need a bigger rewrite, see ANALYSIS_v4.6.md Sec. 5:
#   the n^8 conditioning of the D3/D4 recursion, which leaves the eigenvalues
#   near alpha ~ 0.31 uncertain by about their own separation; F being dragged
#   bodily down to Im alpha ~ -1.15 instead of staying a deformation of the real
#   axis; the radius-7 rolling average applied to F every single iteration.
#
# ---------------------------------------------------------------------------
# Inherited description (v4.6), unchanged in substance:
#
# THE CLOSED LOOP. First version since the abandoned
# briggsv4_ribbon.jl in which the vibrating-ribbon contour is INSIDE the
# spatial<->temporal iteration instead of being painted on afterwards:
#
#     F --(temporal OS)--> omega_F --bounds--> L = RIBBON contour (horizontals
#     at omega_i + risers + arc over omega_f) --(spatial OS, tracked)-->
#     alpha_L_u, alpha_L_l --(potential Phi)--> deforms F --> ...
#
# The ribbon's alpha branches -- including the detour excursions and the sheet
# change across the branch point found near omega ~ 0.25 - 0.05i -- now steer
# F's deformation and the pinch search. omega_f is no longer cosmetic.
#
# Relation to previous versions:
#   - v4.1..v4.5 kept the descent byte-identical to briggsv4.jl (straight L)
#     and drew the ribbon per frame as display only. v4.5's analysis proved the
#     detour changes the branches non-trivially (monodromy ~ 1.9), so the
#     straight-line loop and the ribbon loop are genuinely different problems.
#   - briggsv4_ribbon.jl put the ribbon in the loop first but tracked branches
#     with per-point nearest matching (the v4.1/v4.2 mode-hopping mechanism),
#     feeding noise into Phi every iteration.
#
# v4.6 keeps the descent THEORY of briggsv4.jl unchanged -- same potential
# functions, same descent dynamics, same acceptance rules, same pinch criterion
# -- and changes only contour geometry, tracking mechanism, and the mechanical
# fixes the ribbon geometry forces:
#   1. contour_L() IS hankel_L(omega_i); trial contours are ribbon contours.
#   2. in-loop branches come from the v4.5 tracker (interior seed by the
#      side-of-F rule, side filter, secant predictor, separation guard,
#      adaptive bisection) via ribbon_branches(); contour_alpha_L_conti is dead.
#   3. the omega-level update is evaluated on the HORIZONTAL nodes of the
#      ribbon (the level belongs to them; risers/arc do not descend). The slope
#      term of briggsv4 is kept for structural fidelity; it is identically zero
#      on the constant level vector, as in briggsv4 itself.
#   4. Phi_F endpoint indexing uses end instead of the hardcoded N=100 (the
#      ribbon branches have n_L = 570 points).
#   5. branch_distance is the PAIRWISE minimum (elementwise pairing is
#      meaningless on the ribbon); pinch_tol applies to that.
#   6. per-frame diagnostics: sideclass_u/l, nfallb_u/l, loop_mono_u/l,
#      refinement counts. No descent overlay exists any more -- there is no
#      straight-line run to compare against.
# Same parameters as v4.5: Re=2000, beta=0, num_modes=150, omega_f=0.25,
# r_arc=0.01, lambda=4.0, sigma=3e-5, delta_t=1e-3, zeta=4e-4, 300 iterations.
# No analytical pinch point exists for this case; none is drawn.
# Run from vib_ribbon/src/:  julia briggsv5.1_ribbon.jl
# Kilian Vinzenz Wilhelm
begin
    using Distributed, Plots, BenchmarkTools, FFTW, JSON, Statistics, Printf
    addprocs(7)
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
                eigvals, eigvecs = eigen(A, B)
                return eigvals, eigvecs
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

    # ALPHA contour: F
    @everywhere begin
        # v5.1 CHANGE A: F is a truncation of the WHOLE real alpha axis, so it
        # is centred on 0. N doubles with the range so the spacing is unchanged
        # (d alpha_r = 2.0/599 = 3.34e-3, same as v5's 1.0/299). Keep N EVEN:
        # then alpha_r = 0 is not a grid point, which avoids the degenerate
        # alpha = 0 limit of the temporal problem.
        alpha_r_start = -1.0
        alpha_r_end = 1.0   
        # v5 CHANGE 1: contour resolution. N sets BOTH the F contour and the
        # omega_r grid the ribbon horizontals are built on. 100 -> 300 puts the
        # horizontal spacing at 0.5/299/n_sub = 8.4e-4, below the 1.25e-3 at
        # which the branch continuation converged in the v4.6 post-mortem.
        # This is NOT num_modes -- that stays at 150 (see the header).
        N = 600     
        alpha_r = range(alpha_r_start, alpha_r_end, length=N)
        alpha_i = fill(0.0, N)  
    end
    function contour_F()
        F = [alpha_r[j] + alpha_i[j] * im for j in 1:N]
        return F
    end
    @everywhere F = Vector{Complex{Float64}}[]
    F = contour_F()
    @everywhere function couetteflow_temporal_sing_mode(alpha)
        #eigvals = couetteflow_temporal(alpha)
        eigvals, eigvecs = couetteflow(alpha, nothing, Val(:omega_collocation), Val(:eigen))
        mask = [isfinite(real(eigval)) && isfinite(imag(eigval)) for eigval in eigvals]
        eigvals = eigvals[mask]
        eigvecs = eigvecs[:, mask]
        eigval = eigvals[argmax(imag.(eigvals))]
        eigvec = eigvecs[:,argmax(imag.(eigvals))]
        return eigval, eigvec
    end
    @everywhere function couetteflow_spatial_sing_mode_comparison(
        omega,
        alpha_approximation,
        branch_side::Symbol,
        F,
        normals_F;
        side_tol = 0.0
    )
        eigvals, _ = couetteflow(nothing, omega, Val(:alpha_collocation), Val(:eigen))

        mask = [isfinite(real(eigval)) && isfinite(imag(eigval)) for eigval in eigvals]
        eigvals = eigvals[mask]

        candidates = ComplexF64[]

        for eigval in eigvals
            distances = [abs(eigval - f) for f in F]
            idx_min = argmin(distances)

            f_near = F[idx_min]
            normal = normals_F[idx_min]

            # signed distance relative to F
            proj = real(conj(normal) * (eigval - f_near))

            if branch_side == :upper && proj > side_tol
                push!(candidates, eigval)
            elseif branch_side == :lower && proj < -side_tol
                push!(candidates, eigval)
            end
        end

        if isempty(candidates)
            # fallback: old behavior, but only if side-filter finds nothing
            diffs = abs.(ComplexF64(alpha_approximation) .- ComplexF64.(eigvals))
            return eigvals[argmin(diffs)]
        end

        diffs = abs.(ComplexF64(alpha_approximation) .- candidates)
        return candidates[argmin(diffs)]
    end
    function contour_omega_F(F)
        @everywhere omega_F = Complex{Float64}[]
        omega_F_results = pmap(alpha -> couetteflow_temporal_sing_mode(alpha), F)
        omega_F = [x for (x, _) in omega_F_results]
        return omega_F
    end
    omega_F = contour_omega_F(F)
    # OMEGA contour: L
    @everywhere begin
        # v5.1 CHANGE A: L is a truncation of the whole line Im omega = sigma,
        # so it too is centred on 0. Horizontal spacing stays at
        # 1.0/599/n_sub = 8.35e-4, the value v5 was built around. With N even,
        # omega_r = 0 is not a grid point -- good, because that is the imaginary
        # omega axis where the spectrum is symmetric under alpha -> -conj(alpha).
        omega_r_start = -0.5
        omega_r_end = 0.5     
        omega_i = 0.0
        omega_r = range(omega_r_start, omega_r_end, length=N)
        # ---- vibrating-ribbon marker (the ONLY new input) ----
        # omega_f = real forcing frequency the drawn L-contour detours around.
        # Must lie strictly inside (omega_r_start, omega_r_end). Cosmetic for the
        # descent; it positions the arc in the saved animation frames.
        omega_f = 0.25
        r_arc   = 0.01
    end
    function contour_L()
        # v4.6: the RIBBON contour IS the descent contour -- the loop is closed.
        return hankel_L(omega_i)
    end
    # ---- ribbon contour geometry (v4.3) ----------------------------------
    # The horizontal parts sit ON the descent grid omega_r, so they are exactly the
    # omega points the descent already solved and their alpha values can be taken
    # straight from the descent's tracked branches (no re-classification).
    n_ver = 160          # points per vertical riser
    n_arc = 60           # points on the semicircle over omega_f
    n_bot = 21           # points on the short bottom segment (monodromy test only)
    n_sub = 2            # horizontal refinement: subdivisions per descent interval
    idx_hor_left  = findall(w -> w <= omega_f - r_arc, collect(omega_r))
    idx_hor_right = findall(w -> w >= omega_f + r_arc, collect(omega_r))
    @assert !isempty(idx_hor_left) && !isempty(idx_hor_right) "omega_f +/- r_arc outside the omega_r range"

    # Refine the horizontals while KEEPING the descent nodes as a subset, so the
    # red-vs-grey cross-check is at identical omega with no interpolation.
    function subdivide(xs, m)
        m <= 1 && return collect(float.(xs))
        out = Float64[]
        for i in 1:(length(xs) - 1), k in 0:(m - 1)
            push!(out, xs[i] + (xs[i+1] - xs[i]) * k / m)
        end
        push!(out, float(xs[end]))
        return out
    end
    wl_ref = subdivide(collect(omega_r)[idx_hor_left],  n_sub)
    wr_ref = subdivide(collect(omega_r)[idx_hor_right], n_sub)
    n_hl = length(wl_ref); n_hr = length(wr_ref)
    # positions of the descent nodes inside the full ribbon contour
    map_hl = [(j - 1) * n_sub + 1 for j in 1:length(idx_hor_left)]
    off_hr = n_hl + n_ver + n_arc + n_ver
    map_hr = [off_hr + (j - 1) * n_sub + 1 for j in 1:length(idx_hor_right)]
    i_riser_bot_L = n_hl + 1        # bottom of the left riser,  omega_f - r_arc + i*omega_i
    i_riser_bot_R = off_hr          # bottom of the right riser, omega_f + r_arc + i*omega_i

    riser_left_nodes(omega_i_level)  = ComplexF64.((omega_f - r_arc) .+ im .* collect(LinRange(omega_i_level, 0.0, n_ver)))
    riser_right_nodes(omega_i_level) = ComplexF64.((omega_f + r_arc) .+ im .* collect(LinRange(omega_i_level, 0.0, n_ver)))
    arc_nodes() = ComplexF64.(omega_f .+ r_arc .* cis.(collect(LinRange(pi, 0.0, n_arc))))

    function hankel_L(omega_i_level)
        # NEAT Eq.(24) ribbon contour at horizontal level omega_i_level:
        # horizontals (descent grid, subdivided) + risers + semicircle over omega_f.
        seg1 = ComplexF64.(wl_ref .+ im * omega_i_level)
        seg5 = ComplexF64.(wr_ref .+ im * omega_i_level)
        return vcat(seg1, riser_left_nodes(omega_i_level), arc_nodes(),
                    reverse(riser_right_nodes(omega_i_level)), seg5)
    end
    # short straight segment closing the detour along the bottom, used only for
    # the monodromy test (it lies on the descent line L).
    bottom_nodes(omega_i_level) =
        ComplexF64.(collect(LinRange(omega_f - r_arc, omega_f + r_arc, n_bot)) .+ im * omega_i_level)

    # ---- v5 geometry bookkeeping -------------------------------------------
    # Nodes whose omega does NOT move when the horizontals descend: the top of
    # each riser (pinned to the real axis) and the whole arc. alpha there is a
    # function of omega alone, so it must be frame-independent -- that is what
    # CHANGE 2 checks. The block is contiguous in the contour ordering
    # [seg1 | riser_L | arc | reverse(riser_R) | seg5].
    i_riser_top_L = n_hl + n_ver                       # last node of riser_L
    i_arc_first   = n_hl + n_ver + 1
    i_arc_last    = n_hl + n_ver + n_arc
    i_riser_top_R = i_arc_last + 1                     # first node after the arc
    idx_fixed_omega = collect(i_riser_top_L:i_riser_top_R)
    # apex of the arc: omega = omega_f + i*r_arc, the forcing frequency itself
    i_arc_top = n_hl + n_ver + cld(n_arc, 2)
    # v5.1 CHANGE D: seed by PHYSICAL omega_r, not by an index fraction.
    # v4.6/v5 used start_index = N/4, which on the one-sided grid [0, 0.5] meant
    # omega_r ~ 0.12 -- a place where the branch sits next to F and is well
    # separated from its neighbours. On the symmetric grid [-0.5, 0.5] the same
    # fraction would land at omega_r ~ -0.25 instead.
    seed_omega_r = 0.12
    i_seed_coarse = map_hl[argmin(abs.(collect(omega_r)[idx_hor_left] .- seed_omega_r))]

    # Same ribbon contour with every segment refined by `mult`. The coarse nodes
    # are a SUBSET of the fine ones: subdivide(xs, m) keeps its inputs at stride
    # m, and LinRange(a, b, mult*(n-1)+1) contains LinRange(a, b, n) at stride
    # mult -- for the reversed right riser too. So the comparison in the grid
    # convergence check needs no interpolation. Returns (nodes, coarse->fine map).
    function hankel_L_refined(omega_i_level, mult)
        wl = subdivide(collect(omega_r)[idx_hor_left],  n_sub * mult)
        wr = subdivide(collect(omega_r)[idx_hor_right], n_sub * mult)
        nv = mult * (n_ver - 1) + 1
        na = mult * (n_arc - 1) + 1
        rl = ComplexF64.((omega_f - r_arc) .+ im .* collect(LinRange(omega_i_level, 0.0, nv)))
        rr = ComplexF64.((omega_f + r_arc) .+ im .* collect(LinRange(omega_i_level, 0.0, nv)))
        ar = ComplexF64.(omega_f .+ r_arc .* cis.(collect(LinRange(pi, 0.0, na))))
        nodes = vcat(ComplexF64.(wl .+ im * omega_i_level), rl, ar, reverse(rr),
                     ComplexF64.(wr .+ im * omega_i_level))
        nhl_f = length(wl)
        o1 = nhl_f; o2 = o1 + nv; o3 = o2 + na; o4 = o3 + nv
        map_c = Int[]
        for j in 1:n_hl;  push!(map_c,      (j - 1) * mult + 1); end
        for j in 1:n_ver; push!(map_c, o1 + (j - 1) * mult + 1); end
        for j in 1:n_arc; push!(map_c, o2 + (j - 1) * mult + 1); end
        for j in 1:n_ver; push!(map_c, o3 + (j - 1) * mult + 1); end
        for j in 1:n_hr;  push!(map_c, o4 + (j - 1) * mult + 1); end
        @assert length(map_c) == n_hl + 2 * n_ver + n_arc + n_hr
        @assert maximum(abs.(nodes[map_c] .- hankel_L(omega_i_level))) < 1e-12
        return nodes, map_c
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
    @everywhere function dominant_eigvals(omega, F, normals_F)
        #eigvals = couetteflow_spatial(omega)
        eigvals, eigvecs = couetteflow(nothing, omega, Val(:alpha_collocation), Val(:eigen))
        mask = [isfinite(real(eigval)) && isfinite(imag(eigval)) for eigval in eigvals]
        eigvals = eigvals[mask]
        eigvecs = eigvecs[:, mask]
        eigval_dominant_u = nothing
        eigval_dominant_l = nothing
        eigvec_dominant_u = nothing
        eigvec_dominant_l = nothing
        min_dist_u = Inf
        min_dist_l = Inf
        for (eigval, eigvec) in zip(eigvals, eachcol(eigvecs))
            signed_projections = [real(conj(normal) * (eigval - f)) for (f, normal) in zip(F, normals_F)] # projection
            distances = [abs(eigval - f) for f in F]
            idx_min = argmin(distances)
            proj = signed_projections[idx_min]
            dist_to_line = distances[idx_min]
            if proj > 0.0 && dist_to_line < min_dist_u
                min_dist_u = dist_to_line
                eigval_dominant_u = eigval
                eigvec_dominant_u = copy(eigvec)
            elseif proj < 0.0 && dist_to_line < min_dist_l
                min_dist_l = dist_to_line
                eigval_dominant_l = eigval
                eigvec_dominant_l = copy(eigvec)
            end
        end
        return eigval_dominant_u, eigval_dominant_l, eigvec_dominant_u, eigvec_dominant_l
    end
    function contour_alpha_L_init(L)
        alpha_L_results = pmap(omega -> dominant_eigvals(omega, F, normals_F), L)
        @everywhere alpha_L_u = Complex{Float64}[]
        @everywhere alpha_L_l = Complex{Float64}[]
        for (u, l, _, _) in alpha_L_results
            push!(alpha_L_u, u)
            push!(alpha_L_l, l)
        end
        return alpha_L_u, alpha_L_l
    end
    # v4.6: the real initial branches are computed by ribbon_branches() just
    # before the first frame is written (the tracker is defined further down).
    alpha_L_u = ComplexF64[]
    alpha_L_l = ComplexF64[]
    load_on_workers()

    #######################
    # Eigenvalue tracking #
    #######################
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
                                                                            omega,
                                                                            alpha_tracked,
                                                                            branch_side,
                                                                            F,
                                                                            normals_F
                                                                        )
                (index, alpha_corrected)
            end, args)
            result = tracked[1]
            push!(results, result)
            alpha_current = result[2]
            current_index = result[1]
        end
        return results
    end
    function contour_alpha_L_conti(L)
        N = length(L)
        start_index = floor(Int, N / 4)
        alpha_u_start, alpha_l_start, _, _ = dominant_eigvals(L[start_index], F, normals_F)
        future_u_fwd = @spawn track_branch_pmap(L, start_index, alpha_u_start, +1, :upper)
        future_u_bwd = @spawn track_branch_pmap(L, start_index, alpha_u_start, -1, :upper)

        future_l_fwd = @spawn track_branch_pmap(L, start_index, alpha_l_start, +1, :lower)
        future_l_bwd = @spawn track_branch_pmap(L, start_index, alpha_l_start, -1, :lower)
        results_u_fwd = fetch(future_u_fwd)
        results_u_bwd = fetch(future_u_bwd)
        results_l_fwd = fetch(future_l_fwd)
        results_l_bwd = fetch(future_l_bwd)
        results_u = vcat(results_u_bwd, results_u_fwd[2:end])  # avoid duplication
        results_l = vcat(results_l_bwd, results_l_fwd[2:end])
        sort!(results_u, by = x -> x[1])
        sort!(results_l, by = x -> x[1])
        ordered_alpha_u = [x[2] for x in results_u]
        ordered_alpha_l = [x[2] for x in results_l]
        return ordered_alpha_u, ordered_alpha_l
    end
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
        global omega_r = range(omega_r_start, omega_r_end, length=N)    
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
plot_omega()
plot_alpha()
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
        phi_F = exp(zeta_alpha / (abs(alpha_F - alpha)^s_alpha + epsilon)) - 1.0
        return phi_F
    end
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
        d_d_alpha_r_phi_F = -zeta_alpha * (real(alpha_F) - real(alpha)) * s_alpha * abs(alpha_F - alpha)^(s_alpha - 2) / (abs(alpha_F - alpha)^s_alpha + epsilon)^2.0 * exp(zeta_alpha / (abs(alpha_F - alpha)^s_alpha + epsilon))
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
        d_d_alpha_i_phi_F = -zeta_alpha * (imag(alpha_F) - imag(alpha)) * s_alpha * abs(alpha_F - alpha)^(s_alpha - 2) / (abs(alpha_F - alpha)^s_alpha + epsilon)^2.0 * exp(zeta_alpha / (abs(alpha_F - alpha)^s_alpha + epsilon))
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

            # keep endpoint behavior consistent
            y[1] = y[2]
            y[end] = y[end-1]
        end
        return y
    end
    function complexvec_to_json(vec)
        return [Dict("re" => real(x), "im" => imag(x)) for x in vec]
    end
    filename = "contour_iteration_v5.1_ribbon.json"
end

# ---- ribbon animation frame (v4.3): descent-anchored, adaptively refined ----
@everywhere function spatial_spectrum(omega)
    eigvals, _ = couetteflow(nothing, omega, Val(:alpha_collocation), Val(:eigen))
    return ComplexF64[ev for ev in eigvals if isfinite(real(ev)) && isfinite(imag(ev)) && abs(ev) < 1.0e3]
end

isfin(z) = isfinite(real(z)) && isfinite(imag(z))

# Signed distance of each eigenvalue from the F contour: >0 upper side, <0 lower.
# Computed once per node and shared by both branches.
# v5.1 CHANGE E: allocation-free. This is O(n_evs * N) per contour node and N
# has doubled; the old version built a length-N temporary for every eigenvalue
# at every node. abs2 avoids the sqrt and gives the same argmin.
function side_projections(evs, F, normals_F)
    out = Vector{Float64}(undef, length(evs))
    @inbounds for i in eachindex(evs)
        ev = evs[i]
        k = 1; dbest = Inf
        for j in eachindex(F)
            d = abs2(ev - F[j])
            if d < dbest
                dbest = d; k = j
            end
        end
        out[i] = real(conj(normals_F[k]) * (ev - F[k]))
    end
    return out
end

# Continuation over a PRECOMPUTED set of spectra.
# (B) Candidates are first restricted to the required SIDE OF F -- the rule of
# couetteflow_spatial_sing_mode_comparison() -- and only then is the nearest to
# the predictor taken. If a side is empty the step falls back to nearest-overall,
# exactly as the descent does, and the fallback is counted so it cannot pass
# unnoticed. Secant predictor with the displacement clamped to da_max; a node is
# flagged when the accepted step is too large or the winner is not clearly nearer
# than the runner-up.
# v5 CHANGE 1 (cont.): the step guard has to be sized for the sampling. At the
# v4.6 spacing of 2.53e-3 the typical alpha step on the horizontal was ~3.5e-3,
# so da_max = 0.03 sat ~12x above it and never fired -- the iteration-14 hop was
# accepted with a step of 0.0151. At the v5 spacing of 8.4e-4 the typical step
# is ~1.2e-3, so 0.012 is still ~10x headroom but 2.5x tighter in absolute
# terms. Rejected intervals are bisected and re-solved, so a tighter guard costs
# accuracy nothing -- only time.
const DA_MAX = 0.012
# v5 CHANGE 5: how far a node's alpha may move between two consecutive frames
# before it is worth flagging. Measured on the v5 run: the true motion for
# d omega_i = 4e-3 is <= 1.1e-2, while the nearest rival branch is >= 3.6e-2
# away, so 0.02 sits comfortably between the two.
const DFRAME_MAX = 0.02

# v5.1: counter for nodes where the side-of-F set and the frame-identity set
# have no eigenvalue in common. Global so it does not have to be threaded back
# through track_path / track_outward; reset by track_on_contour.
n_framefb_global = 0

function continue_branch(nodes, spectra, sides, seed, branch_side;
                         da_max = DA_MAX, sep_frac = 0.6, side_tol = 0.0,
                         prev = nothing, dframe_max = DFRAME_MAX)
    n = length(nodes)
    vals = Vector{ComplexF64}(undef, n)
    flag = falses(n)
    n_fallback = 0
    a_prev  = ComplexF64(seed); w_prev  = nodes[1]
    a_prev2 = ComplexF64(seed); w_prev2 = nodes[1]
    have2 = false
    for j in 1:n
        w = nodes[j]
        pred = a_prev
        if have2
            dw  = w - w_prev
            dwp = w_prev - w_prev2
            if abs(dwp) > 1e-14
                disp = (dw / dwp) * (a_prev - a_prev2)
                abs(disp) > da_max && (disp *= da_max / abs(disp))
                pred = a_prev + disp
            end
        end
        evs = spectra[j]; sd = sides[j]
        if isempty(evs) || !isfin(pred)
            vals[j] = a_prev
            flag[j] = true
        else
            keep = Int[]
            for k in eachindex(evs)
                if (branch_side === :upper && sd[k] >  side_tol) ||
                   (branch_side === :lower && sd[k] < -side_tol)
                    push!(keep, k)
                end
            end
            if isempty(keep)                      # side empty -> descent's fallback
                keep = collect(eachindex(evs))
                n_fallback += 1
                flag[j] = true
            end
            # v5.1 CHANGE B: frame-identity FILTER (not predictor). Between
            # frames every node moves by at most |d omega_i| ~ 4e-3 and the
            # branch by <= 1.1e-2, so anything further than dframe_max from
            # prev[j] is a different branch. What survives is then chosen by the
            # along-contour secant, so the curve stays smooth ALONG the contour
            # as well -- v5 used prev as the predictor instead, got frame
            # continuity without along-contour continuity, and tore the curve at
            # the omega* crossing. An empty filter is counted, not silently
            # ignored: a nonzero count means the two continuities disagree,
            # which is exactly what happens at an unresolvable crossing.
            if prev !== nothing && isfin(prev[j])
                keepf = Int[]
                for k in keep
                    abs(evs[k] - prev[j]) <= dframe_max && push!(keepf, k)
                end
                if isempty(keepf)
                    global n_framefb_global += 1
                    flag[j] = true
                else
                    keep = keepf
                end
            end
            dists = [abs(evs[k] - pred) for k in keep]
            k1 = argmin(dists)
            d1 = dists[k1]
            d2 = Inf
            for (kk, dd) in enumerate(dists)
                kk != k1 && dd < d2 && (d2 = dd)
            end
            cand = evs[keep[k1]]
            vals[j] = cand
            (abs(cand - a_prev) > da_max || d1 > sep_frac * d2) && (flag[j] = true)
        end
        a_prev2 = a_prev; w_prev2 = w_prev
        a_prev  = vals[j]; w_prev = w
        have2 = true
    end
    return vals, flag, n_fallback
end

# Track one path, bisecting every rejected interval and re-solving the inserted
# nodes in a single parallel batch. Values are returned at the ORIGINAL nodes only.
function track_path(nodes0, spectra0, sides0, seed, branch_side;
                    max_refine = 4, node_cap = 6000, da_max = DA_MAX, sep_frac = 0.6,
                    prev = nothing, dframe_max = DFRAME_MAX)
    nodes   = copy(nodes0)
    spectra = Vector{Vector{ComplexF64}}(spectra0)
    sides   = Vector{Vector{Float64}}(sides0)
    pv      = prev === nothing ? nothing : ComplexF64.(collect(prev))
    is_orig = trues(length(nodes))
    n_added = 0
    vals, flag, nfb = continue_branch(nodes, spectra, sides, seed, branch_side;
                                      da_max = da_max, sep_frac = sep_frac, prev = pv,
                                      dframe_max = dframe_max)
    for _ in 1:max_refine
        ins_at = Int[]; ins_w = ComplexF64[]
        for j in 2:length(nodes)
            if flag[j]
                push!(ins_at, j)
                push!(ins_w, 0.5 * (nodes[j-1] + nodes[j]))
            end
        end
        (isempty(ins_w) || length(nodes) + length(ins_w) > node_cap) && break
        new_spec  = pmap(spatial_spectrum, ins_w)
        new_sides = [side_projections(sp, F, normals_F) for sp in new_spec]
        n_added += length(ins_w)
        nn = ComplexF64[]; ss = Vector{Vector{ComplexF64}}()
        vv = Vector{Vector{Float64}}(); oo = Bool[]
        pp = pv === nothing ? nothing : ComplexF64[]
        p = 1
        for j in 1:length(nodes)
            while p <= length(ins_at) && ins_at[p] == j
                push!(nn, ins_w[p]); push!(ss, new_spec[p]); push!(vv, new_sides[p]); push!(oo, false)
                # inserted nodes sit midway, so predict midway too
                pv === nothing || push!(pp, 0.5 * (pv[j-1] + pv[j]))
                p += 1
            end
            push!(nn, nodes[j]); push!(ss, spectra[j]); push!(vv, sides[j]); push!(oo, is_orig[j])
            pv === nothing || push!(pp, pv[j])
        end
        nodes = nn; spectra = ss; sides = vv; is_orig = oo; pv = pp
        vals, flag, nfb = continue_branch(nodes, spectra, sides, seed, branch_side;
                                          da_max = da_max, sep_frac = sep_frac, prev = pv,
                                          dframe_max = dframe_max)
    end
    return vals[is_orig], n_added, nfb
end

# (A) March OUTWARD from an interior seed node, both directions -- the protocol of
# contour_alpha_L_conti(), which starts at N/4 rather than at a contour end.
function track_outward(nodes, spectra, sides, i_seed, seed, branch_side; prev = nothing)
    pf = prev === nothing ? nothing : prev[i_seed:end]
    pb = prev === nothing ? nothing : prev[i_seed:-1:1]
    vf, af, nf = track_path(nodes[i_seed:end],    spectra[i_seed:end],    sides[i_seed:end],    seed, branch_side; prev = pf)
    vb, ab, nb = track_path(nodes[i_seed:-1:1],   spectra[i_seed:-1:1],   sides[i_seed:-1:1],   seed, branch_side; prev = pb)
    vals = vcat(reverse(vb[2:end]), vf)          # vb[1] and vf[1] are the same node
    return vals, af + ab, nf + nb
end

# One-shot continuation for the monodromy leg (short, lies on the descent line).
function track_simple(nodes, seed, branch_side; da_max = DA_MAX, sep_frac = 0.6)
    spectra = pmap(spatial_spectrum, nodes)
    sides   = [side_projections(sp, F, normals_F) for sp in spectra]
    vals, _, _ = continue_branch(nodes, spectra, sides, seed, branch_side;
                                 da_max = da_max, sep_frac = sep_frac)
    return vals
end

# The descent's OWN branch rule -- "of the eigenvalues on each side of F, take the
# one closest to F" -- applied to an already-computed spectrum. Identical logic to
# dominant_eigvals(), but it does not re-solve the eigenvalue problem.
function classify_by_F(evs, sd, F)
    best_u = nothing; best_l = nothing
    du = Inf; dl = Inf
    @inbounds for k in eachindex(evs)                 # v5.1 CHANGE E: no temporaries
        ev = evs[k]
        d2 = Inf
        for j in eachindex(F)
            t = abs2(ev - F[j])
            t < d2 && (d2 = t)
        end
        d = sqrt(d2)
        if sd[k] > 0.0 && d < du
            du = d; best_u = ev
        elseif sd[k] < 0.0 && d < dl
            dl = d; best_l = ev
        end
    end
    return best_u, best_l
end

function descent_branch_clean(branch)
    return ComplexF64[x === nothing ? ComplexF64(NaN, NaN) : ComplexF64(x) for x in branch]
end

function seed_or_nearest(cand, evs, F)
    (cand !== nothing && isfin(ComplexF64(cand))) && return ComplexF64(cand)
    isempty(evs) && return ComplexF64(NaN, NaN)
    return ComplexF64(evs[argmin([minimum(abs.(ev .- F)) for ev in evs])])
end

maxabsdiff(a, b) = isempty(a) ? 0.0 : maximum(abs.(ComplexF64.(a) .- ComplexF64.(b)))

# v4.6: the loop's branch computation. Replaces contour_alpha_L_conti -- the
# branch DEFINITION is unchanged (side of F, nearest), but the continuation uses
# the v4.5 machinery: interior seed by the side-of-F rule, side filter, secant
# predictor, separation guard, adaptive bisection. Spectra are cached for the
# per-frame diagnostics.
last_spectra = Vector{Vector{ComplexF64}}()
last_sides   = Vector{Vector{Float64}}()
last_track_diag = (0, 0, 0, 0)
# v5: the body of v4.6's ribbon_branches, with the seed index passed in and the
# module globals NOT touched, so it can also be called on a refined contour by
# the grid-convergence check without clobbering the frame diagnostics.
nframefb_u = 0
nframefb_l = 0
function track_on_contour(Lc, i_seed; prev_u = nothing, prev_l = nothing)
    spectra = pmap(spatial_spectrum, Lc)
    sides   = [side_projections(sp, F, normals_F) for sp in spectra]
    su, sl  = classify_by_F(spectra[i_seed], sides[i_seed], F)
    seed_u  = seed_or_nearest(su, spectra[i_seed], F)
    seed_l  = seed_or_nearest(sl, spectra[i_seed], F)
    global n_framefb_global = 0
    au, add_u, nfb_u = track_outward(Lc, spectra, sides, i_seed, seed_u, :upper; prev = prev_u)
    global nframefb_u = n_framefb_global
    global n_framefb_global = 0
    al, add_l, nfb_l = track_outward(Lc, spectra, sides, i_seed, seed_l, :lower; prev = prev_l)
    global nframefb_l = n_framefb_global
    return au, al, spectra, sides, (add_u, nfb_u, add_l, nfb_l)
end

# v5 CHANGE 5: the previous accepted frame's branches, used as the predictor.
# Empty on frame 1 (falls back to the along-contour secant); ignored on a
# refined contour, where the lengths do not match.
prev_alpha_u = ComplexF64[]
prev_alpha_l = ComplexF64[]
framedev_u = 0.0
framedev_l = 0.0

# The prev-frame predictor is applied ONLY on the horizontals. The risers and
# the arc are left to along-contour continuation on purpose: their omega is
# pinned (exactly, for the arc), so if prev drove them too, alpha there would
# reproduce frame 1 BY CONSTRUCTION and fixedw_u/l could never fail -- the
# "a good check is one that can fail" rule. Left as is, they are reached by
# continuation from the prev-anchored end of the horizontal, so fixedw genuinely
# tests whether the horizontal ended up on the right branch, which is the
# failure mode that broke v4.6 and the first v5 run. prev buys little there
# anyway: riser nodes move by at most d omega_i and the arc not at all.
function mask_to_horizontals(p)
    p === nothing && return nothing
    q = fill(ComplexF64(NaN, NaN), length(p))
    q[1:n_hl] .= p[1:n_hl]
    q[(off_hr + 1):end] .= p[(off_hr + 1):end]
    return q
end

function ribbon_branches(Lc)
    pu = length(prev_alpha_u) == length(Lc) ? mask_to_horizontals(prev_alpha_u) : nothing
    pl = length(prev_alpha_l) == length(Lc) ? mask_to_horizontals(prev_alpha_l) : nothing
    au, al, spectra, sides, diag = track_on_contour(Lc, i_seed_coarse; prev_u = pu, prev_l = pl)
    global last_spectra = spectra
    global last_sides = sides
    global last_track_diag = diag
    return au, al
end

# ---- v5 CHANGE 2: fixed-omega consistency -------------------------------
# alpha at the riser tops and on the arc depends on omega only, and those omega
# do not move, so these values must be bit-identical in every frame. Reference
# is taken from frame 1; fixedw_u/l is the largest deviation from it. Anything
# above fixed_tol means the branch hopped somewhere upstream (this is exactly
# what happened at v4.6 iteration 14, where the arc value moved from
# 0.413445+0.088846i to 0.552821+0.255033i and then stayed wrong for 987 more
# iterations while every other diagnostic looked healthy).
# How far the branches moved since the last accepted frame. INFORMATIONAL
# ONLY -- with prev as the predictor this is small by construction, so it is a
# report, not a test. The independent test is fixedw_u/l below, which compares
# against frame 1 at nodes whose omega never moves.
function frame_deviation!()
    if !isempty(prev_alpha_u) && length(prev_alpha_u) == length(alpha_L_u)
        hor = vcat(collect(1:n_hl), collect((off_hr + 1):length(alpha_L_u)))
        global framedev_u = maximum(abs.(ComplexF64.(alpha_L_u[hor]) .- prev_alpha_u[hor]))
        global framedev_l = maximum(abs.(ComplexF64.(alpha_L_l[hor]) .- prev_alpha_l[hor]))
    else
        global framedev_u = 0.0
        global framedev_l = 0.0
    end
    return framedev_u, framedev_l
end
function commit_frame_branches!()
    global prev_alpha_u = ComplexF64.(alpha_L_u)
    global prev_alpha_l = ComplexF64.(alpha_L_l)
    return nothing
end

# ---- v5.1 CHANGE C: tear detector ---------------------------------------
# The largest along-contour step in the finished branch, AFTER refinement. The
# contour is continuous in omega, so alpha must be continuous along it; a step
# far above the local spacing means the curve is torn between two branches.
# v5 had no such check -- track_path bisects flagged intervals up to max_refine
# times and then silently returns whatever it has, discarding the flags, which
# is how a 0.18-wide discontinuity survived 1000 iterations unreported.
# On the v5 data this fires at frame 15 (0.0347) with frame 14 clean (0.0037).
maxstep_u = 0.0
maxstep_l = 0.0
function tear_check!()
    global maxstep_u = maximum(abs.(diff(ComplexF64.(alpha_L_u))))
    global maxstep_l = maximum(abs.(diff(ComplexF64.(alpha_L_l))))
    return maxstep_u, maxstep_l
end

fixedw_ref_u = ComplexF64[]
fixedw_ref_l = ComplexF64[]
fixedw_u = 0.0
fixedw_l = 0.0
function fixed_omega_check!()
    if isempty(fixedw_ref_u)
        global fixedw_ref_u = ComplexF64.(alpha_L_u[idx_fixed_omega])
        global fixedw_ref_l = ComplexF64.(alpha_L_l[idx_fixed_omega])
        global fixedw_u = 0.0
        global fixedw_l = 0.0
    else
        global fixedw_u = maximum(abs.(ComplexF64.(alpha_L_u[idx_fixed_omega]) .- fixedw_ref_u))
        global fixedw_l = maximum(abs.(ComplexF64.(alpha_L_l[idx_fixed_omega]) .- fixedw_ref_l))
    end
    return fixedw_u, fixedw_l
end

# ---- v5: grid-convergence check ------------------------------------------
# Retrack the SAME contour at grid_check_mult times the sampling and compare at
# the coincident nodes. This is the direct measure of the v4.6 failure: at the
# v4.6 spacing it would have returned ~0.25, at the v5 spacing it should return
# a value comparable to the eigenvalue noise floor (a few 1e-3). Costs roughly
# one extra frame, so it is run every grid_check_every iterations, not always.
grid_check_mult  = 2
grid_check_every = 10
gridconv_u = -1.0     # -1 = not evaluated in this frame (JSON has no NaN)
gridconv_l = -1.0
function grid_convergence_check!()
    nodes_f, map_c = hankel_L_refined(omega_i, grid_check_mult)
    au_f, al_f, _, _, _ = track_on_contour(nodes_f, grid_check_mult * (i_seed_coarse - 1) + 1)
    global gridconv_u = maximum(abs.(ComplexF64.(au_f[map_c]) .- ComplexF64.(alpha_L_u)))
    global gridconv_l = maximum(abs.(ComplexF64.(al_f[map_c]) .- ComplexF64.(alpha_L_l)))
    return gridconv_u, gridconv_l
end

# ---- v5 CHANGE 4: the ribbon observables ---------------------------------
# Number of times Im(alpha) changes sign along a branch. A branch that starts
# above the real alpha axis and ends below it (crossing while L is lowered) is
# a downstream-amplifying wave for x > d; one crossing the other way amplifies
# upstream for x < -d (Goertz Sec. 7.2.4).
# JSON.jl emits a bare NaN/Inf, which is not valid JSON and which MATLAB's
# jsondecode refuses. Diagnostics are mapped to -1 instead.
jsonnum(x) = isfinite(x) ? float(x) : -1.0

function n_axis_crossings(branch)
    n = 0
    for j in 2:length(branch)
        a = imag(ComplexF64(branch[j-1])); b = imag(ComplexF64(branch[j]))
        (isfinite(a) && isfinite(b) && a * b < 0.0) && (n += 1)
    end
    return n
end

function make_ribbon_frame(iter)
    # v4.6: the loop's own alpha_L_u / alpha_L_l ARE the ribbon branches; no
    # separate display pass. Diagnostics reuse the spectra cached by the last
    # ribbon_branches() call. (F is updated after the branches were accepted, so
    # side/mono diagnostics lag the drawn F by one update -- diagnostic only.)
    add_u, nfb_u, add_l, nfb_l = last_track_diag
    scu = 0.0; scl = 0.0
    for i in vcat(map_hl, map_hr)
        cu, cl = classify_by_F(last_spectra[i], last_sides[i], F)
        cu !== nothing && (scu = max(scu, abs(ComplexF64(alpha_L_u[i]) - ComplexF64(cu))))
        cl !== nothing && (scl = max(scl, abs(ComplexF64(alpha_L_l[i]) - ComplexF64(cl))))
    end
    bn = bottom_nodes(omega_i)
    loop_mono_u = abs(track_simple(bn, ComplexF64(alpha_L_u[i_riser_bot_L]), :upper)[end] - ComplexF64(alpha_L_u[i_riser_bot_R]))
    loop_mono_l = abs(track_simple(bn, ComplexF64(alpha_L_l[i_riser_bot_L]), :lower)[end] - ComplexF64(alpha_L_l[i_riser_bot_R]))
    a_u_wf = ComplexF64(alpha_L_u[i_arc_top])
    a_l_wf = ComplexF64(alpha_L_l[i_arc_top])
    ncr_u = n_axis_crossings(alpha_L_u)
    ncr_l = n_axis_crossings(alpha_L_l)
    @printf("        ribbon: side=%.3e/%.3e mono=%.3e/%.3e | refined u=%d l=%d | side-fallbacks u=%d l=%d\n",
            scu, scl, loop_mono_u, loop_mono_l, add_u, add_l, nfb_u, nfb_l)
    @printf("        v5.1  : fixed-omega=%.3e/%.3e | gridconv=%.3e/%.3e | framedev=%.3e/%.3e | maxstep=%.3e/%.3e | frame-fallbacks u=%d l=%d\n",
            fixedw_u, fixedw_l, gridconv_u, gridconv_l, framedev_u, framedev_l,
            maxstep_u, maxstep_l, nframefb_u, nframefb_l)
    @printf("        ribbon: alpha(omega_f)=%+.6f%+.6fi / %+.6f%+.6fi | Im-crossings u=%d l=%d\n",
            real(a_u_wf), imag(a_u_wf), real(a_l_wf), imag(a_l_wf), ncr_u, ncr_l)
    if framedev_u > DFRAME_MAX || framedev_l > DFRAME_MAX
        @printf("        WARN  : framedev %.3e/%.3e exceeds DFRAME_MAX=%.3g -- the branch moved further in one omega_i step than expected (measured max 1.7e-2). Consider a smaller omega step.\n",
                framedev_u, framedev_l, DFRAME_MAX)
    end
    flush(stdout)
    return Dict(
        "iteration" => iter,
        "omega_f"   => omega_f,
        "L"         => complexvec_to_json(L),
        "alpha_L_u" => complexvec_to_json(ComplexF64.(alpha_L_u)),
        "alpha_L_l" => complexvec_to_json(ComplexF64.(alpha_L_l)),
        "F"         => complexvec_to_json(F),
        "omega_F"   => complexvec_to_json(omega_F),
        "sideclass_u" => scu,         "sideclass_l" => scl,
        "nfallb_u"    => nfb_u,       "nfallb_l"    => nfb_l,
        "loop_mono_u" => loop_mono_u, "loop_mono_l" => loop_mono_l,
        # ---- v5 fields ----
        "fixedw_u"    => jsonnum(fixedw_u),   "fixedw_l"   => jsonnum(fixedw_l),
        "gridconv_u"  => jsonnum(gridconv_u), "gridconv_l" => jsonnum(gridconv_l),
        "framedev_u"  => jsonnum(framedev_u), "framedev_l" => jsonnum(framedev_l),
        "maxstep_u"   => jsonnum(maxstep_u),  "maxstep_l"  => jsonnum(maxstep_l),
        "nframefb_u"  => nframefb_u,          "nframefb_l" => nframefb_l,
        "alpha_u_wf"  => Dict("re" => real(a_u_wf), "im" => imag(a_u_wf)),
        "alpha_l_wf"  => Dict("re" => real(a_l_wf), "im" => imag(a_l_wf)),
        "ncross_u"    => ncr_u,       "ncross_l"    => ncr_l,
        "omega_i"     => omega_i,     "n_L"         => length(L),
    )
end

# ---- streaming frame writer -------------------------------------------------
# Up to v4.3 the ENTIRE frame array was re-serialised and rewritten on every
# iteration. That is O(n^2) I/O (~4 GB of writes over 300 frames) and builds a
# >25 MB String each time; if the filesystem ever returns a short write, the file
# is left truncated mid-object and every later read of it fails. Instead the file
# is opened once and each frame is appended in place. The closing "]" is written
# after every frame and then seeked back over, so the file on disk is a COMPLETE,
# VALID JSON array at all times -- if the job is killed or hits a quota, whatever
# was written up to that point still loads.
mutable struct FrameWriter
    io::IOStream
    n::Int
end
function open_frame_writer(fname)
    io = open(fname, "w+")
    write(io, "[")
    p = position(io)
    write(io, "]")
    flush(io)
    seek(io, p)
    return FrameWriter(io, 0)
end
function write_frame!(fw::FrameWriter, frame)
    fw.n > 0 && write(fw.io, ",")
    write(fw.io, JSON.json(frame))
    p = position(fw.io)
    write(fw.io, "]")          # keep the file valid on disk...
    flush(fw.io)
    seek(fw.io, p)             # ...but overwrite it with the next frame
    fw.n += 1
    return fw.n
end
function close_frame_writer!(fw::FrameWriter)
    write(fw.io, "]")
    flush(fw.io)
    close(fw.io)
    return fw.n
end

begin
    iteration_step = 1
    # v4.6: compute the tracked ribbon branches before the first frame.
    alpha_L_u, alpha_L_l = ribbon_branches(L)
    load_on_workers()
    frame_deviation!()            # v5: zero on frame 1 (there is no previous)
    tear_check!()                 # v5.1: baseline along-contour smoothness
    fixed_omega_check!()          # v5: take the fixed-omega reference from frame 1
    grid_convergence_check!()     # v5: and one convergence check straight away
    @printf("v5.1 setup: N=%d, n_L=%d, omega_r in [%.3f, %.3f], alpha_r in [%.3f, %.3f]\n",
            N, length(L), omega_r_start, omega_r_end, alpha_r_start, alpha_r_end)
    @printf("            d(omega_r)=%.3e, d(alpha_r)=%.3e, seed at omega_r=%.5f (index %d)\n",
            (omega_r[2] - omega_r[1]) / n_sub, alpha_r[2] - alpha_r[1],
            real(L[i_seed_coarse]), i_seed_coarse)
    @printf("            DA_MAX=%.3g, DFRAME_MAX=%.3g, tear_tol=%.3g, num_modes=%d\n",
            DA_MAX, DFRAME_MAX, 2.0 * DA_MAX, num_modes)
    frame_writer = open_frame_writer(filename)
    last_frame = make_ribbon_frame(iteration_step)
    write_frame!(frame_writer, last_frame)
    commit_frame_branches!()      # v5: frame 1 becomes the predictor for frame 2
    iteration_step += 1
end
function branch_distance(alpha_L_u, alpha_L_l)
    # v4.6: PAIRWISE minimum. Elementwise pairing is meaningless on the ribbon
    # (the index runs over horizontals + risers + arc, not a common parameter).
    dmin = Inf
    for u in alpha_L_u, l in alpha_L_l
        d = abs(u - l)
        d < dmin && (dmin = d)
    end
    return dmin
end

function contour_L_at(omega_i_trial)
    return hankel_L(omega_i_trial)   # v4.6: trial contours are ribbon contours
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
            "upper/lower overlap: d=%.3e at upper[%d], lower[%d]",
            min_dist, min_i, min_j
        )
    end

    return true, @sprintf(
        "no overlap: min upper/lower distance=%.3e",
        min_dist
    )
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

    nearest_branch = all_branches[idx]
    distance = distances[idx]

    return nearest_branch, distance
end
    
function directional_dt_check(F_old, F_trial, alpha_L_u, alpha_L_l;
                          move_safety = 0.5,
                          global_max_move = 0.02)

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

        # If distance is almost zero, be extremely careful
        if dist < 1e-12
            toward_amounts[j] = abs(move_vec)
            local_allowed_toward[j] = 0.0
            dt_ok = false
            continue
        end

        direction_to_branch = nearest_branch - f_old
        unit_to_branch = direction_to_branch / dist

        # Positive = moving toward branch
        # Negative = moving away from branch
        toward_amount = real(conj(unit_to_branch) * move_vec)

        allowed_toward = move_safety * dist

        toward_amounts[j] = toward_amount
        local_allowed_toward[j] = allowed_toward

        # Only restrict motion TOWARD the branch
        if toward_amount > allowed_toward
            dt_ok = false
        end

        # Still prevent crazy jumps even if moving away
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
# v5 CHANGE 3: 1000 -> 60. In the v4.6 run alpha at the arc took exactly two
# values over 1001 frames (0.413445+0.088846i for iterations 1-13, then the
# hopped value for 14-1001), because the arc sits at fixed omega. There is
# nothing for the ribbon to gain from a long descent.
n_iterations = 60
for k = 1:n_iterations
    #print_iteration_header(k)
    global omega_i, L, alpha_L_u, alpha_L_l, alpha_i, F, omega_F, iteration_step, frame_writer, last_frame

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
        
        #print_block("OMEGA-CONTOUR UPDATE")

        omegaF_max_i = maximum(imag.(omega_F))
        omega_clearance = 2e-7
        omega_lower_bound = omegaF_max_i + omega_clearance


        # ------------------------------------------------------------
        # REPAIR STEP:
        # If the updated omega_F has moved above the current L contour,
        # then the current L is already invalid.
        # In that case, lift L just above omega_F before trying to move downward.
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
        omega_dt_min = 1e-8

        min_valid_omega_candidates = 2 
        max_omega_jump_factor = 40.0

        if !omega_repaired
            for omega_attempt in 1:50
                # v4.6: the level update is evaluated on the HORIZONTAL nodes of
                # the ribbon (the level omega_i belongs to them; risers/arc do not
                # descend). The slope term is kept for structural fidelity to
                # briggsv4; it is identically zero on the constant level vector.
                hor_idx = vcat(collect(1:n_hl), collect((off_hr + 1):length(L)))
                omega_i_vectorization = fill(omega_i, length(hor_idx))
                omega_i_cache = copy(omega_i_vectorization)

                for jj in 2:(length(hor_idx) - 1)
                    j = hor_idx[jj]
                    omega_i_cache[jj] =
                        omega_i_vectorization[jj] +
                        omega_dt * (
                            -lambda
                            + (omega_i_vectorization[jj+1] - omega_i_vectorization[jj-1]) /
                            (real(L[hor_idx[jj+1]]) - real(L[hor_idx[jj-1]])) *
                            d_d_omega_r_Phi_L(L[j])
                            - d_d_omega_i_Phi_L(L[j])
                        )
                end
                
                # acceptance condition for L 

                omega_i_cache = [isfinite(x) && abs(x) < 10.0 ? x : -Inf for x in omega_i_cache]
                greater_candidates = filter(x -> isfinite(x) && x > omega_lower_bound, omega_i_cache)

                #println("omega attempt = ", omega_attempt)
                #println("omega_dt = ", omega_dt)
                #println("max omega_i_cache = ", maximum(omega_i_cache))
                #println("min omega_i_cache = ", minimum(omega_i_cache))
                #println("omega lower bound = ", omega_lower_bound)
                #println("number valid omega candidates = ", length(greater_candidates))

                if length(greater_candidates) >= min_valid_omega_candidates

                    omega_candidate = greater_candidates[argmin(abs.(greater_candidates .- omega_lower_bound))]

                    omega_jump = abs(omega_candidate - omega_i)
                    max_allowed_omega_jump = max_omega_jump_factor * omega_dt
                    if omega_jump <= max_allowed_omega_jump

                        # Build trial L and trial branches BEFORE accepting omega
                        L_trial = contour_L_at(omega_candidate)
                        alpha_L_u_trial, alpha_L_l_trial = ribbon_branches(L_trial)

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
                                k, overlap_reason, omega_dt, 0.5 * omega_dt
                            )

                            omega_dt *= 0.5
                        end

                    else
                        @printf(
                            "[k=%d] omega rejected: jump %.3e > allowed %.3e | reducing dt_omega %.2e -> %.2e\n",
                            k, omega_jump, max_allowed_omega_jump, omega_dt, 0.5 * omega_dt
                        )

                        omega_dt *= 0.5
                    end

                else
                    @printf(
                        "[k=%d] omega rejected: valid=%d < required=%d | reducing dt_omega %.2e -> %.2e\n",
                        k, length(greater_candidates), min_valid_omega_candidates, omega_dt, 0.5 * omega_dt
                    )

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
            #global omega_i = omega_lower_bound
            break
        end
        if omega_repaired
            global L = contour_L()
            load_on_workers()

            @everywhere begin
                normals_F = contour_normals(F)
            end

            global alpha_L_u, alpha_L_l = ribbon_branches(L)
            load_on_workers()
        else
            load_on_workers()

            @everywhere begin
                normals_F = contour_normals(F)
            end
        end
                    #
        alpha_i_cache = copy(alpha_i)

        #print_block("BRANCH DISTANCES")

        d_branch = branch_distance(alpha_L_u, alpha_L_l)
        d_contour = contour_distance(F, alpha_L_u, alpha_L_l)
        # v4.6: pairwise -- F has N=100 points, the ribbon branches have n_L=570
        dist_u = minimum(minimum(abs.(f .- alpha_L_u)) for f in F)
        dist_l = minimum(minimum(abs.(f .- alpha_L_l)) for f in F)
        branch_factor = branch_slowdown_factor(d_branch)
        local_delta_t = delta_t * branch_factor

        #println("dist to upper branch = ", dist_u)
        #println("dist to lower branch = ", dist_l)
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
        # v5 CHANGE 4: pinch_tol is no longer the goal of the run. For the ribbon
        # it is the ABSOLUTE-INSTABILITY ALARM: if the two branches do collide
        # while L is being lowered, and the collision is at Im omega > 0, the
        # harmonic-source argument of Goertz Sec. 7.2.4 does not apply at all and
        # the run must stop and say so.
        pinch_tol = 1e-4
        # v5 CHANGE 2: alpha at fixed-omega nodes must not move between frames.
        fixed_tol = 1e-8
        # v5.1 CHANGE C: and the branch must be continuous ALONG the contour.
        tear_tol = 2.0 * DA_MAX
        fixed_omega_check!()
        tear_check!()
        if k % grid_check_every == 0
            grid_convergence_check!()
        end

        #println("d_branch = ", d_branch, ", d_contour = ", d_contour)
        stop_after_save = false
        stop_reason = ""

        if !isfinite(d_branch) || !isfinite(d_contour)
            @printf(
                "[k=%d] STOP: non-finite distance | d_branch=%.3e | d_contour=%.3e\n",
                k, d_branch, d_contour
            )
            break
        end

        if d_branch < pinch_tol
            stop_after_save = true
            stop_reason = @sprintf(
                "PINCH: branches collided at omega_i=%.6e (%s) -- Briggs saddle; %s",
                omega_i,
                omega_i > 0.0 ? "Im omega > 0" : "Im omega < 0",
                omega_i > 0.0 ?
                    "ABSOLUTE INSTABILITY, the ribbon/harmonic-source argument does not apply" :
                    "not absolutely unstable, ribbon analysis stays valid"
            )
        end

        # v5.1 CHANGE C: a tear means the curve is on two different branches at
        # once. Refinement has already been tried and failed, so stop.
        if maxstep_u > tear_tol || maxstep_l > tear_tol
            stop_after_save = true
            stop_reason = @sprintf(
                "TEAR: branch discontinuous along the contour -- max along-contour step %.3e (upper) / %.3e (lower) exceeds %.3g after refinement. The tracked curve is on two different branches at once; this is what an unresolvable branch crossing looks like (v5 hit one at omega ~ 0.148-0.050i, where the branches close to 0.0027 while the discretisation resolves only 5e-3..1.5e-2). Fix the eigenvalue solver, or indent L around the crossing.",
                maxstep_u, maxstep_l, tear_tol)
        end

        # v5 CHANGE 2: a hop invalidates every later frame, so stop on it.
        if fixedw_u > fixed_tol || fixedw_l > fixed_tol
            stop_after_save = true
            stop_reason = @sprintf(
                "FIXED-OMEGA CHECK FAILED: alpha moved by %.3e (upper) / %.3e (lower) at nodes whose omega never changes -- the branch has hopped upstream. Refine the horizontals (increase N or n_sub) and rerun.",
                fixedw_u, fixedw_l
            )
        end

    #print_block("ALPHA-CONTOUR")

    accepted = false
    alpha_i_cache = copy(alpha_i)
    alpha_i_smooth = copy(alpha_i)

    for attempt in 1:5
        alpha_i_trial = copy(alpha_i)

        for j in 2:(length(alpha_i_trial) - 1)

            alpha_i_r =
                (alpha_i[j+1] - alpha_i[j-1]) /
                (alpha_r[j+1] - alpha_r[j-1])

            alpha_i_rr =
                (alpha_i[j+1] - 2.0 * alpha_i[j] + alpha_i[j-1]) /
                ((alpha_r[j+1] - alpha_r[j]) * (alpha_r[j] - alpha_r[j-1]))

            rhs_j =
                (
                    alpha_i_r * d_d_alpha_r_Phi_F(F[j])
                    - d_d_alpha_i_Phi_F(F[j])
                    + sigma * alpha_i_rr
                ) / (1.0 + alpha_i_r^2)

            # Safety cap: prevents 1e59 blow-up but still allows aggressive motion
            rhs_cap = 100.0
            rhs_j = clamp(rhs_j, -rhs_cap, rhs_cap)

            alpha_i_trial[j] = alpha_i[j] + local_delta_t * rhs_j
        end

        alpha_i_trial[1] = alpha_i_trial[2]
        alpha_i_trial[end] = alpha_i_trial[end-1]

        # strong smoothing of the trial contour
        alpha_i_smooth = rolling_average_filter(alpha_i_trial, 7)

        alpha_i_smooth[1] = alpha_i_smooth[2]
        alpha_i_smooth[end] = alpha_i_smooth[end-1]

        factor = acceptance_factor.(F, Ref(alpha_L_u), Ref(alpha_L_l))

        move_raw = abs.(alpha_i_trial .- alpha_i)
        move_smooth = abs.(alpha_i_smooth .- alpha_i)

        #println("attempt = ", attempt)
        #println("delta_t = ", delta_t)
        #println("max raw move = ", maximum(move_raw))
        #println("max smooth move = ", maximum(move_smooth))
        #println("min factor = ", minimum(factor))
        #println("max factor = ", maximum(factor))

        if all(move_raw .<= factor .* max.(abs.(alpha_i), 1e-12))
            alpha_i_cache = copy(alpha_i_trial)
            accepted = true
            alpha_status = "accepted"
            alpha_attempt_used = attempt
            max_raw_move = maximum(move_raw)
            max_smooth_move = maximum(move_smooth)
            break
        else
            local_delta_t *= 0.5

            alpha_status = "smoothed/reduced"
            alpha_attempt_used = attempt
            max_raw_move = maximum(move_raw)
            max_smooth_move = maximum(move_smooth)
        end
    end
       
        global alpha_i = copy(alpha_i_smooth)
        global F = contour_F()
        load_on_workers()
        global omega_F = contour_omega_F(F)
        @printf(
            "[%04d] jump=%9.3e | dUL=%9.3e | dUF=%9.3e | dLF=%9.3e | dt=%9.3e | ζ=%9.3e | %s/%s\n",
            k,
            omega_jump,
            d_branch,
            dist_u,
            dist_l,
            local_delta_t,
            zeta_common,
            omega_status,
            alpha_status
        )
        flush(stdout)
        load_on_workers()
        frame_deviation!()        # v5: measure BEFORE prev is overwritten
        global last_frame = make_ribbon_frame(iteration_step)
        write_frame!(frame_writer, last_frame)
        commit_frame_branches!()  # v5: this frame becomes the next predictor
        global iteration_step += 1

        if stop_after_save
            @printf(
                "[k=%d] STOP AFTER SAVE: %s | d_branch=%.3e | d_contour=%.3e | omega_i=%.6e\n",
                k, stop_reason, d_branch, d_contour, omega_i
            )
            flush(stdout)
            break
        end
    end    
    #print_iteration_footer(k)
end

n_frames_written = close_frame_writer!(frame_writer)
@printf("\nRun finished: %d frames written to %s\n", n_frames_written, filename)

# ---- v5 CHANGE 4: the ribbon verdict -------------------------------------
# Goertz Sec. 7.2.4. With no absolute instability, L is lowered through the real
# axis with the pole at omega_f kept enclosed by the detour, and the long-time
# response is the residue at omega_f: harmonic in time, with the spatial
# behaviour of alpha_u(omega_f) for x > d and alpha_l(omega_f) for x < -d.
# Downstream growth requires Im alpha_u < 0; upstream growth requires
# Im alpha_l > 0. Neither is a statement about the pinch.
begin
    a_u_wf = ComplexF64(alpha_L_u[i_arc_top])
    a_l_wf = ComplexF64(alpha_L_l[i_arc_top])
    println("\n================ RIBBON RESULT ================")
    @printf("omega_f          = %.6f  (r_arc = %.4f)\n", omega_f, r_arc)
    @printf("final omega_i    = %.6e   after %d iterations\n", omega_i, iteration_step - 1)
    @printf("alpha_u(omega_f) = %+.6f %+.6fi   -> Im %s 0\n",
            real(a_u_wf), imag(a_u_wf), imag(a_u_wf) < 0 ? "<" : ">")
    @printf("alpha_l(omega_f) = %+.6f %+.6fi   -> Im %s 0\n",
            real(a_l_wf), imag(a_l_wf), imag(a_l_wf) > 0 ? ">" : "<")
    @printf("Im-alpha sign changes along the branches: upper %d, lower %d\n",
            n_axis_crossings(alpha_L_u), n_axis_crossings(alpha_L_l))
    println(imag(a_u_wf) < 0.0 ?
        "  => alpha_u has crossed the real alpha axis: spatially AMPLIFYING for x > d" :
        "  => alpha_u has NOT crossed the real alpha axis: no downstream growth")
    println(imag(a_l_wf) > 0.0 ?
        "  => alpha_l has crossed the real alpha axis: spatially AMPLIFYING for x < -d" :
        "  => alpha_l has NOT crossed the real alpha axis: no upstream growth")
    @printf("verification: fixed-omega %.3e / %.3e (tol 1e-8, the independent check)\n", fixedw_u, fixedw_l)
    @printf("              grid convergence %.3e / %.3e, frame deviation %.3e / %.3e (informational)\n",
            gridconv_u, gridconv_l, framedev_u, framedev_l)
    @printf("              max along-contour step %.3e / %.3e (tear_tol %.3g), frame-fallbacks %d / %d\n",
            maxstep_u, maxstep_l, 2.0 * DA_MAX, nframefb_u, nframefb_l)
    @printf("absolute-instability alarm: min pairwise |alpha_u - alpha_l| = %.4e (pinch_tol = 1e-4)\n",
            branch_distance(alpha_L_u, alpha_L_l))
    println("==============================================")
end
flush(stdout)

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
# NOTE: load_step() above is kept for interactive use, but it is NOT called here
# any more. Re-parsing the full (tens of MB) results file at the end of a batch
# run is pointless -- the final frame is already in memory -- and it is what
# crashed the v4.3 job when the file on disk had been truncated.
# It also restored the WRONG arrays in the ribbon scripts: entry["L"] /
# entry["alpha_L_u"] are the 476-point RIBBON contour, not the 100-point descent
# state. The *_descent fields are the ones that define the descent.
iteration_step = last_frame["iteration"]
L         = json_to_complexvec(last_frame["L"])
alpha_L_u = json_to_complexvec(last_frame["alpha_L_u"])
alpha_L_l = json_to_complexvec(last_frame["alpha_L_l"])
F         = json_to_complexvec(last_frame["F"])
omega_F   = json_to_complexvec(last_frame["omega_F"])
println("Loaded iteration: ", iteration_step)
global omega_i = imag(L[1])
global alpha_i = imag.(F)
load_on_workers()
#plot_omega()
#plot_alpha()
function truncate_json!(filename; offset=0)
    json_str = open(filename, "r") do file
        read(file, String)
    end
    data = JSON.parse(json_str)
    n = length(data)
    idx = n + offset
    truncated = data[1:idx]
    open(filename, "w") do file
        write(file, JSON.json(truncated))
    end
    println("Truncated JSON to step with iteration=", truncated[end]["iteration"], " (kept $idx entries).")
end

#truncate_json!("contour_iteration.json"; offset=-1)
load_on_workers()
#plot_alpha()

#F
#plot(test)