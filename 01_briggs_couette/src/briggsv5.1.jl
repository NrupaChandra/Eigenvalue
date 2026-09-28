###############################################################################
#
#   BRIGGS' METHOD  --  absolute/convective instability of plane Couette flow
#   v5.1
#
#   Lower the Bromwich contour L until two alpha-branches pinch the inversion
#   contour F.  F is a gauge: any deformation with no pole crossing gives the
#   same omega_p, so F is free to be shaped and the whole run is
#
#        Im omega_p  =  min over admissible F  of  max over alpha in F  Im omega(alpha)
#
#   omega_i is floored by F's own peak; F descends; repeat.  When L has no room
#   to go lower, the only way forward is for F to deform -- that coupling is the
#   method, and it is unchanged from v4.7.
#
#
#   WHAT THIS IS
#
#   v5's dynamics, kept because the run data says they work, with four fixes and
#   with the accumulated switches, dead paths and v4.1->v5 commentary removed:
#   1082 code lines against v5's 2176.
#
#   An earlier v5.1 also replaced the smoothing and the acceptance rule.  That
#   was wrong and this file reverts it -- see "WHAT WAS TRIED AND REVERTED".
#
#
#   THE FOUR FIXES
#
#   FIX 1  Equilibrate the quadratic pencil before linearizing.  ||A12|| was
#          6.8e12 against ||B11|| = 7.2e5, so below a branch gap of ~3e-3 the
#          returned pair was pure noise -- at |omega-omega_p| = 1e-13 the true
#          gap is 1.4e-6 and the unscaled solve reported 3.4e-3 (2463x).  After
#          the fix it tracks the exact square-root law to 1.000 down to 1.4e-5.
#          det(R P S) = det(R) det(S) det(P), so no root moves.
#
#          THIS IS THE ONE THAT MATTERS.  v4.7's best d_branch was 2.328e-3 and
#          v5's was 2.411e-3 -- five versions of F machinery, 3 % apart, both
#          sitting on a float64 noise floor measured independently at 2.23e-3.
#          Everything else here is protecting a descent that already worked.
#
#   FIX 2  The peak parabola must not extrapolate.  Its excess over the sampled
#          maximum is (A-B)^2/8(A+B) -> max(A,B)/8 on a one-sided triple, and
#          the old |t| > 1 guard was dead code (|t| <= 0.5 always).  omega_i is
#          floored by that number, so an overshoot pushes L back UP: on v5's
#          frames 300-1000 the bound sat above Im(omega_p) on 199 of 701, which
#          is the whole of the reversal that took v5 from +4.28e-8 out to
#          +7.26e-5.  With the guard: 0 frames.
#
#   FIX 3  Floor the smoothing half-width at two grid cells.  v5's width filter
#          is an exact no-op below one cell and switched itself off at iteration
#          200; f_ripple then grew to 7.4e-2 against v4.7's 1e-4..6e-3.  The
#          floor keeps v5's widths exactly while they are wide (15.0 cells at
#          iteration 10, 12.2 at 45, 1.0 at 200) and holds a 1-2-1 filter after
#          that instead of nothing.  SMOOTH_HW_CELLS is the knob.
#
#   FIX 4  ADAPT_D_FLOOR (2.0e-3) and ADAPT_MIN_H (5e-7) were both set to the
#          OLD noise floor.  h = 5e-7 alone caps d at about 1.5e-3, so leaving
#          either would make FIX 1 invisible.
#
#
#   WHAT WAS REMOVED FROM v5  (dead in practice, or replaced)
#
#   - the 50-attempt omega descent loop and the 25-probe bisection.  Measured
#     over the whole v5 run: omega_gap was exactly 0 on 966/1000 frames, the
#     bisection never used more than 1 probe, and omega_i was simply pinned to
#     the peak + 1e-9 every time.  This does that directly, with a short
#     fallback bisection for the frames where it is not admissible.
#   - the zeta running median + rate limiter.  It existed to smooth a noisy
#     d_branch; FIX 1 removes the noise, so zeta is now the plain formula.
#   - the per-node trust region.  It was floored at the global value, so at
#     d = 3.015 a node 1.75e-3 from a branch point could still move 0.1508 --
#     86x its own clearance.  The floor defeated the limit entirely.  This is
#     v4.9's single theta, which is what reached f_ripple = 1.25e-5.
#   - the separate checkpoint file.  The log entry carries omega_r and alpha_r,
#     so a resume needs nothing else.
#   - every on/off switch that reproduced an older version, the v4.4 reference
#     path, the plotting helpers, and the unused filters.
#
#
#   WHAT WAS TRIED AND REVERTED  (the first v5.1, for the record)
#
#   - A ripple-triggered filter INSTEAD of the wide one.  It only removes the 2h
#     mode, so the early run went jagged: at iteration 45, f_ripple = 1.2e-2 and
#     omega_i = -0.150, where v5 was at 4.8e-3 and -0.293.  The wide filter is
#     load-bearing for the first 200 iterations.  Reverted; FIX 3 is the correct
#     version of the same idea.
#   - Strict monotone acceptance ("the peak must fall").  The objective is a max
#     over 100 nodes and is not differentiable; from v5's own peak_move the peak
#     rose on 3.1 % of steps over iterations 1-100, 27.1 % over 200-300 and
#     48.2 % over 300-600.  Refusing those steps left F moving only on an escape
#     clause every eighth iteration.  Reverted to v5's tolerance plus best-trial
#     fallback.
#   - A decaying descent gain.  Because rises are normal it decayed 0.30 -> 0.018
#     in 70 iterations and switched off the one mechanism that had earned
#     anything.  Reverted to a fixed DESC_FRAC.
#   - Scaling the descent term by theta/dt (the CURRENT step) instead of
#     theta/dt0 (the one the attempt loop started at).  That made the descent
#     contribution to the displacement independent of dt, so halving dt shrank
#     only the barrier and every trial kept overshooting.  Over iterations
#     82-111 it left omega_i bouncing +/-5e-4 with zero net drift, needing 3-4
#     attempts per step, where v5 over 80-115 moved -4.8e-2 with 31 steps down,
#     5 up, and attempt 1 accepted on 112 of 141 frames.  Fixed.
#   - A graded F grid, and rejecting trials whose peak is unreadable.  Neither is
#     supported by a measurement yet.  Both removed.
#
#
#   NOTHING HERE KNOWS THE PINCH.  No omega'', no alpha_p, no omega_p, and no
#   constant derived from them.  FIX 1 reads only |A11|,|A12|,|B11|; FIX 2 reads
#   only Im(omega_F); FIX 3 reads only d_branch and the grid spacing.
#
#
#   RUN:      julia briggsv5.1.jl         (set RESUME = false for a clean start)
#   OUTPUT:   contour_iteration_v5.1.json
#
#   WATCH, in the console line, and what would falsify each fix:
#     dUL    the branch gap.  It should fall through 2e-3 instead of sticking
#            near 2.4e-3 the way v4.7 and v5 both did.  If it sticks anyway,
#            the pencil was not the wall and FIX 1 is wrong.
#     w      omega_i.  v5 reached 4.28e-8 from Im(omega_p) at iteration 355 and
#            then drifted out to 7.26e-5.  This should hold near its best.
#     asym   min(A,B)/max(A,B) at the peak node.  A run of values below 0.05
#            means F is going jagged and FIX 3 is too weak -- raise
#            SMOOTH_HW_CELLS toward 4.
#     rip    f_ripple.  v4.7 held 1e-4..6e-3 for 501 iterations; v5 ended at
#            7.4e-2.  Staying in v4.7's band is the target.
#     hw/h   the smoothing width in cells.  Should follow v5's (15.0, 12.2, 8.0,
#            1.0 at iterations 10, 45, 75, 200) and then flatten at 2.0.
#     tol    the allowed rise this step.  Inf for the first 24 iterations, then
#            10 % of the recent descent rate, tightening to 0 as the run stalls.
#
###############################################################################

using Distributed, JSON, Statistics, Printf
addprocs(5)

@everywhere using LinearAlgebra, Statistics

# =============================================================================
# 1.  FLOW AND DISCRETISATION
# =============================================================================
@everywhere begin
    Re        = 2000.0
    beta      = 0.0 + 0.0im
    v_g       = 0.0 + 0.0im
    num_modes = 150
    y_start   = 0.0
    y_end     = 1.0
end

@everywhere begin
    # Chebyshev in COEFFICIENT space: D0 maps coefficients to collocation
    # values, D1..D4 are built by the standard recurrence and then rescaled to
    # y in [y_start, y_end].
    yc = [cos((j - 1) * pi / (num_modes - 1)) for j = 1:num_modes]
    yp = (y_start + y_end) / 2 .- yc * ((y_end - y_start) / 2)

    D0 = zeros(Float64, num_modes, num_modes)
    for j = 1:num_modes
        D0[:, j] .= cos.((j - 1) * acos.(yc))
    end
    D1 = zeros(Float64, num_modes, num_modes)
    D2 = zeros(Float64, num_modes, num_modes)
    D3 = zeros(Float64, num_modes, num_modes)
    D4 = zeros(Float64, num_modes, num_modes)
    D1[:, 2] = D0[:, 1];  D1[:, 3] = 4 * D0[:, 2]
    D2[:, 3] = 4 * D0[:, 1]
    for j = 4:num_modes
        D1[:, j] .= 2 * (j - 1) * D0[:, j - 1] + (j - 1) * D1[:, j - 2] / (j - 3)
        D2[:, j] .= 2 * (j - 1) * D1[:, j - 1] + (j - 1) * D2[:, j - 2] / (j - 3)
        D3[:, j] .= 2 * (j - 1) * D2[:, j - 1] + (j - 1) * D3[:, j - 2] / (j - 3)
        D4[:, j] .= 2 * (j - 1) * D3[:, j - 1] + (j - 1) * D4[:, j - 2] / (j - 3)
    end
    smap = -(y_end - y_start) / 2
    D1 ./= smap^1
    D2 ./= smap^2
    D3 ./= smap^3
    D4 ./= smap^4

    u   = yp                 # plane Couette: u = y
    d2u = 0.0
end

# -----------------------------------------------------------------------------
# FIX 1.  Equilibrate the quadratic pencil P(alpha) = a^2*M + a*C + K before
# linearizing.  In coefficient space T_n'''' reaches n^2(n^2-1)(n^2-4)(n^2-9)/105
# -- about 4e16 at n = 149 after the y-mapping -- so K = -A12 carries 1/Re*D4 at
# 6.8e12 while M = B11 sits at 7.2e5, and the companion form then puts that
# block next to A21 = I.  LAPACK's own balancing cannot repair it: the companion
# structure ties the two block-rows together.
#
# R*P(alpha)*S has exactly the same roots for any nonsingular diagonal R, S.
# Eigenvectors transform as v -> S^{-1} v; nothing here reads them.  No tuning
# constant: gamma = 1, 3.09, 100 and 3088 all give the same answer to 5 digits.
# Cost O(N^2): 1.3 ms against a ~250 ms solve.
# -----------------------------------------------------------------------------
@everywhere function pencil_equilibrate(M, C, K; iters::Int = 12)
    n = size(M, 1)
    W = abs.(M) .+ abs.(C) .+ abs.(K)
    R = ones(Float64, n)
    S = ones(Float64, n)
    for _ in 1:iters
        r = sqrt.(max.(vec(maximum(W; dims = 2)), 1e-300))
        R ./= r
        W ./= r                      # row i /= r[i]
        c = sqrt.(max.(vec(maximum(W; dims = 1)), 1e-300))
        S ./= c
        W ./= transpose(c)           # col j /= c[j]
    end
    return R, S
end

# omega given alpha (temporal problem), or the alpha spectrum given omega
# (spatial problem, the quadratic pencil).
@everywhere function couetteflow_omega(alpha)
    A11 = -im * alpha * (u * ones(ComplexF64, 1, length(u))) .* D2 +
           im * alpha * (u * ones(ComplexF64, 1, length(u))) * (alpha^2 + beta^2) .* D0 +
           im * alpha * (d2u * ones(1, length(u))) .* D0 +
           1 / Re .* D4 - 2 / Re * (alpha^2 + beta^2) .* D2 +
           1 / Re * (alpha^2 + beta^2)^2 .* D0 + alpha * v_g .* D0
    A = [-200im * [D0[1:1, :]; D1[1:1, :]];
         A11[3:num_modes-2, :];
         -200im * [D1[num_modes:num_modes, :]; D0[num_modes:num_modes, :]]]
    B11 = -im .* D2 + im * (alpha^2 + beta^2) .* D0
    B = [[D0[1:1, :]; D1[1:1, :]];
         B11[3:num_modes-2, :];
         [D1[num_modes:num_modes, :]; D0[num_modes:num_modes, :]]]
    return eigvals(A, B)
end

@everywhere function couetteflow_alpha(omega)
    A11 = -2im * omega * D1 - 4 / Re * D3 + 4 / Re * beta^2 * D1 -
           im * (u * ones(ComplexF64, 1, length(u))) .* D2 +
           im * beta^2 * (u * ones(1, length(u))) .* D0 +
           im * (d2u * ones(1, length(u))) .* D0 -
           im * v_g .* D2 + im * v_g * beta^2 .* D0
    A12 = im * omega * D2 - im * omega * beta^2 * D0 + 1 / Re * D4 -
          2 / Re * beta^2 * D2 + 1 / Re * beta^4 * D0

    A11 = [zeros(ComplexF64, 2, num_modes);
           A11[3:num_modes-2, :];
           zeros(ComplexF64, 2, num_modes)]
    A12 = [-200im * [D0[1:1, :]; D1[1:1, :]];
           A12[3:num_modes-2, :];
           -200im * [D1[num_modes:num_modes, :]; D0[num_modes:num_modes, :]]]
    B11 = -4 / Re * D2 - 2im * (u * ones(ComplexF64, 1, length(u))) .* D1 +
           2im * v_g .* D1
    B11 = [zeros(ComplexF64, 2, num_modes);
           B11[3:num_modes-2, :];
           zeros(ComplexF64, 2, num_modes)]

    # ---- FIX 1 -------------------------------------------------------------
    Rd, Sd = pencil_equilibrate(B11, A11, A12)
    A11 = Rd .* A11 .* transpose(Sd)
    A12 = Rd .* A12 .* transpose(Sd)
    B11 = Rd .* B11 .* transpose(Sd)
    # ------------------------------------------------------------------------

    Z = zeros(ComplexF64, num_modes, num_modes)
    Id = Matrix{ComplexF64}(I, num_modes, num_modes)
    A = [A11 A12; Id Z]
    B = [B11 Z;   Z  Id]
    return eigvals(A, B)
end

@everywhere function omega_of_alpha(alpha)
    ev = couetteflow_omega(alpha)
    ev = ev[[isfinite(real(e)) && isfinite(imag(e)) for e in ev]]
    return ev[argmax(imag.(ev))]
end

# =============================================================================
# 2.  THE TWO CONTOURS
#
#   F is a graph over alpha_r:  F[j] = alpha_r[j] + i*alpha_i[j].
#   L is horizontal:            L[j] = omega_r[j] + i*omega_i.
#
#   alpha_r and omega_r are both Vectors, not ranges, because both are refined
#   at run time (omega_r by the adaptive L block, alpha_r by FIX 5).
# =============================================================================
const N_F0 = 100
const N_L0 = 100

global alpha_r = collect(range(0.0, 1.0, length = N_F0))
global alpha_i = zeros(Float64, N_F0)
global omega_r_base = collect(range(0.0, 0.5, length = N_L0))
global omega_r = copy(omega_r_base)
global omega_i = 0.0

global F         = ComplexF64[]
global L         = ComplexF64[]
global omega_F   = ComplexF64[]
global alpha_L_u = ComplexF64[]
global alpha_L_l = ComplexF64[]

@everywhere F         = ComplexF64[]
@everywhere normals_F = ComplexF64[]

contour_F() = ComplexF64[alpha_r[j] + alpha_i[j] * im for j in eachindex(alpha_r)]
contour_L() = ComplexF64[omega_r[j] + omega_i * im   for j in eachindex(omega_r)]
contour_L_at(wi) = ComplexF64[omega_r[j] + wi * im   for j in eachindex(omega_r)]

@everywhere function contour_normals(Fv)
    nrm = ComplexF64[]
    for j in 2:(length(Fv) - 1)
        t = Fv[j+1] - Fv[j-1]
        push!(nrm, im * t / abs(t))
    end
    insert!(nrm, 1, nrm[1])
    push!(nrm, nrm[end])
    return nrm
end

# Ship F and its normals to the workers.  spatial_payload takes both as
# arguments, but they must agree with each other, so they are always pushed
# together and never separately.
function push_F!()
    Fc = deepcopy(F)
    nc = contour_normals(Fc)
    @sync for pid in workers()
        @async remotecall_wait(Core.eval, pid, Main, :(F = $Fc; normals_F = $nc))
    end
    global normals_F = nc
    return nothing
end

contour_omega_F(Fv) = ComplexF64.(pmap(a -> omega_of_alpha(a), Fv))

# =============================================================================
# 3.  BRANCH TRACKING
#
#   spatial_payload does the one eigensolve per omega and classifies each root
#   by the sign of its projection onto the nearest F normal.  The tracking
#   itself is a nearest-neighbour walk done afterwards on the main process, so
#   the expensive part parallelises cleanly.
#
#   TRACK_SIDE_SLACK (v4.7 CHANGE 9) is the continuity guard: F's node spacing
#   sets the resolution of the side test, and once the two branch points are
#   closer than that spacing the sign is decided by a margin far below it.  If
#   honouring the side classification costs SLACK times more movement than the
#   nearest root overall, the classification is what is wrong.
# =============================================================================
const TRACK_SIDE_SLACK = 10.0
const TRACK_MIN_STEP   = 1e-9

global track_overrides = 0
global track_jump_u    = 0.0

@everywhere function spatial_payload(omega, Fv, nrm)
    ev = couetteflow_alpha(omega)
    ev = ev[[isfinite(real(e)) && isfinite(imag(e)) for e in ev]]
    up = ComplexF64[]; up_d = Float64[]
    lo = ComplexF64[]; lo_d = Float64[]
    for e in ev
        dmin = Inf; k = 0
        @inbounds for t in eachindex(Fv)
            dd = abs(e - Fv[t])
            if dd < dmin
                dmin = dd; k = t
            end
        end
        proj = real(conj(nrm[k]) * (e - Fv[k]))
        if proj > 0.0
            push!(up, e); push!(up_d, dmin)
        elseif proj < 0.0
            push!(lo, e); push!(lo_d, dmin)
        end
    end
    return (all = ev, up = up, up_d = up_d, lo = lo, lo_d = lo_d)
end

function dominant_from(pl)
    eu = isempty(pl.up) ? nothing : pl.up[argmin(pl.up_d)]
    el = isempty(pl.lo) ? nothing : pl.lo[argmin(pl.lo_d)]
    return eu, el
end

function select_tracked(pl, alpha_prev, side::Symbol)
    isempty(pl.all) && error("spatial_payload: no finite eigenvalues")
    ap       = ComplexF64(alpha_prev)
    all_best = pl.all[argmin(abs.(ap .- pl.all))]
    cand     = side === :upper ? pl.up : pl.lo
    isempty(cand) && return all_best
    side_best = cand[argmin(abs.(ap .- cand))]
    if abs(ap - side_best) > TRACK_SIDE_SLACK * max(abs(ap - all_best), TRACK_MIN_STEP)
        global track_overrides += 1
        return all_best
    end
    return side_best
end

# Seed in the middle of the grid, then walk outward in both directions.
function track_branches(Lv)
    payloads = pmap(w -> spatial_payload(w, F, normals_F), Lv)
    s = argmin(abs.(real.(Lv) .- omega_r_base[25]))
    au = Vector{ComplexF64}(undef, length(Lv))
    al = Vector{ComplexF64}(undef, length(Lv))
    eu, el = dominant_from(payloads[s])
    (eu === nothing || el === nothing) &&
        error("track_branches: empty branch side at the seed omega = $(Lv[s])")
    au[s] = eu; al[s] = el
    global track_overrides = 0
    for j in (s + 1):length(Lv)
        au[j] = select_tracked(payloads[j], au[j-1], :upper)
        al[j] = select_tracked(payloads[j], al[j-1], :lower)
    end
    for j in (s - 1):-1:1
        au[j] = select_tracked(payloads[j], au[j+1], :upper)
        al[j] = select_tracked(payloads[j], al[j+1], :lower)
    end
    global track_jump_u = length(au) < 2 ? 0.0 :
                          maximum(abs.(au[2:end] .- au[1:end-1]))
    return au, al
end

function track_branches_init(Lv)
    payloads = pmap(w -> spatial_payload(w, F, normals_F), Lv)
    au = Vector{ComplexF64}(undef, length(Lv))
    al = Vector{ComplexF64}(undef, length(Lv))
    for j in eachindex(Lv)
        eu, el = dominant_from(payloads[j])
        (eu === nothing || el === nothing) &&
            error("track_branches_init: empty branch side at omega = $(Lv[j])")
        au[j] = eu; al[j] = el
    end
    return au, al
end

# The two branches must never merge -- that would mean the tracker has latched
# both sides onto the same mode.
function branches_disjoint(au, al; tol = 1e-8)
    dmin = Inf
    for a in au, b in al
        d = abs(a - b)
        d < dmin && (dmin = d)
    end
    return dmin >= tol, @sprintf("min upper/lower distance = %.3e", dmin)
end

branch_distance(au, al) = minimum(abs.(au .- al))
contour_distance(Fv, au, al) =
    minimum(minimum(abs.(f .- vcat(au, al))) for f in Fv)

# =============================================================================
# 4.  THE BARRIER
#
#   Every tracked branch point carries a charge; F is repelled by all of them.
#   phi = exp(zeta / (r^2 + eps)) - 1, integrated along each branch with the
#   local arc-length weight, so clustering L points near the pinch resolves the
#   same charge distribution better rather than changing it.
#
#   zeta sets the repulsion length sqrt(zeta), and it has to shrink with the gap
#   F must thread, or the barrier saturates: at zeta = 4e-4 and d = 1e-3 the
#   exponent is ~1600.  v5 adapted zeta through a running median of d_branch
#   with a 3 %/iteration rate limit, to protect against a noisy d_branch.  FIX 1
#   removes that noise, so the formula is used directly.
#
#   eps_alpha is tied to zeta so the exponent phi_F can ever see is bounded by
#   EXP_ARG_TARGET at every stage of the run (v4.9 CHANGE 10).
# =============================================================================
const S_ALPHA        = 2.0
const ZETA_REF       = 4.0e-4     # the hand-tuned value ...
const ZETA_D_REF     = 3.08e-2    # ... and the d_branch it was healthy at
const ZETA_MIN       = 1.0e-12
const EXP_ARG_TARGET = 10.0

global zeta_alpha    = ZETA_REF
global epsilon_alpha = ZETA_REF / EXP_ARG_TARGET

zeta_for(d) = (isfinite(d) && d > 0) ?
              clamp(ZETA_REF * (d / ZETA_D_REF)^2, ZETA_MIN, ZETA_REF) : ZETA_REF

function set_zeta!(d)
    global zeta_alpha    = zeta_for(d)
    global epsilon_alpha = max(zeta_alpha / EXP_ARG_TARGET, 1e-300)
    return zeta_alpha
end

const EXP_ARG_MAX = 400.0
expc(x) = exp(x < EXP_ARG_MAX ? x : EXP_ARG_MAX)

# Arc-length weights for one branch, computed once per iteration instead of
# once per node per call -- this is the whole of v5's Phi_F/d_d_*_Phi_F
# duplication, which was six near-identical 20-line functions.
function arc_weights(b)
    n = length(b)
    w = Vector{Float64}(undef, n)
    w[1] = abs(b[2] - b[1])
    for j in 2:(n - 1)
        w[j] = abs(0.5 * (b[j+1] - b[j-1]))
    end
    w[n] = abs(b[n] - b[n-1])
    return w
end

global charge_z = ComplexF64[]    # every branch point, both sides
global charge_w = Float64[]       # its arc-length weight

function set_charges!(au, al)
    global charge_z = vcat(au, al)
    global charge_w = vcat(arc_weights(au), arc_weights(al))
    return nothing
end

# Phi and its two gradients at one F point.  One pass over the charges.
function barrier_grad(af)
    gr = 0.0
    gi = 0.0
    @inbounds for k in eachindex(charge_z)
        dz = af - charge_z[k]
        r2 = abs2(dz)
        den = r2 + epsilon_alpha
        e   = expc(zeta_alpha / den)
        # d/dx exp(zeta/(r^2+eps)) = -zeta * 2x / (r^2+eps)^2 * exp(...)
        c = -zeta_alpha * S_ALPHA * e / den^2 * charge_w[k]
        gr += c * real(dz)
        gi += c * imag(dz)
    end
    return gr, gi
end

# =============================================================================
# 5.  FIX 2  --  the peak of Im(omega_F) must not extrapolate
#
#   With A = y_j - y_{j-1} >= 0 and B = y_j - y_{j+1} >= 0 at a discrete max:
#
#       vertex offset  t      = 0.5 (A - B) / (A + B)      -> |t| <= 0.5 ALWAYS
#       excess         y_pk - y_j = (A - B)^2 / (8 (A + B))
#                                 -> max(A,B)/8 when one-sided
#
#   so v5's only guard (|t| > 1) could never fire, and nothing bounded the
#   extrapolation.  omega_i is floored by this number, so an overshoot pushes L
#   back UP: replayed on v5's frames 300-1000 the bound sat ABOVE Im(omega_p) on
#   199 of 701 frames, which is the whole of the post-355 reversal.  With the
#   guard: 0 frames.
#
#   The parabola is right on a smooth F -- v4.6 added it because the sampled
#   maximum genuinely under-reports the continuous one, and it was good to
#   4.2e-8 on the smooth v4.5 frames.  It is only unbounded on a jagged one.
# =============================================================================
const PEAK_ASYM_MIN   = 0.05    # min(A,B)/max(A,B) below this -> no parabola
const PEAK_EXCESS_CAP = 0.5     # and never add more than this * min(A,B)

# (peak value, asymmetry).  asym = NaN means the fit was not attempted.
function peak_imag(wF)
    y = imag.(wF)
    j = argmax(y)
    (j == firstindex(y) || j == lastindex(y)) && return y[j], NaN
    A = y[j] - y[j-1]
    B = y[j] - y[j+1]
    (A < 0 || B < 0 || (A + B) <= 0) && return y[j], NaN
    asym = min(A, B) / max(A, B)
    asym < PEAK_ASYM_MIN && return y[j], asym
    excess = min((A - B)^2 / (8 * (A + B)), PEAK_EXCESS_CAP * min(A, B))
    ypk = y[j] + excess
    return (isfinite(ypk) ? ypk : y[j]), asym
end

# =============================================================================
# 6.  SMOOTHING  --  a triangular window whose width tracks d, with a floor
#
#   This is the one F change that the run data supports, and it is one clamp.
#
#   v4.7 ran a 15-node box filter on EVERY step and held f_ripple between 1e-4
#   and 6e-3 for all 501 iterations.  v5 replaced it with half-width
#   0.05*d_branch, which is an EXACT no-op once it drops below one cell
#   (width_average's window collapses to hi == lo -> out[j] = ys[j]).  That
#   happened at d_branch = 0.202, i.e. iteration 200, and f_ripple then grew
#   monotonically to 7.4e-2.  The first v5.1 removed the wide filter entirely
#   and never got going: at iteration 45 it sat at omega_i = -0.150 with
#   f_ripple = 1.2e-2, where v5 was at -0.293 with 4.8e-3.
#
#   v5's widths, in cells: 15.0 at iteration 10, 12.2 at 45, 8.0 at 75, 1.0 at
#   200, 0.16 at 250, 0.03 at 300.  The floor keeps all of that and stops the
#   collapse.  At exactly 2 cells the triangular window has weights
#   (0.5, 1, 0.5)/2 -- the 1-2-1 filter, transfer (1 + cos kh)/2, which
#   annihilates the 2h mode and is the mildest thing that does.
#
#   SMOOTH_HW_CELLS is the single knob: 7.5 recovers v4.9's box (which flattened
#   F's slope 180x over 1100 iterations), 1.5 is milder, below ~1.2 it is a
#   no-op again.  Watch f_ripple and the F slope and move it if either misbehaves.
# =============================================================================
const SMOOTH_FRAC     = 0.05
const SMOOTH_W_MAX    = 0.1515     # v4.9's 15-node window at h = 1.0101e-2
const SMOOTH_HW_CELLS = 2.0        # floor, in F grid cells

smooth_halfwidth(d) =
    clamp(isfinite(d) ? SMOOTH_FRAC * d : SMOOTH_W_MAX,
          SMOOTH_HW_CELLS * (alpha_r[2] - alpha_r[1]), SMOOTH_W_MAX)

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
                wm = 1.0 - abs(xs[m] - xs[j]) / half_w   # triangular, C0 at the edge
                wm <= 0 && continue
                acc += wm * ys[m]; wsum += wm
            end
            out[j] = wsum > 0 ? acc / wsum : ys[j]
        end
    end
    return out
end

# The chord residual: how far node j sits off the straight line through its two
# neighbours.  This is the 2h content in absolute units.  v4.7 held it in
# 1e-4 .. 6e-3 for its whole run; v5 ended at 7.4e-2.  Nothing acts on it.
function f_ripple_of(Fv)
    n = length(Fv)
    n < 3 && return NaN
    x = real.(Fv); y = imag.(Fv)
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

function linear_interp(xs, ys, xq)
    n = length(xs)
    xq <= xs[1]   && return ys[1]
    xq >= xs[end] && return ys[end]
    j = clamp(searchsortedlast(xs, xq), 1, n - 1)
    t = (xq - xs[j]) / (xs[j+1] - xs[j])
    return ys[j] + t * (ys[j+1] - ys[j])
end

# =============================================================================
# 7.  THE F STEP
#
#     rhs_j = ( y'_j * dPhi/dx - dPhi/dy ) / (1 + y'^2)      barrier
#             - DESC_FRAC * (theta/dt) * descent_j            descent on Im(omega)
#     y_j  += theta * tanh(dt * rhs_j / theta)                smooth limiter
#     y    := implicit sigma diffusion                        tension
#     y    := width_average(y, smooth_halfwidth(d))           section 6
#     accept if no branch crosses F and the peak falls within tolerance;
#     otherwise halve dt and retry, and apply the best of the five regardless.
#
#   THE ACCEPTANCE RULE IS v5's, ON PURPOSE.  The objective is a max over 100
#   nodes -- not differentiable -- and strict monotone descent on it is not
#   achievable with a fixed step direction: a step that lowers the peak at the
#   current argmax raises it somewhere else.  From v5's own peak_move, the peak
#   rose on 3.1 % of steps over iterations 1-100, 27.1 % over 200-300 and 48.2 %
#   over 300-600.  v5 took those steps (tolerance plus a best-trial fallback)
#   and descended.  The first v5.1 refused them and moved only on an escape
#   clause every eighth iteration, which is why it stalled at omega_i = -0.158.
#
#   DESC_FRAC IS FIXED.  The first v5.1 halved a gain whenever the peak rose;
#   because rises are normal, it decayed 0.30 -> 0.018 in 70 iterations and
#   switched off the one mechanism that had earned anything (v5's descent term
#   took omega_i to 4.28e-8 against v4.7's 3.75e-7, an 8.8x gain, by rotating
#   F's slope at the peak to -1.09 against v4.7's -0.06).
#
#   ONE THETA, NOT ONE PER NODE.  v5's per-node trust region was floored at the
#   global value, so at d = 3.015 a node sitting 1.75e-3 from a branch point was
#   still allowed to move 0.1508 -- 86x its own clearance.  The floor defeated
#   the per-node limit entirely, so it is gone; this is v4.9's single theta,
#   which is what reached f_ripple = 1.25e-5.
#
#   The descent direction is exact and free: omega is analytic, so
#   dv/dy = Re(domega/dalpha) by Cauchy-Riemann, and domega/dalpha comes from
#   the chord (omega_F[j+1]-omega_F[j-1])/(F[j+1]-F[j-1]) -- no extra solves.
#   It is weighted by a softmax over the top DESC_NW nodes, so only the
#   neighbourhood of the peak is driven.
#
#   The sigma diffusion is implicit: unconditionally stable, so sigma never
#   constrains dt.
# =============================================================================
const SIGMA              = 3e-5
const DT0                = 1e-3
const MOVE_FRAC          = 0.05    # node displacement cap, as a fraction of d_branch
const GLOBAL_MAX_MOVE    = 0.20
const DESC_FRAC          = 0.30    # share of theta the descent term may use
const DESC_NW            = 9       # nodes defining the softmax temperature
const ALPHA_MAX_ATTEMPTS = 5
const DESCENT_WINDOW     = 24      # iterations in the peak-descent-rate window (even)
const DESCENT_TOL_FRAC   = 0.10    # allowed rise, as a fraction of the recent rate
const CROSS_ESCAPE       = 8       # barren iterations before the crossing veto lifts
const CLEAR_FRAC         = 0.0     # required clearance, as a fraction of d_branch

global delta_t    = DT0
global barren_run = 0
global peak_hist  = Float64[]

function peak_push!(p)
    isfinite(p) || return nothing
    push!(peak_hist, p)
    while length(peak_hist) > DESCENT_WINDOW
        popfirst!(peak_hist)
    end
    return nothing
end

# Allowed rise in the peak for one accepted step.  Inf while the window fills,
# and -> 0 as the run stalls, so the rule tightens itself exactly when the
# random walk is the thing being fought.
function descent_tolerance()
    length(peak_hist) < DESCENT_WINDOW && return Inf
    h = DESCENT_WINDOW >> 1
    older = median(peak_hist[1:h])
    newer = median(peak_hist[h+1:end])
    rate  = (older - newer) / h            # > 0 while descending
    return isfinite(rate) ? DESCENT_TOL_FRAC * max(rate, 0.0) : Inf
end

# omega is analytic, so the chord derivative along F IS the full complex
# derivative -- validated against the analytic omega'' to 9 % over 17 nodes.
# Verbatim v5.  The ends take a ONE-SIDED difference rather than copying their
# neighbour -- the peak of Im(omega_F) sits at an endpoint of F on a good
# fraction of frames (83 of 631 in v5, all during its fast phase), so the
# endpoint derivative is not a formality.
function domega_dalpha(Fv, wF)
    n = length(Fv)
    g = zeros(ComplexF64, n)
    n < 2 && return g
    for j in 2:(n-1)
        dz = Fv[j+1] - Fv[j-1]
        g[j] = dz == 0 ? 0.0 + 0.0im : (wF[j+1] - wF[j-1]) / dz
    end
    dz1 = Fv[2] - Fv[1];    g[1] = dz1 == 0 ? 0.0+0.0im : (wF[2] - wF[1]) / dz1
    dzn = Fv[n] - Fv[n-1];  g[n] = dzn == 0 ? 0.0+0.0im : (wF[n] - wF[n-1]) / dzn
    for j in eachindex(g)
        isfinite(real(g[j])) && isfinite(imag(g[j])) || (g[j] = 0.0 + 0.0im)
    end
    return g
end

# Verbatim v5.  The exponent clamp matters when the top DESC_NW nodes are
# nearly tied -- which they often are -- because without it exp() underflows to
# exactly 0 and the weight vector collapses onto a single node.
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

# Implicit sigma*y'' on a possibly non-uniform grid, with zero-slope ends.
function diffusion_implicit(xs, ys, mscale, dt, sig)
    n = length(ys)
    (n < 3 || dt <= 0 || sig <= 0) && return copy(ys)
    dl = zeros(Float64, n-1); dg = ones(Float64, n); du = zeros(Float64, n-1)
    rhs = copy(ys)
    @inbounds for j in 2:(n - 1)
        hp = xs[j+1] - xs[j]; hm = xs[j] - xs[j-1]
        (hp > 0 && hm > 0) || continue
        a = 2.0 / (hm * (hp + hm))
        c = 2.0 / (hp * (hp + hm))
        f = dt * sig * mscale[j]
        dl[j-1] = -f * a
        du[j]   = -f * c
        dg[j]   = 1.0 + f * (a + c)
    end
    dg[1] = 1.0; du[1] = -1.0; rhs[1] = 0.0
    dg[n] = 1.0; dl[n-1] = -1.0; rhs[n] = 0.0
    out = try
        Tridiagonal(dl, dg, du) \ rhs
    catch
        copy(ys)
    end
    return all(isfinite, out) ? out : copy(ys)
end

# Signed clearance: how far the nearest branch point is on the correct side of
# F, and how many are on the wrong side.  Only pairs whose own gap is below
# gap_bound and whose alpha_r lies inside F's span are tested -- a far-field
# pair outside the span would be compared against a fictitious flat extension
# of F, which is the bug that deadlocked an earlier attempt at this.
function branch_clearance(Fv, au, al; gap_bound = 0.5)
    xs = real.(Fv); ys = imag.(Fv)
    lo = xs[1]; hi = xs[end]
    worst = Inf; ncross = 0
    for k in eachindex(au)
        abs(au[k] - al[k]) > gap_bound && continue
        (lo <= real(au[k]) <= hi) || continue
        (lo <= real(al[k]) <= hi) || continue
        cu = imag(au[k]) - linear_interp(xs, ys, real(au[k]))    # want > 0
        cl = linear_interp(xs, ys, real(al[k])) - imag(al[k])    # want > 0
        c = min(cu, cl)
        c < worst && (worst = c)
        c < 0 && (ncross += 1)
    end
    return worst, ncross
end

# Everything about the step that does NOT depend on dt: the slope of F, the
# barrier gradient at each node, and the descent direction.  All are functions
# of the accepted geometry, so they are computed once per iteration rather than
# once per trial -- the five trials differ only in step length.
function f_predirection()
    n = length(alpha_i)
    yr = zeros(Float64, n)
    gr = zeros(Float64, n)
    gi = zeros(Float64, n)
    for j in 2:(n - 1)
        hp = alpha_r[j+1] - alpha_r[j]
        hm = alpha_r[j]   - alpha_r[j-1]
        yr[j] = (-hp / (hm * (hp + hm))) * alpha_i[j-1] +
                ((hp - hm) / (hp * hm))  * alpha_i[j]   +
                ( hm / (hp * (hp + hm))) * alpha_i[j+1]
        gr[j], gi[j] = barrier_grad(F[j])
    end
    dwda  = domega_dalpha(F, omega_F)
    wdesc = descent_weights(omega_F)
    draw  = [wdesc[j] * real(dwda[j]) for j in 1:n]
    dsc   = maximum(abs, draw)
    return (yr = yr, gr = gr, gi = gi, draw = draw, dsc = dsc)
end

# One trial at step dt.  Returns the filtered trial and the limiter count.
#
# dt0 is the step length the attempt loop STARTED at and never changes; dt is
# the current, possibly halved one.  The descent term is scaled by theta/dt0 so
# that its contribution to the displacement, dt * (theta/dt0) * ..., shrinks
# linearly when the line search backs off -- exactly like the barrier term.
# Scaling it by theta/dt instead makes that contribution dt-independent, so
# halving dt shrinks only the barrier and the trial keeps overshooting.  That
# was the bug that left the first two v5.1 runs bouncing +/-5e-4 about a fixed
# omega_i while v5, over the same iterations, moved -4.8e-2 with 31 steps down
# and 5 up and accepted on attempt 1 on 112 of 141 frames.
function f_trial(dt, dt0, theta, hw, pre)
    n  = length(alpha_i)
    yt = copy(alpha_i)
    nlim = 0
    for j in 2:(n - 1)
        rhs = (pre.yr[j] * pre.gr[j] - pre.gi[j]) / (1.0 + pre.yr[j]^2)
        # Sign: moving y by dy changes Im(omega) by Re(domega/dalpha)*dy, so the
        # descent direction is minus that.
        if pre.dsc > 0
            rhs -= DESC_FRAC * (theta / dt0) * pre.draw[j] / pre.dsc
        end
        isnan(rhs) && (rhs = 0.0)
        ul = dt * rhs / theta
        abs(ul) > 2.0 && (nlim += 1)
        yt[j] = alpha_i[j] + theta * tanh(ul)
    end
    yt[1] = yt[2]; yt[n] = yt[n-1]

    metric = similar(yt)
    @inbounds for j in eachindex(yt)
        jm = max(j - 1, 1); jp = min(j + 1, n)
        dx = alpha_r[jp] - alpha_r[jm]
        sl = dx > 0 ? (yt[jp] - yt[jm]) / dx : 0.0
        metric[j] = 1.0 / (1.0 + sl^2)
    end
    yt = diffusion_implicit(alpha_r, yt, metric, dt, SIGMA)

    yt = width_average(alpha_r, yt, hw)
    yt[1] = yt[2]; yt[n] = yt[n-1]
    return yt, nlim
end

# The full step.  Returns a NamedTuple of everything the log and console want.
function f_step!(d_branch)
    n = length(alpha_i)
    peak_old, _ = peak_imag(omega_F)
    tol = descent_tolerance()

    theta = (isfinite(d_branch) && d_branch > 0) ?
            min(MOVE_FRAC * d_branch, GLOBAL_MAX_MOVE) : GLOBAL_MAX_MOVE
    theta = max(theta, 1e-14)      # tanh(dt*rhs/theta) must never divide by 0
    hw    = smooth_halfwidth(d_branch)

    pre = f_predirection()
    clear_now, ncross_now = branch_clearance(F, alpha_L_u, alpha_L_l)
    clear_req = CLEAR_FRAC * (isfinite(d_branch) ? d_branch : 0.0)
    veto = barren_run < CROSS_ESCAPE

    dt        = delta_t * clamp(d_branch / 0.2, 0.1, 1.0)   # branch slowdown
    dt0       = dt                                          # fixed reference
    best_peak = Inf
    best_y    = Float64[]
    best_w    = ComplexF64[]
    best_at   = 0
    nlim = 0; asym = NaN
    minclr = clear_now; ncross = ncross_now
    status = "none"
    # The one number the last run could not answer: dt was 1.25e-4 on every
    # frame, i.e. three halvings, so the full-size step was always rejected and
    # only 1/8 of it got through.  pm_full is attempt 1's peak_move, so the next
    # log says by how much the full step overshoots and whether it is the
    # barrier or the descent term doing it.
    pm_full = NaN

    for attempt in 1:ALPHA_MAX_ATTEMPTS
        yt, nl = f_trial(dt, dt0, theta, hw, pre)
        if !all(isfinite, yt)
            dt *= 0.5; status = "nonfinite/reduced"; continue
        end
        Ft = ComplexF64[alpha_r[j] + yt[j] * im for j in 1:n]

        # No branch point may end up on the wrong side of F.
        ct, nx = branch_clearance(Ft, alpha_L_u, alpha_L_l)
        if veto && !(ct >= clear_req || ct > clear_now)
            dt *= 0.5; status = "crossing/reduced"; continue
        end

        wt = contour_omega_F(Ft)
        pt, at = peak_imag(wt)
        attempt == 1 && isfinite(pt) && isfinite(peak_old) && (pm_full = pt - peak_old)
        if isfinite(pt) && pt < best_peak
            best_peak = pt; best_y = copy(yt); best_w = copy(wt); best_at = attempt
            nlim = nl; asym = at; minclr = ct; ncross = nx
        end
        if !isfinite(peak_old) || !isfinite(tol) ||
           (isfinite(pt) && pt <= peak_old + tol)
            status = "accepted"
            break
        end
        dt *= 0.5
        status = "descent/reduced"
    end

    # Always take the BEST trial seen, never a rejected one.  best_at is 0 only
    # when every attempt was non-finite or crossing, in which case the geometry
    # is left alone and the next iteration retries with the halved dt.
    if best_at == 0
        global barren_run += 1
        status = "discarded"
    else
        global barren_run = 0
        global alpha_i = copy(best_y)
        global omega_F = copy(best_w)
        global F       = contour_F()
        push_F!()
        peak_push!(best_peak)
    end

    pmove = (best_at == 0 || !isfinite(peak_old)) ? NaN : best_peak - peak_old
    return (status = status, attempts = max(best_at, 1), peak_move = pmove,
            peak_move_full = pm_full,
            n_limited = nlim, peak_asym = asym, min_clear = minclr,
            n_cross = ncross, dt = dt, theta = theta, smooth_hw = hw,
            descent_tol = tol)
end

# =============================================================================
# 8.  THE OMEGA STEP
#
#   omega_i is floored by F's own peak:  omega_i >= peak_imag(omega_F).  That
#   floor IS the coupling -- when L has no room to go lower, the only way down
#   is for F to deform.
#
#   v5 reached that floor through a 50-attempt gradient descent on a second
#   potential plus a 25-probe bisection.  Over the whole 1000-frame run the
#   result was always the same: omega_gap was exactly 0 on 966 frames and the
#   bisection never used more than one probe.  So v5.1 sets omega_i to the floor
#   directly and only bisects when that is not admissible.
# =============================================================================
const OMEGA_CLEARANCE  = 1e-9
const OMEGA_PROBES_MAX = 8
const OMEGA_BISECT_TOL = 1e-14

# Try one omega_i: rebuild L, re-track, check the branches stay separate.
function omega_trial(wi)
    Lt = contour_L_at(wi)
    au, al = track_branches(Lt)
    (all(isfinite, au) && all(isfinite, al)) ||
        return false, "non-finite branch", Lt, au, al
    ok, why = branches_disjoint(au, al)
    return ok, why, Lt, au, al
end

function omega_step!()
    pk, _ = peak_imag(omega_F)
    lb    = pk + OMEGA_CLEARANCE

    if omega_i <= lb
        # F's peak has risen above where L sits: L must come back up.  Nothing
        # to validate -- this is the state the previous iteration certified.
        global omega_i = lb
        global L = contour_L()
        global alpha_L_u, alpha_L_l = track_branches(L)
        return "pinned-up", 0, omega_i - lb
    end

    ok, why, Lt, au, al = omega_trial(lb)
    if ok
        global omega_i = lb
        global L = Lt; global alpha_L_u = au; global alpha_L_l = al
        return "pinned", 1, 0.0
    end

    # lb is not admissible: bisect between it and the height we already hold.
    a_lo = lb; b_hi = omega_i
    bw = omega_i; bL = L; bu = alpha_L_u; bl = alpha_L_l
    probes = 1
    for _ in 2:OMEGA_PROBES_MAX
        probes += 1
        wt = 0.5 * (a_lo + b_hi)
        ok2, _, Lt2, au2, al2 = omega_trial(wt)
        if ok2
            bw = wt; bL = Lt2; bu = au2; bl = al2; b_hi = wt
        else
            a_lo = wt
        end
        (b_hi - a_lo) < OMEGA_BISECT_TOL && break
    end
    if bw < omega_i
        global omega_i = bw
        global L = bL; global alpha_L_u = bu; global alpha_L_l = bl
        return "bisect", probes, omega_i - lb
    end
    @printf("[omega] no admissible height below %.12e (%s)\n", omega_i, why)
    return "stuck", probes, omega_i - lb
end

# =============================================================================
# 9.  ADAPTIVE L REFINEMENT
#
#   d_branch is read at a GRID POINT of L, so once the contour is close the
#   spacing of omega_r is what limits it.  Refine when the winning omega_r has
#   parked and d_min has stopped drifting, and place the new nodes at the vertex
#   of the d^4 parabola rather than by blind subdivision.
#
#   TWO FLOORS WERE RAISED FOR v5.1.  ADAPT_D_FLOOR was 2.0e-3 and ADAPT_MIN_H
#   was 5e-7: both were set to the eigensolver noise floor, which FIX 1 removes.
#   h = 5e-7 alone caps d at 2*sqrt(2*(h/4)/|omega''|) ~ 1.5e-3, so leaving them
#   would have made FIX 1 pointless.  ADAPT_MAX_POINTS still bounds the cost.
# =============================================================================
const ADAPT_ON             = true
const ADAPT_FACTOR         = 4
const ADAPT_HALF_CELLS     = 2
const ADAPT_STALL_WINDOW   = 20
const ADAPT_STALL_TOL      = 1e-2
const ADAPT_WINDOW_CELLS   = 2.0
const ADAPT_MAX_LEVEL      = 20
const ADAPT_MAX_POINTS     = 400
const ADAPT_MIN_H          = 1e-12     # was 5e-7, the old eigensolver floor
const ADAPT_D_FLOOR        = 0.0       # was 2.0e-3, likewise
const ADAPT_MIN_GAIN       = 0.90
const ADAPT_NOGAIN_MAX     = 2
const ADAPT_GRID_SHARE     = 0.90
const ADAPT_GRID_SHARE_MAX = 1.50
# Below this, the horizontal miss explains essentially none of d_branch: the L
# grid is already far finer than anything omega_i can use, and refining it
# again buys nothing.  Removing ADAPT_D_FLOOR and ADAPT_MIN_H took away the
# only stops this loop had, and the SCORING gate above cannot supply one --
# it needs grid_share in [0.90, 1.50] before adapt_nogain can even increment,
# so at grid_share = 0 the round is never scored and ADAPT never disables
# itself.  The v5.1 run reached level 15, N_L = 169, h = 5.2e-7, every round
# reporting "gain 1.000", where v4.7 was still at level 0 and v5 at level 2 at
# the same omega_i.  v5's six refinements happened at grid_share = 0.135,
# 0.082, 0.968, 0.936, 0.573 and 0.178, so 0.05 permits every one of them and
# blocks only the runaway.
const ADAPT_GS_MIN         = 0.05
const ADAPT_CLUSTER_HALF   = 2
const ADAPT_VERTEX_MAXMOVE = 2.0
const ADAPT_FIT_HALF       = 2

global adapt_level  = 0
global adapt_nogain = 0
global adapt_done   = false
global adapt_hist_d = Float64[]
global adapt_hist_w = Float64[]
global adapt_hist_h = Float64[]

function adapt_push!(d, w, h)
    push!(adapt_hist_d, d); push!(adapt_hist_w, w); push!(adapt_hist_h, h)
    while length(adapt_hist_d) > ADAPT_STALL_WINDOW
        popfirst!(adapt_hist_d); popfirst!(adapt_hist_w); popfirst!(adapt_hist_h)
    end
    return nothing
end

adapt_reset!() = (empty!(adapt_hist_d); empty!(adapt_hist_w); empty!(adapt_hist_h); nothing)

# Parked in omega_r (not in index -- once vertex placement has packed a dozen
# nodes into one old cell the winning index hops while barely moving), and
# d_min's drift inside the noise.
function adapt_stalled()
    n = length(adapt_hist_d)
    n < ADAPT_STALL_WINDOW && return false
    href = minimum(adapt_hist_h)
    (isfinite(href) && href > 0) || return false
    (maximum(adapt_hist_w) - minimum(adapt_hist_w)) > ADAPT_WINDOW_CELLS * href && return false
    m = abs(mean(adapt_hist_d))
    m <= 0.0 && return false
    h = n ÷ 2
    drift = abs(mean(adapt_hist_d[(h+1):end]) - mean(adapt_hist_d[1:h]))
    noise = 2.0 * std(adapt_hist_d) / sqrt(n)
    return drift <= ADAPT_STALL_TOL * m + noise
end

function local_h(wr, i)
    n = length(wr)
    n < 2 && return Inf
    i == 1 && return wr[2] - wr[1]
    i == n && return wr[n] - wr[n-1]
    return min(wr[i+1] - wr[i], wr[i] - wr[i-1])
end

# Fit d^4 against omega_r near the argmin.  d ~ 2 sqrt(2|w-w_p|/|w''|) makes d^4
# a parabola in omega_r, so the vertex locates the pinch in omega_r and the
# curvature gives |omega''|.  Returns (vertex, |omega''|).
function pinch_fit(wr, d, i; n = ADAPT_FIT_HALF)
    lo = max(firstindex(d), i - n); hi = min(lastindex(d), i + n)
    (i - lo < 2 || hi - i < 2) && return (NaN, NaN)
    x = wr[lo:hi] .- wr[i]
    y = d[lo:hi] .^ 4
    (all(isfinite, x) && all(isfinite, y)) || return (NaN, NaN)
    xs = maximum(abs, x); ys = maximum(abs, y)
    (isfinite(xs) && isfinite(ys) && xs > 0 && ys > 0) || return (NaN, NaN)
    c = hcat((x ./ xs) .^ 2, x ./ xs, ones(length(x))) \ (y ./ ys)
    (all(isfinite, c) && c[1] > 0) || return (NaN, NaN)
    w_pr = wr[i] - xs * c[2] / (2 * c[1])
    curv = ys * c[1] / xs^2
    return (w_pr, curv > 0 ? 8.0 / sqrt(curv) : NaN)
end

function refine_omega_r(wr, i)
    n = length(wr)
    lo = max(1, i - ADAPT_HALF_CELLS); hi = min(n, i + ADAPT_HALF_CELLS)
    hi <= lo && return copy(wr)
    out = Float64[]
    append!(out, wr[1:(lo-1)])
    for c in lo:(hi-1)
        seg = collect(range(wr[c], wr[c+1], length = ADAPT_FACTOR + 1))
        append!(out, seg[1:(end-1)])
    end
    push!(out, wr[hi])
    append!(out, wr[(hi+1):n])
    return out
end

function insert_cluster(wr, w_c, h)
    (isfinite(w_c) && isfinite(h) && h > 0) || return copy(wr)
    (w_c <= wr[1] || w_c >= wr[end]) && return copy(wr)
    out = collect(wr)
    for j in -ADAPT_CLUSTER_HALF:ADAPT_CLUSTER_HALF
        w = w_c + j * h
        (w <= wr[1] || w >= wr[end]) && continue
        minimum(abs.(out .- w)) < h / 8 && continue
        push!(out, w)
    end
    sort!(out)
    return out
end

function push_omega_r!()
    wr = deepcopy(omega_r)
    @sync for pid in workers()
        @async remotecall_wait(Core.eval, pid, Main, :(omega_r = $wr))
    end
    return nothing
end

# =============================================================================
# 10.  LOG AND RESUME
#
#   One file.  v5 kept a separate checkpoint because the log did not store the
#   adapted omega_r; v5.1 stores omega_r and alpha_r in every entry (about 2 kB
#   against the 24 kB the five complex arrays already cost), so the log alone is
#   enough to restart from.
#
#   SAVE_EVERY exists because the log is rewritten whole on each save and it
#   reaches ~35 MB over 1000 frames -- v5 did that every iteration, which is
#   35 GB of I/O over a run and is what produced the 39 MB write failure.  A
#   resume redoes at most SAVE_EVERY-1 iterations.
# =============================================================================
const LOGFILE     = "contour_iteration_v5.1.json"
const RESUME      = true
const SAVE_EVERY  = 5
const ITER_TARGET = 2000

cvec(v)  = [Dict("re" => real(x), "im" => imag(x)) for x in v]
jnum(x)  = (x isa Real && isfinite(x)) ? x : nothing
uncvec(a) = ComplexF64[complex(x["re"], x["im"]) for x in a]

global log_array = Any[]
global iteration_step = 1

# The previous attempt at a v5.1 died here -- SystemError on opening the .tmp,
# not on the rename -- and took a 700-iteration run with it.  So the temp write
# is guarded too, and any failure at any stage falls back to writing the target
# directly rather than throwing.
# The log writer, and the one rule it now obeys: NEVER truncate the existing
# log.  The previous version fell back to open(path, "w") when the rename kept
# failing, and that rewrite was interrupted at exactly 1 MiB -- iterations 1-34
# survived, everything after was lost.
#
# What actually fails on this machine is Julia's mv unlinking the .tmp
# (IOError: unlink("...json.tmp"): resource busy or locked, EBUSY) -- a virus
# scanner, indexer or cloud-sync agent grabbing the file the moment it is
# written.  A fixed .tmp name makes that permanent, because the held file is
# still there on the next save.  So: a fresh temp name every attempt, and if
# every attempt fails, skip this save and keep the last good log.  The worst
# case is losing SAVE_EVERY iterations of history, never the file.
function atomic_write(path, str; tries = 6)
    for attempt in 1:tries
        tmp = string(path, ".tmp", getpid(), "_", attempt)
        wrote = try
            open(tmp, "w") do io
                write(io, str)
            end
            true
        catch err
            false
        end
        if !wrote
            try; rm(tmp; force = true); catch; end
            sleep(0.25 * attempt)
            continue
        end
        try
            mv(tmp, path; force = true)
            return nothing
        catch err
            try; rm(tmp; force = true); catch; end
            attempt == tries && begin
                @printf("[io] could not replace %s after %d attempts (%s) -- THIS SAVE SKIPPED, the previous log is intact\n",
                        path, tries, sprint(showerror, err))
                flush(stdout)
            end
            sleep(0.25 * attempt)
        end
    end
    return nothing
end

save_log() = atomic_write(LOGFILE, JSON.json(log_array))

# Every entry carries the same key set: MATLAB's jsondecode only builds a struct
# array when the field names match across the whole array.
function log_entry(; iteration, d_branch, d_contour, dist_u, dist_l, h_local,
                     omega_gap, omega_status, omega_probes, alpha_status,
                     attempts, peak_move, peak_move_full, peak_asym, n_limited, f_ripple,
                     min_clear, n_cross, dt_used, theta, smooth_hw,
                     descent_tol, zeta, grid_share, pinch_miss)
    return Dict(
        "iteration"    => iteration,
        "L"            => cvec(L),
        "F"            => cvec(F),
        "alpha_L_u"    => cvec(alpha_L_u),
        "alpha_L_l"    => cvec(alpha_L_l),
        "omega_F"      => cvec(omega_F),
        "omega_r"      => collect(omega_r),
        "alpha_r"      => collect(alpha_r),
        "omega_i"      => omega_i,
        "n_L"          => length(omega_r),
        "n_F"          => length(alpha_r),
        "adapt_level"  => adapt_level,
        "adapt_done"   => adapt_done,
        "delta_t"      => jnum(delta_t),
        "d_branch"     => jnum(d_branch),
        "d_contour"    => jnum(d_contour),
        "dist_u"       => jnum(dist_u),
        "dist_l"       => jnum(dist_l),
        "h_local"      => jnum(h_local),
        "omega_gap"    => jnum(omega_gap),
        "omega_status" => omega_status,
        "omega_probes" => omega_probes,
        "alpha_status" => alpha_status,
        "attempts"     => attempts,
        "peak_move"    => jnum(peak_move),
        "peak_move_full" => jnum(peak_move_full),
        "peak_asym"    => jnum(peak_asym),
        "n_limited"    => n_limited,
        "f_ripple"     => jnum(f_ripple),
        "min_clear"    => jnum(min_clear),
        "n_cross"      => n_cross,
        "dt_used"      => jnum(dt_used),
        "theta"        => jnum(theta),
        "smooth_hw"    => jnum(smooth_hw),
        "descent_tol"  => jnum(descent_tol),
        "zeta_alpha"   => jnum(zeta),
        "grid_share"   => jnum(grid_share),
        "pinch_miss"   => jnum(pinch_miss),
        "track_jump_u" => jnum(track_jump_u),
        "track_overrides" => track_overrides,
    )
end

function try_resume()
    isfile(LOGFILE) || return false
    arr = try
        JSON.parsefile(LOGFILE)
    catch err
        @printf("RESUME: %s is unreadable (%s) -- starting fresh\n",
                LOGFILE, sprint(showerror, err))
        return false
    end
    (arr isa Vector && !isempty(arr)) || return false
    e = arr[end]
    # "descent_tol" and "smooth_hw" only exist in this version's schema, so a
    # log left by the reverted first v5.1 cannot be resumed from -- it would
    # hand back that run's stalled geometry.  It gets archived instead.
    for key in ("iteration", "omega_r", "alpha_r", "F", "omega_i",
                "adapt_level", "adapt_done", "descent_tol", "smooth_hw")
        haskey(e, key) || begin
            @printf("RESUME: last entry has no \"%s\" -- starting fresh\n", key)
            return false
        end
    end
    global log_array      = arr
    global iteration_step = Int(e["iteration"]) + 1
    global omega_r        = Float64.(e["omega_r"])
    global alpha_r        = Float64.(e["alpha_r"])
    global alpha_i        = imag.(uncvec(e["F"]))
    global omega_i        = Float64(e["omega_i"])
    global adapt_level    = Int(e["adapt_level"])
    global adapt_done     = Bool(e["adapt_done"])
    global delta_t        = e["delta_t"] === nothing ? DT0 : Float64(e["delta_t"])
    global F   = contour_F()
    global L   = contour_L()
    push_F!()
    push_omega_r!()
    global omega_F = contour_omega_F(F)
    global alpha_L_u, alpha_L_l = track_branches(L)
    set_charges!(alpha_L_u, alpha_L_l)
    adapt_reset!()
    @printf("RESUMED at iteration %d: omega_i = %.12e, N_L = %d, N_F = %d, adapt level %d\n",
            iteration_step, omega_i, length(omega_r), length(alpha_r), adapt_level)
    flush(stdout)
    return true
end

# =============================================================================
# 11.  INITIALISE
#
# try_resume rebuilds everything from the log, so the fresh-start work below is
# only done when there is nothing to resume from -- it costs 200 eigensolves.
# =============================================================================
resumed = RESUME ? try_resume() : false

if !resumed
    if isfile(LOGFILE)
        stamp = string(round(Int, time()))
        mv(LOGFILE, LOGFILE * ".old." * stamp; force = true)
        println("RESUME off or no usable log: existing log moved to ",
                LOGFILE * ".old." * stamp)
    end
    global log_array = Any[]
    global iteration_step = 1
    global F = contour_F()
    global L = contour_L()
    push_F!()
    push_omega_r!()
    global omega_F = contour_omega_F(F)
    global alpha_L_u, alpha_L_l = track_branches_init(L)
    set_charges!(alpha_L_u, alpha_L_l)
    push!(log_array, log_entry(
        iteration = 1, d_branch = NaN, d_contour = NaN, dist_u = NaN, dist_l = NaN,
        h_local = NaN, omega_gap = NaN, omega_status = "initial", omega_probes = 0,
        alpha_status = "initial", attempts = 0, peak_move = NaN,
        peak_move_full = NaN, peak_asym = NaN,
        n_limited = 0, f_ripple = f_ripple_of(F), min_clear = NaN,
        n_cross = 0, dt_used = NaN, theta = NaN,
        smooth_hw = smooth_halfwidth(NaN), descent_tol = NaN,
        zeta = zeta_alpha, grid_share = NaN, pinch_miss = NaN))
    save_log()
    global iteration_step = 2
end

@printf("\nv5.1 ready.  N_F = %d, N_L = %d, num_modes = %d, Re = %.0f\n",
        length(alpha_r), length(omega_r), num_modes, Re)
@printf("starting at iteration %d, omega_i = %.12e\n\n", iteration_step, omega_i)
flush(stdout)

# =============================================================================
# 12.  MAIN LOOP
# =============================================================================
while iteration_step <= ITER_TARGET
    k = iteration_step

    # ---- omega: press L down onto F's peak --------------------------------
    ostat, oprobes, ogap = omega_step!()
    if ostat == "stuck"
        println("STOP: no admissible omega height.")
        break
    end
    push_F!()
    set_charges!(alpha_L_u, alpha_L_l)

    # ---- geometry ----------------------------------------------------------
    d_vec     = abs.(alpha_L_u .- alpha_L_l)
    i_pinch   = argmin(d_vec)
    d_branch  = d_vec[i_pinch]
    h_local   = local_h(omega_r, i_pinch)
    d_contour = contour_distance(F, alpha_L_u, alpha_L_l)
    dist_u    = minimum(minimum(abs.(f .- alpha_L_u)) for f in F)
    dist_l    = minimum(minimum(abs.(f .- alpha_L_l)) for f in F)

    if !isfinite(d_branch) || !isfinite(d_contour)
        @printf("[k=%d] STOP: non-finite distance (d_branch=%.3e d_contour=%.3e)\n",
                k, d_branch, d_contour)
        break
    end

    # How much of d_branch the horizontal miss in omega_r explains.  ~1 means
    # the L grid is the limiter; << 1 means omega_i still is.
    w_pr_fit, w2_fit = pinch_fit(omega_r, d_vec, i_pinch)
    pinch_miss = abs(omega_r[i_pinch] - w_pr_fit)
    d_horiz    = (isfinite(pinch_miss) && isfinite(w2_fit) && w2_fit > 0) ?
                 2 * sqrt(2 * pinch_miss / w2_fit) : NaN
    grid_share = (isfinite(d_horiz) && d_branch > 0) ? d_horiz / d_branch : NaN

    set_zeta!(d_branch)

    # ---- F: deform it downwards -------------------------------------------
    st = f_step!(d_branch)
    f_rip = f_ripple_of(F)

    @printf("[%04d] dUL=%9.3e dUF=%9.3e dLF=%9.3e | w=%.10f gap=%9.2e %s(%d) | dt=%8.2e z=%8.2e | NL=%3d L%d h=%8.2e gs=%5.2f | hw/h=%6.2f | dpk=%+9.2e(1st%+9.2e a%d) tol=%8.2e asym=%6.3f rip=%8.2e lim=%3d | clr=%+9.2e nx=%2d | %s\n",
            k, d_branch, dist_u, dist_l,
            omega_i, ogap, ostat, oprobes,
            st.dt, zeta_alpha,
            length(omega_r), adapt_level, h_local, grid_share,
            st.smooth_hw / (alpha_r[2] - alpha_r[1]),
            st.peak_move, st.peak_move_full, st.attempts,
            st.descent_tol, st.peak_asym, f_rip, st.n_limited,
            st.min_clear, st.n_cross, st.status)
    flush(stdout)

    push!(log_array, log_entry(
        iteration = k, d_branch = d_branch, d_contour = d_contour,
        dist_u = dist_u, dist_l = dist_l, h_local = h_local,
        omega_gap = ogap, omega_status = ostat, omega_probes = oprobes,
        alpha_status = st.status, attempts = st.attempts,
        peak_move = st.peak_move, peak_move_full = st.peak_move_full,
        peak_asym = st.peak_asym,
        n_limited = st.n_limited, f_ripple = f_rip,
        min_clear = st.min_clear, n_cross = st.n_cross, dt_used = st.dt,
        theta = st.theta, smooth_hw = st.smooth_hw,
        descent_tol = st.descent_tol, zeta = zeta_alpha,
        grid_share = grid_share, pinch_miss = pinch_miss))

    if k % SAVE_EVERY == 0 || k == ITER_TARGET
        save_log()
    end

    # ---- adaptive L refinement --------------------------------------------
    # Last, after the log entry, so every entry describes the grid its numbers
    # were computed on.
    if ADAPT_ON && !adapt_done && st.status != "discarded"
        adapt_push!(d_branch, omega_r[i_pinch], h_local)
        if adapt_stalled() && isfinite(grid_share) && grid_share < ADAPT_GS_MIN
            @printf("[k=%d] adapt: grid_share %.3f < %.2f -- the grid is not the limiter, not refining (h=%.2e, N_L=%d)\n",
                    k, grid_share, ADAPT_GS_MIN, h_local, length(omega_r))
            flush(stdout)
            adapt_reset!()
        elseif adapt_stalled()
            h_new = h_local / ADAPT_FACTOR
            n_add = 2 * ADAPT_HALF_CELLS * (ADAPT_FACTOR - 1)
            if adapt_level >= ADAPT_MAX_LEVEL
                @printf("[k=%d] ADAPT OFF: level cap %d\n", k, ADAPT_MAX_LEVEL)
                global adapt_done = true
            elseif d_branch <= ADAPT_D_FLOOR
                @printf("[k=%d] ADAPT OFF: d_branch %.3e at the floor\n", k, d_branch)
                global adapt_done = true
            elseif h_new < ADAPT_MIN_H
                @printf("[k=%d] ADAPT OFF: spacing floor (h=%.2e)\n", k, h_local)
                global adapt_done = true
            elseif length(omega_r) + n_add > ADAPT_MAX_POINTS
                @printf("[k=%d] ADAPT OFF: point budget (N_L=%d)\n", k, length(omega_r))
                global adapt_done = true
            elseif i_pinch <= 1 || i_pinch >= length(omega_r)
                @printf("[k=%d] adapt: argmin at a grid endpoint -- not refining\n", k)
                adapt_reset!()
            else
                grid_limited = isfinite(grid_share) &&
                               ADAPT_GRID_SHARE <= grid_share <= ADAPT_GRID_SHARE_MAX
                h_cluster = max(h_new, 2 * pinch_miss)
                vertex_ok = isfinite(w_pr_fit) && isfinite(pinch_miss) &&
                            pinch_miss <= ADAPT_VERTEX_MAXMOVE * h_local
                n_before = length(omega_r)
                placement = "subdivide"
                if vertex_ok
                    try_r = insert_cluster(omega_r, w_pr_fit, h_cluster)
                    if length(try_r) > n_before
                        global omega_r = try_r
                        placement = "cluster"
                    else
                        global omega_r = refine_omega_r(omega_r, i_pinch)
                    end
                else
                    global omega_r = refine_omega_r(omega_r, i_pinch)
                end
                global adapt_level += 1
                push_omega_r!()
                global L = contour_L()
                global alpha_L_u, alpha_L_l = track_branches(L)
                set_charges!(alpha_L_u, alpha_L_l)
                adapt_reset!()
                d_after = branch_distance(alpha_L_u, alpha_L_l)
                gain = d_branch > 0 ? d_after / d_branch : NaN
                @printf("[k=%d] ADAPT -> level %d (%s) | h %.3e -> %.3e | N_L %d -> %d | gs=%.2f | d %.4e -> %.4e (gain %.3f)\n",
                        k, adapt_level, placement, h_local,
                        placement == "cluster" ? h_cluster : h_new,
                        n_before, length(omega_r), grid_share, d_branch, d_after, gain)
                flush(stdout)
                # Only score a round when the grid is what is being judged:
                # while omega_i is still the limiter the gain is capped near 1
                # whatever the grid does.
                if grid_limited
                    if isfinite(gain) && gain > ADAPT_MIN_GAIN
                        global adapt_nogain += 1
                    else
                        global adapt_nogain = 0
                    end
                    if adapt_nogain >= ADAPT_NOGAIN_MAX
                        @printf("[k=%d] ADAPT OFF: %d rounds gained under %.0f%% -- the L grid is no longer the limiter\n",
                                k, adapt_nogain, 100 * (1 - ADAPT_MIN_GAIN))
                        global adapt_done = true
                    end
                end
            end
        end
    end

    global iteration_step += 1
end

save_log()
@printf("\nFINISHED at iteration %d.  %d entries in %s\n",
        iteration_step - 1, length(log_array), LOGFILE)
flush(stdout)
