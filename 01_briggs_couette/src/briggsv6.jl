###############################################################################
#
#   BRIGGS v6  --  absolute/convective instability of plane Couette flow
#
#   The method, in four lines.  Im omega_p = min over admissible F of
#   max over alpha in F of Im omega(alpha).  L is horizontal and is floored by
#   F's own peak; F deforms to lower that peak; repeat.  When L has no room to
#   go lower, the only way down is for F to move.  That coupling is the method
#   and it is unchanged since v4.0.
#
#   WHAT IS NEW RELATIVE TO v5.1, AND WHY
#
#   1  The TEMPORAL pencil is equilibrated as well as the spatial one.  A
#      carries D4/Re at ~1e14 against B11 at ~1e9 and the imbalance grows like
#      n^8, so at num_modes = 150 omega(alpha) was good to only ~7e-9.  omega_i
#      is pinned to max Im omega_F, so that error was a floor on the branch gap
#      at d ~ 3.5e-4 no matter what the contour did.  Equilibrated: ~2e-14.
#      Same function, same proof: R*P*S has the same roots for diagonal R, S.
#
#   2  num_modes = 80, not 150.  Unequilibrated the error GREW with N
#      (1.3e-10 at 40, 6.7e-9 at 150, 6.5e-7 at 200); equilibrated it is flat
#      at ~2e-14 for N >= 60 and the solve is several times cheaper, which is
#      what makes the endgame L grid affordable.  Justify it yourself with
#      check_modes() -- it compares N against N and never looks at an answer.
#
#   3  The alpha_r window extends itself.  During the descent the peak of
#      Im omega_F sat at F's RIGHT ENDPOINT, so the ceiling being enforced was
#      the ceiling of a truncated contour, not the real one.  The left end
#      needs nothing: Im omega(-conj(alpha)) = Im omega(alpha) exactly, so
#      alpha_i is even in alpha_r and zero slope at 0 is the symmetry
#      condition, not a guess.
#
#   4  F grades itself.  One rule: the spacing at the peak must be below
#      F_H_FRAC * d_branch.  The fine block is centred on the run's own throat
#      estimate -- the midpoint of the tracked branch pair at argmin|a_u - a_l|.
#
#   NOTHING HERE KNOWS THE PINCH, and the rule is stronger than "no pinch
#   numbers in the source": no constant below was chosen by comparing against a
#   known answer.  Where v5.1's ADAPT_FIT_HALF was picked that way, v6 fits the
#   vertex with a 2-node and a 3-node stencil and uses their disagreement as its
#   own error bar.  compute_pinch2.py is the referee: postprocessing only.
#
#   ONE FILE, APPENDED TO.  The log is JSONL, one line per iteration, opened
#   in append mode.  The last complete line IS the checkpoint, so a crash
#   costs one iteration instead of the run -- which is what happened twice.
#   No separate checkpoint, no atomic rename, no truncation logic.
#
#   RUN:     julia briggsv6.jl          (delete/rename the log for a clean start)
#   OUTPUT:  contour_iteration_v6.jsonl
#   READ:    entries = [JSON.parse(l) for l in eachline(f) if !isempty(strip(l))]
#
#   WATCH, in order.  Each one falsifies a specific change.
#     w      omega_i.  v5's best was 4.28e-8 from Im(omega_p) and it then drifted
#            out to 7.26e-5.  This should go past 4e-8 and stay.  If it sticks
#            near 4e-8, change 1 is wrong.
#     d      the branch gap.  Was pinned near 2.4e-3 in v4.7 and v5 because the
#            solve could not read below it.  It should now fall straight through.
#     NF/hF  F's node count and the spacing at the peak over d.  hF/d should
#            come down to F_H_FRAC and stay; NF should grow in steps.
#     res    3-node vs 5-node disagreement at the peak, as a fraction of the
#            local curvature.  If it stays large while hF/d is at target, the
#            resolution rule is too loose -- lower F_H_FRAC.
#     a      attempts.  Living at MAX_ATTEMPTS means the line search is
#            throttling itself; the tolerance floor is meant to stop that.
#     rip    ripple.  v4.7 held 1e-4..6e-3 for 501 iterations; v5 ended at 7e-2.
#
#   SMOKE TEST BEFORE A LONG RUN:  set ITER_TARGET = 3, run it, confirm three
#   lines appear in the log and the console.  Two patches in this project's
#   history were never parsed by Julia before a long run.
#
###############################################################################

using Distributed, JSON, Statistics, Printf
addprocs(5)
@everywhere using LinearAlgebra, Statistics

# =============================================================================
# 1.  FLOW, DISCRETISATION, AND THE TWO PENCILS
#
#   Chebyshev in coefficient space; D0 maps coefficients to collocation values.
#   U = y so U'' = 0 and every U'' term is dropped.  The -200im boundary rows
#   are the usual trick: they push the four boundary-condition eigenvalues to
#   omega = -200i, far from anything physical.
# =============================================================================
@everywhere begin
    const RE        = 2000.0
    const BETA      = 0.0 + 0.0im
    const VG        = 0.0 + 0.0im
    const NUM_MODES = 80
    const Y0        = 0.0
    const Y1        = 1.0
end

@everywhere begin
    yc = [cos((j - 1) * pi / (NUM_MODES - 1)) for j = 1:NUM_MODES]
    yp = (Y0 + Y1) / 2 .- yc .* ((Y1 - Y0) / 2)

    const D0 = zeros(Float64, NUM_MODES, NUM_MODES)
    const D1 = zeros(Float64, NUM_MODES, NUM_MODES)
    const D2 = zeros(Float64, NUM_MODES, NUM_MODES)
    const D3 = zeros(Float64, NUM_MODES, NUM_MODES)
    const D4 = zeros(Float64, NUM_MODES, NUM_MODES)
    for j = 1:NUM_MODES
        D0[:, j] .= cos.((j - 1) .* acos.(yc))
    end
    D1[:, 2] = D0[:, 1];  D1[:, 3] = 4 * D0[:, 2]
    D2[:, 3] = 4 * D0[:, 1]
    for j = 4:NUM_MODES
        D1[:, j] .= 2 * (j - 1) * D0[:, j - 1] + (j - 1) * D1[:, j - 2] / (j - 3)
        D2[:, j] .= 2 * (j - 1) * D1[:, j - 1] + (j - 1) * D2[:, j - 2] / (j - 3)
        D3[:, j] .= 2 * (j - 1) * D2[:, j - 1] + (j - 1) * D3[:, j - 2] / (j - 3)
        D4[:, j] .= 2 * (j - 1) * D3[:, j - 1] + (j - 1) * D4[:, j - 2] / (j - 3)
    end
    smap = -(Y1 - Y0) / 2
    D1 ./= smap
    D2 ./= smap^2
    D3 ./= smap^3
    D4 ./= smap^4
    const UM = yp * ones(Float64, 1, NUM_MODES)     # UM[i,j] = U(y_i)
end

# Ruiz-style two-sided equilibration.  det(R*P*S) = det(R)det(S)det(P), so no
# root moves; eigenvectors would transform as S^-1 v and nothing here reads
# them.  No tuning constant: the scaling is read off |A|+|B| itself.
@everywhere function equilibrate(mats...; iters::Int = 12)
    n = size(mats[1], 1)
    W = zeros(Float64, n, n)
    for M in mats
        W .+= abs.(M)
    end
    R = ones(Float64, n)
    S = ones(Float64, n)
    for _ in 1:iters
        r = sqrt.(max.(vec(maximum(W; dims = 2)), 1e-300))       # row maxima
        R ./= r
        W ./= r
        c = sqrt.(max.(vec(maximum(W; dims = 1)), 1e-300))       # column maxima
        S ./= c
        W ./= transpose(c)
    end
    return R, S
end

# omega given alpha (temporal).  CHANGE 1: equilibrated.
@everywhere function omega_spectrum(alpha)
    a2 = alpha^2 + BETA^2
    A11 = (-im * alpha) .* UM .* D2 .+ (im * alpha * a2) .* UM .* D0 .+
          (1 / RE) .* D4 .- (2 / RE * a2) .* D2 .+ (1 / RE * a2^2) .* D0 .+
          (alpha * VG) .* D0
    A = [(-200im) .* [D0[1:1, :]; D1[1:1, :]];
         A11[3:NUM_MODES-2, :];
         (-200im) .* [D1[NUM_MODES:NUM_MODES, :]; D0[NUM_MODES:NUM_MODES, :]]]
    B11 = (-im) .* D2 .+ (im * a2) .* D0
    B = [[D0[1:1, :]; D1[1:1, :]];
         B11[3:NUM_MODES-2, :];
         [D1[NUM_MODES:NUM_MODES, :]; D0[NUM_MODES:NUM_MODES, :]]]
    R, S = equilibrate(A, B)
    return eigvals(R .* A .* transpose(S), R .* B .* transpose(S))
end

# the alpha spectrum given omega (spatial, a quadratic pencil linearised in
# companion form).  Equilibrated before linearising -- afterwards is too late,
# because the companion structure ties the two block rows together and LAPACK's
# own balancing cannot undo it.
@everywhere function alpha_spectrum(omega)
    A11 = (-2im * omega) .* D1 .- (4 / RE) .* D3 .+ (4 / RE * BETA^2) .* D1 .-
          im .* UM .* D2 .+ (im * BETA^2) .* UM .* D0 .-
          (im * VG) .* D2 .+ (im * VG * BETA^2) .* D0
    A12 = (im * omega) .* D2 .- (im * omega * BETA^2) .* D0 .+ (1 / RE) .* D4 .-
          (2 / RE * BETA^2) .* D2 .+ (1 / RE * BETA^4) .* D0
    Z2  = zeros(ComplexF64, 2, NUM_MODES)
    A11 = [Z2; A11[3:NUM_MODES-2, :]; Z2]
    A12 = [(-200im) .* [D0[1:1, :]; D1[1:1, :]];
           A12[3:NUM_MODES-2, :];
           (-200im) .* [D1[NUM_MODES:NUM_MODES, :]; D0[NUM_MODES:NUM_MODES, :]]]
    B11 = (-4 / RE) .* D2 .- 2im .* UM .* D1 .+ (2im * VG) .* D1
    B11 = [Z2; B11[3:NUM_MODES-2, :]; Z2]

    R, S = equilibrate(B11, A11, A12)
    A11 = R .* A11 .* transpose(S)
    A12 = R .* A12 .* transpose(S)
    B11 = R .* B11 .* transpose(S)

    Z  = zeros(ComplexF64, NUM_MODES, NUM_MODES)
    Id = Matrix{ComplexF64}(I, NUM_MODES, NUM_MODES)
    return eigvals([A11 A12; Id Z], [B11 Z; Z Id])
end

@everywhere finite_only(ev) =
    ev[[isfinite(real(e)) && isfinite(imag(e)) for e in ev]]

@everywhere function omega_of_alpha(alpha)
    ev = finite_only(omega_spectrum(alpha))
    return ev[argmax(imag.(ev))]
end

# Run this ONCE, offline, to justify NUM_MODES without looking at any answer:
# it compares one discretisation against another, nothing else.  Pick the
# smallest N whose column has collapsed to round-off.  Sample alpha over the
# region F actually occupies, INCLUDING the extended window.
function check_modes(alphas = ComplexF64[0.05 - 1.9im, 0.5 - 2.25im,
                                         1.0 - 2.25im, 1.8 - 2.25im,
                                         0.3 - 3.0im,  0.9 - 3.05im])
    @printf("compare each alpha's omega against the NUM_MODES = %d value in this file\n",
            NUM_MODES)
    for a in alphas
        w = omega_of_alpha(a)
        @printf("  alpha = %+.4f%+.4fi   omega = %+.12f%+.12fi\n",
                real(a), imag(a), real(w), imag(w))
    end
    println("re-run with NUM_MODES set to 40, 60, 100, 150 and compare the columns.")
    return nothing
end

# =============================================================================
# 2.  THE TWO CONTOURS
#
#   F is a graph over alpha_r:  F[j] = alpha_r[j] + i*alpha_i[j].
#   L is horizontal:            L[j] = omega_r[j] + i*omega_i.
#   Both node vectors are refined at run time, so neither is a range.
# =============================================================================
const N_F0         = 100
const N_L0         = 100
const ALPHA_R_MAX  = 4.0      # hard stop on the self-extending window
const L_MAX_POINTS = 400
const F_MAX_POINTS = 220

global alpha_r = collect(range(0.0, 1.0, length = N_F0))
global alpha_i = zeros(Float64, N_F0)
global omega_r = collect(range(0.0, 0.5, length = N_L0))
global omega_i = 0.0

global F = ComplexF64[]
global L = ComplexF64[]
global omega_F = ComplexF64[]
global alpha_L_u = ComplexF64[]
global alpha_L_l = ComplexF64[]

@everywhere F = ComplexF64[]
@everywhere normals_F = ComplexF64[]

contour_F() = ComplexF64[alpha_r[j] + alpha_i[j] * im for j in eachindex(alpha_r)]
contour_L() = ComplexF64[omega_r[j] + omega_i * im for j in eachindex(omega_r)]
contour_L_at(wi) = ComplexF64[omega_r[j] + wi * im for j in eachindex(omega_r)]

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

# F and its normals must agree with each other, so they are always shipped
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

function linear_interp(xs, ys, xq)
    n = length(xs)
    xq <= xs[1]   && return ys[1]
    xq >= xs[end] && return ys[end]
    j = clamp(searchsortedlast(xs, xq), 1, n - 1)
    t = (xq - xs[j]) / (xs[j+1] - xs[j])
    return ys[j] + t * (ys[j+1] - ys[j])
end

# =============================================================================
# 3.  BRANCH TRACKING
#
#   One spatial solve per omega; each root is labelled upper/lower by the sign
#   of its projection onto the nearest F normal, and the branches are then
#   walked outward from the middle as a chain.
#
#   TRACK_SLACK is v4.7's continuity guard and it is not optional: the side
#   test's resolution is F's node spacing, and near the pinch it is being asked
#   to resolve a gap far smaller than that.  Without the guard the upper branch
#   was on the wrong mode on 187 of 501 frames.
# =============================================================================
const TRACK_SLACK = 10.0

global track_overrides = 0
global track_jump_u = 0.0

@everywhere function spatial_payload(omega, Fv, nrm)
    ev = finite_only(alpha_spectrum(omega))
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

function nearest_each_side(pl)
    eu = isempty(pl.up) ? nothing : pl.up[argmin(pl.up_d)]
    el = isempty(pl.lo) ? nothing : pl.lo[argmin(pl.lo_d)]
    return eu, el
end

function select_tracked(pl, alpha_prev, side::Symbol)
    isempty(pl.all) && error("spatial_payload returned no finite eigenvalues")
    ap = ComplexF64(alpha_prev)
    all_best = pl.all[argmin(abs.(ap .- pl.all))]
    cand = side === :upper ? pl.up : pl.lo
    isempty(cand) && return all_best
    side_best = cand[argmin(abs.(ap .- cand))]
    if abs(ap - side_best) > TRACK_SLACK * max(abs(ap - all_best), 1e-9)
        global track_overrides += 1
        return all_best
    end
    return side_best
end

# chain = false is for the very first call, where F is the real axis and every
# node can pick its own nearest pair without a predecessor.
function track_branches(Lv; chain::Bool = true)
    payloads = pmap(w -> spatial_payload(w, F, normals_F), Lv)
    n  = length(Lv)
    au = Vector{ComplexF64}(undef, n)
    al = Vector{ComplexF64}(undef, n)
    global track_overrides = 0
    if !chain
        for j in 1:n
            eu, el = nearest_each_side(payloads[j])
            (eu === nothing || el === nothing) &&
                error("empty branch side at omega = $(Lv[j])")
            au[j] = eu; al[j] = el
        end
    else
        s = max(1, n ÷ 4)
        eu, el = nearest_each_side(payloads[s])
        (eu === nothing || el === nothing) &&
            error("empty branch side at the seed omega = $(Lv[s])")
        au[s] = eu; al[s] = el
        for j in (s+1):n
            au[j] = select_tracked(payloads[j], au[j-1], :upper)
            al[j] = select_tracked(payloads[j], al[j-1], :lower)
        end
        for j in (s-1):-1:1
            au[j] = select_tracked(payloads[j], au[j+1], :upper)
            al[j] = select_tracked(payloads[j], al[j+1], :lower)
        end
    end
    global track_jump_u = n < 2 ? 0.0 : maximum(abs.(au[2:end] .- au[1:end-1]))
    return au, al
end

# The two sides must never be the same mode at the same omega.
branches_ok(au, al) = minimum(abs.(au .- al)) >= 1e-9
contour_distance(Fv, au, al) =
    minimum(minimum(abs.(f .- vcat(au, al))) for f in Fv)

# =============================================================================
# 4.  THE BARRIER
#
#   Each tracked branch point carries a charge with its local arc-length
#   weight; F is repelled by all of them.  zeta must shrink with the gap F has
#   to thread or the exponent grows like 1/d^2, and eps is tied to zeta so the
#   exponent is bounded by EXP_ARG_TARGET at every stage.  ZETA_REF at
#   ZETA_D_REF is the one setting in this file taken from run history rather
#   than derived -- it is the value that was healthy in v4.4/v4.5, and it says
#   nothing about where the pinch is.
#
#   To test whether the barrier is needed at all, set ZETA_REF = 1e-12 and
#   watch nx (crossings): the veto in section 6 enforces the topology on its
#   own, and the barrier's preferred orientation is not the right one.
# =============================================================================
const ZETA_REF        = 4.0e-4
const ZETA_D_REF      = 3.08e-2
const EXP_ARG_TARGET  = 10.0

global zeta_alpha = ZETA_REF
global eps_alpha  = ZETA_REF / EXP_ARG_TARGET
global charge_z   = ComplexF64[]
global charge_w   = Float64[]

function set_zeta!(d)
    z = (isfinite(d) && d > 0) ?
        clamp(ZETA_REF * (d / ZETA_D_REF)^2, 1e-12, ZETA_REF) : ZETA_REF
    global zeta_alpha = z
    global eps_alpha  = max(z / EXP_ARG_TARGET, 1e-300)
    return z
end

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

function set_charges!(au, al)
    global charge_z = vcat(au, al)
    global charge_w = vcat(arc_weights(au), arc_weights(al))
    return nothing
end

function barrier_grad(af)
    gr = 0.0; gi = 0.0
    @inbounds for k in eachindex(charge_z)
        dz  = af - charge_z[k]
        den = abs2(dz) + eps_alpha
        e   = exp(min(zeta_alpha / den, 400.0))
        c   = -2.0 * zeta_alpha * e / den^2 * charge_w[k]
        gr += c * real(dz)
        gi += c * imag(dz)
    end
    return gr, gi
end

# =============================================================================
# 5.  THE PEAK OF Im(omega_F)
#
#   omega_i is floored by this number, so it must not extrapolate.  With
#   A = y_j - y_{j-1} and B = y_j - y_{j+1} at a discrete maximum, the three-
#   point parabola puts the vertex at 0.5(A-B)/(A+B) -- always within half a
#   cell, so v5's |t| > 1 guard was dead code -- and its excess is
#   (A-B)^2/8(A+B), which tends to max(A,B)/8 on a one-sided triple.  On a
#   jagged F that is unbounded and it pushed L back UP, which is what took v5
#   from 4.28e-8 out to 7.26e-5.
#
#   `res` is the resolution diagnostic: how far a 5-node least-squares parabola
#   disagrees with the 3-node one, as a fraction of the local curvature.  It
#   reads only Im(omega_F).
# =============================================================================
const PEAK_ASYM_MIN   = 0.05
const PEAK_EXCESS_CAP = 0.5

function peak_imag(xs, wF)
    y = imag.(wF)
    j = argmax(y)
    (j == firstindex(y) || j == lastindex(y)) && return y[j], NaN, NaN
    A = y[j] - y[j-1]
    B = y[j] - y[j+1]
    (A < 0 || B < 0 || (A + B) <= 0) && return y[j], NaN, NaN
    asym = min(A, B) / max(A, B)
    asym < PEAK_ASYM_MIN && return y[j], asym, NaN
    e3 = min((A - B)^2 / (8 * (A + B)), PEAK_EXCESS_CAP * min(A, B))
    p3 = y[j] + e3
    res = NaN
    if j - 2 >= firstindex(y) && j + 2 <= lastindex(y)
        x  = xs[j-2:j+2] .- xs[j]
        cc = hcat(x .^ 2, x, ones(5)) \ y[j-2:j+2]
        if all(isfinite, cc) && cc[1] < 0
            p5  = cc[3] - cc[2]^2 / (4 * cc[1])
            res = abs(p5 - p3) / max(A, B)
        end
    end
    return (isfinite(p3) ? p3 : y[j]), asym, res
end

# =============================================================================
# 6.  THE F STEP
#
#     rhs  = (y' dPhi/dx - dPhi/dy)/(1 + y'^2)  -  DESC_FRAC (theta/dt0) * descent
#     y   += theta * tanh(dt * rhs / theta)
#     y   := implicit sigma*y''   then   width average
#     accept if no branch crosses F and the peak falls within tolerance;
#     otherwise halve dt, retry, and apply the best of the attempts regardless.
#
#   The acceptance rule is deliberately NOT strict monotone descent.  The
#   objective is a max over nodes, so a step that lowers the peak at the
#   current argmax raises it elsewhere; from v5's own record the peak rose on
#   27-48 % of steps and the run descended anyway.  Refusing those steps was
#   tried and it stalled.
#
#   The tolerance has a FLOOR at a fraction of the measured frame-to-frame
#   jitter.  Without it, tol is proportional to the recent descent rate, which
#   tends to zero as the run slows -- slow, tighter test, smaller step, slower.
#
#   The descent direction is exact and free: omega is analytic, so the chord
#   derivative along F is the full complex derivative, and dv/dy =
#   Re(domega/dalpha).  THE CONTOUR'S ORIENTATION IS AN OUTPUT OF THIS TERM,
#   never an input.  Watch F's slope at the peak; it is also a free check.
# =============================================================================
const SIGMA        = 3e-5
const DT0          = 1e-3
const MOVE_FRAC    = 0.05     # node displacement cap, as a fraction of d_branch
const MAX_MOVE     = 0.20
const DESC_FRAC    = 0.30     # share of theta the descent term may use, FIXED
const DESC_NW      = 9
const MAX_ATTEMPTS = 5
const TOL_WINDOW   = 24       # even
const TOL_FRAC     = 0.10
const TOL_JITTER   = 0.10     # the floor that stops the rule throttling itself
const CROSS_ESCAPE = 8
const SMOOTH_FRAC  = 0.05
const SMOOTH_CELLS = 2.0      # floor, in LOCAL cells at the peak

global delta_t = DT0
global barren_run = 0
global peak_hist = Float64[]

function peak_push!(p)
    isfinite(p) || return nothing
    push!(peak_hist, p)
    while length(peak_hist) > TOL_WINDOW
        popfirst!(peak_hist)
    end
    return nothing
end

function descent_tolerance()
    length(peak_hist) < TOL_WINDOW && return Inf
    h = TOL_WINDOW >> 1
    rate = (median(peak_hist[1:h]) - median(peak_hist[h+1:end])) / h
    jit  = std(diff(peak_hist))
    isfinite(rate) || return Inf
    return max(TOL_FRAC * max(rate, 0.0), TOL_JITTER * (isfinite(jit) ? jit : 0.0))
end

# smoothing half-width in alpha_r, floored at SMOOTH_CELLS local cells at the
# peak.  v4.9's filter was a fixed number of INDICES and flattened F's slope
# 180x over a run; v5's width collapsed to an exact no-op below one cell.
function smooth_halfwidth(d, hpk)
    w = isfinite(d) ? SMOOTH_FRAC * d : 0.1515
    return max(SMOOTH_CELLS * hpk, min(w, 0.1515))
end

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
                wm = 1.0 - abs(xs[m] - xs[j]) / half_w
                wm <= 0 && continue
                acc += wm * ys[m]; wsum += wm
            end
            out[j] = wsum > 0 ? acc / wsum : ys[j]
        end
    end
    return out
end

# deviation of each node from the chord through its neighbours: the 2h content
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

function domega_dalpha(Fv, wF)
    n = length(Fv)
    g = zeros(ComplexF64, n)
    n < 2 && return g
    for j in 2:(n - 1)
        dz = Fv[j+1] - Fv[j-1]
        g[j] = dz == 0 ? 0.0 + 0.0im : (wF[j+1] - wF[j-1]) / dz
    end
    dz1 = Fv[2] - Fv[1];     g[1] = dz1 == 0 ? 0.0 + 0.0im : (wF[2] - wF[1]) / dz1
    dzn = Fv[n] - Fv[n-1];   g[n] = dzn == 0 ? 0.0 + 0.0im : (wF[n] - wF[n-1]) / dzn
    for j in eachindex(g)
        isfinite(real(g[j])) && isfinite(imag(g[j])) || (g[j] = 0.0 + 0.0im)
    end
    return g
end

# softmax over the top DESC_NW nodes, so only the peak's neighbourhood is driven
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

# implicit sigma*y'' on a non-uniform grid, zero slope at both ends.  Implicit
# because F is graded: explicit diffusion would cap dt at h^2/(2 sigma), and at
# the spacings this run reaches that cap is smaller than any useful step.
function diffusion_implicit(xs, ys, mscale, dt, sig)
    n = length(ys)
    (n < 3 || dt <= 0 || sig <= 0) && return copy(ys)
    dl = zeros(Float64, n - 1); dg = ones(Float64, n); du = zeros(Float64, n - 1)
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

# how far the nearest branch point is on the CORRECT side of F, and how many
# are on the wrong side.  Only pairs whose own gap is small and that lie inside
# F's span are tested: a far-field pair outside the span would be compared
# against a fictitious flat extension of F.
function branch_clearance(Fv, au, al; gap_bound = 0.5)
    xs = real.(Fv); ys = imag.(Fv)
    lo = xs[1]; hi = xs[end]
    worst = Inf; ncross = 0
    for k in eachindex(au)
        abs(au[k] - al[k]) > gap_bound && continue
        (lo <= real(au[k]) <= hi) || continue
        (lo <= real(al[k]) <= hi) || continue
        cu = imag(au[k]) - linear_interp(xs, ys, real(au[k]))
        cl = linear_interp(xs, ys, real(al[k])) - imag(al[k])
        c = min(cu, cl)
        c < worst && (worst = c)
        c < 0 && (ncross += 1)
    end
    return worst, ncross
end

# everything about the step that does not depend on dt, computed once per
# iteration rather than once per trial
function f_predirection()
    n = length(alpha_i)
    yr = zeros(Float64, n); gr = zeros(Float64, n); gi = zeros(Float64, n)
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
    return (yr = yr, gr = gr, gi = gi, draw = draw, dsc = maximum(abs, draw))
end

# dt0 is the step the attempt loop STARTED at and never changes.  The descent
# term is scaled by theta/dt0 so its contribution shrinks linearly when the
# line search backs off, exactly like the barrier term; scaling it by theta/dt
# makes that contribution dt-independent and every trial keeps overshooting.
function f_trial(dt, dt0, theta, hw, pre)
    n = length(alpha_i)
    yt = copy(alpha_i)
    nlim = 0
    for j in 2:(n - 1)
        rhs = (pre.yr[j] * pre.gr[j] - pre.gi[j]) / (1.0 + pre.yr[j]^2)
        pre.dsc > 0 && (rhs -= DESC_FRAC * (theta / dt0) * pre.draw[j] / pre.dsc)
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

function f_step!(d_branch, hpk)
    n = length(alpha_i)
    peak_old, _, _ = peak_imag(alpha_r, omega_F)
    tol = descent_tolerance()

    theta = (isfinite(d_branch) && d_branch > 0) ?
            min(MOVE_FRAC * d_branch, MAX_MOVE) : MAX_MOVE
    theta = max(theta, 1e-14)
    hw = smooth_halfwidth(d_branch, hpk)

    pre = f_predirection()
    clear_now, ncross_now = branch_clearance(F, alpha_L_u, alpha_L_l)
    veto = barren_run < CROSS_ESCAPE

    dt  = delta_t * clamp(d_branch / 0.2, 0.1, 1.0)
    dt0 = dt
    best_peak = Inf; best_y = Float64[]; best_w = ComplexF64[]; best_at = 0
    nlim = 0; asym = NaN; res = NaN
    minclr = clear_now; ncross = ncross_now
    status = "none"

    for attempt in 1:MAX_ATTEMPTS
        yt, nl = f_trial(dt, dt0, theta, hw, pre)
        if !all(isfinite, yt)
            dt *= 0.5; status = "nonfinite"; continue
        end
        Ft = ComplexF64[alpha_r[j] + yt[j] * im for j in 1:n]
        ct, nx = branch_clearance(Ft, alpha_L_u, alpha_L_l)
        if veto && !(ct >= 0.0 || ct > clear_now)
            dt *= 0.5; status = "crossing"; continue
        end
        wt = contour_omega_F(Ft)
        pt, at, rs = peak_imag(alpha_r, wt)
        if isfinite(pt) && pt < best_peak
            best_peak = pt; best_y = copy(yt); best_w = copy(wt); best_at = attempt
            nlim = nl; asym = at; res = rs; minclr = ct; ncross = nx
        end
        if !isfinite(peak_old) || !isfinite(tol) ||
           (isfinite(pt) && pt <= peak_old + tol)
            status = "accepted"
            break
        end
        dt *= 0.5
        status = "reduced"
    end

    if best_at == 0
        global barren_run += 1
        status = "discarded"
    else
        global barren_run = 0
        global alpha_i = copy(best_y)
        global omega_F = copy(best_w)
        global F = contour_F()
        push_F!()
        peak_push!(best_peak)
    end

    pmove = (best_at == 0 || !isfinite(peak_old)) ? NaN : best_peak - peak_old
    return (status = status, attempts = max(best_at, 1), peak_move = pmove,
            n_limited = nlim, peak_asym = asym, peak_res = res,
            min_clear = minclr, n_cross = ncross, dt = dt, theta = theta,
            smooth_hw = hw, descent_tol = tol)
end

# =============================================================================
# 7.  THE OMEGA STEP
#
#   omega_i is floored by F's own peak.  Over v5's whole 1000-frame run that
#   floor was reached directly on 966 frames and the bisection never used more
#   than one probe, so v6 sets it directly and bisects only when the floor is
#   not admissible.
# =============================================================================
const OMEGA_CLEARANCE = 1e-9
const OMEGA_PROBES    = 8

function omega_trial(wi)
    Lt = contour_L_at(wi)
    au, al = track_branches(Lt)
    (all(isfinite, au) && all(isfinite, al)) || return false, Lt, au, al
    return branches_ok(au, al), Lt, au, al
end

function omega_step!()
    pk, _, _ = peak_imag(alpha_r, omega_F)
    lb = pk + OMEGA_CLEARANCE

    if omega_i <= lb
        # F's peak has risen above where L sits; L must come back up.  Nothing
        # to validate -- this is the state the previous iteration certified.
        global omega_i = lb
        global L = contour_L()
        global alpha_L_u, alpha_L_l = track_branches(L)
        return "up", omega_i - lb
    end

    ok, Lt, au, al = omega_trial(lb)
    if ok
        global omega_i = lb
        global L = Lt; global alpha_L_u = au; global alpha_L_l = al
        return "pin", 0.0
    end

    a_lo = lb; b_hi = omega_i
    bw = omega_i; bL = L; bu = alpha_L_u; bl = alpha_L_l
    for _ in 2:OMEGA_PROBES
        wt = 0.5 * (a_lo + b_hi)
        ok2, Lt2, au2, al2 = omega_trial(wt)
        if ok2
            bw = wt; bL = Lt2; bu = au2; bl = al2; b_hi = wt
        else
            a_lo = wt
        end
        (b_hi - a_lo) < 1e-14 && break
    end
    if bw < omega_i
        global omega_i = bw
        global L = bL; global alpha_L_u = bu; global alpha_L_l = bl
        return "bis", omega_i - lb
    end
    return "stuck", omega_i - lb
end

# =============================================================================
# 8.  THE THREE GRIDS
#
#   All three read only the run's own state: the measured gap profile, the
#   tracked branch pair, and Im(omega_F).
#
#   L   refine when the winning omega_r has parked, the gap has stopped
#       drifting, and the horizontal miss still explains a real share of the
#       gap.  Nodes go at the vertex of the d^4 parabola, which is a parabola
#       in omega_r because d = 2 sqrt(2|w - w_p|/|w''|).  The vertex is fitted
#       twice, over 2 and over 3 nodes either side, and used only when the two
#       agree -- that is the stencil's own error bar and it replaces v5.1's
#       ADAPT_FIT_HALF, which was the one constant in this project chosen by
#       comparing against a known answer.
#
#   F   one rule: the spacing at the peak must be below F_H_FRAC * d_branch.
#       The fine block is centred on the midpoint of the tracked branch pair at
#       the argmin -- the run's own estimate of where the throat is.
#
#   win extend alpha_r to the right while the peak is near the end or the end
#       value is still rising.  The left end is fixed by symmetry.
# =============================================================================
const STALL_WINDOW = 20
const STALL_TOL    = 1e-2
const WINDOW_CELLS = 2.0
const GS_MIN       = 0.05     # below this the L grid is not the limiter
const CLUSTER_HALF = 2
const F_H_FRAC     = 0.25     # target spacing at the peak, as a share of d
const F_FINE       = 20       # nodes either side of the throat at that spacing
const F_RATIO      = 1.3      # geometric stretch outside the fine block
const F_HOLD       = 10       # iterations the trigger must hold
const WIN_MARGIN   = 2
const WIN_STEP     = 0.25

global hist_d = Float64[]
global hist_w = Float64[]
global hist_h = Float64[]
global f_coarse_run = 0

function hist_push!(d, w, h)
    push!(hist_d, d); push!(hist_w, w); push!(hist_h, h)
    while length(hist_d) > STALL_WINDOW
        popfirst!(hist_d); popfirst!(hist_w); popfirst!(hist_h)
    end
    return nothing
end
hist_clear!() = (empty!(hist_d); empty!(hist_w); empty!(hist_h); nothing)

function stalled()
    n = length(hist_d)
    n < STALL_WINDOW && return false
    href = minimum(hist_h)
    (isfinite(href) && href > 0) || return false
    (maximum(hist_w) - minimum(hist_w)) > WINDOW_CELLS * href && return false
    m = abs(mean(hist_d))
    m <= 0 && return false
    h = n ÷ 2
    drift = abs(mean(hist_d[(h+1):end]) - mean(hist_d[1:h]))
    noise = 2.0 * std(hist_d) / sqrt(n)
    return drift <= STALL_TOL * m + noise
end

local_h(wr, i) = length(wr) < 2 ? Inf :
                 i == 1 ? wr[2] - wr[1] :
                 i == length(wr) ? wr[end] - wr[end-1] :
                 min(wr[i+1] - wr[i], wr[i] - wr[i-1])

# vertex of d^4 against omega_r over n nodes either side
function vertex_fit(wr, d, i, n)
    lo = max(firstindex(d), i - n); hi = min(lastindex(d), i + n)
    (i - lo < 2 || hi - i < 2) && return (NaN, NaN)
    x = wr[lo:hi] .- wr[i]
    y = d[lo:hi] .^ 4
    xs = maximum(abs, x); ys = maximum(abs, y)
    (isfinite(xs) && isfinite(ys) && xs > 0 && ys > 0) || return (NaN, NaN)
    c = hcat((x ./ xs) .^ 2, x ./ xs, ones(length(x))) \ (y ./ ys)
    (all(isfinite, c) && c[1] > 0) || return (NaN, NaN)
    curv = ys * c[1] / xs^2
    return (wr[i] - xs * c[2] / (2 * c[1]), curv > 0 ? 8.0 / sqrt(curv) : NaN)
end

# (vertex, |omega''| estimate, the two stencils' disagreement)
function vertex_estimate(wr, d, i)
    v2, _  = vertex_fit(wr, d, i, 2)
    v3, k3 = vertex_fit(wr, d, i, 3)
    (isfinite(v2) && isfinite(v3)) || return (NaN, NaN, NaN)
    return (0.5 * (v2 + v3), k3, abs(v2 - v3))
end

function insert_cluster!(w_c, h)
    (isfinite(w_c) && isfinite(h) && h > 0) || return false
    (w_c <= omega_r[1] || w_c >= omega_r[end]) && return false
    out = collect(omega_r)
    added = 0
    for j in -CLUSTER_HALF:CLUSTER_HALF
        w = w_c + j * h
        (w <= omega_r[1] || w >= omega_r[end]) && continue
        minimum(abs.(out .- w)) < h / 8 && continue
        push!(out, w); added += 1
    end
    added == 0 && return false
    sort!(out)
    global omega_r = out
    return true
end

# fine block at h0 around c, then geometric to hmax, out to [lo, hi]
function graded_nodes(c, h0, lo, hi, hmax, nmax)
    xs = Float64[c]
    half = nmax ÷ 2                      # each side gets its own budget, so a
    for dir in (1, -1)                   # long right side cannot starve the left
        x = c; h = h0; used = 0
        for k in 1:half
            h = k <= F_FINE ? h0 : min(h * F_RATIO, hmax)
            x += dir * h
            (x <= lo || x >= hi) && break
            push!(xs, x); used += 1
            used >= half && break
        end
    end
    push!(xs, lo); push!(xs, hi)
    sort!(xs)
    out = Float64[xs[1]]
    for x in xs[2:end]
        x - out[end] > h0 / 8 && push!(out, x)
    end
    return out
end

function regrid_F!(c, d)
    lo = alpha_r[1]; hi = alpha_r[end]
    (isfinite(c) && lo < c < hi && isfinite(d) && d > 0) || return false
    h0   = F_H_FRAC * d
    hmax = (hi - lo) / 40
    h0 >= hmax && return false
    xs = graded_nodes(c, h0, lo, hi, hmax, F_MAX_POINTS)
    length(xs) <= length(alpha_r) && return false
    ynew = [linear_interp(alpha_r, alpha_i, x) for x in xs]
    global alpha_r = xs
    global alpha_i = ynew
    return true
end

function extend_window!()
    y = imag.(omega_F)
    n = length(y)
    n < 3 && return false
    j = argmax(y)
    need = (j >= n - WIN_MARGIN) || (y[n] >= y[n-1])
    (need && alpha_r[end] < ALPHA_R_MAX) || return false
    hcoarse = max(alpha_r[end] - alpha_r[end-1], (alpha_r[end] - alpha_r[1]) / 40)
    hi = min(alpha_r[end] + WIN_STEP, ALPHA_R_MAX)
    xs = collect(alpha_r); ys = collect(alpha_i)
    x = xs[end]
    while x + hcoarse < hi && length(xs) < F_MAX_POINTS
        x += hcoarse
        push!(xs, x); push!(ys, ys[end])
    end
    push!(xs, hi); push!(ys, ys[end])
    global alpha_r = xs
    global alpha_i = ys
    return true
end

# =============================================================================
# 9.  LOG AND RESUME  --  one append-only file
#
#   Every line carries the primitives a restart needs (omega_i, omega_r,
#   alpha_r, F) plus the frame's numbers.  Resume reads the LAST PARSEABLE
#   line, so a partial line left by a crash costs one iteration.  The v4.3 and
#   v5.1 crashes both came from rewriting the whole file; this never rewrites.
# =============================================================================
const LOGFILE     = "contour_iteration_v6.jsonl"
const RESUME      = true
const ITER_TARGET = 3000

cvec(v) = [Dict("re" => real(x), "im" => imag(x)) for x in v]
jnum(x) = (x isa Real && isfinite(x)) ? x : nothing
uncvec(a) = ComplexF64[complex(x["re"], x["im"]) for x in a]

global iteration_step = 1

function log_line(d::Dict)
    open(LOGFILE, "a") do io
        println(io, JSON.json(d))
    end
    return nothing
end

function try_resume()
    isfile(LOGFILE) || return false
    last = nothing
    nlines = 0
    for line in eachline(LOGFILE)
        isempty(strip(line)) && continue
        e = try
            JSON.parse(line)
        catch
            nothing                      # partial final line after a crash
        end
        if e !== nothing
            last = e; nlines += 1
        end
    end
    last === nothing && return false
    (get(last, "num_modes", NUM_MODES) == NUM_MODES && get(last, "Re", RE) == RE) ||
        error("$LOGFILE was written with different numerics (num_modes / Re) -- rename it first")
    global iteration_step = Int(last["iteration"]) + 1
    global omega_r = Float64.(last["omega_r"])
    global alpha_r = Float64.(last["alpha_r"])
    global alpha_i = imag.(uncvec(last["F"]))
    global omega_i = Float64(last["omega_i"])
    global F = contour_F()
    global L = contour_L()
    push_F!()
    global omega_F = contour_omega_F(F)
    global alpha_L_u, alpha_L_l = track_branches(L)
    set_charges!(alpha_L_u, alpha_L_l)
    @printf("RESUMED at iteration %d from %d lines: omega_i = %.12e, N_L = %d, N_F = %d\n",
            iteration_step, nlines, omega_i, length(omega_r), length(alpha_r))
    flush(stdout)
    return true
end

# =============================================================================
# 10.  INITIALISE
# =============================================================================
if !(RESUME && try_resume())
    isfile(LOGFILE) && error("$LOGFILE exists; rename or delete it for a clean start")
    global F = contour_F()
    global L = contour_L()
    push_F!()
    global omega_F = contour_omega_F(F)
    global alpha_L_u, alpha_L_l = track_branches(L; chain = false)
    set_charges!(alpha_L_u, alpha_L_l)
    global iteration_step = 1
end

@printf("\nv6 ready.  N_F = %d, N_L = %d, num_modes = %d, Re = %.0f\n",
        length(alpha_r), length(omega_r), NUM_MODES, RE)
@printf("starting at iteration %d, omega_i = %.12e\n\n", iteration_step, omega_i)
flush(stdout)

# =============================================================================
# 11.  MAIN LOOP
#
#   The grids move FIRST, using the previous iteration's diagnostics, so every
#   logged line describes one consistent state and a resume from it is exact.
# =============================================================================
global gs_last = NaN
global vtx_last = NaN
global vspread_last = NaN
global ipk_last = 0

while iteration_step <= ITER_TARGET
    k = iteration_step

    # ---- 1. grids ----------------------------------------------------------
    if k > 1
        moved_F = false
        extend_window!() && (moved_F = true)

        # F: the spacing AT THE THROAT must be below F_H_FRAC * d.  The throat
        # is the midpoint of the tracked branch pair at the argmin -- the run's
        # own estimate, never an input.  Measuring the trigger and centring the
        # fine block at the same place keeps the rule self-consistent.
        if ipk_last > 0 && length(alpha_r) < F_MAX_POINTS
            c  = real(0.5 * (alpha_L_u[ipk_last] + alpha_L_l[ipk_last]))
            d  = abs(alpha_L_u[ipk_last] - alpha_L_l[ipk_last])
            jc = argmin(abs.(alpha_r .- c))
            if isfinite(d) && d > 0 && local_h(alpha_r, jc) > F_H_FRAC * d
                global f_coarse_run += 1
            else
                global f_coarse_run = 0
            end
            if f_coarse_run >= F_HOLD
                global f_coarse_run = 0
                if regrid_F!(c, d)
                    moved_F = true
                    @printf("[k=%d] F regrid: N_F = %d, target h = %.3e at alpha_r = %.6f (d = %.3e)\n",
                            k, length(alpha_r), F_H_FRAC * d, c, d)
                    flush(stdout)
                end
            end
        end

        if moved_F
            global F = contour_F()
            push_F!()
            global omega_F = contour_omega_F(F)
            global alpha_L_u, alpha_L_l = track_branches(L)
            set_charges!(alpha_L_u, alpha_L_l)
            hist_clear!()
        end

        if stalled() && isfinite(gs_last) && gs_last >= GS_MIN &&
           length(omega_r) + 2 * CLUSTER_HALF <= L_MAX_POINTS &&
           isfinite(vtx_last) && isfinite(vspread_last) &&
           vspread_last <= local_h(omega_r, max(ipk_last, 1))
            hcl = max(local_h(omega_r, max(ipk_last, 1)) / 4,
                      2 * abs(omega_r[max(ipk_last, 1)] - vtx_last))
            if insert_cluster!(vtx_last, hcl)
                global L = contour_L()
                global alpha_L_u, alpha_L_l = track_branches(L)
                set_charges!(alpha_L_u, alpha_L_l)
                hist_clear!()
                @printf("[k=%d] L refine: N_L = %d, h = %.3e at omega_r = %.12f (gs = %.2f)\n",
                        k, length(omega_r), hcl, vtx_last, gs_last)
                flush(stdout)
            end
        end
    end

    # ---- 2. omega: press L down onto F's peak -------------------------------
    ostat, ogap = omega_step!()
    if ostat == "stuck"
        println("STOP: no admissible omega height below the current one.")
        break
    end
    set_charges!(alpha_L_u, alpha_L_l)

    # ---- 3. geometry -------------------------------------------------------
    d_vec     = abs.(alpha_L_u .- alpha_L_l)
    i_pinch   = argmin(d_vec)
    d_branch  = d_vec[i_pinch]
    h_local   = local_h(omega_r, i_pinch)
    d_contour = contour_distance(F, alpha_L_u, alpha_L_l)
    if !isfinite(d_branch) || !isfinite(d_contour)
        @printf("[k=%d] STOP: non-finite distance\n", k)
        break
    end

    vtx, w2, vspread = vertex_estimate(omega_r, d_vec, i_pinch)
    miss = isfinite(vtx) ? abs(omega_r[i_pinch] - vtx) : NaN
    gs = (isfinite(miss) && isfinite(w2) && w2 > 0 && d_branch > 0) ?
         2 * sqrt(2 * miss / w2) / d_branch : NaN
    global gs_last = gs
    global vtx_last = vtx; global vspread_last = vspread; global ipk_last = i_pinch

    # spacing at the throat (what the F rule controls) and at the peak (what
    # the smoothing floor uses)
    c_throat = real(0.5 * (alpha_L_u[i_pinch] + alpha_L_l[i_pinch]))
    h_throat = local_h(alpha_r, argmin(abs.(alpha_r .- c_throat)))
    hpk = local_h(alpha_r, argmax(imag.(omega_F)))
    set_zeta!(d_branch)
    hist_push!(d_branch, omega_r[i_pinch], h_local)

    # ---- 4. F: deform it downwards -----------------------------------------
    st = f_step!(d_branch, hpk)
    rip = f_ripple_of(F)

    # F's slope at the peak.  It should settle near Re(w'')/Im(w'') on its own;
    # that is a free check, never an input.
    jp = clamp(argmax(imag.(omega_F)), 2, length(alpha_r) - 1)
    fslope = (alpha_i[jp+1] - alpha_i[jp-1]) / (alpha_r[jp+1] - alpha_r[jp-1])

    # is L still below F's own peak after the step?  it should be.
    pk_after, _, _ = peak_imag(alpha_r, omega_F)
    admissible = omega_i >= pk_after

    @printf("[%04d] d=%9.3e | w=%.10f g=%8.1e %-3s | NL=%3d h=%8.2e gs=%5.2f | NF=%3d hF/d=%6.2f res=%5.2f | dpk=%+9.2e a%d tol=%8.2e | rip=%8.2e lim=%2d clr=%+9.2e nx=%2d | slp=%+6.2f %s %s\n",
            k, d_branch, omega_i, ogap, ostat,
            length(omega_r), h_local, gs,
            length(alpha_r), h_throat / max(d_branch, 1e-300), st.peak_res,
            st.peak_move, st.attempts, st.descent_tol,
            rip, st.n_limited, st.min_clear, st.n_cross,
            fslope, admissible ? "ok" : "LOW", st.status)
    flush(stdout)

    # ---- 5. one line, appended ---------------------------------------------
    log_line(Dict(
        "iteration" => k, "num_modes" => NUM_MODES, "Re" => RE,
        "omega_i" => omega_i, "omega_r" => collect(omega_r),
        "alpha_r" => collect(alpha_r),
        "F" => cvec(F), "L" => cvec(L),
        "alpha_L_u" => cvec(alpha_L_u), "alpha_L_l" => cvec(alpha_L_l),
        "omega_F" => cvec(omega_F),
        "n_F" => length(alpha_r), "n_L" => length(omega_r),
        "d_branch" => jnum(d_branch), "d_contour" => jnum(d_contour),
        "h_local" => jnum(h_local), "h_peak" => jnum(hpk),
        "h_throat" => jnum(h_throat), "throat_r" => jnum(c_throat),
        "grid_share" => jnum(gs), "vertex" => jnum(vtx),
        "vertex_spread" => jnum(vspread), "omega2_fit" => jnum(w2),
        "omega_status" => ostat, "omega_gap" => jnum(ogap),
        "admissible" => admissible,
        "alpha_status" => st.status, "attempts" => st.attempts,
        "peak_move" => jnum(st.peak_move), "peak_asym" => jnum(st.peak_asym),
        "peak_res" => jnum(st.peak_res), "descent_tol" => jnum(st.descent_tol),
        "f_ripple" => jnum(rip), "f_slope" => jnum(fslope),
        "n_limited" => st.n_limited, "n_cross" => st.n_cross,
        "min_clear" => jnum(st.min_clear), "dt_used" => jnum(st.dt),
        "theta" => jnum(st.theta), "smooth_hw" => jnum(st.smooth_hw),
        "zeta_alpha" => jnum(zeta_alpha),
        "track_overrides" => track_overrides, "track_jump_u" => jnum(track_jump_u),
    ))

    global iteration_step += 1
end

@printf("\nFINISHED at iteration %d.  log: %s\n", iteration_step - 1, LOGFILE)
flush(stdout)
