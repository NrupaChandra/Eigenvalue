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
        Re = 5000.0
        beta = 0.0 + 0.0 * im
        num_modes = 150
        start = -1
        terminate = 1
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
        u = 1.0 .- y_colloc_points_new.^2
        d2u = -2.0 .* ones(length(y_colloc_points_new))
    end 
    @everywhere function equilibrate(mats...; iters::Int = 12)
        n = size(mats[1], 1)
        W = zeros(Float64, n, n)
        for M in mats
            W .+= abs.(M)
        end
        R = ones(Float64, n)
        S = ones(Float64, n)
        for _ in 1:iters
            r = sqrt.(max.(vec(maximum(W; dims = 2)), 1e-300))
            R ./= r
            W ./= r
            c = sqrt.(max.(vec(maximum(W; dims = 1)), 1e-300))
            S ./= c
            W ./= transpose(c)
        end
        return R, S
    end

    @everywhere function couetteflow(alpha, omega, mode::Val, mode2::Val)
        setprecision(53) do
            if mode === Val(:omega_collocation)
                A11 = - im * alpha * (u * ones(Complex{Float64}, 1, length(u))) .* D2 + im * alpha * (u * ones(Complex{Float64}, 1, length(u))) * (alpha^2 + beta^2) .* D0 + im * alpha * (d2u * ones(1, length(u))) .* D0 + 1 / Re .* D4 - 2 / Re * (alpha^2 + beta^2) .* D2 + 1 / Re * (alpha^2 + beta^2)^2 .* D0
                A11 = [-200 * im * [D0[1:1, :]; D1[1:1, :]];   A11[3:num_modes - 2, :];    -200 * im * [D1[num_modes:num_modes, :]; D0[num_modes:num_modes, :]]] 
                A = A11  
                B11 = - im .* D2 + im * (alpha^2 + beta^2) .* D0
                B11 = [[D0[1:1, :]; D1[1:1, :]];   B11[3:num_modes - 2, :];    [D1[num_modes:num_modes, :]; D0[num_modes:num_modes, :]]]
                B = B11
                # v7: equilibrate the (A, B) pair.  Same roots, cond(A) 9.7e10 -> 1.9e3.
                Rq, Sq = equilibrate(A, B)
                A = Rq .* A .* transpose(Sq)
                B = Rq .* B .* transpose(Sq)
            elseif mode === Val(:alpha_collocation)
                A11 = -2 * im * omega * D1 - 4 / Re * D3 + 4 / Re * beta^2 * D1 - im * (u * ones(Complex{Float64}, 1, length(u))) .* D2 + im * beta^2 * (u * ones(1, length(u))) .* D0 + im * (d2u * ones(1, length(u))) .* D0
                A12 = im * omega * D2 - im * omega * beta^2 * D0 + 1 / Re * D4 - 2 / Re * beta^2 * D2 + 1 / Re * beta^4 * D0
                A11 .= [zeros(Complex{Float64}, 2, num_modes);                    A11[3:num_modes - 2, :];    zeros(Complex{Float64}, 2, num_modes)]
                A12 .= [-200 * im * [D0[1:1, :]; D1[1:1, :]];   A12[3:num_modes - 2, :];    -200 * im * [D1[num_modes:num_modes, :]; D0[num_modes:num_modes, :]]]
                B11 = - 4 / Re * D2 - 2 * im * (u * ones(Complex{Float64}, 1, length(u))) .* D1
                B11 = [zeros(Complex{Float64}, 2, num_modes);  B11[3:num_modes - 2,:];    zeros(Complex{Float64}, 2, num_modes);]
                Rq, Sq = equilibrate(B11, A11, A12)
                A11 = Rq .* A11 .* transpose(Sq)
                A12 = Rq .* A12 .* transpose(Sq)
                B11 = Rq .* B11 .* transpose(Sq)
                A21 = 1 * Matrix{Complex{Float64}}(I, num_modes, num_modes)
                A22 = zeros(Complex{Float64}, num_modes, num_modes)
                A = [A11 A12; A21 A22]
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
        alpha_r_start = -1.5
        alpha_r_end = 1.5
        N = 300
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
        omega_r_start = -1.5
        omega_r_end = 1.5
        omega_i = 1e-3

        # v4.5: this is the ONLY grid built up front.  v4.4's REFINE_LO /
        # REFINE_HI / PTS_PER_CELL are gone -- the refinement below is placed
        # by the run, not guessed before it.
        omega_r_base = collect(range(omega_r_start, omega_r_end, length=150))
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
    #   up    those on the +normal side of F,  up_d their distance to F,
    #         up_s the arclength along F of their foot point
    #   lo    those on the -normal side of F,  lo_d / lo_s likewise
    #
    # Point-cloud version: the side test is unchanged (sign of the projection
    # onto the nearest node's normal).  up_d / lo_d are now the distance to
    # the F *polyline* rather than to the nearest node, and up_s / lo_s give
    # the position along F, which is the coordinate the cloud frontier is
    # built over (see cloud_frontier).
    # -----------------------------------------------------------------------
    @everywhere arclength(F) = [0.0; cumsum(abs.(diff(F)))]

    # Foot of e on the polyline F, searched on the two segments next to its
    # nearest node k.  Returns (distance, arclength of the foot).  The two end
    # segments are not clamped outwards, so eigenvalues beyond F's ends get
    # s < 0 or s > s_F[end] instead of all collapsing onto the end nodes.
    @everywhere function foot_on_F(e, F, s_F, k)
        n = length(F)
        best_d = Inf
        best_s = s_F[k]
        for i in (k - 1, k)
            (1 <= i && i + 1 <= n) || continue
            a = F[i]
            ab = F[i + 1] - a
            t = real(conj(ab) * (e - a)) / abs2(ab)
            t = clamp(t, i == 1 ? -Inf : 0.0, i + 1 == n ? Inf : 1.0)
            dd = abs(e - (a + t * ab))
            if dd < best_d
                best_d = dd
                best_s = s_F[i] + t * abs(ab)
            end
        end
        return best_d, best_s
    end

    @everywhere function spatial_payload(omega, F, normals_F)
        ev = couetteflow(nothing, omega, Val(:alpha_collocation), Val(:eigvals))
        ev = ev[[isfinite(real(e)) && isfinite(imag(e)) && abs(real(e)) > 1e-2 for e in ev]]
        s_F = arclength(F)
        up = ComplexF64[];  up_d = Float64[];  up_s = Float64[]
        lo = ComplexF64[];  lo_d = Float64[];  lo_s = Float64[]
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
            d_F, s_e = foot_on_F(e, F, s_F, k)
            if proj > 0.0
                push!(up, e);  push!(up_d, d_F);  push!(up_s, s_e)
            elseif proj < 0.0
                push!(lo, e);  push!(lo_d, d_F);  push!(lo_s, s_e)
            end
        end
        return (all = ev, up = up, up_d = up_d, up_s = up_s,
                lo = lo, lo_d = lo_d, lo_s = lo_s)
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

    # -----------------------------------------------------------------------
    # POINT-CLOUD BRANCHES (replaces seed-and-walk tracking).
    #
    # Every omega in L is one independent eigensolve, done by pmap on the
    # workers.  The main process then pools ALL classified eigenvalues of all
    # omegas into two clouds -- cloud_u (+normal side of F) and cloud_l
    # (-normal side) -- drops the ones of unreasonable magnitude, and takes as
    # alpha_L_u / alpha_L_l the frontier of each cloud that faces F.
    #
    # Nothing is tracked, so nothing can jump to a neighbouring mode, and
    # there is no 4/6 seed: whichever mode is closest to F at a given position
    # along F is on the frontier, from whatever omega it came.
    #
    # Frontier.  Each eigenvalue has a foot point on F at arclength s and a
    # distance d to F (spatial_payload).  The s-range of the cloud is cut into
    # bins and the minimum-d eigenvalue of each bin is kept.  A bin that holds
    # no sample of the near branch picks up a farther mode instead, which
    # shows up as a spike (or zig-zag) away from F.  A frontier point is
    # dropped when, within CLOUD_WIN bins on BOTH sides, there are points
    # closer to F by more than the factor 1 + CLOUD_SPIKE_REL; repeated until
    # nothing changes.  Requiring both sides keeps the steep flanks of a dip
    # towards F (e.g. near a pinch), where one side is always farther.
    # The result, ordered by s, is a
    # polyline of actual eigenvalues (CLOUD_SMOOTH = 0) or a moving average
    # of them (CLOUD_SMOOTH > 0).
    #
    # NOTE: the frontier is parametrised by s along F, not by omega, so
    # alpha_L_u[j] no longer belongs to L[j], and the two branches need not
    # have the same length.  Anything that pairs them by index (e.g.
    # branch_distance in briggsv4.4) has to use a set distance instead.
    # -----------------------------------------------------------------------
    const CLOUD_ABS_MAX  = 10.0   # drop |alpha| above this (spurious modes)
    const CLOUD_DIST_MAX = 2.0    # drop eigenvalues farther than this from F
    const CLOUD_NBINS    = 0      # frontier bins along F; 0 -> length(L)
    const CLOUD_WIN      = 3      # spike filter half-width, in bins
    const CLOUD_SPIKE_REL = 0.1   # spike if both sides hold a point with
                                  # d * (1 + CLOUD_SPIKE_REL) < own d
    const CLOUD_SMOOTH   = 0      # moving-average half-width; 0 = raw eigenvalues
    const CLOUD_S_PAD    = 0.0    # keep feet up to this fraction of len(F) beyond
                                  # F's ends; past them the near branch runs out
                                  # and far modes would form the frontier

    global cloud_u = ComplexF64[]   # filtered clouds of the last sweep,
    global cloud_l = ComplexF64[]   # kept for plot_cloud()

    function gather_cloud(payloads, side::Symbol)
        pts = ComplexF64[];  s = Float64[];  d = Float64[]
        S = arclength(F)[end]
        s_min, s_max = -CLOUD_S_PAD * S, (1 + CLOUD_S_PAD) * S
        for pl in payloads
            e, ss, dd = side === :upper ? (pl.up, pl.up_s, pl.up_d) :
                                          (pl.lo, pl.lo_s, pl.lo_d)
            for i in eachindex(e)
                (abs(e[i]) <= CLOUD_ABS_MAX && dd[i] <= CLOUD_DIST_MAX &&
                 s_min <= ss[i] <= s_max) || continue
                push!(pts, e[i]);  push!(s, ss[i]);  push!(d, dd[i])
            end
        end
        return pts, s, d
    end

    function smooth_polyline(z, w)
        (w <= 0 || length(z) < 3) && return z
        out = copy(z)
        for i in 2:(length(z) - 1)
            r = min(w, i - 1, length(z) - i)
            out[i] = sum(@view z[(i - r):(i + r)]) / (2r + 1)
        end
        return out
    end

    function cloud_frontier(pts, s, d; nbins::Int)
        isempty(pts) && error("cloud_frontier: empty cloud after filtering")
        s_lo, s_hi = extrema(s)
        h = (s_hi - s_lo) / nbins
        # minimum-d point of each bin
        best = zeros(Int, nbins)
        for i in eachindex(s)
            b = h > 0 ? clamp(floor(Int, (s[i] - s_lo) / h) + 1, 1, nbins) : 1
            if best[b] == 0 || d[i] < d[best[b]]
                best[b] = i
            end
        end
        idx = filter(!=(0), best)               # ordered by s
        # spike filter: drop points with clearly closer points on both sides
        w = CLOUD_WIN * h
        for _ in 1:length(idx)
            keep = trues(length(idx))
            for a in eachindex(idx)
                mL = Inf
                b = a - 1
                while b >= 1 && s[idx[a]] - s[idx[b]] <= w
                    mL = min(mL, d[idx[b]]);  b -= 1
                end
                mR = Inf
                b = a + 1
                while b <= length(idx) && s[idx[b]] - s[idx[a]] <= w
                    mR = min(mR, d[idx[b]]);  b += 1
                end
                keep[a] = d[idx[a]] <= (1 + CLOUD_SPIKE_REL) * max(mL, mR)
            end
            all(keep) && break
            idx = idx[keep]
        end
        return smooth_polyline(pts[idx], CLOUD_SMOOTH)
    end

    function contour_alpha_L_cloud(L)
        payloads = pmap(w -> spatial_payload(w, F, normals_F), L)
        nb = CLOUD_NBINS > 0 ? CLOUD_NBINS : length(L)
        pu, su, du = gather_cloud(payloads, :upper)
        pl, sl, dl = gather_cloud(payloads, :lower)
        global cloud_u = pu
        global cloud_l = pl
        au = cloud_frontier(pu, su, du; nbins = nb)
        al = cloud_frontier(pl, sl, dl; nbins = nb)
        global track_overrides = 0
        global track_jump_u = length(au) < 2 ? 0.0 : maximum(abs.(diff(au)))
        global track_jump_l = length(al) < 2 ? 0.0 : maximum(abs.(diff(al)))
        return au, al
    end

    function contour_alpha_L_init(L)
        au, al = contour_alpha_L_cloud(L)
        @everywhere alpha_L_u = Complex{Float64}[]
        @everywhere alpha_L_l = Complex{Float64}[]
        return au, al
    end

    # Same name the iteration loop calls.
    contour_alpha_L_conti(L) = contour_alpha_L_cloud(L)

    # -----------------------------------------------------------------------
    # v4.5 seed-and-walk tracking, kept for comparison with the cloud
    # frontier.  One balanced pmap over the whole of L, then the same
    # seed-and-walk selection v4.4 performed, done locally.
    # -----------------------------------------------------------------------
    function contour_alpha_L_walk(L)
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
        idx::Int64 = round(4*length(omega_r_base)/6)
        omega_start_target = omega_r_base[idx]
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
    function plot_cloud()
        scatter(real.(cloud_u), imag.(cloud_u); ms=1.5, msw=0, alpha=0.4, label="cloud u")
        scatter!(real.(cloud_l), imag.(cloud_l); ms=1.5, msw=0, alpha=0.4, label="cloud l")
        plot!(real.(F), imag.(F); color=:black, label="F")
        plot!(real.(alpha_L_u), imag.(alpha_L_u); lw=2, label="alpha_L_u")
        plot!(real.(alpha_L_l), imag.(alpha_L_l); lw=2, label="alpha_L_l")
    end
end

plot_omega()
plot_alpha()
plot_cloud()

L
L = L .- 0.002 * im

