# Kilian Vinzenz Wilhelm

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

    # Same rule as v4.4 couetteflow_spatial_sing_mode_comparison with
    # side_tol = 0: nearest to alpha_prev among the requested side, falling
    # back to nearest among all finite eigenvalues if that side is empty.
    function select_tracked(pl, alpha_prev, side::Symbol)
        cand = side === :upper ? pl.up : pl.lo
        if isempty(cand)
            isempty(pl.all) && error("spatial_payload: no finite eigenvalues")
            return pl.all[argmin(abs.(ComplexF64(alpha_prev) .- pl.all))]
        end
        return cand[argmin(abs.(ComplexF64(alpha_prev) .- cand))]
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

        for j in (s + 1):length(L)
            au[j] = select_tracked(payloads[j], au[j - 1], :upper)
            al[j] = select_tracked(payloads[j], al[j - 1], :lower)
        end
        for j in (s - 1):-1:1
            au[j] = select_tracked(payloads[j], au[j + 1], :upper)
            al[j] = select_tracked(payloads[j], al[j + 1], :lower)
        end
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
        phi_F = expc(zeta_alpha / (abs(alpha_F - alpha)^s_alpha + epsilon)) - 1.0
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
        d_d_alpha_r_phi_F = -zeta_alpha * (real(alpha_F) - real(alpha)) * s_alpha * abs(alpha_F - alpha)^(s_alpha - 2) / (abs(alpha_F - alpha)^s_alpha + epsilon)^2.0 * expc(zeta_alpha / (abs(alpha_F - alpha)^s_alpha + epsilon))
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
        d_d_alpha_i_phi_F = -zeta_alpha * (imag(alpha_F) - imag(alpha)) * s_alpha * abs(alpha_F - alpha)^(s_alpha - 2) / (abs(alpha_F - alpha)^s_alpha + epsilon)^2.0 * expc(zeta_alpha / (abs(alpha_F - alpha)^s_alpha + epsilon))
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
    filename = "contour_iteration_v4.7.json"
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
const ADAPT_MIN_H        = 1e-11
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
# RESUME = false  ->  original behaviour: `filename` is OVERWRITTEN with a
#                    single initial entry and the run starts from scratch.
# RESUME = true   ->  continue the run already stored in `filename`.
#                    v4.5: omega_r and adapt_level are restored from the
#                    stored entry, because the grid is no longer a constant.
# ---------------------------------------------------------------------------
const RESUME = false

begin
    if RESUME && isfile(filename)
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
        adapt_reset_history!()
        push_omega_r()
        load_on_workers()
        @everywhere begin
            normals_F = contour_normals(F)
        end
        @printf("RESUMED from iteration %d: omega_i = %.9e, N_L = %d, adapt_level = %d, %d entries in %s\n",
                resume_entry["iteration"], omega_i, length(omega_r), adapt_level,
                length(resume_array), filename)
        flush(stdout)
    else
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
        )
        current_array = Any[]
        push!(current_array, dict_to_JSON)
        open(filename, "w") do file
            write(file, JSON.json(current_array))
        end
        iteration_step += 1
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
for k = 1:2000
    global omega_i, L, alpha_L_u, alpha_L_l, alpha_i, F, omega_F, iteration_step
    global omega_r, adapt_level, adapt_nogain, adapt_done
    local dict_to_JSON, current_array

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
                    #
        alpha_i_cache = copy(alpha_i)

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
        r_min   = min(dist_u, dist_l)
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
            rhs_cap = 100.0
            rhs_j = isnan(rhs_j) ? 0.0 : clamp(rhs_j, -rhs_cap, rhs_cap)
            alpha_i_trial[j] = alpha_i[j] + local_delta_t * rhs_j
        end
        alpha_i_trial[1] = alpha_i_trial[2]
        alpha_i_trial[end] = alpha_i_trial[end-1]
        alpha_i_smooth = rolling_average_filter(alpha_i_trial, 7)
        alpha_i_smooth[1] = alpha_i_smooth[2]
        alpha_i_smooth[end] = alpha_i_smooth[end-1]
        factor = acceptance_factor.(F, Ref(alpha_L_u), Ref(alpha_L_l))
        move_raw = abs.(alpha_i_trial .- alpha_i)
        move_smooth = abs.(alpha_i_smooth .- alpha_i)
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

        if all(isfinite, alpha_i_smooth)
            global alpha_i = copy(alpha_i_smooth)
        else
            @printf("[k=%d] WARNING: alpha update produced %d non-finite values; step discarded\n",
                    k, count(!isfinite, alpha_i_smooth))
            flush(stdout)
        end
        global F = contour_F()
        load_on_workers()
        global omega_F = contour_omega_F(F)
        @printf(
            "[%04d] jump=%9.3e | dUL=%9.3e | dUF=%9.3e | dLF=%9.3e | dt=%9.3e | z=%9.3e xarg=%8.1f | gap=%9.3e (%s,%d) | NL=%3d L%d h=%8.2e gs=%5.2f miss=%8.2e | %s/%s\n",
            k, omega_jump, d_branch, dist_u, dist_l, local_delta_t, zeta_alpha, exp_arg,
            omega_gap_after, omega_bisect_status, omega_bisect_tries,
            length(omega_r), adapt_level, h_local, grid_share, pinch_miss,
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
        )
        json_str = open(filename, "r") do file
            read(file, String)
        end
        local current_array = JSON.parse(json_str)
        push!(current_array, dict_to_JSON)
        json_str = JSON.json(current_array)
        open(filename, "w") do file
            write(file, json_str)
        end
        global iteration_step += 1

        if stop_after_save
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
            adapt_push!(i_pinch, d_branch, omega_r[i_pinch], h_local)

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
iteration_step, L, alpha_L_u, alpha_L_l, F, omega_F = load_step(filename; offset=0)
println("Loaded iteration: ", iteration_step)
global omega_i = imag(L[1])
global alpha_i = imag.(F)
# v4.5: the grid comes back from the stored L, not from the constant.
global omega_r = real.(L)
push_omega_r()
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
