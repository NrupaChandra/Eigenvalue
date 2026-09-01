# Kilian Vinzenz Wilhelm

# ---------------------------------------------------------------------------
# v4.4  -- close the omega-plane gap between L and omega_F.
#
# PROBLEM in v4.3.  The omega update never aims at omega_F.  It builds one
# candidate height per L point,

# CHANGES
#   1. BISECTION ENDGAME   After the existing candidate step,
#      drive omega_i down towards omega_lower_bound by bisection, accepting the
#      lowest height whose contour_alpha_L_conti + branch_overlap_valid still
#      succeed.  Halves the gap per attempt: 1e-4 -> 1e-8 in ~14 tries.
#      Switch: OMEGA_BISECT (set false to reproduce v4.3 behaviour exactly).
#   2. omega_clearance 2e-7 -> 1e-9, so the target is not the limiter.
#   3. min_useful_omega_jump and omega_dt_min 1e-8 -> 1e-12.  A step below 1e-8
#      was rejected as "useless" and the run gave up, which made converging TO
#      a 1e-8 gap impossible by construction.
#   4. L GRID: 50 points per base cell over omega_r_base[54 .. 62]
#      (dw_r = 1.031e-4, N_L = 484 vs 478).  omega_p is at omega_r = 0.2986491,
#      13 % into cell [60,61], and the candidates that decide each step sit at
#      omega_r in [0.285, 0.288]; base[54..62] = [0.2677, 0.3081] covers both.
#   5. JSON: "omega_F_at_L" saved next to "omega_F".  omega_F is recomputed
#      AFTER the alpha update, so the stored (L, omega_F) pair is one step out
#      of sync and the plotted gap is not the one the algorithm controlled.
#      omega_F_at_L is omega_F as it stood when L was placed.
#   6. Per-iteration log line reports the achieved gap and what limited it.
#
# NOTE the dead term in the candidate loop (kept, see the comment there):
# omega_i_vectorization is constant, so its central difference is identically
# zero and the d_d_omega_r_Phi_L coupling never contributes.
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
        omega_r_start = 0.0
        omega_r_end = 0.5
        omega_i = 0.0

        # Original coarse discretization
        omega_r_base = collect(range(omega_r_start, omega_r_end, length=100))

        # ------------------------------------------------------------------
        # v4.4 local refinement: PTS_PER_CELL points spanning each base cell
        # between omega_r_base[REFINE_LO] and omega_r_base[REFINE_HI].
        #
        #   base spacing                    = 5.0505e-3
        #   base[54] = 0.2676768   base[62] = 0.3080808
        #   refined spacing                 = 5.0505e-3 / 49 = 1.0307e-4
        #   N_L                             = 53 + 393 + 38 = 484
        #
        # omega_p sits at omega_r = 0.2986491 (13 % into base cell [60,61]) and
        # the candidates that decide the omega step sit at omega_r in
        # [0.285, 0.288]; base[54..62] = [0.2677, 0.3081] brackets both.
        # v4.3 was REFINE_LO=47, REFINE_HI=67, PTS_PER_CELL=21 (200 pts over
        # 10 cells), i.e. dw_r = 2.5379e-4.
        # ------------------------------------------------------------------
        REFINE_LO    = 54
        REFINE_HI    = 62
        PTS_PER_CELL = 50      # points spanning one base cell, endpoints included

        # Build the refined block cell by cell, dropping the duplicated left
        # endpoint of every cell after the first.
        omega_refined = collect(range(omega_r_base[REFINE_LO],
                                      omega_r_base[REFINE_LO + 1],
                                      length = PTS_PER_CELL))
        for c in (REFINE_LO + 1):(REFINE_HI - 1)
            cell = collect(range(omega_r_base[c], omega_r_base[c + 1],
                                 length = PTS_PER_CELL))
            append!(omega_refined, cell[2:end])
        end

        # Assemble complete omega discretization.
        # omega_refined already carries base[REFINE_LO] and base[REFINE_HI].
        omega_r = vcat(
            omega_r_base[1:(REFINE_LO - 1)],
            omega_refined,
            omega_r_base[(REFINE_HI + 1):end]
        )

        @assert issorted(omega_r)
        @assert minimum(diff(omega_r)) > 0

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
    alpha_L_u, alpha_L_l = contour_alpha_L_init(L)
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

        # Same physical starting omega as old 100-point grid
        omega_start_target = omega_r_base[25]

        # Find that location on the new refined grid
        start_index = argmin(abs.(real.(L) .- omega_start_target))
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
        global omega_r = vcat(
            omega_r_base[1:(REFINE_LO - 1)],
            omega_refined,
            omega_r_base[(REFINE_HI + 1):end]
        )
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

    # -----------------------------------------------------------------------
    # v4.3a OVERFLOW GUARD  (fixes the LAPACKException(150) crash at k = 122)
    #
    # exp(zeta/(d^s + epsilon)) overflows once d = |alpha_F - alpha_L| < 7.5e-4.
    # Both d_d_alpha_r_Phi_F and d_d_alpha_i_Phi_F carry the same exp factor and
    # enter rhs_j with opposite signs, so the numerator becomes (+Inf) + (-Inf)
    # = NaN, and clamp(NaN, -100, 100) === NaN passes it straight through into
    # alpha_i -> F -> eigen(A, B) -> LAPACK INFO = N = 150.
    #
    # The cap bounds the exponent only.  It starts to engage at d = 1.0e-3; the
    # largest argument reached anywhere in the 121 iterations that ran was 284.9
    # (at d = 1.185e-3), so the whole existing trajectory is unchanged.  Below
    # d = 9.5e-3 rhs_j is already pinned to rhs_cap = 100, so the cap is
    # invisible to the dynamics as well - it only stops the arithmetic from
    # leaving the finite range.
    #
    # The omega-side functions are deliberately NOT capped: their exponent runs
    # to ~1e4 by design, and the resulting NaN is already sanitised downstream
    # by the `isfinite(x) && abs(x) < 10.0` filter on the omega candidates.
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

            # keep endpoint behavior consistent
            y[1] = y[2]
            y[end] = y[end-1]
        end
        return y
    end
    function complexvec_to_json(vec)
        return [Dict("re" => real(x), "im" => imag(x)) for x in vec]
    end
    filename = "contour_iteration_v4.4.json"
end

# ---------------------------------------------------------------------------
# v4.4 omega-descent controls
# ---------------------------------------------------------------------------
const OMEGA_BISECT        = true    # false -> exactly the v4.3 omega update
const OMEGA_GAP_TARGET    = 1e-8    # wanted omega_i - max(imag(omega_F))
const OMEGA_BISECT_ENGAGE = 1e-3    # only bisect once the gap is below this
const OMEGA_BISECT_MAX    = 25      # probes per iteration (each = 1 branch track)
const OMEGA_BISECT_TOL    = 1e-12   # stop when the bracket is this narrow
# COST.  A probe costs one contour_alpha_L_conti, the same call the v4.3 omega
# step already made once or twice per iteration.  When the bound is admissible
# the first probe takes it and the endgame costs ONE extra call per iteration.
# Probes only pile up if there is a wall, and then every such iteration prints
# a "bisect wall ..." line.  25 probes take 1e-3 -> 3e-11, i.e. well past
# OMEGA_GAP_TARGET; if the log shows persistent walls and the run is crawling,
# drop this to ~8 (still 1e-3 -> 4e-6, ~75x better than v4.3) rather than
# raising it.
# Optional extra guard: reject a trial height whose alpha branches come closer
# to F than this.  0.0 = disabled, which is deliberate for the first v4.4 run:
# we want the bisection to find the REAL wall (branch tracking / Phi_F near
# field) rather than an imposed one.  The achieved F-branch distance is logged
# and stored as "d_contour_at_L" every iteration, so if it collapses, set this.
const OMEGA_BISECT_MIN_DFC = 0.0
# OMEGA_BISECT_ENGAGE also caps how far L can fall in one iteration, so the
# early large-scale descent -- where F still has to deform a lot -- runs
# exactly as in v4.3 and only the endgame is aimed.

# ---------------------------------------------------------------------------
# RESUME = false  ->  original behaviour: `filename` is OVERWRITTEN with a
#                    single initial entry and the run starts from scratch.
# RESUME = true   ->  continue the run already stored in `filename`.  Nothing
#                    is overwritten; the loop appends from the last entry on.
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
        load_on_workers()
        @everywhere begin
            normals_F = contour_normals(F)
        end
        @printf("RESUMED from iteration %d: omega_i = %.9e, %d entries in %s\n",
                resume_entry["iteration"], omega_i, length(resume_array), filename)
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
            # same fields as the in-loop entries, so every entry has the same
            # schema.  Before the loop L has not been placed against omega_F
            # yet, so omega_F_at_L is just omega_F.
            "omega_F_at_L" => complexvec_to_json(omega_F),
            "omega_gap" => omega_i - maximum(imag.(omega_F)),
            "omega_bisect_status" => "initial",
            "omega_bisect_tries" => 0,
            # null here: branch_distance / contour_distance are defined below
            # this block, and the distances carry no meaning for the pre-loop
            # state anyway.  Every in-loop entry has numbers.
            "d_branch" => nothing,
            "d_contour" => nothing,
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
            "upper/lower overlap: d=%.3e at upper[%d], lower[%d]",
            min_dist, min_i, min_j
        )
    end

    return true, @sprintf(
        "no overlap: min upper/lower distance=%.3e",
        min_dist
    )
end

# ---------------------------------------------------------------------------
# v4.4: is a trial height for L admissible?
#
# Builds the trial L, tracks both alpha branches along it and applies the same
# acceptance test the v4.3 omega step used (finite branches, no upper/lower
# overlap).  Returns the trial contours as well so a successful probe is not
# recomputed.  This is the expensive call in the bisection -- one
# contour_alpha_L_conti per probe -- which is why the bisection only engages
# once the gap is already small (see OMEGA_BISECT_ENGAGE).
# ---------------------------------------------------------------------------
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
for k = 1:2000
    #print_iteration_header(k)
    global omega_i, L, alpha_L_u, alpha_L_l, alpha_i, F, omega_F, iteration_step
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
        
        #print_block("OMEGA-CONTOUR UPDATE")

        # omega_F as it stands NOW, i.e. the one L is about to be placed
        # against.  Stored in the JSON as "omega_F_at_L"; the "omega_F" field
        # is recomputed after the alpha update and is one step newer.
        omega_F_at_L = copy(omega_F)

        omegaF_max_i = maximum(imag.(omega_F))
        # v4.4: 2e-7 -> 1e-9.  In v4.3 this was never the binding constraint
        # (the candidate spacing at the bound was ~4.5e-4, 1500x larger), but
        # with the bisection below it becomes the actual floor, so it has to
        # sit under OMEGA_GAP_TARGET.
        omega_clearance = 1e-9
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
        # v4.4: 1e-8 -> 1e-12.  Converging TO a 1e-8 gap needs steps smaller
        # than 1e-8; with the old floors the run rejected them as "useless"
        # and then gave up with "No safe omega step found".
        omega_dt_min = 1e-12

        min_valid_omega_candidates = 2
        max_omega_jump_factor = 40.0

        min_useful_omega_jump = 1e-12

        if !omega_repaired
            for omega_attempt in 1:50
                # NOTE (v4.4, unchanged on purpose): omega_i_vectorization is
                # a constant vector and is never written inside this loop --
                # only omega_i_cache is -- so the central difference
                # omega_i_vectorization[j+1] - omega_i_vectorization[j-1] is
                # identically zero and the d_d_omega_r_Phi_L coupling below
                # contributes nothing.  Reconstructing the v4.3 steps with that
                # term dropped reproduces the stored omega_i exactly.  Left in
                # place so v4.4 differs from v4.3 only by the bisection; if L
                # is ever meant to deform rather than translate, this is where
                # it would have to be fixed.
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
                
                # acceptance condition for L 

                omega_i_cache = [isfinite(x) && abs(x) < 10.0 ? x : -Inf for x in omega_i_cache]
                greater_candidates = filter(x -> isfinite(x) && x > omega_lower_bound &&  abs(x - omega_i) > min_useful_omega_jump, omega_i_cache)

                if length(greater_candidates) >= min_valid_omega_candidates

                    # This keeps L close to omega_F, but avoids accepting zero movement.
                    omega_candidate = minimum(greater_candidates)

                    omega_jump = abs(omega_candidate - omega_i)
                    max_allowed_omega_jump = max_omega_jump_factor * omega_dt
                    if omega_jump <= max_allowed_omega_jump

                        # Build trial L and trial branches BEFORE accepting omega
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
                        "[k=%d] omega rejected: valid/useful=%d < required=%d | lower=%.6e | omega_i=%.6e | reducing dt_omega %.2e -> %.2e\n",
                        k,
                        length(greater_candidates),
                        min_valid_omega_candidates,
                        omega_lower_bound,
                        omega_i,
                        omega_dt,
                        0.5 * omega_dt
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

            global alpha_L_u, alpha_L_l = contour_alpha_L_conti(L)
            load_on_workers()
        else
            load_on_workers()

            @everywhere begin
                normals_F = contour_normals(F)
            end
        end

        # -------------------------------------------------------------------
        # v4.4 BISECTION ENDGAME
        #
        # The candidate step above can only land on one of N_L pre-computed
        # heights, and near the bound those sit ~4.5e-4 apart, so it parks a
        # few 1e-4 above omega_F no matter how long it runs (v4.3: 963
        # iterations, best gap 1.24e-5 at k=401, ending at 3.0e-4).
        #
        # Here we aim instead: take the LOWEST omega_i whose branch tracking
        # still succeeds.  omega_lower_bound is the deepest legal height and
        # the current omega_i is known good, so bisect between them.
        #   - probe omega_lower_bound first: if it is admissible we are done in
        #     one branch track, which is the steady-state cost once converged
        #   - otherwise halve the bracket; ~17 probes take 1e-3 -> 1e-8
        # A failed probe is not an error: it means the branch tracker cannot
        # follow L that close to omega_F, and the bracket then reports where
        # that wall is (omega_bisect_status = "wall").
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
                # first probe goes all the way down; after that, bisect
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
                    isempty(omega_wall_reason) ? "-" : omega_wall_reason
                )
                flush(stdout)
            end
        end
                    #
        alpha_i_cache = copy(alpha_i)

        #print_block("BRANCH DISTANCES")

        d_branch = branch_distance(alpha_L_u, alpha_L_l)
        d_contour = contour_distance(F, alpha_L_u, alpha_L_l)
        dist_u = minimum(
            minimum(abs.(f .- alpha_L_u))
            for f in F
        )
        dist_l = minimum(
            minimum(abs.(f .- alpha_L_l))
            for f in F
        )
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
        pinch_tol = 1e-4

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
            stop_reason = "pinch tolerance reached"
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
            # clamp(NaN, ...) === NaN, so NaN has to be caught separately
            rhs_j = isnan(rhs_j) ? 0.0 : clamp(rhs_j, -rhs_cap, rhs_cap)

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
       
        # last line of defence: a non-finite alpha_i would reach eigen(A, B) as
        # a NaN operator.  Discard the step instead of poisoning the contour.
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
            "[%04d] jump=%9.3e | dUL=%9.3e | dUF=%9.3e | dLF=%9.3e | dt=%9.3e | ζ=%9.3e | gap=%9.3e (%s,%d) | %s/%s\n",
            k,
            omega_jump,
            d_branch,
            dist_u,
            dist_l,
            local_delta_t,
            zeta_common,
            omega_gap_after,
            omega_bisect_status,
            omega_bisect_tries,
            omega_status,
            alpha_status
        )
        flush(stdout)
        load_on_workers()
        local dict_to_JSON = Dict(
            "iteration" => iteration_step,
            "L" => complexvec_to_json(L),
            "alpha_L_u" => complexvec_to_json(alpha_L_u),
            "alpha_L_l" => complexvec_to_json(alpha_L_l),
            "F" => complexvec_to_json(F),
            "omega_F" => complexvec_to_json(omega_F),
            # v4.4: omega_F as it stood when L was placed.  "omega_F" above is
            # recomputed after the alpha update, so (L, omega_F) is one step
            # out of sync and the gap read from it is not the controlled one.
            "omega_F_at_L" => complexvec_to_json(omega_F_at_L),
            "omega_gap" => omega_gap_after,
            "omega_bisect_status" => omega_bisect_status,
            "omega_bisect_tries" => omega_bisect_tries,
            # both measured with the post-bisection L and the pre-alpha-update
            # F, i.e. the state at the moment L was placed
            "d_branch" => d_branch,
            "d_contour" => d_contour,
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
                k, stop_reason, d_branch, d_contour, omega_i
            )
            flush(stdout)
            break
        end
    end    
    #print_iteration_footer(k)
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