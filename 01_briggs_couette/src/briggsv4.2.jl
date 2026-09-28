# Kilian Vinzenz Wilhelm
#
# ============================================================================
# briggsv4.2.jl
# ----------------------------------------------------------------------------
# v4.1 with five changes and no modified physics.  See CHANGES_v4.2.md.
#
#   1. L step is monotone downward and clamped to the omega_F ceiling.
#      (v4.1 accepted UPWARD steps whenever every descending candidate fell
#       below the ceiling, giving a +/-4e-3 limit cycle.)
#   2. pinch_gate(): stopping test based on the clearance and the omega_F cap,
#      not on d_branch.
#   3. d_branch is retired as a convergence measure.  It equals
#      2*sqrt(2|omega-omega_p|/|omega''|), so it measures the L grid spacing:
#      its floor here is ~1e-2 and pinch_tol = 1e-4 is unreachable by ~4 orders.
#   4. refine_pinch(): stage-2 local solver for dω/dα = 0, using only the
#      existing temporal operator.
#   5. verify_pinch(): square-root law + the alpha+/alpha- origin test.
#
# Regression target (straight contour, Couette, Re = 2000, beta = 0):
#      alpha_p = 0.572519492 - 3.038860757im
#      omega_p = 0.298649115 - 0.572829695im
#      d2omega/dalpha2 = 0.406438 - 0.138791im
# ============================================================================
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
        omega_r = range(omega_r_start, omega_r_end, length=N)
    end
    function contour_L()
        L = Complex{Float64}[omega_r[j] + omega_i * im for j in 1:N]
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
        d_alpha_u_N = alpha_L_u[N] - alpha_L_u[N-1]
        Phi_F += phi_F(alpha_F, alpha_L_u[N]) * abs(d_alpha_u_N)
        ######
        d_alpha_l_1 = alpha_L_l[2] - alpha_L_l[1]
        Phi_F += phi_F(alpha_F, alpha_L_l[1]) * abs(d_alpha_l_1)
        for j in 2:(length(alpha_L_l) - 1)
            d_alpha_l_j = 0.5 * (alpha_L_l[j+1] - alpha_L_l[j-1])
            Phi_F += phi_F(alpha_F, alpha_L_l[j]) * abs(d_alpha_l_j)
        end
        d_alpha_l_N = alpha_L_l[N] - alpha_L_l[N-1]
        Phi_F += phi_F(alpha_F, alpha_L_l[N]) * abs(d_alpha_l_N)
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
        d_alpha_u_N = alpha_L_u[N] - alpha_L_u[N-1]
        d_d_alpha_r_Phi_F += d_d_alpha_r_phi_F(alpha_F, alpha_L_u[N]) * abs(d_alpha_u_N)
        ######
        d_alpha_l_1 = alpha_L_l[2] - alpha_L_l[1]
        d_d_alpha_r_Phi_F += d_d_alpha_r_phi_F(alpha_F, alpha_L_l[1]) * abs(d_alpha_l_1)
        for j in 2:(length(alpha_L_l) - 1)
            d_alpha_l_j = 0.5 * (alpha_L_l[j+1] - alpha_L_l[j-1])
            d_d_alpha_r_Phi_F += d_d_alpha_r_phi_F(alpha_F, alpha_L_l[j]) * abs(d_alpha_l_j)
        end
        d_alpha_l_N = alpha_L_l[N] - alpha_L_l[N-1]
        d_d_alpha_r_Phi_F += d_d_alpha_r_phi_F(alpha_F, alpha_L_l[N]) * abs(d_alpha_l_N)
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
        d_alpha_u_N = alpha_L_u[N] - alpha_L_u[N-1]
        d_d_alpha_i_Phi_F += d_d_alpha_i_phi_F(alpha_F, alpha_L_u[N]) * abs(d_alpha_u_N)
        ######
        d_alpha_l_1 = alpha_L_l[2] - alpha_L_l[1]
        d_d_alpha_i_Phi_F += d_d_alpha_i_phi_F(alpha_F, alpha_L_l[1]) * abs(d_alpha_l_1)
        for j in 2:(length(alpha_L_l) - 1)
            d_alpha_l_j = 0.5 * (alpha_L_l[j+1] - alpha_L_l[j-1])
            d_d_alpha_i_Phi_F += d_d_alpha_i_phi_F(alpha_F, alpha_L_l[j]) * abs(d_alpha_l_j)
        end
        d_alpha_l_N = alpha_L_l[N] - alpha_L_l[N-1]
        d_d_alpha_i_Phi_F += d_d_alpha_i_phi_F(alpha_F, alpha_L_l[N]) * abs(d_alpha_l_N)
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
    # --- v4.2 ---------------------------------------------------------------
    omega_settle_tol = 1e-10   # clearance below which L is declared to be resting
    gate_clear_tol   = 1e-6    # clearance required by the pinch gate
    gate_hold        = 25      # iterations the ceiling must be steady for
    # ------------------------------------------------------------------------
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
    filename = "contour_iteration_v4.2.json"
    # stage-2 output goes to its OWN file.  Do not append it to the history
    # array: MATLAB's jsondecode needs every element to have the same fields,
    # and plot_contour_deformation.m would get a cell array instead of a struct.
    pinch_filename = "contour_iteration_v4.2_pinch.json"
end

begin
    iteration_step = 1
    dict_to_JSON = Dict(
        "iteration" => iteration_step,
        "L" => complexvec_to_json(L),
        "alpha_L_u" => complexvec_to_json(alpha_L_u),
        "alpha_L_l" => complexvec_to_json(alpha_L_l),
        "F" => complexvec_to_json(F),
        "omega_F" => complexvec_to_json(omega_F),
    )
    current_array = Any[]
    push!(current_array, dict_to_JSON)
    open(filename, "w") do file
        write(file, JSON.json(current_array))
    end
    iteration_step += 1

end
function branch_distance(alpha_L_u, alpha_L_l)
    return minimum(abs.(alpha_L_u .- alpha_L_l))
end

function contour_L_at(omega_i_trial)
    return ComplexF64[omega_r[j] + omega_i_trial * im for j in 1:N]
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
############################################
# v4.2 -- PINCH GATE AND STAGE-2 SOLVER    #
############################################

# ---------------------------------------------------------------------------
# The gate.  Four conditions, none of them involving d_branch:
#   1. the omega_F maximum is interior  (rules out the endpoint case, which
#      occurred in 23 of the last 1302 iterations of the v4.1 run and would
#      otherwise set a bogus ceiling)
#   2. it is a genuine local maximum with a resolved cap (the saddle crossing)
#   3. L is resting on it
#   4. the ceiling has been steady for `hold` iterations
# ---------------------------------------------------------------------------
function pinch_gate(omega_F, omega_i, ceiling_hist; clear_tol = 1e-6, hold = 25)
    Nf = length(omega_F)
    i = argmax(imag.(omega_F))
    interior = 3 <= i <= Nf - 2
    cap = interior &&
          imag(omega_F[i-1]) < imag(omega_F[i]) &&
          imag(omega_F[i]) > imag(omega_F[i+1])
    clearance = omega_i - imag(omega_F[i])
    settled = (0.0 <= clearance) && (clearance < clear_tol)
    steady = false
    if length(ceiling_hist) >= hold
        window = ceiling_hist[end-hold+1:end]
        steady = (maximum(window) - minimum(window)) < clear_tol
    end
    return (interior && cap && settled && steady), i, clearance
end

# ---------------------------------------------------------------------------
# Roots of  c[1] + c[2] s + ... + c[n] s^(n-1)  via the companion matrix.
# (LinearAlgebra only -- no extra dependency.)
# ---------------------------------------------------------------------------
function poly_roots(c::Vector{ComplexF64})
    n = length(c)
    while n > 1 && abs(c[n]) < 1e-30
        n -= 1
    end
    n <= 1 && return ComplexF64[]
    a = c[1:n] ./ c[n]
    C = zeros(ComplexF64, n - 1, n - 1)
    for i in 1:(n-2)
        C[i+1, i] = 1.0
    end
    C[:, end] .= -a[1:n-1]
    return eigvals(C)
end

# ---------------------------------------------------------------------------
# omega(alpha) on ONE branch, selected by continuation from `ref`.
# NOT argmax(imag(...)): that selection is discontinuous and corrupts the fit.
# ---------------------------------------------------------------------------
function omega_branch(alpha, ref)
    eigvals_, _ = couetteflow(ComplexF64(alpha), nothing,
                              Val(:omega_collocation), Val(:eigen))
    mask = [isfinite(real(e)) && isfinite(imag(e)) for e in eigvals_]
    ev = eigvals_[mask]
    return ev[argmin(abs.(ev .- ref))]
end

# ---------------------------------------------------------------------------
# STAGE 2.  Solve dω/dα = 0 by local polynomial fit on a circle.
#
# Why a fit and not a finite-difference Newton: the generalized eigensolve
# carries ~1e-8 noise, so a central difference with h = 2e-5 gives d2ω/dα2
# noise of order 4e-8/h^2 ~ 1e2.  Doing it that way returns |ω''| = 57.6
# instead of 0.429 -- with a convincing-looking convergence history.
# Sampling at R ~ 1e-2..1e-1 keeps the derivatives far above the noise floor.
#
# Returns (alpha_p, omega_p, d2omega_dalpha2, spread, history).
# ---------------------------------------------------------------------------
function refine_pinch(alpha_c, omega_ref; R = 0.05, M = 24, deg = 6, levels = 3)
    ap = ComplexF64(alpha_c)
    wp = ComplexF64(omega_ref)
    w2 = zero(ComplexF64)
    ac = ComplexF64(alpha_c)
    wr = ComplexF64(omega_ref)
    Rk = float(R)
    history = Tuple{Float64,ComplexF64,ComplexF64,Float64,Float64}[]

    for lev in 1:levels
        pts = ComplexF64[ac + Rk * cis(2pi * (m - 1) / M) for m in 1:M]
        vals = ComplexF64[]
        ref = wr
        for z in pts
            ref = omega_branch(z, ref)
            push!(vals, ref)
        end

        # a check that CAN fail: walking the circle must return to the start
        closure = abs(omega_branch(pts[1], vals[end]) - vals[1])
        if closure > 1e-9
            @warn "refine_pinch: continuation did not close" closure R=Rk M=M
        end

        S = ComplexF64[(z - ac) / Rk for z in pts]
        V = ComplexF64[S[m]^k for m in 1:M, k in 0:deg]
        c = V \ vals
        resid = maximum(abs.(V * c .- vals))

        dc  = ComplexF64[k * c[k+1] for k in 1:deg]              # p'(s)
        ddc = ComplexF64[k * (k - 1) * c[k+1] for k in 2:deg]    # p''(s)

        rts = filter(z -> abs(z) < 1.2, poly_roots(dc))
        isempty(rts) && error("refine_pinch: no stationary point of omega(alpha) inside R = $Rk")
        s0 = rts[argmin(abs.(rts))]

        ap = ac + Rk * s0
        wp = evalpoly(s0, c)
        w2 = evalpoly(s0, ddc) / Rk^2
        push!(history, (Rk, ap, wp, resid, closure))

        ac, wr, Rk = ap, wp, Rk / 2
    end

    spread = 0.0
    for i in 1:levels, j in 1:levels
        spread = max(spread, abs(history[i][2] - history[j][2]))
    end
    return ap, wp, w2, spread, history
end

# ---------------------------------------------------------------------------
# Newton on omega(alpha) = omega_target, warm-started from alpha_guess.
# h = 3e-3 keeps the derivative well clear of the eigensolver noise.
# ---------------------------------------------------------------------------
function alpha_from_omega(omega_target, alpha_guess, ref; h = 3e-3, iters = 12)
    z = ComplexF64(alpha_guess)
    r = ComplexF64(ref)
    for _ in 1:iters
        w0 = omega_branch(z, r)
        r = w0
        d1 = (omega_branch(z + h, w0) - omega_branch(z - h, w0)) / (2h)
        abs(d1) < 1e-10 && break
        st = -(w0 - omega_target) / d1
        if abs(st) > 0.35
            st = st * (0.35 / abs(st))
        end
        z = z + st
        abs(st) < 1e-10 && break
    end
    return z, r
end

# ---------------------------------------------------------------------------
# Verification.  Every one of these is capable of failing.
# ---------------------------------------------------------------------------
function verify_pinch(alpha_p, omega_p, w2, F, L, omega_F; verbose = true)
    report = Dict{String,Any}()

    verbose && println("\nsquare-root check at omega = omega_p + i*delta")
    verbose && println("  |a+ - a-| must equal 2|h| + O(h^3);  the midpoint offset is O(h^2)")
    verbose && @printf("  %8s %15s %14s %10s %14s %8s\n",
                       "delta", "2|h| predicted", "gap measured", "rel.err",
                       "|mid-alpha_p|", "/|h|^2")
    sqrt_ok = true
    prev_rel = Inf
    rows = Any[]
    for delta in (1e-2, 1e-3, 1e-4)
        h = sqrt(2im * delta / w2)
        r1, _ = alpha_from_omega(omega_p + 1im * delta, alpha_p + h, omega_p)
        r2, _ = alpha_from_omega(omega_p + 1im * delta, alpha_p - h, omega_p)
        gap = abs(r1 - r2)
        mid = abs(0.5 * (r1 + r2) - alpha_p)
        rel = abs(gap - 2 * abs(h)) / gap
        verbose && @printf("  %8.0e %15.6e %14.6e %10.2e %14.3e %8.3f\n",
                           delta, 2 * abs(h), gap, rel, mid, mid / abs(h)^2)
        rel < prev_rel || (sqrt_ok = false)
        prev_rel = rel
        push!(rows, Dict("delta" => delta, "gap" => gap, "predicted" => 2 * abs(h),
                         "rel_err" => rel, "midpoint_offset" => mid))
    end
    report["sqrt_law"] = rows
    report["sqrt_law_ok"] = sqrt_ok

    verbose && println("\nBriggs origin test: continue both roots as Im(omega) is raised")
    ladder = vcat(exp10.(range(-3, -1, length = 12)),
                  collect(range(0.12, 0.7, length = 14)),
                  collect(range(0.8, 3.5, length = 10)))
    h0 = sqrt(2im * 1e-3 / w2)
    z1 = alpha_p + h0
    z2 = alpha_p - h0
    r1 = omega_p
    r2 = omega_p
    maxstep = 0.0
    crossstep = 0.0
    verbose && @printf("  %13s | %26s | %26s\n", "Im w - Im w_p", "root A", "root B")
    for (n, delta) in enumerate(ladder)
        p1, p2 = z1, z2
        z1, r1 = alpha_from_omega(omega_p + 1im * delta, z1, r1)
        z2, r2 = alpha_from_omega(omega_p + 1im * delta, z2, r2)
        maxstep = max(maxstep, abs(z1 - p1), abs(z2 - p2))
        (imag(z1) > 0) != (imag(p1) > 0) && (crossstep = max(crossstep, abs(z1 - p1)))
        (imag(z2) > 0) != (imag(p2) > 0) && (crossstep = max(crossstep, abs(z2 - p2)))
        if verbose && (n == 1 || n == 12 || n == 20 || n == 26 || n == 30 || n == length(ladder))
            @printf("  %13.3f | %+12.6f%+12.6fim | %+12.6f%+12.6fim\n",
                    delta, real(z1), imag(z1), real(z2), imag(z2))
        end
    end
    origin_ok = (imag(z1) > 0) != (imag(z2) > 0)
    verbose && @printf("  largest step %.2f (informational); step containing a crossing: %.2f %s\n",
                       maxstep, crossstep, crossstep < 1.5 ? "(resolved)" : "(**refine the ladder**)")
    verbose && println("  => roots end on ", origin_ok ? "OPPOSITE" : "the SAME",
                       " sides of Im(alpha)=0 -> ",
                       origin_ok ? "genuine alpha+/alpha- Briggs pinch" :
                                   "**NOT verified as a Briggs pinch**")
    report["origin_test_ok"] = origin_ok
    report["origin_max_step"] = maxstep
    report["origin_crossing_step"] = crossstep
    report["alpha_plus_end"] = Dict("re" => real(z1), "im" => imag(z1))
    report["alpha_minus_end"] = Dict("re" => real(z2), "im" => imag(z2))

    # stage-1 quality -- none of this involves d_branch
    hF = abs(F[2] - F[1])
    dF = minimum(abs.(F .- alpha_p))
    dL = abs(imag(L[1]) - imag(omega_p))
    dC = abs(maximum(imag.(omega_F)) - imag(omega_p))
    report["min_F_to_alpha_p"] = dF
    report["min_F_to_alpha_p_over_hF"] = dF / hF
    report["abs_ImL_minus_Imomega_p"] = dL
    report["abs_maxImomegaF_minus_Imomega_p"] = dC
    if verbose
        println("\nstage-1 quality (no d_branch anywhere in here):")
        @printf("  min|F - alpha_p|                = %.3e   (= %.2f grid cells, h_F = %.3e)\n",
                dF, dF / hF, hF)
        @printf("  |Im L - Im omega_p|             = %.3e\n", dL)
        @printf("  |max Im omega_F - Im omega_p|   = %.3e\n", dC)
    end
    return report
end

# ---------------------------------------------------------------------------
# Cheap self-test of the new machinery.  Runs in well under a second; do this
# before committing to a long run.  Set to false to skip.
# ---------------------------------------------------------------------------
const RUN_V42_SELFTEST = true
if RUN_V42_SELFTEST
    let
        println("---- v4.2 self-test ----")
        # 1. poly_roots against a polynomial with known roots
        true_roots = ComplexF64[0.3 - 1.2im, -2.0 + 0.5im, 1.1 + 0.0im,
                                0.05 - 0.05im, -0.7 - 0.9im]
        cpoly = ComplexF64[1.0]
        for r in true_roots
            cnew = zeros(ComplexF64, length(cpoly) + 1)
            cnew[2:end] .+= cpoly
            cnew[1:end-1] .-= r .* cpoly
            cpoly = cnew
        end
        cpoly .*= (0.37 - 0.9im)
        got = sort(poly_roots(cpoly); by = z -> (real(z), imag(z)))
        exp_ = sort(true_roots; by = z -> (real(z), imag(z)))
        e1 = maximum(abs.(got .- exp_))
        @printf("  poly_roots       max |root error| = %.2e   %s\n",
                e1, e1 < 1e-9 ? "PASS" : "**FAIL**")

        # 2. derivative-coefficient indexing, checked against finite differences
        cc = ComplexF64[2.0 + 1im, -3.0 + 0im, 0.5 + 2im, 1.5 - 1im,
                        0.25 + 0im, -0.1 + 0.3im, 0.04 + 0im]
        dg = length(cc) - 1
        dc  = ComplexF64[k * cc[k+1] for k in 1:dg]
        ddc = ComplexF64[k * (k - 1) * cc[k+1] for k in 2:dg]
        s = 0.31 - 0.17im
        hh = 1e-6
        e2 = abs(evalpoly(s, dc) - (evalpoly(s + hh, cc) - evalpoly(s - hh, cc)) / (2hh))
        e3 = abs(evalpoly(s, ddc) -
                 (evalpoly(s + hh, cc) - 2 * evalpoly(s, cc) + evalpoly(s - hh, cc)) / hh^2)
        @printf("  p' / p'' coeffs  err = %.2e / %.2e   %s\n",
                e2, e3, (e2 < 1e-8 && e3 < 1e-3) ? "PASS" : "**FAIL**")

        # 3. pinch_gate must reject an endpoint-dominated ceiling
        wtest = ComplexF64[0.0 - (1.0 + abs(j - 50) * 0.01)im for j in 1:100]
        ok_a, pk_a, _ = pinch_gate(wtest, imag(wtest[50]) + 5e-7,
                                   fill(imag(wtest[50]), 30))
        wtest[1] = 0.0 + (imag(wtest[50]) + 1e-3)im
        ok_b, pk_b, _ = pinch_gate(wtest, imag(wtest[1]) + 5e-7,
                                   fill(imag(wtest[1]), 30))
        @printf("  pinch_gate       interior peak: fire=%s (want true, peak=%d); endpoint peak: fire=%s (want false)   %s\n",
                ok_a, pk_a, ok_b, (ok_a && !ok_b) ? "PASS" : "**FAIL**")
        println("------------------------")
        flush(stdout)
    end
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
# --- v4.2 state carried across iterations -----------------------------------
global ceiling_hist     = Float64[]      # history of max Im(omega_F)
global pinch_ready      = false          # did the gate fire?
global pinch_peak_index = 0              # argmax Im(omega_F) when it fired
global n_upward_L_steps = 0              # acceptance criterion #2: must stay 0
# ----------------------------------------------------------------------------
for k = 1:2000
    #print_iteration_header(k)
    global omega_i, L, alpha_L_u, alpha_L_l, alpha_i, F, omega_F, iteration_step
    global pinch_ready, pinch_peak_index, n_upward_L_steps
    local dict_to_JSON, current_array

    omega_i_old = omega_i
    omega_status = "none"
    omega_settled = false
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
        # v4.2: was 1e-8.  The step must be allowed to shrink below the
        # clearance, otherwise the loop aborts exactly as it converges.
        omega_dt_min = 1e-12

        min_valid_omega_candidates = 2
        max_omega_jump_factor = 40.0

        if !omega_repaired
            for omega_attempt in 1:50
                omega_i_vectorization = fill(omega_i, N)
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

                # ============================================================
                # v4.2: MONOTONE DESCENT, CLAMPED TO THE CEILING
                #
                # The j-loop above runs 2:(N-1), so omega_i_cache[1] and [N]
                # are never written and stay exactly equal to omega_i.  Those
                # two are v4's zero-length steps; v4.1's min_useful_omega_jump
                # removed them, but then the only survivors once the clearance
                # drops below lambda*omega_dt are candidates the repulsion has
                # pushed UPWARD -- and minimum(...) returned one of those and
                # the code accepted it.  Hence the +/-4e-3 limit cycle.
                #
                # Fix: drop the endpoints explicitly, and require descent.
                # ============================================================
                interior_cache = omega_i_cache[2:end-1]

                greater_candidates = filter(x -> isfinite(x) &&
                                                 x > omega_lower_bound &&
                                                 x < omega_i, interior_cache)

                # can halving omega_dt help?  Only if something tried to
                # descend but landed below the ceiling.  The SIGN of
                # (candidate - omega_i) does not depend on omega_dt, so if
                # nothing descends at all, halving is futile.
                any_descent = any(x -> isfinite(x) && x < omega_i, interior_cache)

                if length(greater_candidates) >= min_valid_omega_candidates

                    # lowest descending candidate = hug the ceiling from above
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
                            # v4.2 acceptance criterion #2: this must never fire
                            if omega_candidate > omega_i
                                n_upward_L_steps += 1
                                @printf("[k=%d] ** BUG: accepted an UPWARD L step %.3e -> %.3e\n",
                                        k, omega_i, omega_candidate)
                            end

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

                elseif any_descent && (omega_i - omega_lower_bound) > omega_settle_tol
                    # descending candidates exist but overshoot the ceiling:
                    # shrink the step.  NEVER substitute an ascending one.
                    omega_dt *= 0.5

                else
                    # L is resting on the ceiling (or the repulsion beats lambda
                    # everywhere).  Hold position, but still refresh the branches
                    # so that Phi_F stays consistent with the F that has moved
                    # since the last accepted step.
                    L_trial = contour_L_at(omega_i)
                    alpha_L_u_trial, alpha_L_l_trial = contour_alpha_L_conti(L_trial)
                    overlap_ok, overlap_reason =
                        branch_overlap_valid(alpha_L_u_trial, alpha_L_l_trial)

                    if overlap_ok
                        global L = L_trial
                        global alpha_L_u = alpha_L_u_trial
                        global alpha_L_l = alpha_L_l_trial

                        omega_settled = true
                        omega_accepted = true          # do NOT break the outer loop
                        omega_status = "settled"
                        omega_attempt_used = omega_attempt
                        omega_valid_count = 0
                        omega_jump = 0.0
                        omega_dt_used = omega_dt
                        break
                    else
                        @printf(
                            "[k=%d] STOP: branches overlap while L rests on the ceiling (%s).\n",
                            k, overlap_reason
                        )
                        @printf("[k=%d]       this is a true coalescence at the working resolution.\n", k)
                        omega_accepted = false
                        break
                    end
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
                    #
        alpha_i_cache = copy(alpha_i)

        #print_block("BRANCH DISTANCES")

        d_branch = branch_distance(alpha_L_u, alpha_L_l)
        d_contour = contour_distance(F, alpha_L_u, alpha_L_l)
        dist_u = minimum(abs.(F .- alpha_L_u))
        dist_l = minimum(abs.(F .- alpha_L_l))
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
        pinch_tol = 1e-4   # v4.2: retained for reference only -- see the stopping
                           # test below.  The convergence decision no longer uses it.

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

        # ================================================================
        # v4.2: the stopping test.
        #
        # d_branch is NOT used.  It equals 2*sqrt(2|omega-omega_p|/|omega''|),
        # so it measures how close the nearest L node happens to land to
        # omega_p, i.e. the L grid spacing -- not convergence.  With
        # N = 100 on omega_r in [0,0.5] its floor is ~1.3e-2, and reaching
        # pinch_tol = 1e-4 would need ~5e8 nodes on L.  `pinch_tol` is kept
        # only as the absolute-instability alarm (Im(omega_p) > 0 in stage 2).
        # ================================================================
        push!(ceiling_hist, omegaF_max_i)
        gate_ok, gate_peak, gate_clearance =
            pinch_gate(omega_F, omega_i, ceiling_hist;
                       clear_tol = gate_clear_tol, hold = gate_hold)

        if gate_ok
            global pinch_ready = true
            global pinch_peak_index = gate_peak
            stop_after_save = true
            stop_reason = "pinch neighbourhood reached"
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
        # v4.2 log line: the quantities that actually converge.
        # (v4.1 printed dUF/dLF as minimum(abs.(F .- alpha_L_u)), which compares
        #  F[j] with alpha_L_u[j] -- two unrelated parameterisations.  Dropped.)
        @printf(
            "[%04d] ImL=%13.9f | clear=%9.3e | peak=%3d | maxImwF=%13.9f | dUL=%9.3e | dt=%9.3e | %s/%s\n",
            k,
            omega_i,
            gate_clearance,
            gate_peak,
            maximum(imag.(omega_F)),
            d_branch,
            local_delta_t,
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

############################################
# v4.2 -- STAGE 2: SOLVE FOR THE PINCH     #
############################################
begin
    println()
    println("====================================================================")
    println("                       STAGE 2: PINCH SOLVER")
    println("====================================================================")
    @printf("upward L steps accepted outside the repair step: %d  (must be 0)\n",
            n_upward_L_steps)

    if !pinch_ready
        println("gate never fired -- stage 2 not run.")
        println("inspect the log: `clear` should fall monotonically and `peak` should")
        println("settle on one interior index.  If `peak` is 1 or N the ceiling is")
        println("being set by an endpoint of F, not by the saddle.")
    else
        i_peak = pinch_peak_index
        # seed: the argmax node of omega_F, averaged with the branch midpoint at
        # the nearest L node (an independent second opinion on the same point)
        j_star   = argmin(abs.(L .- omega_F[i_peak]))
        seed_mid = 0.5 * (alpha_L_u[j_star] + alpha_L_l[j_star])
        seed_F   = F[i_peak]
        alpha_c  = 0.5 * (seed_F + seed_mid)

        @printf("seed (argmax Im omega_F, i=%d) : %+.7f %+.7fim\n",
                i_peak, real(seed_F), imag(seed_F))
        @printf("seed (branch midpoint, L[%d])  : %+.7f %+.7fim   (d_branch there = %.4e)\n",
                j_star, real(seed_mid), imag(seed_mid), abs(alpha_L_u[j_star] - alpha_L_l[j_star]))
        @printf("seeds differ by %.3e\n", abs(seed_F - seed_mid))

        alpha_p, omega_p, w2, spread, hist =
            refine_pinch(alpha_c, omega_F[i_peak]; R = 0.05, M = 24, deg = 6, levels = 3)

        println("\nlocal fit of omega(alpha):")
        for h in hist
            @printf("  R=%8.5f  resid=%8.1e  closure=%8.1e   alpha_p=%+.9f%+.9fim   omega_p=%+.9f%+.9fim\n",
                    h[1], h[4], h[5], real(h[2]), imag(h[2]), real(h[3]), imag(h[3]))
        end
        @printf("  spread over radii: %.2e   %s\n", spread,
                spread < 1e-5 ? "CONVERGED" : "**NOT CONVERGED**")

        println("\n================================================================")
        @printf("  alpha_p       = %+.9f %+.9fim\n", real(alpha_p), imag(alpha_p))
        @printf("  omega_p       = %+.9f %+.9fim\n", real(omega_p), imag(omega_p))
        @printf("  d2w/dalpha2   = %+.6f %+.6fim     |.| = %.6f\n",
                real(w2), imag(w2), abs(w2))
        @printf("  sigma_abs     = Im(omega_p) = %+.9f  ->  %s\n",
                imag(omega_p),
                imag(omega_p) > 0 ? "ABSOLUTELY UNSTABLE" : "absolutely stable")
        println("================================================================")

        report = verify_pinch(alpha_p, omega_p, w2, F, L, omega_F)

        # regression check against the known answer for this configuration
        alpha_p_ref = 0.572519492 - 3.038860757im
        omega_p_ref = 0.298649115 - 0.572829695im
        @printf("\nregression vs the reference pinch: |dalpha| = %.2e   |domega| = %.2e\n",
                abs(alpha_p - alpha_p_ref), abs(omega_p - omega_p_ref))
        println("  (only meaningful for the straight contour, Couette, Re=2000, beta=0,")
        println("   and it is num_modes-sensitive above ~120 -- see CHANGES_v4.2.md §6)")

        report["alpha_p"]   = Dict("re" => real(alpha_p), "im" => imag(alpha_p))
        report["omega_p"]   = Dict("re" => real(omega_p), "im" => imag(omega_p))
        report["d2omega"]   = Dict("re" => real(w2),      "im" => imag(w2))
        report["spread"]    = spread
        report["sigma_abs"] = imag(omega_p)
        report["absolutely_unstable"] = imag(omega_p) > 0
        report["iteration"] = iteration_step - 1
        report["peak_index"] = i_peak
        report["upward_L_steps"] = n_upward_L_steps
        report["num_modes"] = num_modes
        report["Re"] = Re
        report["N"] = N
        open(pinch_filename, "w") do file
            write(file, JSON.json(report))
        end
        println("\nwrote ", pinch_filename)
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