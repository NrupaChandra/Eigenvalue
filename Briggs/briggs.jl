# Kilian Vinzenz Wilhelm
begin
    using Distributed, Plots, BenchmarkTools, FFTW, JSON
    addprocs(32)
    w = workers()
end
begin
    @everywhere using LinearAlgebra
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
    function load_on_workers()
        names = (:L, :F, :omega_F, :alpha_L_u, :alpha_L_l)
        for pid in workers()
            for name in names
                if isdefined(Main, name)
                    data = deepcopy(getfield(Main, name))
                    fetch(@spawnat pid Core.eval(Main, :($name = $data)))
                end
            end
        end
        return nothing
    end
    #=function load_on_workers()
        for (name, data) in [(:L, L), (:F, F), (:omega_F, omega_F), (:alpha_L_u, alpha_L_u), (:alpha_L_l, alpha_L_l)]
            for pid in workers()
                @spawnat pid eval(Main, :($name = $(deepcopy(data))))
                @spawnat pid L, F, omega_F, alpha_L_u, alpha_L_l
            end
        end
        @everywhere L, F, omega_F, alpha_L_u, alpha_L_l
    end=#
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
   #@everywhere F = Vector{Complex{Float64}}[]
    @everywhere F = ComplexF64[]
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
    @everywhere function couetteflow_spatial_sing_mode_comparison(omega, alpha_approximation)
        #eigvals = couetteflow_temporal(alpha)
        eigvals, _ = couetteflow(nothing, omega, Val(:alpha_collocation), Val(:eigen))
        mask = [isfinite(real(eigval)) && isfinite(imag(eigval)) for eigval in eigvals]
        eigvals = eigvals[mask]
        diffs = abs.(Complex{Float64}(alpha_approximation) .- Complex{Float64}.(eigvals))
        index = argmin(diffs)
        eigval = eigvals[index]
        return eigval
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
    #@everywhere L = Vector{Complex{Float64}}[]
    @everywhere L = ComplexF64[]
    L = contour_L()
    #=
    @everywhere begin
        alpha_L_u = Vector{Complex{Float64}}[]
        alpha_L_l = Vector{Complex{Float64}}[]
    end=#
    @everywhere begin
    alpha_L_u = ComplexF64[]
    alpha_L_l = ComplexF64[]
    end
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
            signed_projections = [real(conj(normal) * (eigval - f)) for (f, normal) in zip(F, normals_F)]
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
    @everywhere function track_branch_pmap(L, start_index, alpha_0, direction)
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
                alpha_corrected = couetteflow_spatial_sing_mode_comparison(omega, alpha_tracked)
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
        future_u_fwd = @spawn track_branch_pmap(L, start_index, alpha_u_start, +1)
        future_u_bwd = @spawn track_branch_pmap(L, start_index, alpha_u_start, -1)
        future_l_fwd = @spawn track_branch_pmap(L, start_index, alpha_l_start, +1)
        future_l_bwd = @spawn track_branch_pmap(L, start_index, alpha_l_start, -1)
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
    global zeta_omega = 1e-5 
    #=function adapt_zeta_omega(L, omega_F)
        min_distance = minimum(abs.(L .- omega_F))
        safety = 0.003
        zeta_max = 1.0
        zeta_min = 1e-10
        zeta_alpha = clamp(safety * (min_distance + epsilon), zeta_min, zeta_max)
        return zeta_alpha 
    end=#
    global zeta_alpha = 2e-3
    function adapt_zeta_alpha(F, alpha_L_u, alpha_L_l)
        all_eigenvals = vcat(alpha_L_u, alpha_L_l)
        double_F = vcat(F, F)
        #min_distance = minimum(abs.(F' .- all_eigenvals))
        min_distance = minimum(abs.(double_F .- all_eigenvals))
        safety = 0.008
        zeta_max = 1.0
        zeta_min = 1e-10
        zeta_alpha = clamp(safety * (min_distance + epsilon), zeta_min, zeta_max)
        return zeta_alpha 
    end
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
    lambda = 2e+0
    sigma = 3e-4
    delta_t = 2e-3
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
    function complexvec_to_json(vec)
        return [Dict("re" => real(x), "im" => imag(x)) for x in vec]
    end
    filename = "contour_iteration.json"
end
###
###
###
2+2
2+2
2+2
2+2
###
###
###
###
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
#
for k = 1:200
    # state variables that must update the globals
    global omega_i, L, alpha_L_u, alpha_L_l, alpha_i, F, omega_F, iteration_step, zeta_alpha

    omega_i_vectorization = fill(omega_i, N)
    omega_i_cache = copy(omega_i_vectorization)

    println(zeta_omega)

    for j in 2:(length(omega_i_cache) - 1)
        omega_i_cache[j] = omega_i_vectorization[j] + delta_t * (
            -lambda
            + (omega_i_vectorization[j+1] - omega_i_vectorization[j-1]) /
              (omega_r[j+1] - omega_r[j-1]) * d_d_omega_r_Phi_L(L[j])
            - d_d_omega_i_Phi_L(L[j])
        )
    end

    omega_ref = maximum(imag.(omega_F))
    greater_candidates = filter(x -> x > omega_ref, omega_i_cache)

    if isempty(greater_candidates)
        omega_i = maximum(omega_i_cache)
    else
        omega_i = greater_candidates[argmin(abs.(greater_candidates .- omega_ref))]
    end

    L = contour_L()
    load_on_workers()

    @everywhere begin
        normals_F = contour_normals(F)
    end

    alpha_L_u, alpha_L_l = contour_alpha_L_conti(L)
    load_on_workers()

    alpha_i_cache = copy(alpha_i)
    zeta_alpha = adapt_zeta_alpha(F, alpha_L_u, alpha_L_l)

    accepted = false
    alpha_i_smooth = copy(alpha_i)

    for attempt in 1:100
        alpha_i_trial = copy(alpha_i)

        for j in 2:(length(alpha_i_trial) - 1)
            alpha_i_trial[j] = alpha_i[j] + delta_t * (
                (alpha_i[j+1] - alpha_i[j-1]) / (alpha_r[j+1] - alpha_r[j-1]) * d_d_alpha_r_Phi_F(F[j])
                - d_d_alpha_i_Phi_F(F[j])
                + sigma * (
                    (alpha_i[j+1] - 2.0 * alpha_i[j] + alpha_i[j-1]) /
                    ((alpha_r[j+1] - alpha_r[j]) * (alpha_r[j] - alpha_r[j-1]))
                ) / (
                    1.0 + ((alpha_i[j+1] - alpha_i[j-1]) / (alpha_r[j+1] - alpha_r[j-1]))^2.0
                )
            )
        end

        alpha_i_trial[1] = alpha_i_trial[2]
        alpha_i_trial[end] = alpha_i_trial[end-1]

        alpha_i_smooth = rolling_average_filter(alpha_i_trial, 7)
        factor = acceptance_factor.(F, Ref(alpha_L_u), Ref(alpha_L_l))
        println(factor)

        if all(abs.(alpha_i_trial .- alpha_i) .<= factor .* max.(abs.(alpha_i), 1e-12))
            alpha_i_cache = copy(alpha_i_trial)
            accepted = true
            break
        else
            zeta_alpha *= 0.5
            println("Rejected step, halving zeta_alpha -> ", zeta_alpha)
        end
    end

    println(zeta_alpha)

    alpha_i = copy(alpha_i_smooth)
    F = contour_F()
    load_on_workers()

    omega_F = contour_omega_F(F)
    load_on_workers()

    dict_to_JSON = Dict(
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

    current_array = JSON.parse(json_str)
    push!(current_array, dict_to_JSON)

    json_str = JSON.json(current_array)
    open(filename, "w") do file
        write(file, json_str)
    end

    iteration_step += 1
end
plot_omega()
plot_alpha()
plot_omega_potential()
plot_alpha_potential()

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
iteration_step, L, alpha_L_u, alpha_L_l, F, omega_F = load_step("contour_iteration.json"; offset=0)
println("Loaded iteration: ", iteration_step)
plot_omega()
plot_alpha()
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
### careful, this one truncates, and the new offset will be a different number
###
###
###
truncate_json!("contour_iteration.json"; offset=-1)
###
###
###

#test = spectral_filter(F, 0.05)
test = spectral_filter(alpha_i, 0.05)
load_on_workers()
plot_alpha()

F
plot(test)