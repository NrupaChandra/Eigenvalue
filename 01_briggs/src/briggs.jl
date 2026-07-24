    # Kilian Vinzenz Wilhelm
    begin
        using Distributed, Plots, BenchmarkTools, FFTW, JSON, Statistics
        addprocs(31)
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

        #=function load_on_workers()
            for (name, data) in [(:L, L), (:F, F), (:omega_F, omega_F), (:alpha_L_u, alpha_L_u), (:alpha_L_l, alpha_L_l)]
                for pid in workers()
                    @spawnat pid Core.eval(Main, :($name = $(deepcopy(data)))) #use Core.eval(...) insted of eval(...)
                    #@spawnat pid L, F, omega_F, alpha_L_u, alpha_L_l  # this line causing problems
                end
            end
        # @everywhere L, F, omega_F, alpha_L_u, alpha_L_l #this too, they do not properly sync variables.
        end=#
    
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
        lambda = 5e-1
        sigma = 3e-4
        delta_t = 5e-4
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
    function branch_distance(alpha_L_u, alpha_L_l)
        return minimum(abs.(alpha_L_u .- alpha_L_l))
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
    for k = 1:5000
        print_iteration_header(k)
        global omega_i, L, alpha_L_u, alpha_L_l, alpha_i, F, omega_F, iteration_step
        local dict_to_JSON, current_array
        begin
            #=omega_i_vectorization = fill(omega_i, N)
            omega_i_cache = copy(omega_i_vectorization)
            #global zeta_omega = adapt_zeta_omega(L, omega_F)
            println(zeta_omega)
            for j in 2:(length(omega_i_cache) - 1)
                omega_i_cache[j] =  omega_i_vectorization[j] + delta_t * ( -lambda + (omega_i_vectorization[j+1] - omega_i_vectorization[j-1]) / (omega_r[j+1] - omega_r[j-1]) * d_d_omega_r_Phi_L(L[j]) - d_d_omega_i_Phi_L(L[j]))
            end
            #
            imag(omega_F)[argmax(imag(omega_F))]
            greater_candidates = filter(x -> x > imag(omega_F)[argmax(imag(omega_F))], omega_i_cache)
            global omega_i = greater_candidates[argmin(abs.(greater_candidates .- imag(omega_F)[argmax(imag(omega_F))]))]
            #omega_i_cache[argmax(abs.(omega_i_cache))]
            #omega_i = omega_i_cache[argmax(abs.(omega_i_cache))]
            global L = contour_L()=#
            
            print_block("OMEGA-CONTOUR UPDATE")

            omegaF_max_i = maximum(imag.(omega_F))

            # Clearance between L and the highest point of omega_F.
            # Smaller = more aggressive, closer to pinch.
            omega_clearance = 1e-5
            omega_lower_bound = omegaF_max_i + omega_clearance

            println("zeta_omega = ", zeta_omega)
            println("omegaF_max_i = ", omegaF_max_i)
            println("old omega_i = ", omega_i)
            println("omega lower bound = ", omega_lower_bound)

            # ------------------------------------------------------------
            # REPAIR STEP:
            # If the updated omega_F has moved above the current L contour,
            # then the current L is already invalid.
            # In that case, lift L just above omega_F before trying to move downward.
            # ------------------------------------------------------------
            omega_repaired = false
            omega_accepted = false

            if omega_i <= omega_lower_bound
                println("WARNING: current L is already below admissible omega level.")
                println("This means omega_F has overtaken L after the last F update.")
                println("Repairing L upward and accepting repaired omega contour for this iteration.")

                global omega_i = omega_lower_bound
                global L = contour_L()
                load_on_workers()

                omega_repaired = true
                omega_accepted = true

                # After a repair, reduce delta_t so the next deformation is less aggressive.
                global delta_t = min(delta_t, 2e-4)

                println("repaired omega_i = ", omega_i)
                println("reduced delta_t = ", delta_t)
            end

            omega_dt = delta_t
            omega_dt_min = 1e-8

            # Near the pinch you can later reduce this, but keep 5 for now.
            min_valid_omega_candidates = max(5, Int(round(0.05 * N)))
            max_omega_jump_factor = 6.0

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

                    omega_i_cache = [isfinite(x) && abs(x) < 10.0 ? x : -Inf for x in omega_i_cache]
                    greater_candidates = filter(x -> isfinite(x) && x > omega_lower_bound, omega_i_cache)

                    println("omega attempt = ", omega_attempt)
                    println("omega_dt = ", omega_dt)
                    println("max omega_i_cache = ", maximum(omega_i_cache))
                    println("min omega_i_cache = ", minimum(omega_i_cache))
                    println("omega lower bound = ", omega_lower_bound)
                    println("number valid omega candidates = ", length(greater_candidates))

                    if length(greater_candidates) >= min_valid_omega_candidates

                        omega_candidate =
                            greater_candidates[
                                argmin(abs.(greater_candidates .- omega_lower_bound))
                            ]

                        omega_jump = abs(omega_candidate - omega_i)
                        max_allowed_omega_jump = max_omega_jump_factor * omega_dt

                        if omega_jump <= max_allowed_omega_jump
                            global omega_i = omega_candidate
                            omega_accepted = true

                            println("Accepted omega step.")
                            println("new omega_i = ", omega_i)
                            println("omega jump = ", omega_jump)

                            break
                        else
                            println("Rejected omega step: omega jump too large.")
                            println("omega candidate = ", omega_candidate)
                            println("omega jump = ", omega_jump)
                            println("max allowed omega jump = ", max_allowed_omega_jump)

                            omega_dt *= 0.5
                            println("Reducing omega_dt -> ", omega_dt)
                        end

                    else
                        println("Rejected omega step: too few valid omega candidates.")
                        println("required valid candidates = ", min_valid_omega_candidates)

                        omega_dt *= 0.5
                        println("Reducing omega_dt -> ", omega_dt)
                    end

                    if omega_dt < omega_dt_min
                        println("omega_dt below minimum. No safe omega step found.")
                        omega_accepted = false
                        break
                    end
                end
            end

            #=for omega_attempt in 1:50
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

                omega_i_cache = [isfinite(x) && abs(x) < 10.0 ? x : -Inf for x in omega_i_cache]
                greater_candidates = filter(x -> isfinite(x) && x > omega_lower_bound, omega_i_cache)

                #greater_candidates = filter(x -> x > omega_lower_bound, omega_i_cache)

                println("omega attempt = ", omega_attempt)
                println("omega_dt = ", omega_dt)
                println("max omega_i_cache = ", maximum(omega_i_cache))
                println("min omega_i_cache = ", minimum(omega_i_cache))
                println("omega lower bound = ", omega_lower_bound)
                println("number valid omega candidates = ", length(greater_candidates))

                #=if !isempty(greater_candidates)
                    global omega_i =
                        greater_candidates[
                            argmin(abs.(greater_candidates .- omega_lower_bound))
                        ]

                    omega_accepted = true

                    println("Accepted omega step.")
                    println("new omega_i = ", omega_i)

                    break
                else
                    omega_dt *= 0.5
                    println("Rejected omega step, reducing omega_dt -> ", omega_dt)

                    if omega_dt < omega_dt_min
                        println("omega_dt below minimum. Keeping L just above omega_F.")
                        global omega_i = omega_lower_bound
                        omega_accepted = true
                        break
                    end
                end=#

                if length(greater_candidates) >= min_valid_omega_candidates

                    omega_candidate =
                        greater_candidates[
                            argmin(abs.(greater_candidates .- omega_lower_bound))
                        ]

                    omega_jump = abs(omega_candidate - omega_i)
                    max_allowed_omega_jump = max_omega_jump_factor * omega_dt

                    if omega_jump <= max_allowed_omega_jump
                        global omega_i = omega_candidate
                        omega_accepted = true

                        println("Accepted omega step.")
                        println("new omega_i = ", omega_i)
                        println("omega jump = ", omega_jump)

                        break
                    else
                        println("Rejected omega step: omega jump too large.")
                        println("omega candidate = ", omega_candidate)
                        println("omega jump = ", omega_jump)
                        println("max allowed omega jump = ", max_allowed_omega_jump)

                        omega_dt *= 0.5
                        println("Reducing omega_dt -> ", omega_dt)
                    end

                else
                    println("Rejected omega step: too few valid omega candidates.")
                    println("required valid candidates = ", min_valid_omega_candidates)

                    omega_dt *= 0.5
                    println("Reducing omega_dt -> ", omega_dt)
                end

                if omega_dt < omega_dt_min
                    println("omega_dt below minimum. No safe omega step found.")
                    omega_accepted = false
                    break
                end
            end =#

            if !omega_accepted
                println("No safe omega accepted.")
                #global omega_i = omega_lower_bound
                break
            end

            global L = contour_L()
            load_on_workers()
            @everywhere begin
                normals_F = contour_normals(F)
            end
   
            global alpha_L_u, alpha_L_l = contour_alpha_L_conti(L)     
            load_on_workers()
            #
            alpha_i_cache = copy(alpha_i)
            global zeta_alpha = adapt_zeta_alpha(F, alpha_L_u, alpha_L_l)

            print_block("BRANCH DISTANCES")

            d_branch = branch_distance(alpha_L_u, alpha_L_l)
            d_contour = contour_distance(F, alpha_L_u, alpha_L_l)
            dist_u = minimum(abs.(F .- alpha_L_u))
            dist_l = minimum(abs.(F .- alpha_L_l))

            println("dist to upper branch = ", dist_u)
            println("dist to lower branch = ", dist_l)
            dt_min = 1e-10
            dt_max = 2e-3

            dt_start = min(delta_t, dt_max)

            #=move_safety = 0.5
            local_dist = local_contour_distances(F, alpha_L_u, alpha_L_l)
            local_move_tol = move_safety .* max.(local_dist, 1e-8)=#

            if d_contour > 0.01
                move_safety = 0.5
            else
                move_safety = 0.25
            end
            global_max_move = 0.02

            pinch_tol = 1e-4

            println("d_branch = ", d_branch, ", d_contour = ", d_contour)

            if d_branch < pinch_tol
                println("Pinch tolerance reached. Stopping.")
                break
            end

            print_block("ALPHA-CONTOUR ADAPTIVE STEP")

            accepted = false
            local_delta_t = dt_start
            for attempt in 1:100
                println("adaptive attempt =", attempt)
                println("local_delta_t = ", local_delta_t)
                flush(stdout)
                # 1. Compute RHS/update direction first
                rhs_alpha = zeros(Float64, length(alpha_i))

                for j in 2:(length(alpha_i) - 1)

                    alpha_i_r = (alpha_i[j+1] - alpha_i[j-1]) /
                                (alpha_r[j+1] - alpha_r[j-1])

                    alpha_i_rr = (alpha_i[j+1] - 2.0 * alpha_i[j] + alpha_i[j-1]) /
                                ((alpha_r[j+1] - alpha_r[j]) * (alpha_r[j] - alpha_r[j-1]))

                    rhs_alpha[j] =
                        (
                            alpha_i_r * d_d_alpha_r_Phi_F(F[j])
                            - d_d_alpha_i_Phi_F(F[j])
                            + sigma * alpha_i_rr
                        ) / (1.0 + alpha_i_r^2)
                end

                # boundary consistency for RHS
                rhs_alpha[1] = rhs_alpha[2]
                rhs_alpha[end] = rhs_alpha[end-1]

            # 2. Smooth the RHS, not the contour itself.

            if d_contour > 0.02
                rhs_smooth = smooth_alpha_curve(rhs_alpha; radius=3, passes=2)
                smooth_weight = 1.0
            elseif d_contour > 0.006
                rhs_smooth = smooth_alpha_curve(rhs_alpha; radius=2, passes=1)
                smooth_weight = 0.7
            else
                rhs_smooth = smooth_alpha_curve(rhs_alpha; radius=1, passes=1)
                smooth_weight = 0.3
            end

            rhs_alpha = (1.0 - smooth_weight) .* rhs_alpha .+ smooth_weight .* rhs_smooth

            rhs_alpha[1] = rhs_alpha[2]
            rhs_alpha[end] = rhs_alpha[end-1]

            # 3. Build trial contour.
            # Now the entire movement is proportional to local_delta_t.
            alpha_i_trial = alpha_i .+ local_delta_t .* rhs_alpha

            alpha_i_trial[1] = alpha_i_trial[2]
            alpha_i_trial[end] = alpha_i_trial[end-1]
            
                #=4. Check movement
                move = abs.(alpha_i_trial .- alpha_i)

                dt_ok = all(move .<= local_move_tol)
                zeta_ok = true

                println("attempt = ", attempt)
                println("local_delta_t = ", local_delta_t)
                println("max move = ", maximum(move))
                println("min local_move_tol = ", minimum(local_move_tol))
                println("max move/tol = ", maximum(move ./ local_move_tol))
                println("dt_ok = ", dt_ok) =#

                # 4. Check movement directionally
                F_trial = [alpha_r[j] + alpha_i_trial[j] * im for j in 1:N]

                dt_ok, toward_amounts, local_allowed_toward, full_moves =
                    directional_dt_check(
                        F,
                        F_trial,
                        alpha_L_u,
                        alpha_L_l;
                        move_safety = move_safety,
                        global_max_move = global_max_move
                    )

                zeta_ok = true

                # For diagnostics
                positive_toward = max.(toward_amounts, 0.0)
                ratio_toward = positive_toward ./ max.(local_allowed_toward, 1e-12)

                println("attempt = ", attempt)
                println("local_delta_t = ", local_delta_t)
                println("max full move = ", maximum(full_moves))
                println("max toward branch = ", maximum(positive_toward))
                println("min allowed toward = ", minimum(local_allowed_toward))
                println("max toward/allowed = ", maximum(ratio_toward))
                println("global_max_move = ", global_max_move)
                println("dt_ok = ", dt_ok)

                #=if zeta_ok && dt_ok
                    alpha_i_cache = copy(alpha_i_trial)
                    accepted = true

                    global delta_t = min(1.2 * local_delta_t, dt_max)

                    println("Accepted step with local_delta_t = ", local_delta_t, ", zeta_alpha = ", zeta_alpha)
                    println("Next delta_t = ", delta_t)
                    println("max move = ", maximum(move))
                    println("min local_move_tol = ", minimum(local_move_tol))
                    break=#

                if zeta_ok && dt_ok
                    alpha_i_cache = copy(alpha_i_trial)
                    accepted = true

                    global delta_t = min(1.2 * local_delta_t, dt_max)

                    println("Accepted step with local_delta_t = ", local_delta_t, ", zeta_alpha = ", zeta_alpha)
                    println("Next delta_t = ", delta_t)
                    println("max full move = ", maximum(full_moves))
                    println("max toward branch = ", maximum(positive_toward))
                    println("min allowed toward = ", minimum(local_allowed_toward))
                    println("max toward/allowed = ", maximum(ratio_toward))
                    break
                else
                    local_delta_t *= 0.5

                    if local_delta_t < dt_min
                        println("delta_t below minimum. Rejecting this iteration.")
                        break
                    end
                end
            end
            #=for attempt in 1:100
                alpha_i_trial = copy(alpha_i)
                for j in 2:(length(alpha_i_trial) - 1)
                    alpha_i_trial[j] = alpha_i[j] + local_delta_t * ((alpha_i[j+1] - alpha_i[j-1]) / (alpha_r[j+1] - alpha_r[j-1]) * d_d_alpha_r_Phi_F(F[j]) - d_d_alpha_i_Phi_F(F[j]) + sigma * ((alpha_i[j+1] - 2.0 * alpha_i[j] + alpha_i[j-1]) / ((alpha_r[j+1] - alpha_r[j]) * (alpha_r[j] - alpha_r[j-1]) ))  / (1.0 + ((alpha_i[j+1] - alpha_i[j-1]) / (alpha_r[j+1] - alpha_r[j-1]))^2.0))
                end
                alpha_i_trial[1] = alpha_i_trial[2]
                alpha_i_trial[end] = alpha_i_trial[end-1]
                #global alpha_i_smooth = spectral_filter(alpha_i_trial, 0.05)
                #global alpha_i_smooth = rolling_average_filter(alpha_i_trial, 7)
                #alpha_i_trial = rolling_average_filter(alpha_i_trial, 7)
                alpha_i_trial = rolling_average_filter(alpha_i_trial, 3)  #reduction of over smoothing 
                #factor = acceptance_factor.(F, Ref(alpha_L_u), Ref(alpha_L_l))  

                #zeta_ok = all(abs.(alpha_i_trial .- alpha_i) .<= factor .* max.(abs.(alpha_i), 1e-12)) # already existing condition implemented by kilian
                #dt_ok = maximum(abs.(alpha_i_trial .- alpha_i)) <= dt_tol # condition 1 for adaptive time stepping 
                #move = abs.(alpha_i_trial .- alpha_i)
                #dt_ok = all(move .<= local_move_tol)
                #println(factor)
                move = abs.(alpha_i_trial .- alpha_i)
                dt_ok = all(move .<= local_move_tol)
                zeta_ok = true   # disable old relative condition for now

                println("max move = ", maximum(move))
                println("min local_move_tol = ", minimum(local_move_tol))
                println("dt_ok = ", dt_ok)

            #= if zeta_ok && dt_ok
                    alpha_i_cache = copy(alpha_i_trial)
                    accepted = true
                    println("Accepted step with local_delta_t = ", local_delta_t, ", zeta_alpha = ", zeta_alpha)
                    println("Final accepted dt: ", local_delta_t)  #debug line for delta T
                    println("max move = ", maximum(abs.(alpha_i_trial .- alpha_i)))
                    break=#
                if zeta_ok && dt_ok
                    alpha_i_cache = copy(alpha_i_trial)
                    accepted = true

                    global delta_t = min(1.2 * local_delta_t, dt_max)

                    println("Accepted step with local_delta_t = ", local_delta_t, ", zeta_alpha = ", zeta_alpha)
                    println("Next delta_t = ", delta_t)
                    println("max move = ", maximum(move))
                    println("min local_move_tol = ", minimum(local_move_tol))
                    break    
                else
                    if !zeta_ok
                        global zeta_alpha *= 0.5
                    end

                    local_delta_t *= 0.5

                    if local_delta_t < dt_min
                        println("delta_t below minimum. Rejecting this iteration.")
                        break
                    end
                end

                #=if all(abs.(alpha_i_trial .- alpha_i) .<= factor .* max.(abs.(alpha_i), 1e-12))
                    alpha_i_cache = copy(alpha_i_trial)
                    accepted = true
                    #break
                else
                    global zeta_alpha *= 0.5
                    println("Rejected step, halving zeta_alpha -> ", zeta_alpha)

                end=#
            end=#

            if !accepted
                println("No acceptable step found. Stopping simulation instead of repeating same frozen state.")
                break
            end
            #=if !accepted
                println("No acceptable step found. Keeping old alpha_i.")
                alpha_i_cache = copy(alpha_i)
            end=#
            println(zeta_alpha)
            # ZERO DERIVATIVE BOUNDARY CONDITION #
            #alpha_i_cache[1] = alpha_i_cache[2]
            #alpha_i_cache[N] = alpha_i_cache[N-1]
            #alpha_i_smooth = spectral_filter(alpha_i_cache, 0.05)
            #global alpha_i = copy(alpha_i_smooth)
        
            global alpha_i = copy(alpha_i_cache) #Because alpha_i_smooth is just the last trial, and that last trial may have been rejected.
            global F = contour_F()
            load_on_workers()
            global omega_F = contour_omega_F(F)
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
        end    
        print_iteration_footer(k)
    end
    #plot_omega()
    #plot_alpha()
    #plot_omega_potential()
    #plot_alpha_potential()

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
    ### careful, this one truncates, and the new offset will be a different number
    ###
    ###
    ###
    truncate_json!("contour_iteration.json"; offset=-1)
    ###
    ###
    ###

    #test = spectral_filter(F, 0.05)
    load_on_workers()
    #plot_alpha()

    #F
    #plot(test)