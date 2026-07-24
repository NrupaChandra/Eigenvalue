    # Kilian Vinzenz Wilhelm
    begin
        using Distributed, Plots, BenchmarkTools, FFTW, JSON, Printf, Statistics
        addprocs(7)
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
        #@everywhere F = Vector{Complex{Float64}}[]
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
        @everywhere function couetteflow_spatial_pair_comparison(
            omega,
            alpha_u_guess,
            alpha_l_guess,
            F,
            normals_F;
            side_tol = 1e-10
        )
            eigvals, _ =
                couetteflow(
                    nothing,
                    omega,
                    Val(:alpha_collocation),
                    Val(:eigen)
                )

            mask = [
                isfinite(real(eigval)) &&
                isfinite(imag(eigval))
                for eigval in eigvals
            ]

            eigvals = ComplexF64.(eigvals[mask])

            if length(eigvals) < 2
                error("Fewer than two finite spatial eigenvalues found at omega = $omega")
            end

            # Track branch identity only by continuity from the previous point.
            # F and normals_F stay in the signature for compatibility, but they
            # no longer decide which eigenvalue is allowed to be selected.
            distances_u = abs.(eigvals .- ComplexF64(alpha_u_guess))
            index_u = argmin(distances_u)

            distances_l = abs.(eigvals .- ComplexF64(alpha_l_guess))
            distances_l[index_u] = Inf  # keep the two tracked roots distinct
            index_l = argmin(distances_l)

            return eigvals[index_u], eigvals[index_l]
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
            omega_r_end = 1.0     
            omega_i = 0.0
            omega_r = range(omega_r_start, omega_r_end, length=N)
        end
        function contour_L()
            omega_f = 0.5
            r = 0.01

            # Fixed number of points in every segment
            n_hor_left  = 110
            n_ver_left  = 25
            n_arc       = 40
            n_ver_right = 25
            n_hor_right = 110

            # Horizontal left
            seg1 = ComplexF64.(
                collect(LinRange(
                    omega_r_start,
                    omega_f - r,
                    n_hor_left
                )) .+ im * omega_i
            )

            # Vertical upward
            seg2 = ComplexF64.(
                (omega_f - r) .+
                im .* collect(LinRange(
                    omega_i,
                    0.0,
                    n_ver_left
                ))
            )

            # Upper semicircle
            theta = collect(LinRange(pi, 0.0, n_arc))

            seg3 = ComplexF64.(
                omega_f .+ r .* cis.(theta)
            )

            # Vertical downward
            seg4 = ComplexF64.(
                (omega_f + r) .+
                im .* collect(LinRange(
                    0.0,
                    omega_i,
                    n_ver_right
                ))
            )

            # Horizontal right
            seg5 = ComplexF64.(
                collect(LinRange(
                    omega_f + r,
                    omega_r_end,
                    n_hor_right
                )) .+ im * omega_i
            )

            return vcat(
                seg1[1:end-1],
                seg2[1:end-1],
                seg3[1:end-1],
                seg4[1:end-1],
                seg5
            )
        end
        #@everywhere L = Vector{Complex{Float64}}[]
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
        global alpha_L_u, alpha_L_l = contour_alpha_L_init(L)
        load_on_workers()

    #######################
    # Eigenvalue tracking #
    #######################

    function contour_alpha_L_conti(
        L,
        old_L,
        old_alpha_L_u,
        old_alpha_L_l
    )
        nL = length(L)

        if isempty(old_L) ||
        length(old_alpha_L_u) != nL ||
        length(old_alpha_L_l) != nL

            return contour_alpha_L_init(L)
        end

        # Start on the smooth left horizontal part of L
        start_index =
            clamp(
                floor(Int, nL / 4),
                1,
                nL
            )

        new_alpha_u =
            Vector{ComplexF64}(undef, nL)

        new_alpha_l =
            Vector{ComplexF64}(undef, nL)

        # Correct the starting eigenvalues
        alpha_u_start,
        alpha_l_start =
            couetteflow_spatial_pair_comparison(
                L[start_index],
                old_alpha_L_u[start_index],
                old_alpha_L_l[start_index],
                F,
                normals_F
            )

        new_alpha_u[start_index] =
            alpha_u_start

        new_alpha_l[start_index] =
            alpha_l_start

        # Track forward using the neighbouring eigenvalue as the guess
        alpha_u_previous =
            alpha_u_start

        alpha_l_previous =
            alpha_l_start

        if start_index < nL
            for j in (start_index + 1):nL
                alpha_u_previous,
                alpha_l_previous =
                    couetteflow_spatial_pair_comparison(
                        L[j],
                        alpha_u_previous,
                        alpha_l_previous,
                        F,
                        normals_F
                    )

                new_alpha_u[j] =
                    alpha_u_previous

                new_alpha_l[j] =
                    alpha_l_previous
            end
        end

        # Track backward from the same starting point
        alpha_u_previous =
            alpha_u_start

        alpha_l_previous =
            alpha_l_start

        if start_index > 1
            for j in (start_index - 1):-1:1
                alpha_u_previous,
                alpha_l_previous =
                    couetteflow_spatial_pair_comparison(
                        L[j],
                        alpha_u_previous,
                        alpha_l_previous,
                        F,
                        normals_F
                    )

                new_alpha_u[j] =
                    alpha_u_previous

                new_alpha_l[j] =
                    alpha_l_previous
            end
        end

        return new_alpha_u, new_alpha_l
    end
    end
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

            min_distance =
                minimum(
                    abs.(
                        reshape(F, :, 1) .-
                        reshape(all_eigenvals, 1, :)
                    )
                )

            # Keep the largest exponential argument approximately bounded:
            # zeta_alpha / (d^s_alpha + epsilon) <= target_exponent.
            target_exponent = 0.20
            zeta_max = 2e-3
            zeta_min = 1e-12

            return clamp(
                target_exponent *
                (min_distance^s_alpha + epsilon),
                zeta_min,
                zeta_max
            )
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
        lambda = 5e-1
        sigma = 3e-4
        # Fixed-size contour relaxation with simple backtracking only when
        # the proposed F contour is non-finite or actually intersects a branch.
        delta_t_alpha = 5e-4
        delta_t_alpha_min = 1e-8
        delta_t_alpha_max = 2e-3
        n_alpha_substeps = 5
        max_alpha_backtracks = 12

        delta_t_omega = 2e-3
        delta_t_omega_min = 1e-8
        delta_t_omega_max = 2e-3

        # This guards only the L update. It no longer blocks the F relaxation.
        branch_jump_tol = 2e-2
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
        filename = "contour_iteration_v3.1.json"
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
        return minimum(
            abs.(
                reshape(alpha_L_u, :, 1) .-
                reshape(alpha_L_l, 1, :)
            )
        )
    end

    function local_contour_distances(F, alpha_L_u, alpha_L_l)
        all_branches = vcat(alpha_L_u, alpha_L_l)

        return [
            minimum(abs.(f .- all_branches))
            for f in F
        ]
    end

    function contour_distance(F, alpha_L_u, alpha_L_l)
        return minimum(
            local_contour_distances(F, alpha_L_u, alpha_L_l)
        )
    end
    function maximum_branch_change(
        new_branch,
        old_branch
    )
        if length(new_branch) != length(old_branch)
            error(
                "New and old branches have different lengths."
            )
        end

        changes = abs.(new_branch .- old_branch)

        return maximum(changes)
    end
    function branch_side_valid(
        F_trial,
        alpha_L_u,
        alpha_L_l;
        side_tol = 1e-10
    )
        normals_trial = contour_normals(F_trial)

        upper_min_projection = Inf
        lower_max_projection = -Inf

        for alpha in alpha_L_u
            distances = abs.(alpha .- F_trial)
            idx_min = argmin(distances)

            projection =
                real(
                    conj(normals_trial[idx_min]) *
                    (alpha - F_trial[idx_min])
                )

            upper_min_projection =
                min(upper_min_projection, projection)
        end

        for alpha in alpha_L_l
            distances = abs.(alpha .- F_trial)
            idx_min = argmin(distances)

            projection =
                real(
                    conj(normals_trial[idx_min]) *
                    (alpha - F_trial[idx_min])
                )

            lower_max_projection =
                max(lower_max_projection, projection)
        end

        valid =
            upper_min_projection >= -side_tol &&
            lower_max_projection <= side_tol

        return valid,
               upper_min_projection,
               lower_max_projection
    end

    # A genuine topology failure means that the proposed F polyline
    # intersects one of the tracked spatial-branch polylines.
    cross2(a::Complex, b::Complex) = real(a) * imag(b) - imag(a) * real(b)

    function segments_intersect(z1, z2, w1, w2; tol = 1e-12)
        r = z2 - z1
        s = w2 - w1
        denominator = cross2(r, s)

        # Parallel or nearly parallel segments are not counted as a crossing.
        if abs(denominator) <= tol
            return false
        end

        t = cross2(w1 - z1, s) / denominator
        u = cross2(w1 - z1, r) / denominator

        return (-tol <= t <= 1.0 + tol) &&
               (-tol <= u <= 1.0 + tol)
    end

    function contour_crosses_branch(F_trial, branch)
        for i in 1:(length(F_trial) - 1)
            for j in 1:(length(branch) - 1)
                if segments_intersect(
                    F_trial[i],
                    F_trial[i + 1],
                    branch[j],
                    branch[j + 1]
                )
                    return true
                end
            end
        end
        return false
    end

    for k = 1:400
        global omega_i, L, alpha_L_u, alpha_L_l, alpha_i, F, omega_F
        global iteration_step, zeta_alpha, delta_t_omega, delta_t_alpha
        local dict_to_JSON, current_array, json_str

        ########################
        # 1. Try moving L once #
        ########################

        old_omega_i = omega_i
        old_L = copy(L)
        old_alpha_L_u = copy(alpha_L_u)
        old_alpha_L_l = copy(alpha_L_l)

        omegaF_max_i = maximum(imag.(omega_F))
        omega_clearance = 2e-7
        omega_lower_bound = omegaF_max_i + omega_clearance

        omega_i_cache = [
            omega_i +
            delta_t_omega * (
                -lambda -
                d_d_omega_i_Phi_L(
                    ComplexF64(omega_r_value, omega_i)
                )
            )
            for omega_r_value in omega_r
        ]

        greater_candidates =
            filter(
                x -> isfinite(x) && x > omega_lower_bound,
                omega_i_cache
            )

        if omega_i <= omega_lower_bound
            omega_i = omega_lower_bound
        elseif !isempty(greater_candidates)
            omega_i =
                greater_candidates[
                    argmin(abs.(greater_candidates .- omega_lower_bound))
                ]
        end

        L = contour_L()
        load_on_workers()

        @everywhere begin
            normals_F = contour_normals(F)
        end

        omega_update_accepted = false

        try
            candidate_alpha_L_u, candidate_alpha_L_l =
                contour_alpha_L_conti(
                    L,
                    old_L,
                    old_alpha_L_u,
                    old_alpha_L_l
                )

            jump_u = maximum_branch_change(
                candidate_alpha_L_u,
                old_alpha_L_u
            )

            jump_l = maximum_branch_change(
                candidate_alpha_L_l,
                old_alpha_L_l
            )

            if max(jump_u, jump_l) > branch_jump_tol
                @printf(
                    "[k=%04d] L HELD: large branch change | upper=%.3e | lower=%.3e | F will still relax\n",
                    k,
                    jump_u,
                    jump_l
                )

                omega_i = old_omega_i
                L = old_L
                alpha_L_u = old_alpha_L_u
                alpha_L_l = old_alpha_L_l
                delta_t_omega =
                    max(0.5 * delta_t_omega, delta_t_omega_min)
            else
                alpha_L_u = candidate_alpha_L_u
                alpha_L_l = candidate_alpha_L_l
                omega_update_accepted = true
            end
        catch err
            @printf(
                "[k=%04d] L HELD: branch tracking error; F will still relax\n",
                k
            )
            showerror(stdout, err)
            println()

            omega_i = old_omega_i
            L = old_L
            alpha_L_u = old_alpha_L_u
            alpha_L_l = old_alpha_L_l
            delta_t_omega =
                max(0.5 * delta_t_omega, delta_t_omega_min)
        end

        load_on_workers()

        ############################################
        # 2. Always relax F while the branches hold #
        ############################################

        zeta_alpha =
            adapt_zeta_alpha(
                F,
                alpha_L_u,
                alpha_L_l
            )

        total_alpha_move = 0.0
        accepted_substeps = 0

        for substep in 1:n_alpha_substeps
            local_delta_t = delta_t_alpha
            accepted_alpha = false

            for attempt in 1:max_alpha_backtracks
                raw_increment = zeros(Float64, length(alpha_i))

                for j in 2:(length(alpha_i) - 1)
                    alpha_slope =
                        (alpha_i[j + 1] - alpha_i[j - 1]) /
                        (alpha_r[j + 1] - alpha_r[j - 1])

                    alpha_curv =
                        (alpha_i[j + 1] - 2.0 * alpha_i[j] + alpha_i[j - 1]) /
                        (
                            (alpha_r[j + 1] - alpha_r[j]) *
                            (alpha_r[j] - alpha_r[j - 1])
                        )

                    rhs_j =
                        (
                            alpha_slope *
                            d_d_alpha_r_Phi_F(F[j])
                            -
                            d_d_alpha_i_Phi_F(F[j])
                            +
                            sigma * alpha_curv
                        ) /
                        (1.0 + alpha_slope^2)

                    # No arbitrary displacement clamp.
                    raw_increment[j] = local_delta_t * rhs_j
                end

                raw_increment[1] = raw_increment[2]
                raw_increment[end] = raw_increment[end - 1]

                smooth_increment =
                    rolling_average_filter(raw_increment, 7)

                alpha_i_trial = alpha_i .+ smooth_increment
                alpha_i_trial[1] = alpha_i_trial[2]
                alpha_i_trial[end] = alpha_i_trial[end - 1]

                finite_trial = all(isfinite, alpha_i_trial)

                F_trial =
                    ComplexF64[
                        alpha_r[j] + alpha_i_trial[j] * im
                        for j in eachindex(alpha_r)
                    ]

                crosses_upper =
                    finite_trial &&
                    contour_crosses_branch(F_trial, alpha_L_u)

                crosses_lower =
                    finite_trial &&
                    contour_crosses_branch(F_trial, alpha_L_l)

                if finite_trial && !crosses_upper && !crosses_lower
                    substep_move = maximum(abs.(alpha_i_trial .- alpha_i))
                    total_alpha_move = max(total_alpha_move, substep_move)

                    alpha_i = copy(alpha_i_trial)
                    F = copy(F_trial)
                    accepted_alpha = true
                    accepted_substeps += 1

                    # Slowly recover the step after a clean update.
                    delta_t_alpha =
                        min(1.05 * local_delta_t, delta_t_alpha_max)
                    break
                end

                # Backtrack only for a real intersection or non-finite values.
                local_delta_t *= 0.5

                @printf(
                    "       F backtrack | substep=%d | attempt=%d | finite=%s | crossU=%s | crossL=%s | new dt_alpha=%.3e\n",
                    substep,
                    attempt,
                    string(finite_trial),
                    string(crosses_upper),
                    string(crosses_lower),
                    local_delta_t
                )

                if local_delta_t < delta_t_alpha_min
                    break
                end
            end

            if !accepted_alpha
                delta_t_alpha =
                    max(local_delta_t, delta_t_alpha_min)
                println("       F held for the remaining substeps.")
                break
            end
        end

        # The spatial eigenvalues depend on L and omega, not on F itself.
        # Therefore, do not reclassify/retrack them merely because F moved.
        load_on_workers()
        @everywhere begin
            normals_F = contour_normals(F)
        end

        old_omegaF_max_i = omegaF_max_i
        omega_F = contour_omega_F(F)
        load_on_workers()

        if omega_update_accepted
            delta_t_omega =
                min(1.1 * delta_t_omega, delta_t_omega_max)
        end

        new_omegaF_max_i = maximum(imag.(omega_F))

        d_branch = branch_distance(alpha_L_u, alpha_L_l)
        d_contour = contour_distance(F, alpha_L_u, alpha_L_l)

        dist_u = minimum([
            minimum(abs.(f .- alpha_L_u))
            for f in F
        ])

        dist_l = minimum([
            minimum(abs.(f .- alpha_L_l))
            for f in F
        ])

        @printf(
            "[k=%04d] dUL=%9.3e | dFU=%9.3e | dFL=%9.3e | dF=%9.3e\n",
            k,
            d_branch,
            dist_u,
            dist_l,
            d_contour
        )

        @printf(
            "       F max substep move=%.3e | accepted substeps=%d/%d | dt_alpha=%.3e | dt_omega=%.3e | omega_F max imag: %.6e -> %.6e\n",
            total_alpha_move,
            accepted_substeps,
            n_alpha_substeps,
            delta_t_alpha,
            delta_t_omega,
            old_omegaF_max_i,
            new_omegaF_max_i
        )

        flush(stdout)

        ################
        # 3. Save JSON #
        ################

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

        open(filename, "w") do file
            write(file, JSON.json(current_array))
        end

        iteration_step += 1
    end
    #plot_omega()
    #plot_alpha()
    plot_omega_potential()
    plot_alpha_potential()

    begin 
        plot(L)
        plot!(omega_F)
    end