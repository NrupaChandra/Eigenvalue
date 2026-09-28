# Kilian Vinzenz Wilhelm
2+2
begin
    using Distributed, Plots, BenchmarkTools, FFTW, JSON, Statistics, Printf
    addprocs(7)
    w = workers()
end
@everywhere begin
    C1 = -1.0 + 0.96im
    C2 = 0.4 + 0.8im
    C3 = 1.0 + 0.0im

    alpha_pinch_exact = -C2 / (2*C3)
    omega_pinch_exact = C1 - C2^2 / (4*C3)

    println("Exact pinch alpha = ", alpha_pinch_exact)
    println("Exact pinch omega = ", omega_pinch_exact)
end

@everywhere function algebraic_omega(alpha)
    return C1 + C2*alpha + C3*alpha^2
end
@everywhere function find_valid_start_index(L, F, normals_F)
    N = length(L)
    preferred = floor(Int, N / 2)

    candidate_indices = sort(collect(1:N), by = i -> abs(i - preferred))

    for idx in candidate_indices
        alpha_u_start, alpha_l_start, _, _ = dominant_eigvals(L[idx], F, normals_F)

        if alpha_u_start !== nothing && alpha_l_start !== nothing
            return idx, alpha_u_start, alpha_l_start
        end
    end

    error("No valid start index found: could not classify upper and lower roots.")
end

@everywhere function algebraic_alpha_roots(omega)
    disc = sqrt(C2^2 - 4*C3*(C1 - omega))
    alpha1 = (-C2 + disc) / (2*C3)
    alpha2 = (-C2 - disc) / (2*C3)
    return ComplexF64[alpha1, alpha2]
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
    # ALPHA contour: F starts on the real alpha line
    @everywhere begin
        alpha_r_start = -1.0
        alpha_r_end = 0.5
        N = 141

        alpha_r = range(alpha_r_start, alpha_r_end, length=N)
        #alpha_i = fill(imag(alpha_pinch_exact), N)
        alpha_i = fill(0.0, N)
    end
    function contour_F()
        F = [alpha_r[j] + alpha_i[j] * im for j in 1:N]
        return F
    end
    @everywhere F = Vector{Complex{Float64}}[]
    F = contour_F()
    @everywhere function couetteflow_temporal_sing_mode(alpha)
        omega = algebraic_omega(alpha)
        return omega, nothing
    end
    @everywhere function couetteflow_spatial_sing_mode_comparison(
        omega,
        alpha_approximation,
        branch_side::Symbol,
        F,
        normals_F;
        side_tol = 0.0
    )
        eigvals = algebraic_alpha_roots(omega)

        candidates = ComplexF64[]

        for eigval in eigvals
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
            diffs = abs.(ComplexF64(alpha_approximation) .- eigvals)
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
        omega_r_start = -1.2
        omega_r_end = -0.5
        omega_r = range(omega_r_start, omega_r_end, length=N)
    end

    omega_clearance_initial = 0.1
    omega_i = maximum(imag.(omega_F)) + omega_clearance_initial
    @everywhere omega_i = $omega_i
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
        eigvals = algebraic_alpha_roots(omega)

        eigval_dominant_u = nothing
        eigval_dominant_l = nothing

        min_dist_u = Inf
        min_dist_l = Inf

        for eigval in eigvals
            distances = [abs(eigval - f) for f in F]
            idx_min = argmin(distances)

            f_near = F[idx_min]
            normal = normals_F[idx_min]

            proj = real(conj(normal) * (eigval - f_near))
            dist_to_line = distances[idx_min]

            if proj > 0.0 && dist_to_line < min_dist_u
                min_dist_u = dist_to_line
                eigval_dominant_u = eigval
            elseif proj < 0.0 && dist_to_line < min_dist_l
                min_dist_l = dist_to_line
                eigval_dominant_l = eigval
            end
        end

        return eigval_dominant_u, eigval_dominant_l, nothing, nothing
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
        start_index, alpha_u_start, alpha_l_start =
            find_valid_start_index(L, F, normals_F)
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

        omega_clearance_initial = 0.1
        global omega_i = maximum(imag.(omega_F)) + omega_clearance_initial

        global omega_r = range(omega_r_start, omega_r_end, length=N)
        global L = contour_L()

        load_on_workers()

        @everywhere begin
            global normals_F = contour_normals(F)
        end

        global alpha_L_u, alpha_L_l = contour_alpha_L_init(L)
        load_on_workers()

        global iteration_step = 1
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
    filename = "contour_iteration_algebraic_testv6.json"
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
function isfinite_complex(z)
    return z isa Complex && isfinite(real(z)) && isfinite(imag(z))
end

function finite_complex_vec(v)
    return ComplexF64[z for z in v if isfinite_complex(z)]
end

function min_distance_between(A, B)
    A_f = finite_complex_vec(A)
    B_f = finite_complex_vec(B)

    if isempty(A_f) || isempty(B_f)
        return Inf
    end

    dmin = Inf
    for a in A_f
        for b in B_f
            dmin = min(dmin, abs(a - b))
        end
    end

    return dmin
end
function branch_distance(alpha_L_u, alpha_L_l)
    gaps = Float64[]

    for i in eachindex(alpha_L_u, alpha_L_l)
        u = alpha_L_u[i]
        l = alpha_L_l[i]

        if isfinite_complex(u) && isfinite_complex(l)
            push!(gaps, abs(u - l))
        end
    end

    if isempty(gaps)
        return Inf
    end

    return minimum(gaps)
end

function contour_L_at(omega_i_trial)
    return ComplexF64[omega_r[j] + omega_i_trial * im for j in 1:N]
end

function branch_overlap_valid(alpha_u, alpha_l; overlap_tol = 1e-8)

    alpha_u_f = finite_complex_vec(alpha_u)
    alpha_l_f = finite_complex_vec(alpha_l)

    if isempty(alpha_u_f) || isempty(alpha_l_f)
        return false, "upper/lower branch contains no finite values"
    end

    min_dist = Inf
    min_i = 0
    min_j = 0

    for i in eachindex(alpha_u_f)
        for j in eachindex(alpha_l_f)
            d = abs(alpha_u_f[i] - alpha_l_f[j])

            if isfinite(d) && d < min_dist
                min_dist = d
                min_i = i
                min_j = j
            end
        end
    end

    if min_dist < overlap_tol
        return false, @sprintf(
            "upper/lower overlap: d=%.3e at finite upper[%d], finite lower[%d]",
            min_dist, min_i, min_j
        )
    end

    return true, @sprintf(
        "no overlap: min upper/lower distance=%.3e",
        min_dist
    )
end
function local_contour_distances(F, alpha_L_u, alpha_L_l)
    all_branches = finite_complex_vec(vcat(alpha_L_u, alpha_L_l))

    if isempty(all_branches)
        return fill(Inf, length(F))
    end

    distances = Float64[]

    for f in F
        if isfinite_complex(f)
            push!(distances, minimum(abs.(f .- all_branches)))
        else
            push!(distances, Inf)
        end
    end

    return distances
end

function contour_distance(F, alpha_L_u, alpha_L_l)
    ds = local_contour_distances(F, alpha_L_u, alpha_L_l)
    ds_finite = filter(isfinite, ds)

    if isempty(ds_finite)
        return Inf
    end

    return minimum(ds_finite)
end

function nearest_branch_info(f, alpha_L_u, alpha_L_l)
    all_branches = finite_complex_vec(vcat(alpha_L_u, alpha_L_l))

    if isempty(all_branches) || !isfinite_complex(f)
        return 0.0 + 0.0im, Inf
    end

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
    d_safe = 0.05      # start slowing below this
    min_factor = 0.2  # never go below  10% speed

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
function estimate_pinch_point(alpha_L_u, alpha_L_l, L)

    best_gap = Inf
    best_idx = 0

    for i in eachindex(alpha_L_u, alpha_L_l, L)
        u = alpha_L_u[i]
        l = alpha_L_l[i]
        omega = L[i]

        if isfinite_complex(u) && isfinite_complex(l) && isfinite_complex(omega)
            gap = abs(u - l)

            if gap < best_gap
                best_gap = gap
                best_idx = i
            end
        end
    end

    if best_idx == 0
        return NaN + NaN*im, NaN + NaN*im, Inf, 0
    end

    alpha_pinch_num = 0.5 * (alpha_L_u[best_idx] + alpha_L_l[best_idx])
    omega_pinch_num = L[best_idx]

    return alpha_pinch_num, omega_pinch_num, best_gap, best_idx
end

blocked_count = 0

for k = 1:1500
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

        min_valid_omega_candidates = 1
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

                # Do not allow unchanged endpoints to determine the new horizontal L level
                omega_i_candidates = omega_i_cache[2:end-1]

                omega_i_candidates = [
                    isfinite(x) && abs(x) < 10.0 ? x : -Inf
                    for x in omega_i_candidates
                ]

                greater_candidates = filter(
                    x -> isfinite(x) && x > omega_lower_bound && abs(x - omega_i) > 1e-14,
                    omega_i_candidates
                )

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

        omega_blocked = false

        if !omega_accepted
            println("No safe omega accepted; keeping L fixed and updating F only.")

            omega_status = "blocked"
            omega_blocked = true
            omega_jump = 0.0
            omega_dt_used = 0.0

            # Keep the old omega contour and existing branches
            global omega_i = omega_i_old
            global L = contour_L()

            # Do NOT recompute alpha_L_u / alpha_L_l here.
            # They already correspond to the fixed L from the previous accepted step.
            load_on_workers()

            @everywhere begin
                normals_F = contour_normals(F)
            end
        end
        if omega_repaired
            global L = contour_L()
            load_on_workers()

            @everywhere begin
                normals_F = contour_normals(F)
            end

            global alpha_L_u, alpha_L_l = contour_alpha_L_conti(L)
            load_on_workers()

        elseif omega_accepted && !omega_blocked
            load_on_workers()

            @everywhere begin
                normals_F = contour_normals(F)
            end
        end

        if omega_status == "blocked"
            global blocked_count += 1
        else
            global blocked_count = 0
        end

        if blocked_count > 30
            println("STOP: omega remained blocked for 30 consecutive iterations.")
            break
        end

        alpha_i_cache = copy(alpha_i)

        #print_block("BRANCH DISTANCES")
        n_bad_F = count(z -> !isfinite_complex(z), F)
        n_bad_u = count(z -> !isfinite_complex(z), alpha_L_u)
        n_bad_l = count(z -> !isfinite_complex(z), alpha_L_l)
        n_bad_omegaF = count(z -> !isfinite_complex(z), omega_F)
        if n_bad_F + n_bad_u + n_bad_l + n_bad_omegaF > 0
            @printf(
                "[k=%d] WARNING: non-finite values | F=%d | upper=%d | lower=%d | omega_F=%d\n",
                k, n_bad_F, n_bad_u, n_bad_l, n_bad_omegaF
            )
        end

        d_branch = branch_distance(alpha_L_u, alpha_L_l)
        d_contour = contour_distance(F, alpha_L_u, alpha_L_l)
        dist_u = min_distance_between(F, alpha_L_u)
        dist_l = min_distance_between(F, alpha_L_l)
        branch_factor = branch_slowdown_factor(d_branch)
        local_delta_t = delta_t * branch_factor

        # If L is blocked, only let F creep slowly.
        # This prevents the F-only step from blowing up.
        if omega_status == "blocked"
            local_delta_t = min(local_delta_t, 1e-5)
        end

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
        pinch_tol = 1e-2

        #println("d_branch = ", d_branch, ", d_contour = ", d_contour)
        stop_after_save = false
        stop_reason = ""

    if !isfinite(d_branch)
        @printf(
            "[k=%d] STOP: non-finite branch distance | d_branch=%.3e\n",
            k, d_branch
        )
        break
    end

    if !isfinite(d_contour)
        @printf(
            "[k=%d] WARNING: non-finite contour distance | d_contour=%.3e; using conservative move_safety\n",
            k, d_contour
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

    F_old = copy(F)

    for attempt in 1:8
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
            rhs_j = clamp(rhs_j, -rhs_cap, rhs_cap)

            alpha_i_trial[j] = alpha_i[j] + local_delta_t * rhs_j
        end

        alpha_i_trial[1] = alpha_i_trial[2]
        alpha_i_trial[end] = alpha_i_trial[end-1]

        # Smooth the actual contour that will be tested and applied
        alpha_i_smooth = rolling_average_filter(alpha_i_trial, 3)

        alpha_i_smooth[1] = alpha_i_smooth[2]
        alpha_i_smooth[end] = alpha_i_smooth[end-1]

        if any(!isfinite, alpha_i_smooth)
            local_delta_t *= 0.5
            alpha_status = "reduced-nonfinite"
            alpha_attempt_used = attempt
            continue
        end
        F_trial = ComplexF64[
            alpha_r[j] + alpha_i_smooth[j] * im for j in 1:N
        ]
        if any(z -> !isfinite(real(z)) || !isfinite(imag(z)), F_trial)
            local_delta_t *= 0.5
            alpha_status = "reduced-nonfinite-F"
            alpha_attempt_used = attempt
            continue
        end
        dt_ok, toward_amounts, local_allowed_toward, full_moves =
            directional_dt_check(
                F_old,
                F_trial,
                alpha_L_u,
                alpha_L_l;
                move_safety = move_safety,
                global_max_move = global_max_move
            )

        move_raw = abs.(alpha_i_trial .- alpha_i)
        move_smooth = abs.(alpha_i_smooth .- alpha_i)

        max_raw_move = maximum(move_raw)
        max_smooth_move = maximum(move_smooth)

        if dt_ok
            alpha_i_cache = copy(alpha_i_smooth)
            accepted = true
            alpha_status = "accepted"
            alpha_attempt_used = attempt
            break
        else
            local_delta_t *= 0.5
            alpha_status = "reduced"
            alpha_attempt_used = attempt
        end

        if local_delta_t < dt_min
            alpha_status = "rejected"
            break
        end
    end

    if accepted
        global alpha_i = copy(alpha_i_cache)
    else
        # Important: do not apply a rejected alpha update
        alpha_status = "rejected"
        global alpha_i = copy(alpha_i)
    end
        global F = contour_F()
        load_on_workers()
        global omega_F = contour_omega_F(F)

        if any(z -> !isfinite(real(z)) || !isfinite(imag(z)), F) || any(z -> !isfinite(real(z)) || !isfinite(imag(z)), omega_F)

            @printf("[k=%d] STOP: F or omega_F became non-finite after alpha update.\n", k)
            break
        end
        alpha_pinch_num, omega_pinch_num, branch_gap, pinch_idx = estimate_pinch_point(alpha_L_u, alpha_L_l, L)

        pinch_alpha_error = abs(alpha_pinch_num - alpha_pinch_exact)
        pinch_omega_error = abs(omega_pinch_num - omega_pinch_exact)
        @printf(
            "[%04d] jump=%9.3e | dUL=%9.3e | dUF=%9.3e | dLF=%9.3e | dt=%9.3e | err_alpha_pinch=%9.3e | err_omega_pinch=%9.3e | %s/%s\n",
            k,
            omega_jump,
            d_branch,
            dist_u,
            dist_l,
            local_delta_t,
            pinch_alpha_error,
            pinch_omega_error,
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
            alpha_pinch_num, omega_pinch_num, branch_gap, pinch_idx =
                estimate_pinch_point(alpha_L_u, alpha_L_l, L)

            println()
            println("====================================================================")
            println("                         PINCH POINT FOUND")
            println("====================================================================")

            @printf("Stopped at iteration:             %d\n", k)
            @printf("Stop reason:                      %s\n", stop_reason)
            @printf("Closest branch index:             %d\n", pinch_idx)

            @printf("Numerical alpha pinch:            %.12e %+.12ei\n",
                    real(alpha_pinch_num), imag(alpha_pinch_num))

            @printf("Corresponding omega value:        %.12e %+.12ei\n",
                    real(omega_pinch_num), imag(omega_pinch_num))

            @printf("Upper branch at pinch index:      %.12e %+.12ei\n",
                    real(alpha_L_u[pinch_idx]), imag(alpha_L_u[pinch_idx]))

            @printf("Lower branch at pinch index:      %.12e %+.12ei\n",
                    real(alpha_L_l[pinch_idx]), imag(alpha_L_l[pinch_idx]))

            @printf("Upper-lower branch gap:           %.12e\n", branch_gap)

            @printf("Exact alpha pinch:                %.12e %+.12ei\n",
                    real(alpha_pinch_exact), imag(alpha_pinch_exact))

            @printf("Exact omega pinch:                %.12e %+.12ei\n",
                    real(omega_pinch_exact), imag(omega_pinch_exact))

            @printf("Alpha pinch error:                %.12e\n",
                    abs(alpha_pinch_num - alpha_pinch_exact))

            @printf("Omega pinch error:                %.12e\n",
                    abs(omega_pinch_num - omega_pinch_exact))

            println("====================================================================")
            flush(stdout)

            break
        end
    end    
    #print_iteration_footer(k)
end

alpha_pinch_num, omega_pinch_num, branch_gap, pinch_idx =
    estimate_pinch_point(alpha_L_u, alpha_L_l, L)

println()
println("====================================================================")
println("                         FINAL PINCH ESTIMATE")
println("====================================================================")

@printf("Closest branch index:             %d\n", pinch_idx)

@printf("Numerical alpha pinch:            %.12e %+.12ei\n",
        real(alpha_pinch_num), imag(alpha_pinch_num))

@printf("Corresponding omega value:        %.12e %+.12ei\n",
        real(omega_pinch_num), imag(omega_pinch_num))

@printf("Upper-lower branch gap:           %.12e\n", branch_gap)

@printf("Exact alpha pinch:                %.12e %+.12ei\n",
        real(alpha_pinch_exact), imag(alpha_pinch_exact))

@printf("Exact omega pinch:                %.12e %+.12ei\n",
        real(omega_pinch_exact), imag(omega_pinch_exact))

@printf("Alpha pinch error:                %.12e\n",
        abs(alpha_pinch_num - alpha_pinch_exact))

@printf("Omega pinch error:                %.12e\n",
        abs(omega_pinch_num - omega_pinch_exact))

println("====================================================================")
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
#iteration_step, L, alpha_L_u, alpha_L_l, F, omega_F = load_step(filename; offset=0)
#println("Loaded iteration: ", iteration_step)
#global omega_i = imag(L[1])
#global alpha_i = imag.(F)
#load_on_workers()
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