# briggsv4.5_ribbon.jl -- same DESCENT as v4.1-v4.4 (byte-identical to briggsv4.jl
# apart from file I/O and the iteration cap). NO theory changed. Nothing copied.
#
# What v4.4 got wrong (diagnosed from contour_iteration_v4.4_ribbon.json):
#   the tracking itself was clean -- max step 0.017 against a 0.03 guard, zero
#   guard trips, zero refinements, 100% of points on the upper side of F at every
#   frame -- but the branch it tracked DRIFTED AWAY FROM F:
#       median |alpha - F|   it=11: 0.035   it=61: 0.786   it=301: 2.079   (red)
#                            it=11: 0.035   it=61: 0.114   it=301: 0.612   (descent)
#   Briggs needs the branches ADJACENT to F -- the ones that collide at the pinch.
#   A branch 2 units away is a valid eigenvalue branch and completely irrelevant.
#
#   Cause 1, the seed. v4.4 seeded at the contour's left END, omega_r = 0. At
#   purely imaginary omega the spectrum is symmetric under alpha -> -conj(alpha),
#   so the seed is choosing between MIRROR TWINS (verified to 7 digits: red
#   +0.7836-1.9284i vs descent -0.7836-1.9284i at it=301). The tiebreak is
#   "closer to F", and the margin was 2.6% at it=21 -- a coin flip. It flipped at
#   iteration 21 and every later frame faithfully followed the wrong twin.
#   Cause 2, the corrector. It took the nearest eigenvalue to the predictor with
#   no side-of-F constraint, so nothing kept it adjacent to F once it had left.
#
# v4.5 = fixes A + B, both of which apply the EXISTING descent rule; no value is
# ever copied from the descent and every alpha is computed from this script's own
# freshly solved spectra:
#   A. Seed at the interior node the descent seeds at (start_index = N/4,
#      omega_r ~ 0.121) instead of the contour end, then march outward in both
#      directions -- exactly the protocol of contour_alpha_L_conti(). At that node
#      the branch sits next to F and there is no mirror twin to confuse it with.
#   B. The corrector now filters candidates by side of F before taking the nearest
#      to the predictor, i.e. the rule of couetteflow_spatial_sing_mode_comparison().
#      As there, if no eigenvalue lies on the required side the step falls back to
#      the nearest overall -- and every such fallback is COUNTED and reported, so
#      the filter can never silently swallow a branch that genuinely leaves F.
#
# Diagnostics per frame (stored in the JSON and printed):
#   xcheck_u/l    tracked ribbon branch vs the descent's tracked branch (red vs
#                 grey), at identical omega. NOT forced to zero -- if the two
#                 disagree that is a result, not a bug.
#   sideclass_u/l tracked branch vs the side-of-F classification evaluated
#                 independently at the same point, from the same spectrum.
#   nfallb_u/l    number of steps where the side filter found no candidate and
#                 fell back to nearest-overall (see B). Should stay at or near 0.
#   loop_mono_u/l MONODROMY of the closed detour: up the left riser, around the
#                 arc and down the right riser, vs straight along the bottom.
#                 Nonzero => the loop encloses a branch point of alpha(omega).
# No analytical pinch point exists for this case; none is drawn.
# Run from vib_ribbon/src/:  julia briggsv4.5_ribbon.jl
# Kilian Vinzenz Wilhelm
begin
    using Distributed, Plots, BenchmarkTools, FFTW, JSON, Statistics, Printf
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
        # ---- vibrating-ribbon marker (the ONLY new input) ----
        # omega_f = real forcing frequency the drawn L-contour detours around.
        # Must lie strictly inside (omega_r_start, omega_r_end). Cosmetic for the
        # descent; it positions the arc in the saved animation frames.
        omega_f = 0.25
        r_arc   = 0.01
    end
    function contour_L()
        # DESCENT contour: straight horizontal line -- identical to briggsv4.
        L = Complex{Float64}[omega_r[j] + omega_i * im for j in 1:N]
        return L 
    end
    # ---- ribbon contour geometry (v4.3) ----------------------------------
    # The horizontal parts sit ON the descent grid omega_r, so they are exactly the
    # omega points the descent already solved and their alpha values can be taken
    # straight from the descent's tracked branches (no re-classification).
    n_ver = 160          # points per vertical riser
    n_arc = 60           # points on the semicircle over omega_f
    n_bot = 21           # points on the short bottom segment (monodromy test only)
    n_sub = 2            # horizontal refinement: subdivisions per descent interval
    idx_hor_left  = findall(w -> w <= omega_f - r_arc, collect(omega_r))
    idx_hor_right = findall(w -> w >= omega_f + r_arc, collect(omega_r))
    @assert !isempty(idx_hor_left) && !isempty(idx_hor_right) "omega_f +/- r_arc outside the omega_r range"

    # Refine the horizontals while KEEPING the descent nodes as a subset, so the
    # red-vs-grey cross-check is at identical omega with no interpolation.
    function subdivide(xs, m)
        m <= 1 && return collect(float.(xs))
        out = Float64[]
        for i in 1:(length(xs) - 1), k in 0:(m - 1)
            push!(out, xs[i] + (xs[i+1] - xs[i]) * k / m)
        end
        push!(out, float(xs[end]))
        return out
    end
    wl_ref = subdivide(collect(omega_r)[idx_hor_left],  n_sub)
    wr_ref = subdivide(collect(omega_r)[idx_hor_right], n_sub)
    n_hl = length(wl_ref); n_hr = length(wr_ref)
    # positions of the descent nodes inside the full ribbon contour
    map_hl = [(j - 1) * n_sub + 1 for j in 1:length(idx_hor_left)]
    off_hr = n_hl + n_ver + n_arc + n_ver
    map_hr = [off_hr + (j - 1) * n_sub + 1 for j in 1:length(idx_hor_right)]
    i_riser_bot_L = n_hl + 1        # bottom of the left riser,  omega_f - r_arc + i*omega_i
    i_riser_bot_R = off_hr          # bottom of the right riser, omega_f + r_arc + i*omega_i

    riser_left_nodes(omega_i_level)  = ComplexF64.((omega_f - r_arc) .+ im .* collect(LinRange(omega_i_level, 0.0, n_ver)))
    riser_right_nodes(omega_i_level) = ComplexF64.((omega_f + r_arc) .+ im .* collect(LinRange(omega_i_level, 0.0, n_ver)))
    arc_nodes() = ComplexF64.(omega_f .+ r_arc .* cis.(collect(LinRange(pi, 0.0, n_arc))))

    function hankel_L(omega_i_level)
        # NEAT Eq.(24) ribbon contour at horizontal level omega_i_level:
        # horizontals (descent grid, subdivided) + risers + semicircle over omega_f.
        seg1 = ComplexF64.(wl_ref .+ im * omega_i_level)
        seg5 = ComplexF64.(wr_ref .+ im * omega_i_level)
        return vcat(seg1, riser_left_nodes(omega_i_level), arc_nodes(),
                    reverse(riser_right_nodes(omega_i_level)), seg5)
    end
    # short straight segment closing the detour along the bottom, used only for
    # the monodromy test (it lies on the descent line L).
    bottom_nodes(omega_i_level) =
        ComplexF64.(collect(LinRange(omega_f - r_arc, omega_f + r_arc, n_bot)) .+ im * omega_i_level)
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
    filename = "contour_iteration_v4.5_ribbon.json"
end

# ---- ribbon animation frame (v4.3): descent-anchored, adaptively refined ----
@everywhere function spatial_spectrum(omega)
    eigvals, _ = couetteflow(nothing, omega, Val(:alpha_collocation), Val(:eigen))
    return ComplexF64[ev for ev in eigvals if isfinite(real(ev)) && isfinite(imag(ev)) && abs(ev) < 1.0e3]
end

isfin(z) = isfinite(real(z)) && isfinite(imag(z))

# Signed distance of each eigenvalue from the F contour: >0 upper side, <0 lower.
# Computed once per node and shared by both branches.
function side_projections(evs, F, normals_F)
    out = Vector{Float64}(undef, length(evs))
    for (i, ev) in enumerate(evs)
        k = argmin(abs.(ev .- F))
        out[i] = real(conj(normals_F[k]) * (ev - F[k]))
    end
    return out
end

# Continuation over a PRECOMPUTED set of spectra.
# (B) Candidates are first restricted to the required SIDE OF F -- the rule of
# couetteflow_spatial_sing_mode_comparison() -- and only then is the nearest to
# the predictor taken. If a side is empty the step falls back to nearest-overall,
# exactly as the descent does, and the fallback is counted so it cannot pass
# unnoticed. Secant predictor with the displacement clamped to da_max; a node is
# flagged when the accepted step is too large or the winner is not clearly nearer
# than the runner-up.
function continue_branch(nodes, spectra, sides, seed, branch_side;
                         da_max = 0.03, sep_frac = 0.6, side_tol = 0.0)
    n = length(nodes)
    vals = Vector{ComplexF64}(undef, n)
    flag = falses(n)
    n_fallback = 0
    a_prev  = ComplexF64(seed); w_prev  = nodes[1]
    a_prev2 = ComplexF64(seed); w_prev2 = nodes[1]
    have2 = false
    for j in 1:n
        w = nodes[j]
        pred = a_prev
        if have2
            dw  = w - w_prev
            dwp = w_prev - w_prev2
            if abs(dwp) > 1e-14
                disp = (dw / dwp) * (a_prev - a_prev2)
                abs(disp) > da_max && (disp *= da_max / abs(disp))
                pred = a_prev + disp
            end
        end
        evs = spectra[j]; sd = sides[j]
        if isempty(evs) || !isfin(pred)
            vals[j] = a_prev
            flag[j] = true
        else
            keep = Int[]
            for k in eachindex(evs)
                if (branch_side === :upper && sd[k] >  side_tol) ||
                   (branch_side === :lower && sd[k] < -side_tol)
                    push!(keep, k)
                end
            end
            if isempty(keep)                      # side empty -> descent's fallback
                keep = collect(eachindex(evs))
                n_fallback += 1
                flag[j] = true
            end
            dists = [abs(evs[k] - pred) for k in keep]
            k1 = argmin(dists)
            d1 = dists[k1]
            d2 = Inf
            for (kk, dd) in enumerate(dists)
                kk != k1 && dd < d2 && (d2 = dd)
            end
            cand = evs[keep[k1]]
            vals[j] = cand
            (abs(cand - a_prev) > da_max || d1 > sep_frac * d2) && (flag[j] = true)
        end
        a_prev2 = a_prev; w_prev2 = w_prev
        a_prev  = vals[j]; w_prev = w
        have2 = true
    end
    return vals, flag, n_fallback
end

# Track one path, bisecting every rejected interval and re-solving the inserted
# nodes in a single parallel batch. Values are returned at the ORIGINAL nodes only.
function track_path(nodes0, spectra0, sides0, seed, branch_side;
                    max_refine = 4, node_cap = 2000, da_max = 0.03, sep_frac = 0.6)
    nodes   = copy(nodes0)
    spectra = Vector{Vector{ComplexF64}}(spectra0)
    sides   = Vector{Vector{Float64}}(sides0)
    is_orig = trues(length(nodes))
    n_added = 0
    vals, flag, nfb = continue_branch(nodes, spectra, sides, seed, branch_side;
                                      da_max = da_max, sep_frac = sep_frac)
    for _ in 1:max_refine
        ins_at = Int[]; ins_w = ComplexF64[]
        for j in 2:length(nodes)
            if flag[j]
                push!(ins_at, j)
                push!(ins_w, 0.5 * (nodes[j-1] + nodes[j]))
            end
        end
        (isempty(ins_w) || length(nodes) + length(ins_w) > node_cap) && break
        new_spec  = pmap(spatial_spectrum, ins_w)
        new_sides = [side_projections(sp, F, normals_F) for sp in new_spec]
        n_added += length(ins_w)
        nn = ComplexF64[]; ss = Vector{Vector{ComplexF64}}()
        vv = Vector{Vector{Float64}}(); oo = Bool[]
        p = 1
        for j in 1:length(nodes)
            while p <= length(ins_at) && ins_at[p] == j
                push!(nn, ins_w[p]); push!(ss, new_spec[p]); push!(vv, new_sides[p]); push!(oo, false); p += 1
            end
            push!(nn, nodes[j]); push!(ss, spectra[j]); push!(vv, sides[j]); push!(oo, is_orig[j])
        end
        nodes = nn; spectra = ss; sides = vv; is_orig = oo
        vals, flag, nfb = continue_branch(nodes, spectra, sides, seed, branch_side;
                                          da_max = da_max, sep_frac = sep_frac)
    end
    return vals[is_orig], n_added, nfb
end

# (A) March OUTWARD from an interior seed node, both directions -- the protocol of
# contour_alpha_L_conti(), which starts at N/4 rather than at a contour end.
function track_outward(nodes, spectra, sides, i_seed, seed, branch_side)
    vf, af, nf = track_path(nodes[i_seed:end],    spectra[i_seed:end],    sides[i_seed:end],    seed, branch_side)
    vb, ab, nb = track_path(nodes[i_seed:-1:1],   spectra[i_seed:-1:1],   sides[i_seed:-1:1],   seed, branch_side)
    vals = vcat(reverse(vb[2:end]), vf)          # vb[1] and vf[1] are the same node
    return vals, af + ab, nf + nb
end

# One-shot continuation for the monodromy leg (short, lies on the descent line).
function track_simple(nodes, seed, branch_side; da_max = 0.03, sep_frac = 0.6)
    spectra = pmap(spatial_spectrum, nodes)
    sides   = [side_projections(sp, F, normals_F) for sp in spectra]
    vals, _, _ = continue_branch(nodes, spectra, sides, seed, branch_side;
                                 da_max = da_max, sep_frac = sep_frac)
    return vals
end

# The descent's OWN branch rule -- "of the eigenvalues on each side of F, take the
# one closest to F" -- applied to an already-computed spectrum. Identical logic to
# dominant_eigvals(), but it does not re-solve the eigenvalue problem.
function classify_by_F(evs, sd, F)
    best_u = nothing; best_l = nothing
    du = Inf; dl = Inf
    for (k, ev) in enumerate(evs)
        d = minimum(abs.(ev .- F))
        if sd[k] > 0.0 && d < du
            du = d; best_u = ev
        elseif sd[k] < 0.0 && d < dl
            dl = d; best_l = ev
        end
    end
    return best_u, best_l
end

function descent_branch_clean(branch)
    return ComplexF64[x === nothing ? ComplexF64(NaN, NaN) : ComplexF64(x) for x in branch]
end

function seed_or_nearest(cand, evs, F)
    (cand !== nothing && isfin(ComplexF64(cand))) && return ComplexF64(cand)
    isempty(evs) && return ComplexF64(NaN, NaN)
    return ComplexF64(evs[argmin([minimum(abs.(ev .- F)) for ev in evs])])
end

maxabsdiff(a, b) = isempty(a) ? 0.0 : maximum(abs.(ComplexF64.(a) .- ComplexF64.(b)))

function make_ribbon_frame(iter)
    Lh = hankel_L(omega_i)
    spectra = pmap(spatial_spectrum, Lh)
    sides   = [side_projections(sp, F, normals_F) for sp in spectra]

    # (A) seed at the node the DESCENT seeds at: start_index = N/4 on the omega_r
    # grid, which the ribbon horizontals contain exactly. The rule is recomputed
    # here from this frame's own spectrum -- the descent's value is not used.
    i_seed = map_hl[max(1, min(length(map_hl), floor(Int, N / 4)))]
    su, sl = classify_by_F(spectra[i_seed], sides[i_seed], F)
    seed_u = seed_or_nearest(su, spectra[i_seed], F)
    seed_l = seed_or_nearest(sl, spectra[i_seed], F)

    au, add_u, nfb_u = track_outward(Lh, spectra, sides, i_seed, seed_u, :upper)
    al, add_l, nfb_l = track_outward(Lh, spectra, sides, i_seed, seed_l, :lower)
    @assert length(au) == length(Lh) "ribbon branch/contour length mismatch"

    ud = descent_branch_clean(alpha_L_u)
    ld = descent_branch_clean(alpha_L_l)

    # cross-check: red vs grey at identical omega, NOT forced to agree
    hidx = vcat(map_hl, map_hr); didx = vcat(idx_hor_left, idx_hor_right)
    xcheck_u = maxabsdiff(au[hidx], ud[didx])
    xcheck_l = maxabsdiff(al[hidx], ld[didx])

    # continuation vs side-of-F classification, from the same spectrum
    scu = 0.0; scl = 0.0
    for i in hidx
        cu, cl = classify_by_F(spectra[i], sides[i], F)
        cu !== nothing && (scu = max(scu, abs(au[i] - ComplexF64(cu))))
        cl !== nothing && (scl = max(scl, abs(al[i] - ComplexF64(cl))))
    end

    # monodromy of the closed detour
    bn = bottom_nodes(omega_i)
    loop_mono_u = abs(track_simple(bn, au[i_riser_bot_L], :upper)[end] - au[i_riser_bot_R])
    loop_mono_l = abs(track_simple(bn, al[i_riser_bot_L], :lower)[end] - al[i_riser_bot_R])

    @printf("        ribbon: xcheck=%.3e/%.3e side=%.3e/%.3e mono=%.3e/%.3e | refined u=%d l=%d | side-fallbacks u=%d l=%d\n",
            xcheck_u, xcheck_l, scu, scl, loop_mono_u, loop_mono_l, add_u, add_l, nfb_u, nfb_l)
    flush(stdout)

    return Dict(
        "iteration" => iter,
        "omega_f"   => omega_f,
        "L"         => complexvec_to_json(Lh),
        "alpha_L_u" => complexvec_to_json(au),
        "alpha_L_l" => complexvec_to_json(al),
        "F"         => complexvec_to_json(F),
        "omega_F"   => complexvec_to_json(omega_F),
        "L_descent"         => complexvec_to_json(L),
        "alpha_L_u_descent" => complexvec_to_json(ud),
        "alpha_L_l_descent" => complexvec_to_json(ld),
        "xcheck_u"    => xcheck_u,    "xcheck_l"    => xcheck_l,
        "sideclass_u" => scu,         "sideclass_l" => scl,
        "nfallb_u"    => nfb_u,       "nfallb_l"    => nfb_l,
        "loop_mono_u" => loop_mono_u, "loop_mono_l" => loop_mono_l,
    )
end

# ---- streaming frame writer -------------------------------------------------
# Up to v4.3 the ENTIRE frame array was re-serialised and rewritten on every
# iteration. That is O(n^2) I/O (~4 GB of writes over 300 frames) and builds a
# >25 MB String each time; if the filesystem ever returns a short write, the file
# is left truncated mid-object and every later read of it fails. Instead the file
# is opened once and each frame is appended in place. The closing "]" is written
# after every frame and then seeked back over, so the file on disk is a COMPLETE,
# VALID JSON array at all times -- if the job is killed or hits a quota, whatever
# was written up to that point still loads.
mutable struct FrameWriter
    io::IOStream
    n::Int
end
function open_frame_writer(fname)
    io = open(fname, "w+")
    write(io, "[")
    p = position(io)
    write(io, "]")
    flush(io)
    seek(io, p)
    return FrameWriter(io, 0)
end
function write_frame!(fw::FrameWriter, frame)
    fw.n > 0 && write(fw.io, ",")
    write(fw.io, JSON.json(frame))
    p = position(fw.io)
    write(fw.io, "]")          # keep the file valid on disk...
    flush(fw.io)
    seek(fw.io, p)             # ...but overwrite it with the next frame
    fw.n += 1
    return fw.n
end
function close_frame_writer!(fw::FrameWriter)
    write(fw.io, "]")
    flush(fw.io)
    close(fw.io)
    return fw.n
end

begin
    iteration_step = 1
    frame_writer = open_frame_writer(filename)
    last_frame = make_ribbon_frame(iteration_step)
    write_frame!(frame_writer, last_frame)
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
for k = 1:300
    #print_iteration_header(k)
    global omega_i, L, alpha_L_u, alpha_L_l, alpha_i, F, omega_F, iteration_step, frame_writer, last_frame

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
                greater_candidates = filter(x -> isfinite(x) && x > omega_lower_bound, omega_i_cache)

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
        @printf(
            "[%04d] jump=%9.3e | dUL=%9.3e | dUF=%9.3e | dLF=%9.3e | dt=%9.3e | ζ=%9.3e | %s/%s\n",
            k,
            omega_jump,
            d_branch,
            dist_u,
            dist_l,
            local_delta_t,
            zeta_common,
            omega_status,
            alpha_status
        )
        flush(stdout)
        load_on_workers()
        global last_frame = make_ribbon_frame(iteration_step)
        write_frame!(frame_writer, last_frame)
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

n_frames_written = close_frame_writer!(frame_writer)
@printf("\nRun finished: %d frames written to %s\n", n_frames_written, filename)
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
# NOTE: load_step() above is kept for interactive use, but it is NOT called here
# any more. Re-parsing the full (tens of MB) results file at the end of a batch
# run is pointless -- the final frame is already in memory -- and it is what
# crashed the v4.3 job when the file on disk had been truncated.
# It also restored the WRONG arrays in the ribbon scripts: entry["L"] /
# entry["alpha_L_u"] are the 476-point RIBBON contour, not the 100-point descent
# state. The *_descent fields are the ones that define the descent.
iteration_step = last_frame["iteration"]
L         = json_to_complexvec(last_frame["L_descent"])
alpha_L_u = json_to_complexvec(last_frame["alpha_L_u_descent"])
alpha_L_l = json_to_complexvec(last_frame["alpha_L_l_descent"])
F         = json_to_complexvec(last_frame["F"])
omega_F   = json_to_complexvec(last_frame["omega_F"])
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