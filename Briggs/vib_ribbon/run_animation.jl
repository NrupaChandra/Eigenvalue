# run_animation.jl
# ---------------------------------------------------------------------------
# Produces output in the SAME per-iteration JSON schema as the old Briggs runs
# (contour_iteration*.json) so the existing MATLAB animator
# (../plot_contour_deformation.m) can replay it.
#
# It runs the causal descent as a k-loop: the temporal contour L is a horizontal
# line that is lowered in Im(omega) step by step (k = 1 .. nsteps). At every step
# the two spatial branches alpha_L_u / alpha_L_l are computed over the WHOLE L
# line and one frame is appended to the JSON file, exactly like the old code.
#
# The spatial contour F stays on the real alpha axis (for Couette no spatial root
# approaches it, so F need not deform); omega_F is the dominant temporal
# eigenvalue over F and is therefore constant across frames. Only L and the
# branches move -- that is the causal descent, animated.
#
# Schema per frame (identical to briggsv4.jl / contour_test_v3.1.jl):
#   { "iteration", "L", "alpha_L_u", "alpha_L_l", "F", "omega_F" }
# each complex vector stored as [ {"re":..,"im":..}, ... ].
#
# Usage (from this directory):
#   julia run_animation.jl
# Then in MATLAB point plot_contour_deformation.m at the produced file.
# ---------------------------------------------------------------------------

include("ribbon_core.jl")     # brings in models.jl
using JSON, Printf

# ---- run parameters (edit here) -------------------------------------------
const RE          = 2000.0
const NUM_MODES   = 150            # fine resolution for the solve
const NUM_COARSE  = 100            # coarse resolution for the spurious filter
const OMEGA_F     = 0.5            # ribbon forcing frequency (read alpha+- here afterwards)

const ALPHA_R0    = 0.0            # F : real alpha axis, [ALPHA_R0, ALPHA_R1]
const ALPHA_R1    = 1.0
const N_F         = 100            # points on F

const OMEGA_R0    = 0.0            # L : horizontal line over [OMEGA_R0, OMEGA_R1]
const OMEGA_R1    = 1.0
const M_L         = 100            # points on L  (old runs used up to 296)

const SIGMA0      = 0.3            # start height Im(omega)
const NSTEPS      = 200            # number of descent frames  (old loop was 1..500)
const OUTFILE     = "contour_iteration_ribbon.json"
# ---------------------------------------------------------------------------

complexvec_to_json(v) = [Dict("re" => real(x), "im" => imag(x)) for x in v]

frame_dict(k, L, au, al, F, omega_F) = Dict(
    "iteration" => k,
    "L"          => complexvec_to_json(L),
    "alpha_L_u"  => complexvec_to_json(au),
    "alpha_L_l"  => complexvec_to_json(al),
    "F"          => complexvec_to_json(F),
    "omega_F"    => complexvec_to_json(omega_F),
)

# horizontal L line at height sigma
line_L(sigma) = ComplexF64[wr + im * sigma for wr in range(OMEGA_R0, OMEGA_R1, length = M_L)]

# nearest-to-axis upper / lower spatial root at one omega (F = real axis)
function leaders_at(model, omega, coarse)
    roots, _ = physical_spatial_roots(model, omega; model_coarse = coarse)
    up = filter(a -> imag(a) > 0, roots)
    lo = filter(a -> imag(a) < 0, roots)
    au = isempty(up) ? ComplexF64(NaN, NaN) : up[argmin(imag.(up))]
    al = isempty(lo) ? ComplexF64(NaN, NaN) : lo[argmax(imag.(lo))]
    return au, al
end

# replace any NaN entries by the nearest finite neighbour (keeps JSON MATLAB-safe)
function backfill!(v)
    good = findall(x -> isfinite(real(x)) && isfinite(imag(x)), v)
    isempty(good) && error("a branch is entirely non-finite; adjust omega range / sigma0")
    for j in eachindex(v)
        if !(isfinite(real(v[j])) && isfinite(imag(v[j])))
            j_near = good[argmin(abs.(good .- j))]
            v[j] = v[j_near]
        end
    end
    return v
end

# first frame: classify independently at every point of the L line
function branches_init(model, L, coarse)
    au = ComplexF64[]; al = ComplexF64[]
    for w in L
        a, b = leaders_at(model, w, coarse)
        push!(au, a); push!(al, b)
    end
    return backfill!(au), backfill!(al)
end

# later frames: track each point from the previous frame by continuity
function branches_track(model, L, prev_au, prev_al, coarse)
    au = similar(prev_au); al = similar(prev_al)
    for (j, w) in enumerate(L)
        roots, _ = physical_spatial_roots(model, w; model_coarse = coarse)
        if isempty(roots)
            au[j] = prev_au[j]; al[j] = prev_al[j]; continue
        end
        du = abs.(roots .- prev_au[j]); iu = argmin(du); au[j] = roots[iu]
        dl = abs.(roots .- prev_al[j]); dl[iu] = Inf
        al[j] = length(roots) == 1 ? roots[iu] : roots[argmin(dl)]
    end
    return au, al
end

function main()
    model  = CouetteModel(; Re = RE, num_modes = NUM_MODES)
    coarse = CouetteModel(; Re = RE, num_modes = NUM_COARSE)

    # F : real alpha axis (fixed)
    F = ComplexF64[a + 0im for a in range(ALPHA_R0, ALPHA_R1, length = N_F)]

    # omega_F : dominant temporal eigenvalue over F (constant, F fixed)
    omega_F = ComplexF64[]
    for a in F
        ev = temporal_spectrum(model, a)
        push!(omega_F, ev[argmax(imag.(ev))])
    end

    sig = collect(range(SIGMA0, 0.0, length = NSTEPS))

    L = line_L(sig[1])
    au, al = branches_init(model, L, coarse)

    frames = Any[frame_dict(1, L, au, al, F, omega_F)]
    @printf("k=%4d  sigma=%.5f  min|au-al|=%.3e\n", 1, sig[1], minimum(abs.(au .- al)))

    for k in 2:NSTEPS
        L = line_L(sig[k])
        au, al = branches_track(model, L, au, al, coarse)
        push!(frames, frame_dict(k, L, au, al, F, omega_F))
        if k % 10 == 0 || k == NSTEPS
            @printf("k=%4d  sigma=%.5f  min|au-al|=%.3e\n", k, sig[k], minimum(abs.(au .- al)))
        end
        flush(stdout)
    end

    open(OUTFILE, "w") do io
        write(io, JSON.json(frames))
    end
    @printf("\nwrote %d frames -> %s  (old schema; feed to plot_contour_deformation.m)\n",
            length(frames), OUTFILE)

    # read the causal branches at the forcing frequency for reference
    a_up, a_lo = leaders_at(model, OMEGA_F + 0.0im, coarse)
    @printf("at omega_f=%.3f :  alpha+ = %.5f%+.5fi   alpha- = %.5f%+.5fi\n",
            OMEGA_F, real(a_up), imag(a_up), real(a_lo), imag(a_lo))
end

isinteractive() || main()
