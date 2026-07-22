# ribbon_core.jl
# ---------------------------------------------------------------------------
# Model-agnostic Briggs logic for the vibrating-ribbon signalling problem.
#
# Ribbon = interior time-periodic point forcing f' = delta(x) delta(y-y') exp(-i w_f t)
# (NEAT paper Sec. III.B). The homogeneous OSE eigenvalue problems are unchanged;
# the ribbon only fixes omega to a REAL forcing frequency w_f and asks for the
# CAUSAL spatial branches alpha_plus(w_f) (feeds x>0) and alpha_minus(w_f) (feeds x<0).
#
# Prototype 1 computes only those branches by a fixed-w_f causal descent:
#   start at omega = w_f + i*sigma0 (sigma0 above all temporal singularities),
#   classify roots by side of the real F_alpha axis (sign of Im alpha),
#   lower Im omega toward 0 while tracking the two leaders by continuity,
#   stop and report if a genuine upper/lower pinch obstructs the descent,
#   otherwise return the causal roots at omega = w_f + i*0^+.
#
# Depends only on the model interface in models.jl:
#   spatial_spectrum(model, omega) -> (alphas, vecs)
# Serial. No contour relaxation, no potentials.
# ---------------------------------------------------------------------------

include("models.jl")
using Printf

# --- spectrum hygiene ------------------------------------------------------
# Keep physical roots: inside |alpha| < cutoff, and (if a coarser model is
# supplied) agreeing between the two resolutions to a relative tolerance.
# Analytic models pass model_coarse = nothing and cutoff = Inf.
function physical_spatial_roots(model, omega;
                                model_coarse = nothing,
                                cutoff::Float64 = 30.0,
                                rtol::Float64 = 1e-4)
    roots, vecs = spatial_spectrum(model, omega)
    keep = abs.(roots) .< cutoff
    roots = roots[keep]
    vecs = vecs === nothing ? nothing : vecs[:, keep]

    if model_coarse !== nothing
        coarse, _ = spatial_spectrum(model_coarse, omega)
        coarse = coarse[abs.(coarse) .< cutoff]
        agree = Bool[]
        for a in roots
            d = isempty(coarse) ? Inf : minimum(abs.(a .- coarse))
            push!(agree, d <= rtol * (1 + abs(a)))
        end
        roots = roots[agree]
        vecs = vecs === nothing ? nothing : vecs[:, agree]
    end
    return roots, vecs
end

# --- causal classification at the top of the descent -----------------------
# At omega = w_f + i*sigma0 the spatial contour F_alpha is the real axis, so
# the causal side is simply the sign of Im(alpha). Leaders = least-damped root
# on each side (closest to the real axis).
function classify_leaders(roots)
    upper = filter(a -> imag(a) > 0, roots)
    lower = filter(a -> imag(a) < 0, roots)
    isempty(upper) && error("no upper (Im>0) spatial root at sigma0; raise sigma0")
    isempty(lower) && error("no lower (Im<0) spatial root at sigma0; raise sigma0")
    alpha_plus  = upper[argmin(imag.(upper))]   # feeds x>0 (downstream)
    alpha_minus = lower[argmax(imag.(lower))]   # feeds x<0 (upstream)
    return alpha_plus, alpha_minus, length(upper), length(lower)
end

# nearest root to a guess (optionally excluding one index)
function _nearest(roots, guess; exclude::Int = 0)
    d = abs.(roots .- guess)
    exclude != 0 && (d[exclude] = Inf)
    i = argmin(d)
    return roots[i], i
end

# --- fixed-w_f causal descent ----------------------------------------------
# Returns a NamedTuple with the causal leaders and diagnostics.
function descend(model, w_f;
                 sigma0::Float64 = 0.3,
                 nsteps::Int = 20,
                 model_coarse = nothing,
                 pinch_tol::Float64 = 1e-3,
                 pinch_im_floor::Float64 = 1e-4,
                 jump_warn::Float64 = 0.2,
                 seed = nothing,
                 verbose::Bool = false)

    sig_grid = collect(range(sigma0, 0.0, length = nsteps + 1))

    roots, _ = physical_spatial_roots(model, w_f + im * sigma0; model_coarse = model_coarse)
    if seed === nothing
        ap, am, nu, nl = classify_leaders(roots)
    else
        ap, _ = _nearest(roots, seed[1])
        am, _ = _nearest(roots, seed[2])
        nu = count(a -> imag(a) > 0, roots); nl = count(a -> imag(a) < 0, roots)
    end

    min_gap   = Inf
    min_clear = minimum(abs.(imag.(roots)))
    trace = NamedTuple[]
    status = "OK"; omega_stop = w_f + 0.0im

    for s in sig_grid[2:end]
        omega = w_f + im * s
        rts, _ = physical_spatial_roots(model, omega; model_coarse = model_coarse)

        ap_new, ip = _nearest(rts, ap)
        am_new, _  = _nearest(rts, am; exclude = ip)

        jump = max(abs(ap_new - ap), abs(am_new - am))
        # gap of the upper leader to any root currently below it (opposite family)
        below = filter(a -> imag(a) < imag(ap_new), rts)
        gap = isempty(below) ? Inf : minimum(abs.(ap_new .- below))
        min_gap   = min(min_gap, gap)
        min_clear = min(min_clear, minimum(abs.(imag.(rts))))

        push!(trace, (sigma = s, alpha_plus = ap_new, alpha_minus = am_new,
                      gap = gap, jump = jump, n_up = count(a -> imag(a) > 0, rts),
                      n_lo = count(a -> imag(a) < 0, rts)))
        verbose && @printf("  sigma=%.4f  a+=%.5f%+.5fi  a-=%.5f%+.5fi  gap=%.3e  jump=%.3e\n",
                           s, real(ap_new), imag(ap_new), real(am_new), imag(am_new), gap, jump)

        ap, am = ap_new, am_new

        if gap < pinch_tol && s > pinch_im_floor
            status = "PINCH"; omega_stop = omega
            break
        end
    end

    # causal well-posedness check at the axis: leaders must still straddle it.
    straddles = imag(ap) > 0 && imag(am) < 0
    if status == "OK" && !straddles
        status = "NO_STRADDLE"   # a+ crossed the axis -> not convectively well-posed here
    end

    return (alpha_plus = ap, alpha_minus = am, status = status,
            omega_stop = omega_stop, min_gap = min_gap, min_clear = min_clear,
            n_up = nu, n_lo = nl, trace = trace)
end

# --- outer sweep over real forcing frequency -------------------------------
# Warm-starts each w_f from the previous leaders; every `check_every`-th
# frequency an independent cold start verifies branch identity.
function sweep(model, w_f_grid;
               sigma0::Float64 = 0.3, nsteps::Int = 20,
               model_coarse = nothing, check_every::Int = 10,
               tol_check::Float64 = 1e-6, verbose::Bool = true)

    results = NamedTuple[]
    prev = nothing
    for (k, w_f) in enumerate(w_f_grid)
        r = descend(model, w_f; sigma0 = sigma0, nsteps = nsteps,
                    model_coarse = model_coarse, seed = prev)

        if k % check_every == 0
            rc = descend(model, w_f; sigma0 = sigma0, nsteps = nsteps,
                         model_coarse = model_coarse, seed = nothing)
            dcheck = max(abs(r.alpha_plus - rc.alpha_plus),
                         abs(r.alpha_minus - rc.alpha_minus))
            dcheck > tol_check && @warn "warm/cold mismatch" w_f dcheck
        end

        push!(results, (w_f = w_f, alpha_plus = r.alpha_plus,
                        alpha_minus = r.alpha_minus, status = r.status,
                        min_gap = r.min_gap, min_clear = r.min_clear))
        verbose && @printf("w_f=%+.4f  a+=%.5f%+.5fi  a-=%.5f%+.5fi  %s  clear=%.3e\n",
                           w_f, real(r.alpha_plus), imag(r.alpha_plus),
                           real(r.alpha_minus), imag(r.alpha_minus), r.status, r.min_clear)
        prev = (r.alpha_plus, r.alpha_minus)
    end
    return results
end
