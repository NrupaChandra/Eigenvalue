# run_classify.jl
# ---------------------------------------------------------------------------
# Driver for ONE forcing frequency. Runs the causal descent on a chosen model
# and prints the causal spatial branches alpha_plus / alpha_minus.
#
# Usage (from the vib_ribbon directory):
#     julia run_classify.jl                 # runs all three built-in checks
#
# Edit the `main` calls at the bottom to try other (model, w_f, sigma0).
# ---------------------------------------------------------------------------

include("ribbon_core.jl")
using Printf

function report(title, model, w_f; sigma0, nsteps, model_coarse = nothing,
                expect = nothing)
    println("\n==================================================================")
    println(" $title")
    println(" model: ", model.name, "   w_f = ", w_f, "   sigma0 = ", sigma0)
    println("==================================================================")

    r = descend(model, w_f; sigma0 = sigma0, nsteps = nsteps,
                model_coarse = model_coarse, verbose = true)

    @printf("\n status        : %s\n", r.status)
    @printf(" alpha_plus    : %.10f %+.10fi\n", real(r.alpha_plus), imag(r.alpha_plus))
    @printf(" alpha_minus   : %.10f %+.10fi\n", real(r.alpha_minus), imag(r.alpha_minus))
    @printf(" min gap (a+/-): %.3e\n", r.min_gap)
    @printf(" min axis clear: %.3e   (n_up=%d n_lo=%d)\n", r.min_clear, r.n_up, r.n_lo)
    r.status == "PINCH" && @printf(" pinch at omega: %.6f %+.6fi\n",
                                   real(r.omega_stop), imag(r.omega_stop))

    if expect !== nothing
        ep, em = expect
        @printf(" exact a+      : %.10f %+.10fi   |err|=%.3e\n",
                real(ep), imag(ep), abs(r.alpha_plus - ep))
        @printf(" exact a-      : %.10f %+.10fi   |err|=%.3e\n",
                real(em), imag(em), abs(r.alpha_minus - em))
    end
    return r
end

function main()
    # -- Case S: stable analytic. Descent should COMPLETE (status OK) and the
    #    leaders should match the closed form to ~1e-10 for any sigma0/nsteps.
    mS = analytic_stable()
    wS = 0.5
    report("CASE S (stable analytic) -- expect OK, exact match",
           mS, wS; sigma0 = 0.6, nsteps = 30, expect = exact_roots(mS, wS + 0.0im))

    # -- Case U: absolutely unstable analytic. Descending the line Re w = Re(w0)
    #    must HIT the branch point at omega0 = -0.88 + 0.8i (status PINCH),
    #    where a+ and a- coincide at alpha0 = -0.2 - 0.4i.
    #    NOTE sigma0 must be high enough that the two roots still STRADDLE the
    #    real axis at the top (for this fixture they straddle only above
    #    sigma ~ 1.3; sigma0 = 3.0 is safe). This is what "above all relevant
    #    temporal singularities" means in practice for an unstable case.
    mU = analytic_unstable()
    a0, w0 = saddle(mU)
    report("CASE U (unstable analytic) -- expect PINCH at omega0",
           mU, real(w0); sigma0 = 3.0, nsteps = 60)
    @printf(" (reference: omega0 = %.4f %+.4fi , alpha0 = %.4f %+.4fi)\n",
            real(w0), imag(w0), real(a0), imag(a0))

    # -- Couette: two resolutions for the spurious filter. Expect OK, with
    #    least-damped a+ ~ 0.7496 + 0.1016i and a- ~ 0.059 - 3.275i at w_f=0.5.
    println("\n(building Couette operators N=100 and N=150 ...)")
    mC   = CouetteModel(; Re = 2000.0, num_modes = 150)
    mCc  = CouetteModel(; Re = 2000.0, num_modes = 100)
    report("COUETTE Re=2000 beta=0 -- expect OK, causal branches",
           mC, 0.5; sigma0 = 0.3, nsteps = 20, model_coarse = mCc)
end

isinteractive() || main()
