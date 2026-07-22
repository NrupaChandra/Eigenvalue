# run_sweep.jl
# ---------------------------------------------------------------------------
# Driver for the outer sweep over real forcing frequency w_f.
# Warm-starts each frequency from the previous leaders and cold-checks every
# `check_every`-th point. Saves one JSON record per w_f to results/.
#
# Usage (from the vib_ribbon directory):
#     julia run_sweep.jl
# ---------------------------------------------------------------------------

include("ribbon_core.jl")
using JSON, Printf

function save_sweep(results, filename)
    mkpath(dirname(filename))   # create results/ if it does not exist yet
    arr = [Dict("w_f" => r.w_f,
                "alpha_plus"  => Dict("re" => real(r.alpha_plus),  "im" => imag(r.alpha_plus)),
                "alpha_minus" => Dict("re" => real(r.alpha_minus), "im" => imag(r.alpha_minus)),
                "status" => r.status, "min_gap" => r.min_gap, "min_clear" => r.min_clear)
           for r in results]
    open(filename, "w") do io
        write(io, JSON.json(arr))
    end
    println("saved ", length(arr), " records -> ", filename)
end

function main()
    # -- analytic stable sweep: closed-form check at every frequency ---------
    mS = analytic_stable()
    grid = collect(range(-0.9, 0.9, length = 37))
    println("== Case S sweep (stable analytic) ==")
    resS = sweep(mS, grid; sigma0 = 0.6, nsteps = 25, check_every = 5)
    # exact-match assertion
    worst = 0.0
    for r in resS
        ep, em = exact_roots(mS, r.w_f + 0.0im)
        worst = max(worst, abs(r.alpha_plus - ep), abs(r.alpha_minus - em))
    end
    @printf("Case S: worst |error| vs exact over sweep = %.3e (expect < 1e-10)\n", worst)
    save_sweep(resS, joinpath("results", "sweep_analytic_stable.json"))

    # -- Couette sweep -------------------------------------------------------
    println("\n(building Couette operators N=100 and N=150 ...)")
    mC  = CouetteModel(; Re = 2000.0, num_modes = 150)
    mCc = CouetteModel(; Re = 2000.0, num_modes = 100)
    gridC = collect(range(0.1, 0.9, length = 17))
    println("== Couette sweep Re=2000 beta=0 ==")
    resC = sweep(mC, gridC; sigma0 = 0.3, nsteps = 20, model_coarse = mCc, check_every = 4)
    save_sweep(resC, joinpath("results", "sweep_couette_Re2000.json"))
end

isinteractive() || main()
