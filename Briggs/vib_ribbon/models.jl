# models.jl
# ---------------------------------------------------------------------------
# The PHYSICS layer for the vibrating-ribbon prototype.
#
# A "model" provides the dispersion relation D(alpha, omega) = 0 through two
# functions the Briggs core depends on:
#
#     spatial_spectrum(model, omega)  -> (alphas, vecs)   # roots in alpha at fixed omega
#     temporal_spectrum(model, alpha) -> omegas           # roots in omega at fixed alpha
#
# Analytic models additionally provide exact_roots(model, omega) for validation.
#
# Nothing in this file knows about causal descent, contours, or sweeps. That is
# the whole point: ribbon_core.jl is written against the interface below, so the
# Couette model and the analytic test models are interchangeable.
#
# Serial only. No Distributed, no globals leaking between models.
# ---------------------------------------------------------------------------

using LinearAlgebra

abstract type RibbonModel end

# ===========================================================================
# 1. Analytic quadratic test model
#    Convention taken from Briggs/briggs_test.jl:
#        D(alpha, omega) = omega - (C1 + C2*alpha + C3*alpha^2) = 0
#    so  omega(alpha) = C1 + C2*alpha + C3*alpha^2
#        alpha_pm(omega) = (-C2 +/- sqrt(C2^2 - 4 C3 (C1 - omega))) / (2 C3)
#        saddle:  alpha0 = -C2/(2 C3),  omega0 = C1 - C2^2/(4 C3)
# ===========================================================================
struct AnalyticModel <: RibbonModel
    C1::ComplexF64
    C2::ComplexF64
    C3::ComplexF64
    name::String
end

# The two validation fixtures (constants independently verified in SPEC.md).
#   Case U: absolutely unstable  -> Im(omega0) = +0.80  -> descent must PINCH.
#   Case S: stable               -> Im(omega0) = -1.04  -> descent completes.
analytic_unstable() = AnalyticModel(-1.0 + 0.96im, 0.4 + 0.8im, 1.0 + 0.0im,
                                    "analytic-unstable (Case U)")
analytic_stable()   = AnalyticModel(-1.0 - 0.96im, 0.4 + 0.8im, 1.0 - 0.5im,
                                    "analytic-stable (Case S)")

function spatial_spectrum(m::AnalyticModel, omega::Number)
    disc = sqrt(m.C2^2 - 4.0 * m.C3 * (m.C1 - omega))
    a1 = (-m.C2 + disc) / (2.0 * m.C3)
    a2 = (-m.C2 - disc) / (2.0 * m.C3)
    return ComplexF64[a1, a2], nothing        # analytic model has no eigenvectors
end

function temporal_spectrum(m::AnalyticModel, alpha::Number)
    return ComplexF64[m.C1 + m.C2 * alpha + m.C3 * alpha^2]
end

# Exact causal roots, ordered (alpha_plus, alpha_minus) by Im (upper first).
function exact_roots(m::AnalyticModel, omega::Number)
    roots, _ = spatial_spectrum(m, omega)
    a1, a2 = roots[1], roots[2]
    return imag(a1) >= imag(a2) ? (a1, a2) : (a2, a1)
end

saddle(m::AnalyticModel) = (-m.C2 / (2 * m.C3), m.C1 - m.C2^2 / (4 * m.C3))

# ===========================================================================
# 2. Couette model
#    Orr-Sommerfeld for plane Couette flow U(y) = y on [0, 1], beta = 0.
#    Matrix construction copied verbatim (structure) from Briggs/briggsv4.jl
#    lines 12-89. Spatial problem is the exact quadratic companion pencil of
#    the quartic OSE (linearised via v = phi*exp(-alpha y); verified elsewhere).
# ===========================================================================
struct CouetteOps
    num_modes::Int
    Re::Float64
    beta::ComplexF64
    D0::Matrix{Float64}
    D1::Matrix{Float64}
    D2::Matrix{Float64}
    D3::Matrix{Float64}
    D4::Matrix{Float64}
    u::Vector{Float64}
    v_g::ComplexF64
    d2u::Float64
end

function build_couette_ops(num_modes::Int; Re::Float64 = 2000.0,
                           beta::ComplexF64 = 0.0 + 0.0im,
                           v_g::ComplexF64 = 0.0 + 0.0im,
                           start::Float64 = 0.0, terminate::Float64 = 1.0)
    y = [cos((j - 1) * pi / (num_modes - 1)) for j = 1:num_modes]
    ynew = ((start + terminate) / 2) .- y * ((terminate - start) / 2)

    D0 = zeros(Float64, num_modes, num_modes)
    for j = 1:num_modes
        D0[:, j] .= cos.((j - 1) * acos.(y))
    end

    D1 = zeros(Float64, num_modes, num_modes)
    D2 = zeros(Float64, num_modes, num_modes)
    D3 = zeros(Float64, num_modes, num_modes)
    D4 = zeros(Float64, num_modes, num_modes)
    D1[:, 1:3] .= [zeros(num_modes) D0[:, 1] 4 * D0[:, 2]]
    D2[:, 1:3] .= [zeros(num_modes) zeros(num_modes) 4 * D0[:, 1]]
    for j = 4:num_modes
        D1[:, j] .= 2 * (j - 1) * D0[:, j-1] + (j - 1) * D1[:, j-2] / (j - 3)
        D2[:, j] .= 2 * (j - 1) * D1[:, j-1] + (j - 1) * D2[:, j-2] / (j - 3)
        D3[:, j] .= 2 * (j - 1) * D2[:, j-1] + (j - 1) * D3[:, j-2] / (j - 3)
        D4[:, j] .= 2 * (j - 1) * D3[:, j-1] + (j - 1) * D4[:, j-2] / (j - 3)
    end
    sc = -(terminate - start) / 2
    D1 ./= sc^1
    D2 ./= sc^2
    D3 ./= sc^3
    D4 ./= sc^4

    return CouetteOps(num_modes, Re, beta, D0, D1, D2, D3, D4, ynew, v_g, 0.0)
end

struct CouetteModel <: RibbonModel
    Re::Float64
    beta::ComplexF64
    num_modes::Int
    ops::CouetteOps
    name::String
end

function CouetteModel(; Re::Float64 = 2000.0, beta::ComplexF64 = 0.0 + 0.0im,
                      num_modes::Int = 150)
    ops = build_couette_ops(num_modes; Re = Re, beta = beta)
    return CouetteModel(Re, beta, num_modes, ops, "couette Re=$Re N=$num_modes")
end

# helper: diagonal(u) acting from the left, matching (u*ones(1,N)).*M in briggsv4
_scale_rows(u, M) = u .* M

# Temporal generalised eigenproblem A v = omega B v (size N), fixed alpha.
function _couette_temporal_matrices(ops::CouetteOps, alpha::Number)
    N = ops.num_modes; Re = ops.Re; beta = ops.beta; u = ops.u; d2u = ops.d2u
    v_g = ops.v_g
    D0, D1, D2, D4 = ops.D0, ops.D1, ops.D2, ops.D4
    k2 = alpha^2 + beta^2

    A = -im * alpha * _scale_rows(u, D2) .+ im * alpha * k2 * _scale_rows(u, D0) .+
        im * alpha * d2u * D0 .+ (1 / Re) * D4 .- (2 / Re) * k2 * D2 .+
        (1 / Re) * k2^2 * D0 .+ alpha * v_g * D0
    A = ComplexF64.(A)
    A[1:2, :]     .= -200im * [D0[1:1, :]; D1[1:1, :]]
    A[N-1:N, :]   .= -200im * [D1[N:N, :]; D0[N:N, :]]

    B = -im * D2 .+ im * k2 * D0
    B = ComplexF64.(B)
    B[1:2, :]     .= [D0[1:1, :]; D1[1:1, :]]
    B[N-1:N, :]   .= [D1[N:N, :]; D0[N:N, :]]
    return A, B
end

# Spatial generalised eigenproblem (companion form, size 2N), fixed omega.
function _couette_spatial_matrices(ops::CouetteOps, omega::Number)
    N = ops.num_modes; Re = ops.Re; beta = ops.beta; u = ops.u; d2u = ops.d2u
    v_g = ops.v_g
    D0, D1, D2, D3, D4 = ops.D0, ops.D1, ops.D2, ops.D3, ops.D4

    A11 = -2im * omega * D1 .- (4 / Re) * D3 .+ (4 / Re) * beta^2 * D1 .-
          im * _scale_rows(u, D2) .+ im * beta^2 * _scale_rows(u, D0) .+
          im * d2u * D0 .- im * v_g * D2 .+ im * v_g * beta^2 * D0
    A12 = im * omega * D2 .- im * omega * beta^2 * D0 .+ (1 / Re) * D4 .-
          (2 / Re) * beta^2 * D2 .+ (1 / Re) * beta^4 * D0
    A11 = ComplexF64.(A11); A12 = ComplexF64.(A12)
    A11[1:2, :] .= 0; A11[N-1:N, :] .= 0
    A12[1:2, :]   .= -200im * [D0[1:1, :]; D1[1:1, :]]
    A12[N-1:N, :] .= -200im * [D1[N:N, :]; D0[N:N, :]]

    A21 = Matrix{ComplexF64}(I, N, N)
    A22 = zeros(ComplexF64, N, N)
    A = [A11 A12; A21 A22]

    B11 = -(4 / Re) * D2 .- 2im * _scale_rows(u, D1) .+ 2im * v_g * D1
    B11 = ComplexF64.(B11)
    B11[1:2, :] .= 0; B11[N-1:N, :] .= 0
    B = [B11 zeros(ComplexF64, N, N); zeros(ComplexF64, N, N) Matrix{ComplexF64}(I, N, N)]
    return A, B
end

_finite(vals) = [isfinite(real(v)) && isfinite(imag(v)) for v in vals]

function spatial_spectrum(m::CouetteModel, omega::Number)
    A, B = _couette_spatial_matrices(m.ops, omega)
    F = eigen(A, B)
    mask = _finite(F.values)
    # eigenvector top N block is the physical v-tilde amplitude (companion form)
    return ComplexF64.(F.values[mask]), F.vectors[1:m.num_modes, mask]
end

function temporal_spectrum(m::CouetteModel, alpha::Number)
    A, B = _couette_temporal_matrices(m.ops, alpha)
    vals = eigen(A, B).values
    return ComplexF64.(vals[_finite(vals)])
end
