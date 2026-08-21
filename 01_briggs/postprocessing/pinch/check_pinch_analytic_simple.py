"""
check_pinch_exact_analytic.py

The two analytical Python functions below are copied directly from the
original Julia/PyCall source:

    fnum_python(...)
    fnum_derivatives_python(...)

They are intentionally left unchanged.

This script only adds a small wrapper at the bottom which:
  1. reads alpha_p and omega_p from the pinch JSON,
  2. evaluates the analytical dispersion relation D(alpha_p, omega_p),
  3. evaluates the analytical group velocity dω/dα,
  4. prints the results.

For a pinch / saddle point we expect:

    D(alpha_p, omega_p) ≈ 0
    dω/dα               ≈ 0
"""

import json
import mpmath


# =====================================================================
# USER SETTINGS
# =====================================================================

PINCH_JSON = r"C:\Git\Eigenvalue\01_briggs\results\json\contour_iteration_v4.2_pinch_simple.json"

PRECISION_DECIMAL_POINTS = 40
DEFAULT_RE = 2000.0


# =====================================================================
# EXACT ANALYTICAL DISPERSION-RELATION FUNCTION FROM ORIGINAL SCRIPT
# =====================================================================

def fnum_python(a, omega, re, precision_decimal_points_py):
    mpmath.mp.dps = precision_decimal_points_py

    def integrand_ipa0(y):
        return mpmath.exp(a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_ina2(y):
        return mpmath.exp(-a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_ipa2(y):
        return mpmath.exp(a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_ina0(y):
        return mpmath.exp(-a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    attempts = 1
    while attempts <= 15:
        try:
            ipa0, error_ipa0, *_ = mpmath.quad(lambda y: integrand_ipa0(y), [0, 1], full_output=True, error=True)
            ina2, error_ina2, *_ = mpmath.quad(lambda y: integrand_ina2(y), [0, 1], full_output=True, error=True)
            ipa2, error_ipa2, *_ = mpmath.quad(lambda y: integrand_ipa2(y), [0, 1], full_output=True, error=True)
            ina0, error_ina0, *_ = mpmath.quad(lambda y: integrand_ina0(y), [0, 1], full_output=True, error=True)

            max_error = max(error_ipa0, error_ina2, error_ipa2, error_ina0)

            disprel = (-1) / (2 * a) * (ipa0 * ina2 - ipa2 * ina0)

            if max_error > mpmath.power(10, -precision_decimal_points_py):
                precision_decimal_points_py *= 2
                mpmath.mp.dps = precision_decimal_points_py
                attempts += 1
            else:
                return disprel, error_ipa0, error_ina2, error_ipa2, error_ina0

        except Exception as e:
            print(f"Exception occurred: {e}")
            precision_decimal_points_py *= 2
            mpmath.mp.dps = precision_decimal_points_py
            attempts += 1
            continue

    disprel = mpmath.nan
    error_ipa0 = mpmath.nan
    error_ina2 = mpmath.nan
    error_ina0 = mpmath.nan
    error_ipa2 = mpmath.nan

    return disprel, error_ipa0, error_ina2, error_ipa2, error_ina0


# =====================================================================
# EXACT ANALYTICAL DERIVATIVE FUNCTION FROM ORIGINAL SCRIPT
# =====================================================================

def fnum_derivatives_python(a, omega, re, precision_decimal_points_py):
    mpmath.mp.dps = precision_decimal_points_py

    def integrand_ipa0(y):
        return mpmath.exp(a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_ina2(y):
        return mpmath.exp(-a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_ipa2(y):
        return mpmath.exp(a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_ina0(y):
        return mpmath.exp(-a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_d_ipa0_d_omega(y):
        return mpmath.exp(a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1)

    def integrand_d_ina2_d_omega(y):
        return mpmath.exp(-a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1)

    def integrand_d_ipa2_d_omega(y):
        return mpmath.exp(a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1)

    def integrand_d_ina0_d_omega(y):
        return mpmath.exp(-a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1)

    def integrand_d_ipa0_d_a(y):
        return y * mpmath.exp(a * y) * mpmath.airyai(1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (
                    omega * re - a * y * re + 1j * mpmath.power(a, 2))) + mpmath.exp(a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1) * (1j * a * y * re + 2j * omega * re + 4 * mpmath.power(a, 2)) / (
                    3 * a * mpmath.power((-1j * a * re), (2 / 3)))

    def integrand_d_ina2_d_a(y):
        return (-y) * mpmath.exp(-a * y) * mpmath.airybi(1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (
                    omega * re - a * y * re + 1j * mpmath.power(a, 2))) + mpmath.exp(-a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1) * (1j * a * y * re + 2j * omega * re + 4 * mpmath.power(a, 2)) / (
                    3 * a * mpmath.power((-1j * a * re), (2 / 3)))

    def integrand_d_ipa2_d_a(y):
        return y * mpmath.exp(a * y) * mpmath.airybi(1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (
                    omega * re - a * y * re + 1j * mpmath.power(a, 2))) + mpmath.exp(a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1) * (1j * a * y * re + 2j * omega * re + 4 * mpmath.power(a, 2)) / (
                    3 * a * mpmath.power((-1j * a * re), (2 / 3)))

    def integrand_d_ina0_d_a(y):
        return (-y) * mpmath.exp(-a * y) * mpmath.airyai(1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (
                    omega * re - a * y * re + 1j * mpmath.power(a, 2))) + mpmath.exp(-a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), (1 / 3)) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1) * (1j * a * y * re + 2j * omega * re + 4 * mpmath.power(a, 2)) / (
                    3 * a * mpmath.power((-1j * a * re), (2 / 3)))

    attempts = 1
    while attempts <= 5:
        try:
            ipa0, error_ipa0, *_ = mpmath.quad(lambda y: integrand_ipa0(y), [0, 1], full_output=True, error=True)
            ina2, error_ina2, *_ = mpmath.quad(lambda y: integrand_ina2(y), [0, 1], full_output=True, error=True)
            ipa2, error_ipa2, *_ = mpmath.quad(lambda y: integrand_ipa2(y), [0, 1], full_output=True, error=True)
            ina0, error_ina0, *_ = mpmath.quad(lambda y: integrand_ina0(y), [0, 1], full_output=True, error=True)
            d_ipa0_d_omega, error_d_ipa0_d_omega, *_ = mpmath.quad(lambda y: integrand_d_ipa0_d_omega(y), [0, 1],
                                                                   full_output=True, error=True)
            d_ina2_d_omega, error_d_ina2_d_omega, *_ = mpmath.quad(lambda y: integrand_d_ina2_d_omega(y), [0, 1],
                                                                   full_output=True, error=True)
            d_ipa2_d_omega, error_d_ipa2_d_omega, *_ = mpmath.quad(lambda y: integrand_d_ipa2_d_omega(y), [0, 1],
                                                                   full_output=True, error=True)
            d_ina0_d_omega, error_d_ina0_d_omega, *_ = mpmath.quad(lambda y: integrand_d_ina0_d_omega(y), [0, 1],
                                                                   full_output=True, error=True)
            d_ipa0_d_a, error_d_ipa0_d_a, *_ = mpmath.quad(lambda y: integrand_d_ipa0_d_a(y), [0, 1], full_output=True,
                                                           error=True)
            d_ina2_d_a, error_d_ina2_d_a, *_ = mpmath.quad(lambda y: integrand_d_ina2_d_a(y), [0, 1], full_output=True,
                                                           error=True)
            d_ipa2_d_a, error_d_ipa2_d_a, *_ = mpmath.quad(lambda y: integrand_d_ipa2_d_a(y), [0, 1], full_output=True,
                                                           error=True)
            d_ina0_d_a, error_d_ina0_d_a, *_ = mpmath.quad(lambda y: integrand_d_ina0_d_a(y), [0, 1], full_output=True,
                                                           error=True)

            max_error = max(error_ipa0, error_ina2, error_ipa2, error_ina0, error_d_ipa0_d_omega, error_d_ina2_d_omega,
                            error_d_ipa2_d_omega, error_d_ina0_d_omega, error_d_ipa0_d_a, error_d_ina2_d_a,
                            error_d_ipa2_d_a, error_d_ina0_d_a)

            d_disprel_d_omega = (-1) / (2 * a) * mpmath.power((-1j * a * re), (1 / 3)) / a * (
                        d_ipa0_d_omega * ina2 + ipa0 * d_ina2_d_omega - (d_ipa2_d_omega * ina0 + ipa2 * d_ina0_d_omega))
            d_disprel_d_a = 1 / (2 * mpmath.power(a, 2)) * (ipa0 * ina2 - ipa2 * ina0) - 1 / (2 * a) * (
                        d_ipa0_d_a * ina2 + ipa0 * d_ina2_d_a - (d_ipa2_d_a * ina0 + ipa2 * d_ina0_d_a))

            d_omega_d_alpha = -d_disprel_d_a / d_disprel_d_omega

            if max_error > mpmath.power(10, -precision_decimal_points_py):
                precision_decimal_points_py *= 2
                mpmath.mp.dps = precision_decimal_points_py
                attempts += 1
            else:
                return (d_omega_d_alpha, error_ipa0, error_ina2, error_ipa2, error_ina0, error_d_ipa0_d_omega,
                        error_d_ina2_d_omega, error_d_ipa2_d_omega, error_d_ina0_d_omega, error_d_ipa0_d_a,
                        error_d_ina2_d_a, error_d_ipa2_d_a, error_d_ina0_d_a)

        except Exception as e:
            print(f"Exception occurred: {e}")
            precision_decimal_points_py *= 2
            mpmath.mp.dps = precision_decimal_points_py
            attempts += 1
            continue

    d_omega_d_alpha = mpmath.nan
    error_ipa0 = mpmath.nan
    error_ina2 = mpmath.nan
    error_ina0 = mpmath.nan
    error_ipa2 = mpmath.nan
    error_d_ipa0_d_omega = mpmath.nan
    error_d_ina2_d_omega = mpmath.nan
    error_d_ipa2_d_omega = mpmath.nan
    error_d_ina0_d_omega = mpmath.nan
    error_d_ipa0_d_a = mpmath.nan
    error_d_ina2_d_a = mpmath.nan
    error_d_ipa2_d_a = mpmath.nan
    error_d_ina0_d_a = mpmath.nan

    return (d_omega_d_alpha, error_ipa0, error_ina2, error_ipa2, error_ina0, error_d_ipa0_d_omega, error_d_ina2_d_omega,
            error_d_ipa2_d_omega, error_d_ina0_d_omega, error_d_ipa0_d_a, error_d_ina2_d_a, error_d_ipa2_d_a,
            error_d_ina0_d_a)


# =====================================================================
# EVERYTHING BELOW THIS LINE IS ONLY THE NEW VALIDATION WRAPPER
# =====================================================================

with open(PINCH_JSON, "r") as f:
    pinch = json.load(f)

Re = mpmath.mpf(repr(pinch.get("Re", DEFAULT_RE)))

alpha_p = mpmath.mpc(
    repr(pinch["alpha_p"]["re"]),
    repr(pinch["alpha_p"]["im"])
)

omega_p = mpmath.mpc(
    repr(pinch["omega_p"]["re"]),
    repr(pinch["omega_p"]["im"])
)


print("=" * 72)
print("ANALYTIC PINCH CHECK — ORIGINAL FUNCTIONS UNCHANGED")
print("=" * 72)

print("\nInput pinch:")
print("alpha_p =", alpha_p)
print("omega_p =", omega_p)
print("Re      =", Re)

print("\nEvaluating original analytical dispersion relation...")
D_result = fnum_python(
    alpha_p,
    omega_p,
    Re,
    PRECISION_DECIMAL_POINTS
)

D_value = D_result[0]

print("Evaluating original analytical derivative function...")
derivative_result = fnum_derivatives_python(
    alpha_p,
    omega_p,
    Re,
    PRECISION_DECIMAL_POINTS
)

domega_dalpha = derivative_result[0]


print("\nRESULTS")
print("-" * 72)

print("D(alpha_p, omega_p) =")
print(D_value)
print()
print("|D| =")
print(abs(D_value))

print("\ndomega/dalpha =")
print(domega_dalpha)
print()
print("|domega/dalpha| =")
print(abs(domega_dalpha))

