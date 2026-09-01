"""
check_pinch_analytic_simple.py

Reads alpha_p / omega_p and the seed from a compute_pinch*.py output JSON and
evaluates the analytical dispersion relation D and the group velocity
domega/dalpha at both points.  At the pinch both should be ~0; at the seed
both should be large.

fnum_python() and fnum_derivatives_python() are the original Julia/PyCall
functions, with the cube-root exponents changed from Python floats to exact
mpmath thirds.
"""

import json
import mpmath

# =====================================================================
# USER SETTINGS
# =====================================================================

PINCH_JSON = r"C:\git\Eigenvalue\01_briggs\results\json\contour_iteration_v4.4_pinch_analytic.json"

PRECISION_DECIMAL_POINTS = 30



# =====================================================================
# ANALYTICAL DISPERSION-RELATION FUNCTION FROM ORIGINAL SCRIPT
# =====================================================================

def fnum_python(a, omega, re, precision_decimal_points_py):
    mpmath.mp.dps = precision_decimal_points_py

    # exact thirds, built at the current working precision
    THIRD = mpmath.mpf(1) / 3

    def integrand_ipa0(y):
        return mpmath.exp(a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_ina2(y):
        return mpmath.exp(-a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_ipa2(y):
        return mpmath.exp(a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_ina0(y):
        return mpmath.exp(-a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    attempts = 1
    while attempts <= 15:
        # the integrands close over the NAME, so rebuild at the current mp.dps
        THIRD = mpmath.mpf(1) / 3
        try:
            ipa0, error_ipa0, *_ = mpmath.quad(lambda y: integrand_ipa0(y), [0, 1], full_output=True, error=True)
            ina2, error_ina2, *_ = mpmath.quad(lambda y: integrand_ina2(y), [0, 1], full_output=True, error=True)
            ipa2, error_ipa2, *_ = mpmath.quad(lambda y: integrand_ipa2(y), [0, 1], full_output=True, error=True)
            ina0, error_ina0, *_ = mpmath.quad(lambda y: integrand_ina0(y), [0, 1], full_output=True, error=True)

            max_error = max(error_ipa0, error_ina2, error_ipa2, error_ina0)

            disprel = (-1) / (2 * a) * (ipa0 * ina2 - ipa2 * ina0)

            # NOTE: absolute error vs 10^-dps on integrals of size ~1e13; left
            # as-is, and in practice this branch never fires
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
# ANALYTICAL DERIVATIVE FUNCTION FROM ORIGINAL SCRIPT
# =====================================================================

def fnum_derivatives_python(a, omega, re, precision_decimal_points_py):
    mpmath.mp.dps = precision_decimal_points_py

    # exact thirds, built at the current working precision
    THIRD = mpmath.mpf(1) / 3
    TWOTHIRD = mpmath.mpf(2) / 3

    def integrand_ipa0(y):
        return mpmath.exp(a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_ina2(y):
        return mpmath.exp(-a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_ipa2(y):
        return mpmath.exp(a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_ina0(y):
        return mpmath.exp(-a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)))

    def integrand_d_ipa0_d_omega(y):
        return mpmath.exp(a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1)

    def integrand_d_ina2_d_omega(y):
        return mpmath.exp(-a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1)

    def integrand_d_ipa2_d_omega(y):
        return mpmath.exp(a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1)

    def integrand_d_ina0_d_omega(y):
        return mpmath.exp(-a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1)

    def integrand_d_ipa0_d_a(y):
        return y * mpmath.exp(a * y) * mpmath.airyai(1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (
                    omega * re - a * y * re + 1j * mpmath.power(a, 2))) + mpmath.exp(a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1) * (1j * a * y * re + 2j * omega * re + 4 * mpmath.power(a, 2)) / (
                    3 * a * mpmath.power((-1j * a * re), TWOTHIRD))

    def integrand_d_ina2_d_a(y):
        return (-y) * mpmath.exp(-a * y) * mpmath.airybi(1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (
                    omega * re - a * y * re + 1j * mpmath.power(a, 2))) + mpmath.exp(-a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1) * (1j * a * y * re + 2j * omega * re + 4 * mpmath.power(a, 2)) / (
                    3 * a * mpmath.power((-1j * a * re), TWOTHIRD))

    def integrand_d_ipa2_d_a(y):
        return y * mpmath.exp(a * y) * mpmath.airybi(1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (
                    omega * re - a * y * re + 1j * mpmath.power(a, 2))) + mpmath.exp(a * y) * mpmath.airybi(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1) * (1j * a * y * re + 2j * omega * re + 4 * mpmath.power(a, 2)) / (
                    3 * a * mpmath.power((-1j * a * re), TWOTHIRD))

    def integrand_d_ina0_d_a(y):
        return (-y) * mpmath.exp(-a * y) * mpmath.airyai(1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (
                    omega * re - a * y * re + 1j * mpmath.power(a, 2))) + mpmath.exp(-a * y) * mpmath.airyai(
            1 / (a * re) * mpmath.power((-1j * a * re), THIRD) * (omega * re - a * y * re + 1j * mpmath.power(a, 2)),
            derivative=1) * (1j * a * y * re + 2j * omega * re + 4 * mpmath.power(a, 2)) / (
                    3 * a * mpmath.power((-1j * a * re), TWOTHIRD))

    attempts = 1
    while attempts <= 5:
        # the integrands close over the NAMES, so rebuild at the current mp.dps
        THIRD = mpmath.mpf(1) / 3
        TWOTHIRD = mpmath.mpf(2) / 3
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

            d_disprel_d_omega = (-1) / (2 * a) * (
                mpmath.power((-1j * a * re), THIRD)
            ) / a * (
                        d_ipa0_d_omega * ina2 + ipa0 * d_ina2_d_omega - (d_ipa2_d_omega * ina0 + ipa2 * d_ina0_d_omega))
            d_disprel_d_a = 1 / (2 * mpmath.power(a, 2)) * (ipa0 * ina2 - ipa2 * ina0) - 1 / (2 * a) * (
                        d_ipa0_d_a * ina2 + ipa0 * d_ina2_d_a - (d_ipa2_d_a * ina0 + ipa2 * d_ina0_d_a))

            d_omega_d_alpha = -d_disprel_d_a / d_disprel_d_omega

            # NOTE: same absolute-vs-relative tolerance question as in fnum_python
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
# WRAPPER
# =====================================================================

def read_point(entry):
    """Full-precision strings if the JSON has them (re_str/im_str), floats otherwise."""
    if "re_str" in entry and "im_str" in entry:
        return mpmath.mpc(entry["re_str"], entry["im_str"]), "full precision"
    return mpmath.mpc(repr(entry["re"]), repr(entry["im"])), "double precision"


def check_point(label, expectation, alpha, omega, Re, note=""):
    """Evaluate D and domega/dalpha at ONE point and print them."""
    print("\n" + "-" * 72)
    print(label)
    print("-" * 72)
    print("alpha  =", alpha)
    print("omega  =", omega)
    print("expect :", expectation)
    if note:
        print("note   :", note)

    print("\nEvaluating original analytical dispersion relation...")
    D_value = fnum_python(alpha, omega, Re, PRECISION_DECIMAL_POINTS)[0]

    print("Evaluating original analytical derivative function...")
    domega_dalpha = fnum_derivatives_python(alpha, omega, Re, PRECISION_DECIMAL_POINTS)[0]

    # dD/domega, by central difference -- needed only to give |D| a scale
    h = mpmath.mpf(10) ** (-6)
    Dp = fnum_python(alpha, omega + h, Re, PRECISION_DECIMAL_POINTS)[0]
    Dm = fnum_python(alpha, omega - h, Re, PRECISION_DECIMAL_POINTS)[0]
    dD_domega = (Dp - Dm) / (2 * h)

    print("\nD(alpha, omega) =")
    print(D_value)
    print()
    print("|D| =")
    print(abs(D_value))

    print("\ndomega/dalpha =")
    print(domega_dalpha)
    print()
    print("|domega/dalpha| =")
    print(abs(domega_dalpha))

    # uncomment to report |D| as a distance to the surface, in omega
    # print("\n|dD/domega| =")
    # print(abs(dD_domega))
    # print("\n|D| / |dD/domega|   (distance to the surface, in omega) =")
    # print(abs(D_value) / abs(dD_domega))

    return D_value, domega_dalpha, dD_domega


with open(PINCH_JSON, "r") as f:
    pinch = json.load(f)

# both compute_pinch.py and compute_pinch2.py put the final point last
result = pinch["level_results"][-1]

# set the precision BEFORE parsing -- mpmath.mpc() parses at the current mp.dps
mpmath.mp.dps = PRECISION_DECIMAL_POINTS

Re = mpmath.mpf(repr(pinch["settings"]["Re"]))
alpha_p, source = read_point(result["alpha_p"])
omega_p, _ = read_point(result["omega_p"])

# the seed -- the point Newton STARTED from
alpha_seed, seed_source = read_point(pinch["seed"])

if "omega_seed" in pinch:
    omega_seed, _ = read_point(pinch["omega_seed"])
    seed_note = "omega seed read from the JSON"
else:
    omega_seed = omega_p
    seed_note = ("no omega_seed in the JSON -- omega_p substituted, so this "
                 "row varies alpha only")

print("=" * 72)
print("ANALYTIC PINCH CHECK - ORIGINAL FUNCTIONS, EXACT-THIRDS FIX APPLIED")
print("=" * 72)

print("\nfile        :", PINCH_JSON)
print("Re          =", Re)
print("dps         =", PRECISION_DECIMAL_POINTS)
print("pinch read  :", source)
print("seed read   :", seed_source)
print("seed is", mpmath.nstr(abs(alpha_seed - alpha_p), 6),
      "away from the pinch in alpha")

print("\nRESULTS")
# floors at dps 60: |D| ~ 1.9e-33, |domega/dalpha| ~ 8e-47; anything printed
# below those is noise.  They apply to the PINCH row only.

check_point("SEED  -- the point Newton started from",
            "both LARGE; this is not the pinch",
            alpha_seed, omega_seed, Re, seed_note)

D_value, domega_dalpha, dD_domega = check_point(
    "PINCH -- the converged point",
    "both ~ 0, down to the arithmetic floor above",
    alpha_p, omega_p, Re)
