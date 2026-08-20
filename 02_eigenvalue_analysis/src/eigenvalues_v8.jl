# Kilian Vinzenz Wilhelm

## Packages
# Parallel computing, managing resources
using Pkg

Pkg.update()
#Pkg.add("PyCall")
#Pkg.add("Distributed")
#Pkg.add("LinearAlgebra")
#Pkg.add("Plots")
Pkg.add("JSON")
using JSON
using Base.Filesystem
using Distributed
addprocs(7)
w = workers()
#@everywhere using Pkg
@everywhere using PyCall, Distributed, LinearAlgebra
masterpath::String = "\\\\C:\\Users\\gnrup\\OneDrive\\Desktop\\hiwi\\FDY"

@everywhere begin
    ## Define precision levels
    global precision_decimal_points_pos::Int64 = 20
    global precision_decimal_points_neg::Int64 = 150
    global precision_decimal_points_bigfloat::BigFloat = 0.0  # Initialize the variable
    ## Precision
    # ____________________
    # BigFloat:
    # precision_decimal_points::Int64 = 8
    setprecision(BigFloat, precision_decimal_points_neg; base=10)
    # precision_decimal_points_bigfloat::BigFloat = big(precision_decimal_points)
    # in coordination with the precision of Python of mp.dps = 128
    # ____________________
    ## Input: Parameters
    # ____________________ #
    Re::BigFloat = big"2800.0"   # Reynolds number
    omega::BigFloat = big"0.25" # frequency
    beta::BigFloat = big"0.0"    # spanwise wavenumber, for the 2D case always 0
    # ____________________ #
    ## Differentiation matrices (former Dmat.m)
    # ____________________ #
    num_modes::Int64 = 150   # number of collocation points --> result in 2*num_modes modes in EV spectrum
    start::Int64 = 0
    terminate::Int64 = 1
    # ____________________ #
    function diff_mat(num_modes, start, terminate)
        num_modes::Int64 = round(abs(num_modes))
        y_colloc_points::Vector{BigFloat} = [cos((i - 1) * big(pi) / (num_modes - 1)) for i = 1:num_modes]
        global y_colloc_points_new = ((start + terminate) / 2) .- y_colloc_points * ((terminate - start) / 2)
        D0_static::Matrix{BigFloat} = zeros(BigFloat, num_modes, num_modes) # manual rotation
        for i = 1:num_modes 
            D0_static[:, i] .= cos.((i - 1) * acos.(y_colloc_points))
        end
        D1_static::Matrix{BigFloat} = [zeros(num_modes, 1)    D0_static[:, 1]         4 * D0_static[:, 2]]
        D2_static::Matrix{BigFloat} = [zeros(num_modes, 1)    zeros(num_modes, 1)     4 * D0_static[:, 1]]
        D3_static::Matrix{BigFloat} = [zeros(num_modes, 1)    zeros(num_modes, 1)     zeros(num_modes, 1)]
        D4_static::Matrix{BigFloat} = [zeros(num_modes, 1)    zeros(num_modes, 1)     zeros(num_modes, 1)]
        D1_static_V2::Matrix{BigFloat} = zeros(BigFloat, num_modes, num_modes)
        D2_static_V2::Matrix{BigFloat} = zeros(BigFloat, num_modes, num_modes)
        D3_static_V2::Matrix{BigFloat} = zeros(BigFloat, num_modes, num_modes)
        D4_static_V2::Matrix{BigFloat} = zeros(BigFloat, num_modes, num_modes)
        D0_static_V2::Matrix{BigFloat} = D0_static
        D1_static_V2[:, 1:3] .= D1_static
        D2_static_V2[:, 1:3] .= D2_static
        D3_static_V2[:, 1:3] .= D3_static
        D4_static_V2[:, 1:3] .= D4_static
        for i = 4:num_modes
            D1_static_V2[:, i] .= 2 * (i - 1) * D0_static_V2[:, i - 1] + (i - 1) * D1_static_V2[:, i - 2] / (i - 3)   
            D2_static_V2[:, i] .= 2 * (i - 1) * D1_static_V2[:, i - 1] + (i - 1) * D2_static_V2[:, i - 2] / (i - 3)
            D3_static_V2[:, i] .= 2 * (i - 1) * D2_static_V2[:, i - 1] + (i - 1) * D3_static_V2[:, i - 2] / (i - 3)
            D4_static_V2[:, i] .= 2 * (i - 1) * D3_static_V2[:, i - 1] + (i - 1) * D4_static_V2[:, i - 2] / (i - 3)
        end
        D1_static_V3::Matrix{BigFloat} = D1_static_V2 / (-(terminate - start) / 2)^1
        D2_static_V3::Matrix{BigFloat} = D2_static_V2 / (-(terminate - start) / 2)^2
        D3_static_V3::Matrix{BigFloat} = D3_static_V2 / (-(terminate - start) / 2)^3
        D4_static_V3::Matrix{BigFloat} = D4_static_V2 / (-(terminate - start) / 2)^4
        global D0 = D0_static
        global D1 = D1_static_V3
        global D2 = D2_static_V3
        global D3 = D3_static_V3
        global D4 = D4_static_V3
    end
    diff_mat(num_modes, start, terminate)

    ## Flow type
    # Couette flow OS equation
    function couetteflow(num_modes, y_colloc_points_new, omega, beta, Re, D0, D1, D2, D3, D4)   
        u::Vector{BigFloat} = y_colloc_points_new # mean velocity
        d2u::BigFloat = 0
        A11::Matrix{Complex{BigFloat}} = 2 * im * omega * D1 + 4 / Re * D3 - 4 / Re * beta^2 * D1 + im * (u * ones(Complex{BigFloat}, 1, length(u))) .* D2 - im * beta^2 * (u * ones(1, length(u))) .* D0 - im * (d2u * ones(1, length(u))) .* D0
        A12::Matrix{Complex{BigFloat}} = -im * omega * D2 + im * omega * beta^2 * D0 - 1 / Re * D4 + 2 / Re * beta^2 * D2 - 1 / Re * beta^4 * D0
        A11 .= [zeros(Complex{BigFloat}, 2, num_modes);                    A11[3:num_modes - 2, :];    zeros(Complex{BigFloat}, 2, num_modes)]
        A12 .= [-200 * im * [D0[1:1, :]; D1[1:1, :]];   A12[3:num_modes - 2, :];    -200 * im * [D1[num_modes:num_modes, :]; D0[num_modes:num_modes, :]]]
        A21::Matrix{Complex{BigFloat}} = 1 * Matrix{Complex{BigFloat}}(I, num_modes, num_modes)
        A22::Matrix{Complex{BigFloat}} = zeros(Complex{BigFloat}, num_modes, num_modes)  
        global A = [A11 A12; A21 A22]   
        B11::Matrix{Complex{BigFloat}} = 4 / Re * D2 + 2 * im * (u * ones(Complex{BigFloat}, 1, length(u))) .* D1
        B11 = [zeros(Complex{BigFloat}, 2, num_modes);  B11[3:num_modes - 2,:];    zeros(Complex{BigFloat}, 2, num_modes);]
        B12::Matrix{Complex{BigFloat}} = zeros(Complex{BigFloat}, num_modes, num_modes)
        B21::Matrix{Complex{BigFloat}} = zeros(Complex{BigFloat}, num_modes, num_modes)
        B22::Matrix{Complex{BigFloat}} = 1 * Matrix{Complex{BigFloat}}(I, num_modes, num_modes)
        global B = [B11 B12; B21 B22] 
    end
    # ____________________ #
    couetteflow(num_modes, y_colloc_points_new, omega, beta, Re, D0, D1, D2, D3, D4)
    # ____________________ #

    ## Eigenvalue spectrum
    # Eigenvalues via Julia
    A_normal = convert(Matrix{Complex{Float64}}, A)
    B_normal = convert(Matrix{Complex{Float64}}, B)
    lambda_AB::Vector{Complex{BigFloat}} = big(1) .* eigen(A_normal, B_normal).values
    lambda_BA::Vector{Complex{BigFloat}} = big(1) ./ eigen(B_normal, A_normal).values
    #_________________
    #data = lambda_AB
    data::Vector{Complex{BigFloat}} = lambda_BA
    
    # ____________________
    ## Dispersion relation
    # python function:
    py"""
    import mpmath

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
        while attempts <= 5:
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
    
    """
    # pre-allocation function Fnum
    # direct python variables are not pre-allocated due to problems
    disprel::Complex{BigFloat} = 0.0 + 0.0 * im
    error_ipa0::BigFloat = 0.0
    error_ina2::BigFloat = 0.0
    error_ina0::BigFloat = 0.0
    error_ipa2::BigFloat = 0.0
    max_error_integral::BigFloat = 0.0
    attempts::Int64 = 1
    max_attempts::Int64 = 100
    #___________________

    # Function to set precision based on the sign of the imaginary part
    function set_precision_imag(a)
        if imag(a) > 0.0
            #setprecision(BigFloat, precision_decimal_points_pos; base=10)
            global precision_decimal_points_bigfloat = big(precision_decimal_points_pos)
        elseif imag(a) < 0.0
            #setprecision(BigFloat, precision_decimal_points_neg; base=10)
            global precision_decimal_points_bigfloat = big(precision_decimal_points_neg)
        end
    end
    #_____________________
    
    # function Fnum in julia, precision for positive and negative imaginary can be set here.
    function Fnum(a, omega)
        max_attempts::Int64 = 100
        attempts = 1
        while attempts < max_attempts
            try
                set_precision_imag(a) 
                # Call the Python function with the appropriate precision
                @inline disprel_py, error_ipa0_py, error_ina2_py, error_ipa2_py, error_ina0_py = py"fnum_python"(a, omega, Re, precision_decimal_points_bigfloat)
                # Convert the results from Python back to Julia
                disprel = disprel_py
                error_ipa0 = error_ipa0_py
                error_ina2 = error_ina2_py
                error_ipa2 = error_ipa2_py
                error_ina0 = error_ina0_py
                max_error_integral = max(abs(error_ipa0), abs(error_ina2), abs(error_ipa2), abs(error_ina0))
                return disprel, max_error_integral
            catch ex
                println("An error occurred when calling the Python function:")
                println(ex)
                a += big"1.0" * (big"10.0"^(-precision_decimal_points_bigfloat + big"5.0") + big"10.0"^(-precision_decimal_points_bigfloat + big"5.0") * im)
                attempts += 1
            end    
        end
    
        disprel = NaN
        error_ipa0 = NaN
        error_ina2 = NaN
        error_ipa2 = NaN
        error_ina0 = NaN
        max_error_integral = max(abs(error_ipa0), abs(error_ina2), abs(error_ipa2), abs(error_ina0))
        return disprel, max_error_integral
    end
    # pre-filtering for iteration using Muller's method / generalized secant method
    global check_length_data::Int64 = length(data)
    global i::Int64 = 1
    while i <= check_length_data
        # ____________________
        if isnan(real(data[i])) == true|| isnan(imag(data[i])) == true || isinf(abs(real(data[i]))) == true  || isinf(abs(imag(data[i]))) == true || abs(real(data[i])) > 5 || abs(imag(data[i])) > 20
        # ____________________
            deleteat!(data, i)
            global check_length_data = length(data) 
            global i -= 1
        end
        global i += 1
    end
end # everywhere end


# initial residual array
ires::Vector{Complex{BigFloat}} = zeros(Complex{BigFloat}, length(data))
max_error_ires::Vector{BigFloat} = zeros(BigFloat, length(data))
results = pmap(x -> Fnum(x, omega), data)
for i in eachindex(data)
    ires[i], max_error_ires[i] = results[i]
    println("$i -th mode, initial residuum: $(ires[i])")
end

#______________________

## Iteration using Muller's method / Generalized Secant Method
@everywhere begin
    # Function to get error threshold based on the sign of the imaginary part
    function get_error_threshold(a)
        if imag(a) > 0.0
            return big(10.0^(-precision_decimal_points_pos))
        elseif imag(a) < 0.0
            return big(10.0^(-precision_decimal_points_neg))
        else
            return big(10.0^(-precision_decimal_points_neg))  # Default error threshold
        end
    end
    # termination condition
    # Adjust error threshold based on the sign of the imaginary part
    
    #e::BigFloat = big(10.0^(-precision_decimal_points_neg + 5))
    # 3 point factor iterative: factor method
    cn::BigFloat = 1.0 - 0.0 * big(10.0^(-precision_decimal_points_pos + 2))
    cn_m1::BigFloat = 1.0 - 1.0 * big(10.0^(-precision_decimal_points_pos + 2))
    cn_m2::BigFloat = 1.0 - 2.0 * big(10.0^(-precision_decimal_points_pos + 2))
    # pre-allocation 
    global data_cond_skip::Vector{Bool} = fill(false, (length(data))) # in Matlab data(ii, 2)
    data_cond_skip_individual_vector::Vector{Bool} = fill(false, (length(data)))
    # data_cond_skip_individual_value_future is not pre-allocated due to problems
    data_cond_skip_individual_value_fetch::Complex{BigFloat} = 0.0 + 0.0 * im
    global data_disprel_it_mode::Vector{Complex{BigFloat}} = zeros(Complex{BigFloat}, length(data)) # in Matlab data(ii, 3)
    data_disprel_it_mode_individual_vector::Vector{Complex{BigFloat}} = zeros(Complex{BigFloat}, length(data))
    # data_disprel_it_mode_individual_value_future is not pre-allocated due to problems
    data_disprel_it_mode_individual_value_fetch::Complex{BigFloat} = 0.0 + 0.0 * im
    data_disprel_ini_input::Vector{Complex{BigFloat}} = zeros(Complex{BigFloat}, length(data)) # in Matlab data(ii, 4)
    data_disprel_it_input::Vector{Complex{BigFloat}} = zeros(Complex{BigFloat}, length(data)) # in in Matlab data(ii, 5)
    # pre-allocation auxiliary variables
    x0::Complex{BigFloat} = 0.0 + 0.0 * im
    xn_m2::Complex{BigFloat} = 0.0 + 0.0 * im
    xn_m1::Complex{BigFloat} = 0.0 + 0.0 * im
    xn::Complex{BigFloat} = 0.0 + 0.0 * im
    it::Int64 = 0
    xn_lit::Complex{BigFloat} = 0.0 + 0.0 * im
    res_lit::Complex{BigFloat} = 0.0 + 0.0 * im
    q_it::Complex{BigFloat} = 0.0 + 0.0 * im
    disprel_xn::Complex{BigFloat} = 0.0 + 0.0 * im
    disprel_xn_m1::Complex{BigFloat} = 0.0 + 0.0 * im
    disprel_xn_m2::Complex{BigFloat} = 0.0 + 0.0 * im
    A_it::Complex{BigFloat} = 0.0 + 0.0 * im
    B_it::Complex{BigFloat} = 0.0 + 0.0 * im
    C_it::Complex{BigFloat} = 0.0 + 0.0 * im
    xp::Complex{BigFloat} = 0.0 + 0.0 * im
    xm::Complex{BigFloat} = 0.0 + 0.0 * im
    ap::BigFloat = 0.0
    am::BigFloat = 0.0
    xn_p1::Complex{BigFloat} = 0.0 + 0.0 * im
    temp_res::Complex{BigFloat} = 0.0 + 0.0 * im
    temp_max_error_integral::BigFloat = 0.0
end
#
# loop over all entries




@everywhere begin
    function iteration(a, omega)
        ires, _ = Fnum(a, omega)
        e::BigFloat = 0.0
        e = get_error_threshold(a)
        # conditional steps if NaN in either real or imaginary part,
        # or if +-Inf in either real or imaginary part,
        # or exact zero in both real and imaginary part
        if isnan(real(ires)) == true || isnan(imag(ires)) == true || isinf(abs(real(ires))) == true || isinf(abs(imag(ires))) == true
            data_cond_skip_individual = true # in Matlab data(ii, 2)
            #continue
        else
            x0 = a
        end
        xn_m2 = cn_m2 * x0
        xn_m1 = cn_m1 * x0
        xn = cn * x0
        # iterator
        it = 3
        # pre-allocation
        xn_lit = 0.0 + 0.0 * im
        res_lit = 0.0 + 0.0 * im
        println("n \t xn \t \t \t f(xn)")
        println("0 \t $(xn_m2) \n $(Fnum(xn_m2, omega)[1])")
        println("1 \t $(xn_m1) \n $(Fnum(xn_m1, omega)[1])")
        println("2 \t $(xn) \n $(ires)")
        data_disprel_ini_input = ires
        # initializing to avoid function call in condition
        data_disprel_it_input = data_disprel_ini_input
        # iteration
        while abs(data_disprel_it_input) > e
            time_control = @time begin
                # 3-point term
                q_it = (xn - xn_m1) / (xn_m1 - xn_m2)
                # definitions: A_it, B_it, C_it; in Matlab A, B, C
                disprel_xn, _ = Fnum(xn, omega)
                disprel_xn_m1, _ = Fnum(xn_m1, omega)
                disprel_xn_m2, _ = Fnum(xn_m2, omega)
                A_it = q_it * disprel_xn - q_it * (1 + q_it) * disprel_xn_m1 + q_it^2 * disprel_xn_m2
                B_it = (2 * q_it + 1) * disprel_xn - (1 + q_it)^2 * disprel_xn_m1 + q_it^2 * disprel_xn_m2
                C_it = (1 + q_it) * disprel_xn
                # result comparison: two solutions (+ / -) due to square root
                xp = xn - (xn - xn_m1) * ((2 * C_it) / (B_it + sqrt(B_it^2 - 4 * A_it * C_it)))  # plus
                xm = xn - (xn - xn_m1) * ((2 * C_it) / (B_it - sqrt(B_it^2 - 4 * A_it * C_it)))  # minus
                # resulting residuals
                ap = abs(Fnum(xp, omega)[1]) # plus
                am = abs(Fnum(xm, omega)[1]) # minus
                # assignment for next iteration # x_(n+1)
                xn_p1 = 0.0 + 0.0 * im
                if am < ap  
                    xn_p1 = xm  
                else
                    xn_p1 = xp
                end
            end # time_control end
            # temporary residual value: assign (conditional requirement)
            temp_res, temp_max_error_integral = Fnum(xn_p1, omega)
            println("$it \t $xn_p1 \n $temp_res \n $temp_max_error_integral")
            # assignment for next iteration
            xn_m2 = xn_m1
            xn_m1 = xn
            xn = xn_p1
            # iterator
            it += 1
            #
            # identified conditionals
            #
            # disprel with iterated input
            data_disprel_it_input = temp_res
            # djlp = 0.025  # decimal jump limit parameter, not used anymore
            if isnan(real(xn)) == true || isnan(imag(xn)) == true || isnan(real(data_disprel_it_input)) == true || isnan(imag(data_disprel_it_input)) == true
                data_cond_skip_individual = true # in Matlab data(ii, 2)
                break
            end
            # ____________________
            if it >= 101
            # ____________________
            data_cond_skip_individual = true # in Matlab data(ii, 2)
                break
            end
            # output
            data_cond_skip_individual = false # in Matlab data(ii, 2)  # REITERABLE
            # assign current residual for conditional
            xn_lit = xn
            res_lit = data_disprel_it_input
        end # iteration end
        # iterated wavenumber mode 
        data_disprel_it_mode = xn # in Matlab data(ii, 3)
        # iteration count visualization
        println("$(it - 1) iterations")
        # numerical value comparison
        x0  # initial mode
        xn  # iterated mode
        f_x0 = data_disprel_ini_input # initial residual
        f_xn = data_disprel_it_input # iterated residual
        println("----------------------")
        return data_disprel_it_mode, data_cond_skip_individual # loop end
    end
end

data_disprel_it_mode::Vector{Complex{BigFloat}} = zeros(Complex{BigFloat}, length(data))
data_cond_skip::Vector{Bool} = fill(false, (length(data)))
results = pmap(x -> iteration(x, omega), data)
for i in eachindex(data)
    data_disprel_it_mode[i], data_cond_skip[i] = results[i]
    println("$i -th mode, iterated mode: $(data_disprel_it_mode[i])")
end
# Copy successfull modes
@everywhere begin
    data_disprel_it_mode_sucess::Vector{Complex{BigFloat}} = zeros(length(data_disprel_it_mode))
    # copy successfull modes
    for i in eachindex(data_disprel_it_mode_sucess)
        if data_cond_skip[i] == false || isnan(data_disprel_it_mode[i]) == true || isnan(data_disprel_it_mode[i]) == true
            data_disprel_it_mode_sucess[i] = data_disprel_it_mode[i]
        end
    end
    # filter out null entries
    filter!(x -> x != 0.0 + 0.0 * im, data_disprel_it_mode_sucess)
end


# Group velocity
@everywhere begin
    # python function
    py"""
    import mpmath


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
    
    """
    # pre-allocation function Fnum_derivatives
    # direct python variables are not pre-allocated due to problems
    d_omega_d_alpha::Complex{BigFloat} = 0.0 + 0.0 * im
    error_d_ipa0_d_omega::BigFloat = 0.0 
    error_d_ina2_d_omega::BigFloat = 0.0
    error_d_ipa2_d_omega::BigFloat = 0.0
    error_d_ina0_d_omega::BigFloat = 0.0
    error_d_ipa0_d_a::BigFloat = 0.0
    error_d_ina2_d_a::BigFloat = 0.0
    error_d_ipa2_d_a::BigFloat = 0.0
    error_d_ina0_d_a::BigFloat = 0.0
    # function Fnum in julia
    function Fnum_derivatives(a, omega)
        attempts = 1
        while attempts < max_attempts
            try
                @inline d_omega_d_alpha_py, error_ipa0_py, error_ina2_py, error_ipa2_py, error_ina0_py, error_d_ipa0_d_omega_py, error_d_ina2_d_omega_py, error_d_ipa2_d_omega_py, error_d_ina0_d_omega_py, error_d_ipa0_d_a_py, error_d_ina2_d_a_py, error_d_ipa2_d_a_py, error_d_ina0_d_a_py,  = py"fnum_derivatives_python"(a, omega, Re, precision_decimal_points_bigfloat)
                d_omega_d_alpha = d_omega_d_alpha_py
                error_ipa0 = error_ipa0_py
                error_ina2 = error_ina2_py
                error_ipa2 = error_ipa2_py
                error_ina0 = error_ina0_py
                error_d_ipa0_d_omega = error_d_ipa0_d_omega_py
                error_d_ina2_d_omega = error_d_ina2_d_omega_py
                error_d_ipa2_d_omega = error_d_ipa2_d_omega_py
                error_d_ina0_d_omega = error_d_ina0_d_omega_py
                error_d_ipa0_d_a = error_d_ipa0_d_a_py
                error_d_ina2_d_a = error_d_ina2_d_a_py
                error_d_ipa2_d_a = error_d_ipa2_d_a_py
                error_d_ina0_d_a = error_d_ina0_d_a_py
                max_error_integral = max(abs(error_ipa0), abs(error_ina2), abs(error_ipa2), abs(error_ina0), abs(error_d_ipa0_d_omega), abs(error_d_ina2_d_omega), abs(error_d_ipa2_d_omega), abs(error_d_ina0_d_omega), abs(error_d_ipa0_d_a), abs(error_d_ina2_d_a), abs(error_d_ipa2_d_a), abs(error_d_ina0_d_a))
                return d_omega_d_alpha, max_error_integral
            catch ex
                println("An error occurred when calling the Python function:")
                println(ex)
                a += big"1.0" * (big"10.0"^(-precision_decimal_points_bigfloat + big"5.0") + big"10.0"^(-precision_decimal_points_bigfloat + big"5.0") * im) 
                attempts += 1
            end    
        end
        d_omega_d_alpha = NaN
        error_ipa0 = NaN
        error_ina2 = NaN
        error_ipa2 = NaN
        error_ina0 = NaN
        error_d_ipa0_d_omega = NaN
        error_d_ina2_d_omega = NaN
        error_d_ipa2_d_omega = NaN
        error_d_ina0_d_omega = NaN
        error_d_ipa0_d_a = NaN
        error_d_ina2_d_a = NaN
        error_d_ipa2_d_a = NaN
        error_d_ina0_d_a = NaN
        max_error_integral = max(abs(error_ipa0), abs(error_ina2), abs(error_ipa2), abs(error_ina0), abs(error_d_ipa0_d_omega), abs(error_d_ina2_d_omega), abs(error_d_ipa2_d_omega), abs(error_d_ina0_d_omega), abs(error_d_ipa0_d_a), abs(error_d_ina2_d_a), abs(error_d_ipa2_d_a), abs(error_d_ina0_d_a))
        return d_omega_d_alpha, max_error_integral
    end 
end


# Group velocity
groupvelocity::Vector{Complex{BigFloat}} = zeros(Complex{BigFloat}, length(data_disprel_it_mode_sucess))
max_error_group::Vector{BigFloat} = zeros(BigFloat, length(data))
results = pmap(x -> Fnum_derivatives(x, omega), data_disprel_it_mode_sucess)
for i in eachindex(data_disprel_it_mode_sucess)
    groupvelocity[i], max_error_group[i] = results[i]
    println("$i -th mode, group velocity: $(groupvelocity[i])")
end



# Check condition on main
instab_condition_2d_value::Vector{BigFloat} = zeros(BigFloat, length(data_disprel_it_mode_sucess))
instab_condition_3d_value::Vector{BigFloat} = zeros(BigFloat, length(data_disprel_it_mode_sucess))
for i in eachindex(data_disprel_it_mode_sucess)
    # 2D
    instab_condition_2d_value[i] = real(groupvelocity)[i] * imag(data_disprel_it_mode_sucess)[i] 
    # 3D
    instab_condition_3d_value[i] = real(groupvelocity)[i] * imag(data_disprel_it_mode_sucess)[i] + imag(groupvelocity)[i] * real(data_disprel_it_mode_sucess)[i]
end
    


# output JSON
dict_ev_iterated::Vector{Dict} = fill(Dict(), length(data_disprel_it_mode_sucess))
for i in eachindex(data_disprel_it_mode_sucess)
    dict_ev_iterated[i] = Dict(
        "value" => data_disprel_it_mode_sucess[i],
        "group_velocity" => groupvelocity[i],
        "value_instab_2d" => instab_condition_2d_value[i],
        "value_instab_3d" => instab_condition_3d_value[i]
    )
end
dict_to_JSON = Dict(
    "reynolds" => Re,
    "omega" => omega,
    "ev_uniterated" => data,
    "ev_iterated" => dict_ev_iterated
)
# file locking
const LOCK_FILE_pre = "file.lock"
const LOCK_FILE = masterpath * LOCK_FILE_pre
function lock_file()
    open(LOCK_FILE, "w") do file
        write(file, "locked")
    end
end
function unlock_file()
    rm(LOCK_FILE)
end
function is_locked()
    return isfile(LOCK_FILE)
end
filename_pre = "OS_spatial_ev_master.json"
filename = masterpath * filename_pre
# Read the existing JSON file
current_array = []
while is_locked()
    sleep(1)  # Wait and retry if the file is locked
end
lock_file()  # Acquire the lock
if isfile(filename) && filesize(filename) > 0
    open(filename, "r") do file
        global json_str = read(file, String)
    end
    current_array = JSON.parse(json_str)
end
push!(current_array, dict_to_JSON)
json_str = JSON.json(current_array)
open(filename, "w") do file
    write(file, json_str)
end
unlock_file()


