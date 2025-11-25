# ------------------------------------------------------------------------------------------
#                            SHIFT AND INJECTION PATHS
# ------------------------------------------------------------------------------------------
shift_current(a, b, c, p::ParamsHF, ωlist, part; kws...) =
    nonlinear_terms(a, b, c, p, ωlist, :SHIFT, part; kws...)   

shift_current(a, b, c, p::ParamsHF, of, ωlist, part; kws...) =
    nonlinear_terms(a, b, c, p, of, ωlist, :SHIFT, part; kws...)   

injection_current(a, b, c, p::ParamsHF, ωlist, part; kws...) =
    nonlinear_terms(a, b, c, p, ωlist, :INJECTION, part; kws...)    

injection_current(a, b, c, p::ParamsHF, of, ωlist, part; kws...) =
    nonlinear_terms(a, b, c, p, of, ωlist, :INJECTION, part; kws...)    
# ---------------------------------------------------------------------------------------------
#                                   INTEGRATION
# ---------------------------------------------------------------------------------------------
function nonlinear_terms(a, b, c, p::ParamsHF, ωlist, current_type::Symbol, part::Symbol; 
        η = 0.1, evals = 100, kws...)
    function _integral_nonlinear!(v, ν; kws...)
        v .= integral_nonlinear(ωlist, a, b, c, ParamsHF(p, ν = ν), η, evals, current_type, part; kws...)
        println("$(ν) valley")
    end
    println("Number of iterations: ", evals, "| η = ", η)
    vals = zeros(Float64, length(ωlist))
    vals2 = similar(vals)
    _integral_nonlinear!(vals, 1; kws...)
    _integral_nonlinear!(vals2, -1, kws...)
    return ωlist, vals, vals2
 end 

 function nonlinear_terms(a, b, c, p::ParamsHF, of, ωlist, current_type::Symbol, part::Symbol; 
    η = 0.1, evals = 100, kws...)
    println("Number of iterations: ", evals, "| η = ", η)
    vals = integral_nonlinear(ωlist, a, b, c, p, of, η, evals, current_type, part; kws...)
    vals2 = zeros(Float64, length(ωlist))
    return ωlist, vals, vals2
end 

function integral_nonlinear(ωlist::Array, a, b, c, p::ParamsHF, η, evals, current_type, part; kws...)
    M, xmin, xmax = int_boundaries(p)
    check_current(current_type); check_part(part)
    integral_current = (current_type == :SHIFT ? shift_current_ω : injection_current_ω) #CD is not implemented, yet
    integrand(q) = real(integral_current(part, ωlist, a, b, c, q, p, η))
    angstroms_to_nm = 1/10; spin_dof = 2
    cnst = spin_dof * C * ħ_ev_s * angstroms_to_nm
    return bz_integration(integrand, p, ωlist, evals) * cnst # the prefactor is the jacobian     
end

function integral_nonlinear(ωlist::Array, a, b, c, p::ParamsHF, of, η, evals, current_type, part; kws...)
    M, xmin, xmax = int_boundaries(p)
    check_current(current_type); check_part(part)
    integral_current = (current_type == :SHIFT ? shift_current_ω : injection_current_ω) #CD is not implemented, yet
    integrand(q) = real(integral_current(part, ωlist, a, b, c, q, p, of, η))
    angstroms_to_nm = 1/10; spin_dof = 2
    cnst = spin_dof * C * ħ_ev_s * angstroms_to_nm
    return bz_integration(integrand, p, ωlist, evals) * cnst # the prefactor is the jacobian     
end
# ---------------------------------------------------------------------------------------------
#                                 AUX FUNCTIONS. TESTS
# ---------------------------------------------------------------------------------------------
function check_current(current_type)
    try
        if current_type == :SHIFT || current_type == :INJECTION
            nothing
        else
            throw(ArgumentError("Variable is neither :SHIFT nor :INJECTION"))
        end
    catch e
        return "Error: $(e)"
    end
end

function check_part(part)
    try
        if part == :REAL || part == :IMAG
            nothing 
        else
            throw(ArgumentError("Variable is neither :REAL nor :IMAG"))
        end
    catch e
        return "Error: $(e)"
    end
end
