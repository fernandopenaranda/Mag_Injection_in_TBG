include("non_linear_responses.jl")

# ------------------------------------------------------------------------------------------
# σab_inter_linear
# ------------------------------------------------------------------------------------------

"""
    `σab_inter_linear(a, b, p, ωlist; save = false, η = 10^-3, evals = 6400) `
returns the interband linear optical conductivity for the BM model using:
`σab^inter = -C sum_{m!=n} ωnm (rnm^a rmn^b/(ωmn -ω + iη) * fnm)`
 - `σab(ω)` defined as `J_a = σab(ω) J_b` is computed over a list of frequencies
`ω` in `ωlist`. 
 -`a` and `b` are the in-plane (:x, :y) directions.
 -`p::ParamsBM` is a struct that contains the system information.  
 - `η` defines the energy broadening of the Lorentzian peaks.
 - `evals` sets the mesh size in the BZ integration
 - `save =  true` to automatically save the output
 - `method = :fast` (default) efficently evaluates the integral computing just once the mesh 
 for all frequencies. `method = :slow` is slightly more stable at similar `evals` although
 it is much slower.

"""
function σab_inter_linear(a, b, p, ωlist; grid = :none, monolayer = :false, part = :real, kws...)
    if grid == :none
        if monolayer == true
            println("Calculation for a single graphene monolayer")
            if  return_kwargs(kws, :mass_term) != :layered
                throw(ArgumentError("Mass term must be :layered" ))
            else nothing end
            p = ParamsBM(p, tAA = 0., t3 = 0., mass = 1000) 
        else nothing end    
        return σab_inter_linear_fast(a, b, p, ωlist, part; kws...)
    else 
        return σab_inter_linear_uniformgrid_mat(a, b, p, ωlist, part; kws...)
    end
end



function σab_inter_linear_fast(a, b, p, ωlist, part; save = false, η = 10^-3,
        evals = 10000, kws...)
    vals = zeros(Float64, length(ωlist))
    errs = similar(vals)
    vals, errs = integral_linear_substitution_fast(ωlist, a, b, p, η, evals, part; kws...)
    if isa(ωlist, Array)
        ifelse(save == true , 
            saveconductivity(p, ωlist, linearunits(vals), errs, η, evals, a, b, part;
                kws...), nothing)
    else nothing end
    return vals
 end

 """
 In this function I compute the integral of the linear optical conductivity
 over the Brilluoin zone (BZ), to do so I rotate the axis by a change of variables
 x = v + u + M[1]; y = u - v.
 It performs an adaptive integral in k space doubling the k mesh by adding more points in the
 fast changing regions of the BZ. This is done untill prescribed tolerance is achieved. 
 It uses hcubature techniques.
 """
function integral_linear_substitution_fast(ωlist::Array, a, b, p, η, evals, part; kws...)
    println("Array")
    println(evals)
    M, xmin, xmax = int_boundaries(p)
    real_or_imag = ifelse(part == :real, real, imag)
    dha, dhb = dhs(a, b, p; kw...)
    integrand_fast(q) = 
        real_or_imag(σab_inter_linear_ω_fast(ωlist, dha, dhb, q, p, η; kws...))
    val, err = hcubature(length(ωlist), (x,v) -> v[:] = 
        integrand_fast(x), xmin, [xmax[1]/2, xmax[2]];
        error_norm = Cubature.INDIVIDUAL, reltol = 1e-5, abstol = 0, maxevals = evals) 
    bz_surface  =  (1/(2pi*p.a0))^2 
    return bz_surface * val, err 
end

 function linearunits(vals)    
    σmono = 1/4 # units of e^2/hbar resolved in spin and valley. Each cone contributes with a 1/16 e^2/hbar
    valleys = 2
    spins = 2
    return vals * valleys * spins/ σmono
 end



"""
    `integral_linear_substitution_fast(ω::Float64, dha::SparseMatrixCSC, dhb::SparseMatrixCSC,
        p, η, evals, part; method = :Cubature, kws...)`
computes the integral at frequency ω. 
    Two adaptive integration methods are implemented:
    method == :Cubature (default) uses adaptive quadrature techniques in 2d
    method == :MonteCarlo uses adaptive MonteCarlo techniques based on VegaMC algorithm
    method == :Splines uses a spline interpolation for the integral

"""
function integral_linear_substitution_fast(ω::Float64, dha::SparseMatrixCSC, dhb::SparseMatrixCSC,
    p, η, evals, part; method = :Cubature, kws...)        
    M, xmin, xmax = int_boundaries(p)
    real_or_imag = ifelse(part == :real, real, imag)
    integrand_fast(q) = real_or_imag(σab_inter_linear_ω_fast(ω, dha, dhb, q, p, η; kws...))
    # method with hcubature    
    if method == :Cubature
        val, err = hcubature(integrand_fast, xmin, [xmax[1]/2, xmax[2]], reltol = 1e-4, abstol = 0, maxevals = evals) 
    elseif method == :MonteCarlo 
    # method with MCIntegration adaptive vegas
        xy = Continuous([(0.0, 1.0), (0.0, 1.0)])
        res = integrate(((x, y), c)-> integrand_fast([x[1],y[1]]); var = xy, solver = :vegasmc, neval = evals)
        println(res)    
        val= res[1][1]
        err = res[1][2]
    else 
        nothing
    end
    bz_surface  =  (1/(2*pi*p.a0))^2
    return bz_surface * val, err
end

function integral_linear_substitution_fast(ω::Float64, a::Symbol, b::Symbol, p, η, evals, part; method = :Cubature, kws...)        
    M, xmin, xmax = int_boundaries(p)
    real_or_imag = ifelse(part == :real, real, imag)
    dha, dhb = dhs(a, b, p; kw...)
    integrand_fast(q) = real_or_imag(σab_inter_linear_ω_fast(ω, dha, dhb, q, p, η; kws...))
    # method with hcubature    
    if method == :Cubature
        val, err = hcubature(integrand_fast, xmin, [xmax[1]/2, xmax[2]], reltol = 1e-8, abstol = 0, maxevals = evals) 
    else 
    # method with MCIntegration adaptive vegas
        xy = Continuous([(0.0, 1.0), (0.0, 1.0)])
        res = integrate(((x, y), c)-> integrand_fast([x[1],y[1]]); var = xy, solver = :vegasmc, neval = evals)
        println(res)    
        val= res[1][1]
        err = res[1][2]
    end
    bz_surface = (1/(2*pi*p.a0))^2
    return bz_surface * val, err # In units of sigma_mono = e^2/(4hbar)        
end

function dhs(a, b, p; kw...)                                                                       
    unbounded = (:x, :y)
    if !(a ∈ unbounded && b ∈ unbounded)
        throw(ArgumentError(
            "Only linear photocurrents with unbounded cartesian indices are implemented"))
    end
    dha = d_bistritzer_hamiltonian(p, a; kw...)
    if a == b
        return dha, dha
    else
        dhb = d_bistritzer_hamiltonian(p, b; kw...)
        return dha, dhb
    end
end

function σab_inter_linear_ω_fast(ωlist::Array, dha, dhb, q, p, η; kws...)
    h = bistritzer_hamiltonian(p, q; kws...)
    ϵs, ψs = eigen(Matrix(h))
    mat = linear_integrand_fast(ϵs, ψs, p, dha, dhb)
    return [π .* sum_nondiag(mat .* lorentz(ϵs, ω, η) .* -ω) for ω in ωlist] 
end

function σab_inter_linear_ω_fast(ω::Float64, dha, dhb, q, p, η; kws...)
    h = bistritzer_hamiltonian(p, q; kws...)
    ϵs, ψs = eigen(Matrix(h))
    mat = linear_integrand_fast(ϵs, ψs, p, dha, dhb)
    return π .* sum_nondiag(mat .* lorentz(ϵs, ω, η) .* -ω)
end

function linear_integrand_fast(ϵs, ψs, p, dha, dhb)
    unbounded = (:x, :y)
    dim = length(ϵs)
    return f(ϵs, p.mu) .* (r(ϵs, ψs, dha) .* t(r(ϵs, ψs, dhb))) 
end

