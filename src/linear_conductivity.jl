# ------------------------------------------------------------------------------------------
# σab_inter_linear
# ------------------------------------------------------------------------------------------
""" ωlist in meV 
ωlist, conds = σab_inter_linear(:x, :x, p, collect(0:1:20), evals = 4000);
    plot_linear_conductivity(ωlist, conds, :x,:x, part = :real) """
function σab_inter_linear(a, b, p::ParamsHF, ωlist; save = false, η = 0.1, evals = 10000, part = :real, kws...)
    conds = zeros(Float64, length(ωlist))
    println(evals)
    half_dim = length(ωlist)÷2
    conds[1:half_dim] .= integral_linear(ωlist[1:half_dim], a, b, p, η, evals, part)
    conds[half_dim+1:length(ωlist)] .= integral_linear(ωlist[half_dim+1:length(ωlist)], a, b, p, η, evals, part)
    return ωlist, conds
end

function σab_inter_linear_ω(ωlist::Array, dha, dhb, q, p, η)
    h = hf_hamiltonian(p::ParamsHF, q)
    ϵs, ψs = eigen(Matrix(h))
    mat = linear_integrand(ϵs, ψs, p, dha, dhb)
    return [π .* sum_nondiag(mat .* lorentz(ϵs, ω, η) .* -ω) for ω in ωlist] 
end

linear_integrand(ϵs, ψs, p, dha, dhb) = 
    f(ϵs, p.μ, 0.0) .* (r(ϵs, ψs, dha) .* t(r(ϵs, ψs, dhb)))

function integral_linear(ωlist::Array, a, b, p, η, evals, part)
    # println("Evals: ", evals)
    M, xmin, xmax = int_boundaries(p)
    real_or_imag = ifelse(part == :real, real, imag)
    dha, dhb = dhs(a, b, p)
    integrand(q) = real_or_imag(σab_inter_linear_ω(ωlist, dha, dhb, q, p, η))
    # integrand(q) = hf_jdos_ω(ωlist, p, q, η)
    return bz_integration(integrand, p, ωlist, evals)
end

# function bz_integration(f, p, ωlist, evals) #BUENA UNCOMMENT COMO EN EL PAPER
#     M, xmin, xmax = int_boundaries(p)
#     val, err = hcubature(length(ωlist), (x,v) -> v[:] = f(x), 
#         [-xmax[1]/4, xmin[2]], [xmax[1]/2 - xmax[1]/4, xmax[2]]; reltol = 1e-5, abstol=0, maxevals=evals);
#     bz_surface  =  (1/(2pi*a0))^2 
#     return bz_surface .* val
# end

# function bz_integration(f, p, evals)
#     M, xmin, xmax = int_boundaries(p)
#     val, err = hcubature(length(ωlist), (x,v) -> v[:] = f(x), 
#         [0, xmin[2]], [xmax[1]/2, xmax[2]]; reltol = 1e-5, abstol=0, maxevals=evals);
#     bz_surface  =  (1/(2pi*a0))^2 
#     return bz_surface .* val
# end

function bz_integration(f, p, ωlist, evals) 
    M, xmin, xmax = int_boundaries(p)
    val, err = Cubature.hcubature(length(ωlist), (x,v) -> v[:] = f(x), 
        [0, xmin[2]], [xmax[1]/2, xmax[2]]; reltol = 1e-5, abstol=0, maxevals=evals);
    bz_surface  =  (1/(2pi*a0))^2 
    return bz_surface .* val
end


function bz_grid_integration(f, p, ωlist, evals)
    M, xmin, xmax = int_boundaries(p)
    xm = -xmax[1]/4 
    ym =  xmin[2] 
    xp = xmax[1]/2 - xmax[1]/4 
    yp =  xmax[2]
    vals = zeros(length(ωlist))
    for o in 1:length( ωlist )
        for i in 0:evals
            for j in 0:evals
                qx = xm + i/evals * (xp-xm)
                qy = ym + i/evals * (yp-ym)
                vals .+= f([qx,qy], ωlist)
                println(vals)
                if isnan(vals[2])
                    G1, G2 = bravais_vectors(p)
                    K1 =  κ(G1, G2)
                    println([qx,qy], " vals ", vals)
                else nothing end
            end
        end
    end
    bz_surface  =  (1/(2pi*a0))^2 
    return bz_surface .* vals
end


# ----------------------------------------------------------------------
#                               JDOS
# ----------------------------------------------------------------------

function hf_jdos(p::ParamsHF, ωlist; η = 0.1, evals = 10000, kws...)
    jdos = zeros(Float64, length(ωlist))
    half_dim = length(ωlist)÷2
    jdos[1:half_dim] .= integral_jdos(ωlist[1:half_dim], p, η, evals)
    jdos[half_dim+1:length(ωlist)] .= integral_jdos(ωlist[half_dim+1:length(ωlist)], p, η, evals)
    return ωlist, jdos
end

function integral_jdos(ωlist::Array, p, η, evals)
    integrand(q) = hf_jdos_ω(ωlist, p, q, η)
    return bz_integration(integrand, p, ωlist, evals)
end

function hf_jdos_ω(ωlist::Array, p::ParamsHF, q, η)
    h = hf_hamiltonian(p::ParamsHF, q)
    ϵs, ψs = eigen(Matrix(h))
    return [sum_nondiag(-f(ϵs, 0, 0) .* lorentz(ϵs, ω, η)) for ω in ωlist]   
end

function linearunits(vals)    
    σmono = 1/4 # units of e^2/hbar resolved in spin and valley. Each cone contributes with a 1/16 e^2/hbar
    valleys = 2
    spins = 2
    return vals * valleys * spins/ σmono
end

flatten(v) = collect(Base.Iterators.flatten(v))

#=

DOCUMENTER
returns the interband linear optical conductivity for the BM model using:
"""
`σab_inter_linear(a, b, p, ωlist; save = false, η = 10^-3, evals = 6400) `
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

"""
In this function I compute the integral of the linear optical conductivity
over the Brilluoin zone (BZ), to do so I rotate the axis by a change of variables
x = v + u + M[1]; y = u - v.
It performs an adaptive integral in k space doubling the k mesh by adding more points in the
fast changing regions of the BZ. This is done untill prescribed tolerance is achieved. 
It uses hcubature techniques.
"""
=#

