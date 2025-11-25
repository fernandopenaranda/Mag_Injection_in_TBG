"    # this is the ∂tJ = ηijk thus 1/γ does not appear as a prefactorm
# Both Lorentzian should have the same size and the sum must go over n>m.
# Note however that this could give rise to finite conductivity at ω = 0 due to the
# broadening in the Lorentzians. "
function injection_current_ω(part, ωlist, a, b, c, q, p, η; T = 0) 
    h = hf_hamiltonian(p, q)
    ϵs, ψs = eigen(Matrix(h))
    mat = injection_current_integrand(part, ϵs, ψs, p, q, a, b, c, T, 0.0) 
    return [ħ_ev_s * π/2 * sum_lowerdiag((mat) .* (-lorentz(ϵs, ω, η) .- 0lorentz(ϵs, -ω, η)))  
        for ω in ωlist]  # Units: Å^3 /eV # check if you need a 1000 prefactor
end

function injection_current_ω(part, ωlist, a, b, c, q, p, of, η; T = 0)
    whichham = (p.twovalleystwospins == false ? hf_twovalleyshamiltonian : hf_valley_spin_hamiltonian)
    h = whichham(p, q, of) 
    ϵs, ψs = eigen(Matrix(h))
    mat = injection_current_integrand(part, ϵs, ψs, p, q, a, b, c, T, 0, of)
    return [ ħ_ev_s * π/2 * sum_lowerdiag(mat .* (-lorentz(ϵs, ω, η) .- lorentz(ϵs, -ω, η))) for ω in ωlist]  # Units: Å^3 /eV # # check if you need a 1000 prefactor
end

injection_current_integrand(part, ϵs, ψs, p, q, a, b, c, T, doping) = 
    f(ϵs, doping, T) .* Imat_injection(part, a, b, c, ϵs, ψs, p, q)

injection_current_integrand(part, ϵs, ψs, p, q, a, b, c, T, doping, of) =  
     f!(Imat_injection(part, a, b, c, ϵs, ψs, p, of, q), ϵs, doping, T)

""" returns the symmetrized product sum """
function Imat_injection(part, a, b, c, ϵs, ψs, p, q; current = :SHIFT)
    bounded_dirs = findall(x -> x == :z, (a, b, c))
    if !isempty(bounded_dirs)
        throw(ArgumentError("We cannot compute out out-of-plane currents since the z 
          direction is bounded but also since we lost teh z label we cannot compute oblique incidence"))
    else
        if p.λ == 0.0 #path with delta-like localization. slightly faster
            # println("Delta-like localization method")
            return Imat_injection_uuu(part, a, b, c, ϵs, ψs, p, current)
        else # path with exponentially localized wannier orbs
            # println("Exponential localization method")
        
            return Imat_injection_uuu(part, a, b, c, ϵs, ψs, p, q, current)
        end
    end
end


function Imat_injection(part, a, b, c, ϵs, ψs, p, of, q; current = :SHIFT)
    bounded_dirs = findall(x -> x == :z, (a, b, c))
    if !isempty(bounded_dirs)
        throw(ArgumentError("We cannot compute out out-of-plane currents since the z 
          direction is bounded but also since we lost teh z label we cannot compute oblique incidence"))
    else
        if p.λ == 0.0 #path with delta-like localization. slightly faster
            # println("Delta-like localization method")
            return Imat_injection_uuu(part, a, b, c, ϵs, ψs, p, of, current)
        else # path with exponentially localized wannier orbs
            # println("Exponential localization method")
            return Imat_injection_uuu(part, a, b, c, ϵs, ψs, p, of, q, current)
        end
    end
end

"""
 a, b, c correspond to unbounded directions
 """
function Imat_injection_uuu(part, a, b, c, ϵs, ψs, p, current)
    dha = dhf_hamiltonian(p, a)
    Δa = Δ(ψs, dha)
    inj_func = (part == :REAL ? inj_real : inj_imag)
    if a == b == c 
        ra = r(ϵs, ψs, dha)
        return inj_func(Δa, ra, ra)
    elseif b == c 
        rb = r(ϵs, ψs, dhf_hamiltonian(p, b)) 
        return inj_func(Δa, rb, rb)
    else
        rb = r(ϵs, ψs, dhf_hamiltonian(p, b))
        rc = r(ϵs, ψs, dhf_hamiltonian(p, c))
        return inj_func(Δa, rc, rb)
    end
end

function Imat_injection_uuu(part, a, b, c, ϵs, ψs, p, q::Array, current)
    dha = dhf_hamiltonian(p, q, a)
    Δa = Δ(ψs, dha)
    inj_func = (part == :REAL ? inj_real : inj_imag)
    if a == b == c 
        ra = r(ϵs, ψs, dha)
        return inj_func(Δa, ra, ra)
    elseif b == c 
        rb = r(ϵs, ψs, dhf_hamiltonian(p, q, b)) 
        return inj_func(Δa, rb, rb)
    else
        rb = r(ϵs, ψs, dhf_hamiltonian(p, q, b))
        rc = r(ϵs, ψs, dhf_hamiltonian(p, q, c))
        return inj_func(Δa, rc, rb)
    end
end


function Imat_injection_uuu(part, a, b, c, ϵs, ψs, p, of::SparseMatrixCSC, current)
    dha = dhf_hamiltonian(p, a, of)
    Δa = Δ(ψs, dha)
    inj_func = (part == :REAL ? inj_real : inj_imag)
    if a == b == c 
        ra = r(ϵs, ψs, dha)
        return inj_func(Δa, ra, ra)
    elseif b == c 
        rb = r(ϵs, ψs, dhf_hamiltonian(p, b, of)) 
        return inj_func(Δa, rb, rb)
    else
        rb = r(ϵs, ψs, dhf_hamiltonian(p, b, of))
        rc = r(ϵs, ψs, dhf_hamiltonian(p, c, of))
        return inj_func(Δa, rc, rb)
    end
end

function Imat_injection_uuu(part, a, b, c, ϵs, ψs, p, of, q, current) 
    dha = dhf_hamiltonian(p, q, a, of)
    Δa = Δ(ψs, dha)
    inj_func = (part == :REAL ? inj_real : inj_imag)
    if a == b == c 
        ra = r(ϵs, ψs, dha)
        return inj_func(Δa, ra, ra)
    elseif b == c 
        rb = r(ϵs, ψs, dhf_hamiltonian(p, q, b, of)) 
        return inj_func(Δa, rb, rb)
    else
        rb = r(ϵs, ψs, dhf_hamiltonian(p, q, b, of))
        rc = r(ϵs, ψs, dhf_hamiltonian(p, q, c, of))
        return inj_func(Δa, rc, rb)
    end
end


inj_real(Δ, r1, r2) = Δ .* t(real(r1 .* t(r2)))
inj_imag(Δ, r1, r2) = Δ .* t(imag.(r1 .* t(r2)))

#EXTRA


#########
# TWO VALLEY MODEL SPINFUL MODEL

##########



# vsfreq = true,
# if  vsfreq
# else # route to variation with μ
#     integrand_arr = spzeros(ComplexF64,length(list))
#     ω = 0.005 # CHANGE BY HAND
#     for (i, μ) in enumerate(list)
#         h = bistritzer_hamiltonian(ParamsBM(p, mu = μ), q; kws...)
#         ϵs, ψs = eigen(Matrix(h))
#         mat = injection_current_integrand_noTRS(ϵs, ψs, ParamsBM(p, mu = μ), q, a, b, c, T, 0.0; kws...) 
#         integrand_arr[i] = sum_lowerdiag((mat) .* (-lorentz(ϵs, ω, η) .- lorentz(ϵs, -ω, η)))
#     end
#         return (ħ_ev_s * π/2) .* integrand_arr
# end


#= DOCS

"""
`injection_current_\omega (ωlist, a, b, c, q, p, η; jdos = false, kws...)`
 BPGE with chosen current (else) at all ω in 
`ωlist` for a given momentum `q` which will be integrated over the BZ.
    See `σab_inter_linear(a, b, p, ωlist; kws...)`.
    returns Imag(η_abc)
 where η_abc = ∂tJ """
 =#