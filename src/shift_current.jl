"""
    `shift_current_ω(ωlist, a, b, c, q, p, η; jdos = false, kws...)`
returns the JDOS (if `jdos = true`) or the BPGE with chosen current (else) at all ω in 
`ωlist` for a given momentum `q` which will be integrated over the BZ.
    See `σab_inter_linear(a, b, p, ωlist; kws...)`.
"""
function shift_current_ω(part, ωlist, a, b, c, q, p, η;  omega = 0.005, T = 0)
    h = hf_hamiltonian(p, q) 
    ϵs, ψs = eigen(Matrix(h))
    mat = shift_current_integrand(part, ϵs, ψs, p, q, a, b, c, T, 0)
    return [- π/2 * sum_lowerdiag(mat .* (-lorentz(ϵs, ω, η) .- lorentz(ϵs, -ω, η))) for ω in ωlist] * 1000  # Units: Å^3 /eV # tje 1000 is required because lorentz units are meV not evs # what should be
end

shift_current_integrand(part, ϵs, ψs, p, q, a, b, c, T, doping) =  
    f(ϵs, doping, T) .* Imat_shift(part, a, b, c, ϵs, ψs, p, q)

""" returns the symmetrized product sum """
function Imat_shift(part, a, b, c, ϵs, ψs, p, q)
    bounded_dirs = findall(x -> x ==:z, (a,b,c))
    if !isempty(bounded_dirs)
        throw(ArgumentError("We cannot compute out out-of-plane currents since the z 
          direction is bounded but also since we lost teh z label we cannot compute oblique incidence"))
    else 
        if p.λ == 0.0 #path with delta-like localization. slightly faster
            return Imat_shift_uuu(part, a, b, c, ϵs, ψs, p)
        else # path with exponentially localized wannier orbs
            return Imat_shift_uuu(part, a, b, c, ϵs, ψs, p, q)
        end
    end
end

"""
 a, b, c correspond to unbounded directions
 """
function Imat_shift_uuu(part, a, b, c, ϵs, ψs, p)
    omega, ra, Δa = omega_r_Δ(ϵs, ψs, dhf_hamiltonian(p, a))
    r_cov = similar(ra)
    if a == b == c
        r_cov .= r_covariant(omega, ra, Δa, ra, Δa)
        return whichBPGE(part, r_cov, ra, r_cov ,ra)
    else 
        _, rb, Δb = omega_r_Δ(ϵs, ψs, dhf_hamiltonian(p, b))
        _, rc, Δc = omega_r_Δ(ϵs, ψs, dhf_hamiltonian(p, c))
        r_cov_ca = r_covariant(omega, rc, Δc, ra, Δa)
        r_cov_ba = r_covariant(omega, rb, Δb, ra, Δa)
        return whichBPGE(part, r_cov_ca, rb, r_cov_ba, rc)
    end
end

function Imat_shift_uuu(part, a, b, c, ϵs, ψs, p, q::Vector{Float64})
    omega, ra, Δa = omega_r_Δ(ϵs, ψs, dhf_hamiltonian(p, q, a))
    r_cov = similar(ra)
    if a == b == c
        r_cov .= r_covariant(omega, ra, Δa, ra, Δa)
        return whichBPGE(part, r_cov, ra, r_cov ,ra)
    else 
        _, rb, Δb = omega_r_Δ(ϵs, ψs, dhf_hamiltonian(p, q, b))
        _, rc, Δc = omega_r_Δ(ϵs, ψs, dhf_hamiltonian(p, q, c))
        r_cov_ca = r_covariant(omega, rc, Δc, ra, Δa)
        r_cov_ba = r_covariant(omega, rb, Δb, ra, Δa)
        return whichBPGE(part, r_cov_ca, rb, r_cov_ba, rc)
    end
end

#########
# TWO VALLEY MODEL SPINFUL MODEL

##########

function shift_current_ω(part, ωlist, a, b, c, q, p, of, η;  omega = 0.005, T = 0)
    whichham = (p.twovalleystwospins == false ? hf_twovalleyshamiltonian : hf_valley_spin_hamiltonian)
    h = whichham(p, q, of) 
    ϵs, ψs = eigen(Matrix(h))
    return [abs(sum_lowerdiag(η ./ ((ϵs .- ω).^2 .+ η .^2))) for ω in ωlist]   # DOS comment or uncomment
end

shift_current_integrand(part, ϵs, ψs, p, q, a, b, c, T, doping, of) =  
    f!(Imat_shift(part, a, b, c, ϵs, ψs, p, of, q), ϵs, doping, T)


function Imat_shift(part, a, b, c, ϵs, ψs, p, of, q)
    bounded_dirs = findall(x -> x ==:z, (a,b,c))
    if !isempty(bounded_dirs)
        throw(ArgumentError("We cannot compute out out-of-plane currents since the z 
          direction is bounded but also since we lost teh z label we cannot compute oblique incidence"))
    else 
        if p.λ == 0.0 #path with delta-like localization. slightly faster
            # println("Delta-like localization method")
            return Imat_shift_uuu(part, a, b, c, ϵs, ψs, p, of)
        else # path with exponentially localized wannier orbs
            # println("Exponential localization method")
            return Imat_shift_uuu(part, a, b, c, ϵs, ψs, p, of, q)
        end
    end
end

function Imat_shift_uuu(part, a, b, c, ϵs, ψs, p, of::Union{SparseMatrixCSC, Matrix})
    omega, ra, Δa = omega_r_Δ(ϵs, ψs, dhf_hamiltonian(p, a, of))
    r_cov = zeros(ComplexF64, size(ra,1), size(ra,1))
    if a == b == c
        r_cov .= r_covariant(omega, ra, Δa, ra, Δa)
        return whichBPGE(part, r_cov, ra, r_cov ,ra)
    else 
        _, rb, Δb = omega_r_Δ(ϵs, ψs, dhf_hamiltonian(p, b, of))
        _, rc, Δc = omega_r_Δ(ϵs, ψs, dhf_hamiltonian(p, c, of))
        r_cov_ca = r_covariant(omega, rc, Δc, ra, Δa)
        r_cov_ba = r_covariant(omega, rb, Δb, ra, Δa)
        return whichBPGE(part, r_cov_ca, rb, r_cov_ba, rc)
    end
end

function Imat_shift_uuu(part, a, b, c, ϵs, ψs, p, of, q)
    omega, ra, Δa = omega_r_Δ(ϵs, ψs, dhf_hamiltonian(p, q, a, of))
    r_cov = similar(ra)
    if a == b == c
        r_cov .= r_covariant(omega, ra, Δa, ra, Δa)
        return whichBPGE(part, r_cov, ra, r_cov ,ra)
    else 
        _, rb, Δb = omega_r_Δ(ϵs, ψs, dhf_hamiltonian(p, q, b, of))
        _, rc, Δc = omega_r_Δ(ϵs, ψs, dhf_hamiltonian(p, q, c, of))
        r_cov_ca = r_covariant(omega, rc, Δc, ra, Δa)
        r_cov_ba = r_covariant(omega, rb, Δb, ra, Δa)
        return whichBPGE(part, r_cov_ca, rb, r_cov_ba, rc)
    end
end