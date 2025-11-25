#_________________________________________________________________________________________
# position operators
#_________________________________________________________________________________________

r(ϵs, ψs, dh) = -1im .* vel(ψs, dh) ./ Ω(ϵs)                                    # Units: Å  
function zop(dim)                                                                           
    dim2 = div(dim, 2)  
    return  3.3/2 .* vcat(ones(dim2), -1 .* ones(dim2))       # 3.3 Å = interlayer distance                                            
end

function r_covariant(omega, ra, Δa, rb, Δb)                                    # Units: Å^2
    raa_cov = zeros(ComplexF64, size(ra,1), size(ra,1))#similar(ra)
    dim = size(omega, 1) 
    for n in 1:dim
        for m in 1:dim
            if m != n && n > m 
                aux = 0
                for p in 1:dim
                    if p != n && p != m
                        aux += 1im * (omega[n,p]*ra[n,p]*rb[p,m] - omega[p,m]*rb[n,p]*ra[p,m])
                        # Units: eV * Å * Å  # WARN DIVISIONS BY DIFFERENT OMEGAnm
                    else
                        aux += 0
                    end
                end
                raa_cov[n,m] = -((ra[n,m] * Δb[n,m] + rb[n,m] * Δa[n,m]) + aux)/omega[n,m]
                # Units: (Å * Å * eV) / eV
            else 
                nothing
            end
        end
    end
    return raa_cov
end

function z_covariant(zmat, rmat)
    dim = size(zmat, 1)
    mat = zeros(ComplexF64, size(rmat,1),size(rmat,1))
    @inbounds begin
        for n in 1:dim
            for m in 1:dim
                if m != n && n > m 
                    for p in 1:dim
                        if  p == n 
                            val = zmat[n,p] * rmat[p,m] 
                        elseif p == m 
                            val = -rmat[n,p] * zmat[p, m]
                        else
                            val = zmat[n,p] * rmat[p,m] - rmat[n,p] * zmat[p, m]
                        end
                        mat[n,m] += val    
                    end
                else nothing end
            end
        end
    end
    return -1im .* mat
end

#_________________________________________________________________________________________
# observables
#_________________________________________________________________________________________

lorentz(ϵs, ω, η) = 1/π .* η ./ (δ(ϵs, ω).^2 .+ η^2) 
δ(ϵs, ω) = ω .- Ω(ϵs)  # the minus sign comes from delta(e_f-e_i-omega); ωnm = -ωmn
# this definition is weak against degeneracies but it is highly efficient (benchmark)
Ω(ϵs) = ϵs .- t(ϵs) +  1e-7 * Diagonal(ϵs) .+ 1e-35# the last term is just to avoid infinities it'll cancel out later
# Δ(ψs, dh) = vel(ψs, dh) - t(vel(ψs, dh))
vel(ψs, dh) = ψs' * dh * ψs # It really is hbar v  # Units: meV Å  

function omega_r_Δ(ϵs, ψs, dh)
    vmat = vel(ψs, dh)
    mat = similar(vmat)
    omega = Ω(ϵs)
    for i in 1:size(vmat, 1)
        for j in 1:size(vmat, 1)
            mat[i, j] = vmat[i, i] - vmat[j, j]
        end
    end
    return omega, -1im .* vmat ./ omega, mat 
end

function Δ(ψs, dh)
    vmat = vel(ψs, dh)
    mat = similar(vmat)
    for i in 1:size(vmat, 1)
        for j in 1:size(vmat, 1)
            mat[i, j] = vmat[i, i] - vmat[j, j]
        end
    end
    return mat
end

"""
fn - fm where f are the Fermi Dirac distribution in the limit of T = 0.
I set p.
"""

function f!(mat, ϵs, μ, T, tol = 1e10) 
    replace!(mat, NaN => tol, NaN + NaN*im => tol *(1+im), 1im * NaN => tol*1im )
    mat .*= f(ϵs, μ, T)
end

f(ϵs, μ, T) = [fn(ϵs[i], μ, T) - fn(ϵs[j], μ, T) for i in 1:length(ϵs), j in 1:length(ϵs)]

fn(ϵn,μ::Float64) = ifelse(ϵn < μ, 1.0, 0.0)
function fn(ϵn, μ, T)
    if T == 0
        return ifelse(ϵn < μ, 1.0, 0.0)
    else
        return 1/(exp((ϵn - μ)/(k_B * T)) + 1)
    end
end
    
#-------------------------------------------------------------------------------------------
# Aux operations
#-------------------------------------------------------------------------------------------
function whichBPGE(part, r_cov, r, r_cov2, r2)
    if part == :REAL
        return symmetrized_imagsum(r_cov, r, r_cov2, r2) 
    else
        return antisymmetrized_realsum(r_cov, r, r_cov2, r2)
    end
end

function whichBPGE(part, r_cov, r)
    if part == :REAL
        return symmetrized_imagsum(r_cov, r)
    else 
        return antysymmetrized_realsum(r_cov, r)
    end
end


symmetrized_imagsum(mat1, mat2)  =  imag(mat1 .* t(mat2) .+ mat2 .* t(mat1))
symmetrized_imagsum(mat1, mat2, mat3, mat4) = imag((mat1 .* t(mat2)) .+ (mat3 .* t(mat4))) 

function antisymmetrized_realsum(mat1, mat2) 
    auxmat = similar(mat1)
    dim = size(mat1, 1)
    for n in 1:dim
        for m in 1:dim
            auxmat[n,m] = real(mat1[n,m] * mat2[m,n] - mat2[n,m] * mat1[m,n])
        end
    end
    return auxmat
end

antisymmetrized_realsum(mat1, mat2, mat3, mat4) = real(mat1 .* t(mat2) .- mat3 .* t(mat4))

t(mat) = transpose(mat)
sum_nondiag(mat) = sum(mat) - sum(Diagonal(mat))
sum_lowerdiag(mat::Array) =  sum_lowerdiag(mat, size(mat,1))
sum_lowerdiag(mat, dim) = sum(mat .* (tril(ones(dim, dim), 0) - Diagonal(ones(dim,dim))))

function return_kwargs(kws, name::Symbol)
    kwargs = NamedTuple(kws)
    # println(kwargs)
    if haskey(kwargs, name)
        return kwargs[name]
    else  nothing
    end
end

#_________________________________________________________________________________________
# Derivatives and geometry
#_________________________________________________________________________________________

function dhs(a, b, p)                                                                       
    unbounded = (:x, :y)
    if !(a ∈ unbounded && b ∈ unbounded)
        throw(ArgumentError("Only linear photocurrents with unbounded cartesian indices are implemented"))
    end
    dha =  dhf_hamiltonian(p, a)
    if a == b
        return dha, dha
    else
        dhb =  dhf_hamiltonian(p, b)
        return dha, dhb
    end
end

function dhs(a, b, k, p)                                                                       
    unbounded = (:x, :y)
    if !(a ∈ unbounded && b ∈ unbounded)
        throw(ArgumentError("Only linear photocurrents with unbounded cartesian indices are implemented"))
    end
    dha =  dhf_hamiltonian(p, k, a)
    if a == b
        return dha, dha
    else
        dhb =  dhf_hamiltonian(p, k, b)
        return dha, dhb
    end
end


"""
boundaries of the integration function for the HF model (same as for the BM model)
"""
function int_boundaries(p::ParamsHF)
    G1, G2 = bravais_vectors(p)
    M = m(G1, G2)
    v_length = norm(G2 - G1)
    h_length = norm(G2 + G1)
    κ1 = κ(G1, G2)
    return M, [0, -v_length/2], [2M[1], v_length/2 ]
end