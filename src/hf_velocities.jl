""" k derivative of the HF model.
If λ is finite we will have a k dependency in the derivative of the Hamiltonian.
Taking the limit λ->0 makes the derivative k independent. 
We implement here the two methods:
    `dhf(p, a)` for the `λ = 0` limit.
    `dhf(p, k, a)` for λ ≠ 0.
Note that all the effects of a finite λ enter as an order 2 perturbation in the velocity operator.
I think it is safe to ignore it, but be aware that h should go always with λ = 0 to be consistent. """
function dhf_hamiltonian(p::ParamsHF, a::Symbol)
    indices = gs_indices(p)
    ham_dim = ham_matrix_size(indices)
    dhf_mat = spzeros(ComplexF64, ham_dim, ham_dim)
    d_conduction_electrons!(dhf_mat, p, indices, a)
    d_conduction_localized_coupling!(dhf_mat, p, indices, a)
    dhf_mat
end

""" k-independent k derivative of the Hamiltonian w.r.t. the two spins two valley model"""
function dhf_hamiltonian(p, a, of::Union{SparseMatrixCSC, Matrix})
    dim = size(of,1)
    dh_pv = dhf_hamiltonian(ParamsHF(p, ν = 1), a)
    dh_nv = dhf_hamiltonian(ParamsHF(p, ν = -1), a)
    if dim == 4
        return [dh_pv 0I; 0I dh_nv]
    elseif dim == 8
        aux_dh = [dh_pv 0I; 0I dh_nv]
        return [aux_dh 0I;0I aux_dh] # two valleys two spins
    end
end

function dhf_hamiltonian(p::ParamsHF, k, a::Symbol)
    indices = gs_indices(p)
    ham_dim = ham_matrix_size(indices)
    dhf_mat = spzeros(ComplexF64, ham_dim, ham_dim)
    d_conduction_electrons!(dhf_mat, p, indices, a)
    d_conduction_localized_coupling!(dhf_mat, p, indices, k, a)
    d_substrate!(dhf_mat, p, k, a)
    return dhf_mat
end

""" k-independent k derivative of the Hamiltonian w.r.t. the two spins two valley model"""
function dhf_hamiltonian(p, k, a::Symbol, of)
    dim = size(of,1)
    dh_pv = dhf_hamiltonian(ParamsHF(p, ν = 1), k, a)
    dh_nv = dhf_hamiltonian(ParamsHF(p, ν = -1), k, a)
    if dim == 4
        return [dh_pv 0I; 0I dh_nv]
    elseif dim == 8
        aux_dh = [dh_pv 0I; 0I dh_nv]
        return [aux_dh 0I;0I aux_dh] # two valleys two spins
    end
end

"""hamiltonian of the conduction electrons in the basis spanned by g in Gmesh"""
function d_conduction_electrons!(mat, p, indices, a)    
    σvec = SA[p.ν.*σx, σy]
    g = zeros(Float64, 2)
    ham_dim = size(mat, 1)
    cc_mat = zeros(ComplexF64, 4, 4)
    g1, g2 = bravais_vectors(p)
    for i in 3:4:ham_dim-3
        dhc!(cc_mat, p, a)
        mat[i:i+3, i:i+3] .= cc_mat
    end
end

function dhc!(mat, p, a)
    dhc12 = (a == :x ?  dhc12x : dhc12y)
    mat[1:2,3:4] .= conj.(dhc12(p))
    mat[3:4,1:2] .= dhc12(p)
    mat[3:4,3:4] .= 0σx
end

dhc12x(p) = p.v*p.ν/a0 .* σ0   
dhc12y(p) = -1im*p.v/a0 .* σz

function d_conduction_localized_coupling!(mat, p, indices, a)
    dhcf! = (a == :x ? dhcfx! : dhcfy!)
    σvec = SA[p.ν.*σx, σy] 
    g = zeros(Float64, 2)
    ham_dim = size(mat, 1)
    cf_mat = zeros(ComplexF64, 2, 2)
    g1, g2 = bravais_vectors(p)
    for i in 3:4:ham_dim
        dhcf!(cf_mat, p)
        mat[i:i+1, 1:2] .= cf_mat # the 2 is to keep the ff box in the upper left of the ham matrix
        mat[1:2, i:i+1] .= cf_mat'
    end
    mat
end

function d_conduction_localized_coupling!(mat, p, indices, k, a)
    dhcf! = (a == :x ? dhcfx! : dhcfy!)
    σvec = SA[p.ν.*σx, σy] 
    g = zeros(Float64, 2)
    ham_dim = size(mat, 1)
    cf_mat = zeros(ComplexF64, 2, 2)
    cf_gammas_mat = zeros(ComplexF64, 2, 2)
    g1, g2 = bravais_vectors(p)
    for i in 3:4:ham_dim
        gs_vector!(g, g1, g2, indices[Int(i-2+3)÷4])
        dhcf!(cf_mat, k+g, p)
        mat[i:i+1, 1:2] .= cf_mat # the 2 is to keep the ff box in the upper left of the ham matrix
        mat[1:2, i:i+1] .= cf_mat'
        mat[i+2:i+3, 1:2] .= cf_gammas_mat # the 2 is to keep the ff box in the upper left of the ham matrix
        mat[1:2, i+2:i+3] .= cf_gammas_mat'
    end
    mat
end

function d_substrate!(mat, p, k, a)
    mat[1:2,1:2] .= p.sigmazlayerz * d_non_local(p, k, a) .* σz + 
        p.layerz * d_non_local(p, k, a) .* σ0
end

function dhcfx!(mat, p)
    mat .= p.vp*p.ν .* σx
end

function dhcfy!(mat, p)
    mat .= p.vp .* σy
end

function dhgamma1_2_fx!(mat, p)
    mat .= p.vpp*p.ν .* σx
end

function dhgamma1_2_fy!(mat, p)
    mat .= -p.vpp .* σy
end

function dhcfx!(mat, k, p)
    mat .= ((p.γ.* σ0 + (p.vp*k[2]/a0).*σy) 
        .* (-k[1]/a0 * p.λ^2) + p.vp*p.ν*(1-(k[1]/a0)^2 *p.λ^2).*σx) .* exp(-(norm(k/a0)^2 * p.λ^2)/2)
end

function dhcfy!(mat, k, p)
    mat .= ((p.γ.* σ0 + (p.vp*p.ν*k[1]/a0).* σx) 
        .* (-k[2]/a0* p.λ^2) + p.vp*(1-(k[2]/a0)^2 *p.λ^2).*σy) .* exp(-(norm(k/a0)^2 * p.λ^2)/2)
end