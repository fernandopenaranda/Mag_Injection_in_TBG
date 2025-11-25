"""
hamiltonian structure: [[ff] [fc] ;[cf] [cc]]
[ff] = 0σ0 localized electrons (0 hybridization) possibly add a mu
dim([cc]) = 4length(indices).
H = μ I f⁺ₖfₖ + ∑ₖ h(k+G) c⁺_{k+g}+gc_{k+g} + ∑ₖ ∑_g h(k+g) fₖc_{k+g} + h.c.
the k is restricted to the MBZ. The unbounded k of the conduction electrons
is regularized into a mesh of reciprocal vectors G assuming that the cutoff is
infinite which leads to periodicity. 
See notes.
"""
function bare_hf_hamiltonian(p::ParamsHF, k)
    indices = gs_indices(p)
    ham_dim = ham_matrix_size(indices)
    hf_mat = spzeros(ComplexF64, ham_dim, ham_dim)
    f_electrons!(hf_mat, p, k)
    conduction_electrons!(hf_mat, p, indices, k)
    conduction_localized_coupling!(hf_mat, p, indices, k)
end

function hf_hamiltonian(p::ParamsHF,k) # with all modifiers
    indices = gs_indices(p)
    ham_dim = ham_matrix_size(indices)
    hf_mat = spzeros(ComplexF64, ham_dim, ham_dim)
    f_electrons!(hf_mat, p, k)
    conduction_electrons!(hf_mat, p, indices, k)
    conduction_localized_coupling!(hf_mat, p, indices, k)
    hf_hartree_mod!(hf_mat, p)
    hf_vafek_mod!(hf_mat, p, indices)
    if p.VP == true
        hf_valleypolarized!(hf_mat, p, indices)
    end
    hf_mat
end

""" energy of the f electrons.
The second term is the sigma_z mass (coming from hBN where sigma is the sublattice).
In the HF proyection it reads as a σz × -σz × σz in the f c c subspace
where σ is now the orbital d.o.f."""
function f_electrons!(mat, p, k)
    mat[1:2,1:2] .+= -p.μ .* σ0 + p.sigmaz .* σz + p.sigmazlayerz .* σz * non_local(p, k) + 
        p.layerz .* σ0 * non_local(p, k)
end
"""hamiltonian of the conduction electrons in the basis spanned by g in Gmesh"""
function conduction_electrons!(mat, p, indices, k)
    σvec = SA[p.ν.*σx, σy]
    g = zeros(Float64, 2)
    ham_dim = size(mat, 1)
    cc_mat = zeros(ComplexF64, 4, 4)
    g1, g2 = bravais_vectors(p)
    for i in 3:4:ham_dim-3
        gs_vector!(g, g1, g2, indices[Int(i-2+3)÷4])
        hc!(cc_mat, k+g, p)
        mat[i:i+3, i:i+3] .= cc_mat
    end
end

function hc!(mat, k, p)
    mat[1:2,1:2] .= -p.μ .* σ0 -  p.sigmaz .* σz
    mat[1:2,3:4] .= conj.(hc12(k, p))
    mat[3:4,1:2] .= hc12(k, p)
    mat[3:4,3:4] .= p.M .* σx - p.μ .* σ0 + p.sigmaz .* σz
end

hc12(k, p) = p.v .* (p.ν*k[1]/a0 .* σ0 .- 1im*k[2]/a0 .* σz)
# hc12(k, p) = p.v .* (p.ν*k[1] .* σ0 .- 1im*k[2] .* σz)

function conduction_localized_coupling!(mat, p, indices, k)
    σvec = SA[p.ν.*σx, σy] 
    g = zeros(Float64, 2)
    ham_dim = size(mat, 1)
    cf_mat = zeros(ComplexF64, 2, 2)
    cf_gammas_mat = zeros(ComplexF64, 2, 2)
    g1, g2 = bravais_vectors(p)
    for i in 3:4:ham_dim
        gs_vector!(g, g1, g2, indices[Int(i-2+3)÷4])
        hcf!(cf_mat, k+g, p)
        hgamma1_2_f!(cf_gammas_mat, k+g, p)
        mat[i:i+1, 1:2] .= cf_mat
        mat[1:2, i:i+1] .= cf_mat'
        mat[i+2:i+3, 1:2] .= cf_gammas_mat
        mat[1:2, i+2:i+3] .= cf_gammas_mat'
    end
    mat
end

function hcf!(mat, k, p)
    mat .= (p.γ .* σ0 + p.vp .* (p.ν * k[1]/a0 .* σx .+ k[2]/a0 .* σy)
        # ) .* exp(-(norm(k/a0)^2 * p.λ^2)/2) +  p.c2ybreakingmass .* σ0
        # ) .* exp(-(norm(k/a0)^2 * p.λ^2)/2) +  1im * p.c2ybreakingmass .* σ0 # opposite contributions two valleys
        ) .* exp(-(norm(k/a0)^2 * p.λ^2)/2) + p.ν * 1im * p.c2ybreakingmass .* σ0 #M25 breaks C2y same ocntributions (A1, -, +) and the p.ν is required to anticommute with c2y two valleys
end

""" v_⋆⋆ term """
function hgamma1_2_f!(mat, k, p)
    mat .= p.vpp * (p.ν * k[1]/a0 .* σx .- k[2]/a0 .* σy) .* exp(-(norm(k/a0)^2 * p.λ^2)/2)
end

# function hcf!(mat, k, p)
#     mat .= (p.γ .* σ0 + p.vp .* (p.ν * k[1] .* σx .+ k[2] .* σy)
#         # ) .* exp(-(norm(k/a0)^2 * p.λ^2)/2) +  p.c2ybreakingmass .* σ0
#         # ) .* exp(-(norm(k/a0)^2 * p.λ^2)/2) +  1im * p.c2ybreakingmass .* σ0 # opposite contributions two valleys
#         ) .* exp(-(norm(k)^2 * p.λ^2)/2) + p.ν * 1im * p.c2ybreakingmass .* σ0 #M25 breaks C2y same ocntributions (A1, -, +) and the p.ν is required to anticommute with c2y two valleys
# end




################ 
# Modifiers
############### 
function hf_hartree_mod!(mat, p) # CHECK THIS
    mat[1:2,1:2] .+= sign(p.nf-1)* (p.nf-1)^2 * p.U1/2 .* σ0 # it should be this
end

function hf_mass_mod!(mat, p)
    mat[1:2,1:2] .+= p.sigmaz .* σz #mass that breaks C2x (A2,-,-) preserves C2y
end

function hf_vafek_mod!(mat, p, indices)
    g = zeros(Float64, 2)
    ham_dim = size(mat, 1)
    cc_mat = zeros(ComplexF64, 4, 4)
    g1, g2 = bravais_vectors(p)
    cc_mat[1:2,1:2] = p.vafek .* +σ0
    cc_mat[3:4,3:4] = p.vafek .* -σ0
    for i in 3:4:ham_dim-3
        gs_vector!(g, g1, g2, indices[Int(i-2+3)÷4])
        mat[i:i+3, i:i+3] .+= cc_mat
    end
end

#= 
Valley polarized GS at ν = 0
Made out of Hu, Hj. Hw = Hv = 0
Eu = Cnst and Ej= EV=Ew = 0
I do not introduce selfconsistency for the moment (this is what they call the one-shot GS)
    ## nf, μ, n = selfconsistency(p, 0; tol = 1e-3, α = 0.4) I do not implement here self consistency,
    # i do the one-shot thing
    =#
    
function hf_valleypolarized!(mat, p, indices; val = 1) # acts on the single valley model no self consistency
    # density density local (no hartree) Note that it is finite because at charge neutrality one valley is above and the other is below.
    mat[1:2,1:2] .+= -val*p.ν * p.U1/2 .* σ0
    # exchange interaction in the conduction bands corresponding to the gamma_1+gamma_2 sector
    g = zeros(Float64, 2)
    ham_dim = size(mat, 1)
    g1, g2 = bravais_vectors(p)
    cc_mat = val * p.ν * p.J .* σ0
    for i in 3:4:ham_dim-3
        gs_vector!(g, g1, g2, indices[Int(i-2+3)÷4])
        mat[i+2:i+3, i+2:i+3] .+= cc_mat
    end
end

"""
the two valleys are explicitly included. 
there is still a spins degeneracy assumed.
of is the density matrix of the f electrons, only required for the correlated gs

hf_plotbands(ParamsHF(p, sigmaz=0, μ =0, vafek= 0, J = -8, U1 = 10, VP = true, twovalleys = true, KIVC=false))
"""
function hf_twovalleyshamiltonian(p, k, of = [1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0])
    h = [bare_hf_hamiltonian(ParamsHF(p, ν = 1), k) 0I; 0I   bare_hf_hamiltonian(ParamsHF(p, ν = -1), k)]
    if p.KIVC == true
        kivc!(h, p)
    end
    if p.VP == true
        # nothing
        # vp_ν0!(h, p) # gs VP at ν = 0 without self consistency
        # of = Of(p, p.μ)
        # of = [1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0] # guess neutrality
        # of = [0.5 0 0 0; 0 0.5 0 0; 0 0 0 0; 0 0 0 0] # -1 filling

        vp_sc_α!(h, of, p)

    end
    return h 
end

"""
this generates the non selfconsistent gs at zero filling just for numerical tests
    Ex: 
hf_plotbands(ParamsHF(p, sigmaz=0, μ =0, vafek= 0, J = -8, U1 = 35, VP = true, twovalleys = true))
"""
function vp_ν0!(mat, p)
    indices = gs_indices(p)
    halfdim = size(mat,1) ÷ 2
    aux_mat = spzeros(ComplexF64, halfdim, halfdim)
    hf_valleypolarized!(aux_mat, p, indices, val = 1)
    mat[1:halfdim, 1:halfdim] .+= aux_mat
    aux_mat .*= 0
    hf_valleypolarized!(aux_mat, p, indices, val = -1)
    mat[1+halfdim:end, 1+halfdim:end] .+= aux_mat
end

function vp_sp_ν0!(mat, p)
    indices = gs_indices(p)
    quartdim = size(mat,1) ÷ 4
    aux_mat = spzeros(ComplexF64, quartdim, quartdim)
    # spin up
    hf_valleypolarized!(aux_mat, p, indices, val = 1)
    mat[1:quartdim, 1:quartdim] .+= aux_mat
    aux_mat .*= 0
    hf_valleypolarized!(aux_mat, p, indices, val = -1)
    mat[1+quartdim:2quartdim, 1+quartdim:2quartdim] .+= aux_mat
    aux_mat .*= 0
    # spin down
    hf_valleypolarized!(aux_mat, p, indices, val = 1)
    mat[2quartdim+1:3quartdim, 2quartdim+1:3quartdim] .+= aux_mat
    aux_mat .*= 0
    hf_valleypolarized!(aux_mat, p, indices, val = -1)
    mat[3quartdim+1:end, 3quartdim+1:end] .+= aux_mat
end

"""
this generates the selfconsistent gs at arbitrary filling with n ∈ [-2,2] 
hf_plotbands(ParamsHF(p, sigmaz=0, μ =0, vafek= 0, J = -8, U1 = 35, VP = true, twovalleys = true))
we only consider terms in U1 and J, we impose that Of is diagonal in the valley dof. 
α respects the nomenclature of PHYSICAL REVIEW B 103, 035427 (2021) and it is the VP phase discribed by Song et al in the paper.
"""
function vp_sc_α!(mat, of, p)
    indices = gs_indices(p)
    halfdim = size(mat,1) ÷ 2
    mat[1:2,1:2] .+=  p.U1 .* ( (trace(of)-2+1/2).*1I - of[1:2,1:2]) #-2p.U1*1I # positive valley
    mat[1+halfdim:2+halfdim,1+halfdim:2+halfdim] .+=  p.U1 .* ( (trace(of)-2+1/2) .* 1I - of[3:4,3:4]) #-2p.U1*1I # negative valley
    # exchange interaction in the conduction bands corresponding to the gamma_1+gamma_2 sector
    g = zeros(Float64, 2)
    ham_dim = size(mat, 1)
    g1, g2 = bravais_vectors(p)
    cc_mat_pv = p.J .* (2of[1:2,1:2] - 1I)
    cc_mat_mv = p.J .* (2of[3:4,3:4] - 1I)
    for i in 3:4:halfdim-3
        gs_vector!(g, g1, g2, indices[Int(i-2+3)÷4])
        mat[i+2:i+3, i+2:i+3] .+= cc_mat_pv 
        mat[i+2+halfdim:i+3+halfdim, i+2+halfdim:i+3+halfdim] .+= cc_mat_mv

    end
end


"""
generates the HF model with all degrees of freedom included
"""

function hf_valley_spin_hamiltonian(p::ParamsHF, k, of)
    hdim = ham_matrix_size(p)
    h = spzeros(ComplexF64, hdim, hdim)
    hf_valley_spin_hamiltonian!(h, p, k, of)
end

function hf_valley_spin_hamiltonian!(h, p::ParamsHF, k, of)
    hdim = size(h,1)
    hp = bare_hf_hamiltonian(ParamsHF(p, ν = 1), k)
    hn = bare_hf_hamiltonian(ParamsHF(p, ν = -1), k)
    tauz_pert = 0 *1I # REMOVE this is a small perturbation required to select a manifold with the same valley polarization. And thus a non jumping inj current. Possibly J fixed this problem
    h[1:hdim÷4,1:hdim÷4] .= hp + tauz_pert
    h[1+hdim÷4:2hdim÷4,1+hdim÷4:2hdim÷4]   .= hn - tauz_pert
    h[1+2hdim÷4:3hdim÷4,1+2hdim÷4:3hdim÷4] .= hp - tauz_pert

    # h[1+2hdim÷4:3hdim÷4,1+2hdim÷4:3hdim÷4] .= hp + 0tauz_pert
    h[1+3hdim÷4:4hdim÷4,1+3hdim÷4:4hdim÷4] .= hn - tauz_pert
    if p.VP == true # this key is not the best name to CHANGE
        vp_sp_sc!(h, of, p) # correlations of typep U and J, perturbation type 
    end
    return h
end

"""
Correlation of type U and J only. U and J only involve f electrons *this is a reasonable approximation
in the Mean field decoupling.
"""
function vp_sp_sc!(mat, of, p::ParamsHF)
    indices = gs_indices(p)
    quartdim = size(mat,1) ÷ 4

    mat[1:2, 1:2] .+=  p.U1 .* ((trace(of) - 3.5) .* diagm(ones(2)) - of[1:2,1:2])                                         # positive valley spin up
    mat[1+quartdim: 2+quartdim, 1+quartdim:2+quartdim]    .+= p.U1 .* ((trace(of) - 3.5) .* diagm(ones(2)) - of[3:4,3:4])  # negative valley spin up
    mat[2quartdim+1:2quartdim+2, 2quartdim+1:2quartdim+2] .+= p.U1 .* ((trace(of) - 3.5) .* diagm(ones(2)) - of[5:6,5:6])  # positive valley spin down
    mat[3quartdim+1:3quartdim+2, 3quartdim+1:3quartdim+2] .+= p.U1 .* ((trace(of) - 3.5) .* diagm(ones(2)) - of[7:8,7:8])  # negative valley spin down
    
    mat[1:2, 1:2] .+= 6p.U2 .* (trace(of)-4) .* diagm(ones(2))# positive valley spin up
    mat[1 + quartdim: 2 + quartdim, 1+quartdim:2+quartdim] .+= 6p.U2 .* (trace(of)-4) .*  diagm(ones(2)) # negative valley spin up
    mat[2quartdim+1:2quartdim+2, 2quartdim+1:2quartdim+2]  .+= 6p.U2 .* (trace(of)-4) .* diagm(ones(2))# positive valley spin down
    mat[3quartdim+1:3quartdim+2, 3quartdim+1:3quartdim+2]  .+= 6p.U2 .* (trace(of)-4) .* diagm(ones(2))# negative valley spin down
    #exchange interaction in the conduction bands corresponding to the gamma_1+gamma_2 sector
    cc_mat_pv_su = -2 * p.J .* (of[1:2,1:2] - 1/2*1I) # positive valley spin up
    cc_mat_mv_su = -2 * p.J .* (of[3:4,3:4] - 1/2*1I) # negative valley spin up
    cc_mat_pv_sd = -2 * p.J .* (of[5:6,5:6] - 1/2*1I) # positive valley spin down
    cc_mat_mv_sd = -2 * p.J .* (of[7:8,7:8] - 1/2*1I) # negative valley spin down

    for i in 3:4:quartdim-3
        mat[i+2:i+3, i+2:i+3]                                         .+= cc_mat_pv_su 
        mat[i+2+quartdim:i+3+quartdim, i+2+quartdim:i+3+quartdim]     .+= cc_mat_mv_su
        mat[i+2+2quartdim:i+3+2quartdim, i+2+2quartdim:i+3+2quartdim] .+= cc_mat_pv_sd
        mat[i+2+3quartdim:i+3+3quartdim, i+2+3quartdim:i+3+3quartdim] .+= cc_mat_mv_sd
    end
end
    #diagonal terms as in the paper

    # as in Fernando. It cannot be since H_U won't distinguish local f interactions between same orbitals

    # mat[1:2, 1:2] .+=  2*p.U1 * (of[1:2,1:2] - 4 .* diagm(ones(2))) 
    # mat[1+quartdim:2+quartdim, 1+quartdim:2+quartdim] .+=  2*p.U1 * (of[1:2,1:2] - 4 .* diagm(ones(2))) 
    # mat[2quartdim+1:2quartdim+2, 2quartdim+1:2quartdim+2] .+=  2*p.U1 * (of[1:2,1:2] - 4 .* diagm(ones(2))) 
    # mat[3quartdim+1:3quartdim+2, 3quartdim+1:3quartdim+2] .+=  2*p.U1 * (of[1:2,1:2] - 4 .* diagm(ones(2))) 

# function vp_sp_sc_β!(mat, of, p)
#     vp_sp_sc_α!(mat, of, p)
#     quartmat = size(mat,1) ÷ 4
#     quartid = diagm(ones(quartmat))
#     par = 6
#     mat .+= [-quartid .* 2par 0I 0I 0I; #sup vK # first occupied
#              0I quartid .* 1par  0I 0I; #sup vK'
#              0I 0I -quartid .* 1*par 0I; #sdown vK #fisrt occupied
#              0I 0I 0I quartid .* 0par]  #sdown vK'
# end


# function vp_sp_sc_β_fsector!(mat, of, p)
#     vp_sp_sc_α!(mat, of, p)
#     quartmat = size(mat,1) ÷ 4
#     quartid = diagm(vcat([1,1],zeros(quartmat-2)))
#     par = 6
#     mat .+= [-quartid .* 2par 0I 0I 0I; #sup vK # first occupied
#              0I quartid .* 1par  0I 0I; #sup vK'
#              0I 0I -quartid .* 1*par 0I; #sdown vK #fisrt occupied
#              0I 0I 0I quartid .* 0par]  #sdown vK'
# end

trace(mat) = sum(diag(mat))

# particle hole representation in the momentum basis
# -σyτzK
function isphsym(h)
    phs = zeros(ComplexF64, size(h,1), size(h,1))
    quartdim = size(phs,1)÷4
    cdim = quartdim - 2 # substract the dimension of the f electrons
    phs .*= 0.0 
    # f sector
    phs[1:2, 1:2] .+= -σy                                          # positive valley spin up
    phs[1+quartdim: 2+quartdim, 1+quartdim:2+quartdim]    .+= σy # negative valley spin up
    phs[2quartdim+1:2quartdim+2, 2quartdim+1:2quartdim+2] .+= -σy # positive valley spin down
    phs[3quartdim+1:3quartdim+2, 3quartdim+1:3quartdim+2] .+= σy # negative valley spin down
    # c sector
    for i in 1:2:quartdim-2
        phs[i+2:i+3, i+2:i+3]                                         .+= σy
        phs[i+2+quartdim:i+3+quartdim, i+2+quartdim:i+3+quartdim]     .+= -σy
        phs[i+2+2quartdim:i+3+2quartdim, i+2+2quartdim:i+3+2quartdim] .+= σy
        phs[i+2+3quartdim:i+3+3quartdim, i+2+3quartdim:i+3+3quartdim] .+= -σy
    end
    # for j in 1:4:cdim÷4
    # phs[3:quartdim,3:quartdim] .+=  σy_blocks(cdim÷2) ./2
    # phs[3+quartdim:2quartdim, 3+quartdim:2quartdim] .+= -σy_blocks(cdim÷2) ./2
    # phs[2quartdim+3:3quartdim,2quartdim+3:3quartdim] .+= σy_blocks(cdim÷2) ./2
    # phs[3quartdim+3:4quartdim,3quartdim+3:4quartdim] .+= -σy_blocks(cdim÷2) ./2
    # end
    sparse(phs)
end

function σy_blocks(n::Int)
    # Create a block diagonal matrix with n σ_y matrices
    return sparse(kron(Diagonal(ones(n)), σy))
end