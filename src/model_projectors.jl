using Cubature
using LinearAlgebra
using Cubature, QuadGK

const α1 = 0.7882# here \theta is independent of M so we take one value equal to 1.05
const α2 = 0.6155 # these are approximate values.
λ1(p) =  0.1791 * a0/(2pi*p.θ/(2*360))
λ2(p) =  0.1850 * a0/(2pi*p.θ/(2*360))

P_proj(p::ParamsHF, k) = P_proj(paramsBM(1.05, 3,  p.ν), k) 

"""projector of the Hamiltonian in the BM basis -could be the Hamiltonian- into the 6 lowest bands"""
Q_H_Q(p; k = 0.) = Q_X_Q(p, bistritzer_hamiltonian(p,k), k)

"""projector of a generic observable in the BM basis -could be the Hamiltonian- into the 6 lowest bands (2 f electrons + 4 c electrons).
returns a matrix of the size of hbm whose matrix elements are
∑u_α,G,n (k) * conj(u_β,G',n (k)). n is the band index restricted to the lowest 6 bands. """
function P_operator(p::ParamsBM, k)
    _, evecs = eigen(Matrix(bistritzer_hamiltonian(p,k)))
    half_dim = size(evecs,1)÷2
    u_ns = [evecs[:,i] for i in half_dim-2:half_dim+3]          # 6 lowest energy bands of the BM model
    dim_klat = Int(1 + 3p.nmax * (1 + p.nmax))                  # dimension of the plane-wave basis in the BM model
    dim_ham = 2 * 2 * dim_klat                                  # σ ∘ k_space ∘ τ 
    P_proj = zeros(ComplexF64, dim_ham, dim_ham)
    for n in 1:6                                                # six lowest energy bands
        P_proj .+= u_ns[n] * u_ns'[n]
    end
    P_proj
end 

""" Q proyector operator to the two Wannier functions of the local orbitals.
Eq. S61 in terms of the fourier transforms of the Wannier orbitals. Note that this only works for θ = 1.05°.
This is the operator for a given valley `p.ν`"""
function Q_operator(p::ParamsBM, k) #works
      Gs = Gmesh(Gindices(p), rec_vecs(p)) # vectors of the k-star.
      dim_klat = Int(1 + 3p.nmax * (1 + p.nmax)) # dimension of the plane-wave basis in the BM model
      dim_ham = 2 * 2 * dim_klat# σ ∘ k_space ∘ τ 
      hd = dim_ham ÷ 2
      Q_proj = zeros(ComplexF64, dim_ham, dim_ham)
      for (i, gi) in enumerate(Gs) # dimension of the k_space
        for (j, gj) in enumerate(Gs)
            Q_proj[2i-1:2i, 2j-1:2j] =  vv_an_sublat(p, 1, 1, p.ν, k, Gs[i], Gs[j])             # top-top
            Q_proj[2i-1+hd:2i+hd, 2j-1+hd:2j+hd] = vv_an_sublat(p, 2, 2, p.ν, k, gi, gj)        # bot-bot
            Q_proj[2i-1:2i, 2j-1+hd:2j+hd] =  vv_an_sublat(p, 1, 2, p.ν, k, gi, gj)             # top-bot
            Q_proj[2i-1+hd:2i+hd, 2j-1:2j] =  vv_an_sublat(p, 2, 1, p.ν, k, gi, gj)             # bot-top
          end
      end
      Q_proj
end


""" Proyection of an observable in the continuum basis in the c-c space
        ⟨c,k|O|c',k⟩ = ∑_{Qβ,Q'β'}  u_{Qβ,c}^*  u_{Q'β',c'} O_{Qβ,Q'β'}
    with O_{Qβ, Q'β'} = ⟨Q,β,k|O|Q'β'k⟩  matrix element of the operator O in the continuum basis.
    β is the sublattice. c = {1,2,3,4} c-electron dof of the Heavy Fermion model"""
function Oc(p, obs, k)
    P_op = P_operator(p, k)
    Q_op = Q_operator(p, k)
    F = eigen(P_op-Q_op) # P-Q is the proyector into the 4 c orbitals
    sorted_indices = sortperm(abs.(F.values), rev=true)[1:4] # select the 4 eigenvectors with eigenvalue exactly 1 there should be 4 corresponding to the c-electrons
    # end   
    u_tilde = [F.vectors[:, s] for s in sorted_indices]
    c_space = zeros(ComplexF64, 4, 4)
    hd = size(obs,1) ÷ 2
    uuop(c,cp) = (conj.(u_tilde[c])' * transpose(u_tilde[cp])) .* obs
    for c in 1:4
        for cp in 1:4
            aux_mat = uuop(c,cp)
            for (i, gi) in enumerate(Gs) # dimension of the k_space
                for (j, gj) in enumerate(Gs)
                    c_space[c,cp] += sum(aux_mat[2i-1:2i, 2j-1:2j]+aux_mat[2i-1+hd:2i+hd, 2j-1+hd:2j+hd]
                                        +aux_mat[2i-1+hd:2i+hd, 2j-1:2j]+aux_mat[2i-1:2i, 2j-1+hd:2j+hd])
                end
            end
            
        end
    end
    return c_space
end

"""
proyection of an obs O in the continuum basis into the subspace of f orbitals (2x2 matrix)
β, βp are sublattice indices of the continuum 
α, αp are the f orbitals 
conj(v_an(p, li, η, k, Q(p, gi, li), β, α)) * OBS * v_an(p, lj, η, k, Q(p, gj, lj), βp, αp)
"""
function Of(p::ParamsBM, obs, k)
    Gs = Gmesh(Gindices(p), rec_vecs(p)) # vectors of the k-star
    dim_klat = Int(1 + 3p.nmax * (1 + p.nmax)) # dimension of the plane-wave basis in the BM model
    dim_ham = 2 * 2 * dim_klat # σ ∘ k_space ∘ τ 
    hd = dim_ham ÷ 2
    f_space = zeros(ComplexF64, 2, 2)
    for (i, gi) in enumerate(Gs) # dimension of the k_space
        for (j, gj) in enumerate(Gs)
            # f_space[α,αp] .+= conj(v_an(p, 1, η, k, Q(p, gi, 1), β, α)) * obs[2i-1:2i-1, 2j-1:2j-1] * v_an(p, 1, η, k, Q(p, gj, 1), βp, αp) +
            vec_obs = [obs[2i-1, 2j-1], obs[2i, 2j], obs[2i, 2j-1], obs[2i-1, 2j], obs[2i-1+hd, 2j-1+hd], obs[2i+hd, 2j+hd], obs[2i+hd, 2j-1+hd], obs[2i-1+hd, 2j+hd],
                obs[2i-1, 2j-1+hd], obs[2i, 2j+hd], obs[2i, 2j-1+hd], obs[2i-1, 2j+hd],  obs[2i-1+hd, 2j-1], obs[2i+hd, 2j], obs[2i+hd, 2j-1], obs[2i-1+hd, 2j]]
            f_space[1,1] += transpose(f_ααp(p, k, gi, gj, p.ν, 1, 1)) * vec_obs
            f_space[2,2] += transpose(f_ααp(p, k, gi, gj, p.ν, 2, 2)) * vec_obs
            f_space[1,2] += transpose(f_ααp(p, k, gi, gj, p.ν, 1, 2)) * vec_obs
            f_space[2,1] += transpose(f_ααp(p, k, gi, gj, p.ν, 2, 1)) * vec_obs
        end
    end
    return f_space
end

f_ααp(p, k, gi, gj, η, α, αp) = vcat(f_ααp(p, 1, 1, k, gi, gj, η, α, αp), 
                                     f_ααp(p, 2, 2, k, gi, gj, η, α, αp), 
                                     f_ααp(p, 1, 2, k, gi, gj, η, α, αp),
                                     f_ααp(p, 2, 1, k, gi, gj, η, α, αp))

f_ααp(p, l, lp, k, gi, gj, η, α, αp) = 
    [conj(v_an(p, l, η, k, Q(p, gi, l), 1, α))  * v_an(p, lp, η, k, Q(p, gj, lp), 1, αp),
     conj(v_an(p, l, η, k, Q(p, gi, l), 2, α))  * v_an(p, lp, η, k, Q(p, gj, lp), 2, αp),
     conj(v_an(p, l, η, k, Q(p, gi, l), 2, α))  * v_an(p, lp, η, k, Q(p, gj, lp), 1, αp),
     conj(v_an(p, l, η, k, Q(p, gi, l), 1, α))  * v_an(p, lp, η, k, Q(p, gj, lp), 2, αp)]


Q(p::ParamsBM, gi, l) = Q(gi, l, rec_vecs(p))
Q(rcp::Rec_vecs, gi, l) = gi - ifelse(l == 1, rcp.κ1, rcp.κ2) # Q vectors for the top and bottom sublattices.

vv_an_sublat(p, li, lj, η, k, gi, gj) = [vv_an(p, li, lj, η, k, Q(p, gi, li), Q(p, gj,lj), 1, 1)  vv_an(p, li, lj, η, k, Q(p, gi, li), Q(p, gj, lj), 1, 2);
                                     vv_an(p, li, lj,  η,  k, Q(p, gi,li), Q(p, gj,lj), 2, 1) vv_an(p, li, lj, η, k, Q(p, gi,li), Q(p, gj,lj), 2, 2)]

vv_an(p, l, lp, η, k, Q, Qp, α, β) = conj(v_an(p, l, η, k, Q, α, 1)) * v_an(p, lp, η, k, Qp, β, 1) + 
                                     conj(v_an(p, l, η, k, Q, α, 2)) * v_an(p, lp, η, k, Qp, β, 2)


""" Intregral of the ωij in Eq S45. Supplementary Heavy fermion model, 
l is the layer index which determines which G's can be passed as arguments.
Remember that in the continuum model are two sets of momentum transfers corresponding to the G-lattices in the top and bottom layer
Test: e.g. v11_analytical(pbm, 1,1,[0.,0.],[0.,0.]) ===  v(pbm, [0,0.], 1, 1, ω11).

/ sqrt(3√3/2) * (a0/(2pi*1.05/(2*360))) is the unitcell area. Note that the parameters only work for θ = 1.05° otherwise you'll need
to change α1/2 λ1/2, ...
"""
function v_an(p, l, η, k, Q, α, β)
    if α == β
        if α == 1
            v11_analytical(p, l, η, k, Q)
        elseif α == 2
            v22_analytical(p, l, η, k, Q)
        end
    else
        if α == 1
            v12_analytical(p, l, η, k, Q)
        elseif α == 2
            v21_analytical(p, l, η, k, Q)
        end
    end
end

v11_analytical(p, l, η, k, Q) = α1 * λ1(p) *  √(2π) * exp(1im* π/4 * l * η) * exp(-1/2 * λ1(p)^2 * norm(k-Q)^2) / (sqrt(3√3/2) * (a0/(2pi*1.05/(2*360))))
v22_analytical(p, l, η, k, Q) = α1 * λ1(p) *  √(2π) * exp(-1im* π/4 * l * η) * exp(-1/2 * λ1(p)^2 * norm(k-Q)^2) /  (sqrt(3√3/2) * (a0/(2pi*1.05/(2*360))))
v21_analytical(p, l, η, k, Q) = -l * α2 * λ2(p)^2 * √(2π) * exp(1im * π/4 * l * η) * exp(-1/2 * λ2(p)^2 * norm(k-Q)^2) * (η*(k-Q)[2] -1im * (k-Q)[1]) /  (sqrt(3√3/2) * (a0/(2pi*1.05/(2*360))))
v12_analytical(p, l, η, k, Q) =  l * α2 * λ2(p)^2 * √(2π) * exp(-1im * π/4 * l * η) * exp(-1/2 * λ2(p)^2 * norm(k-Q)^2) * (-η*(k-Q)[2] -1im * (k-Q)[1]) /  (sqrt(3√3/2) * (a0/(2pi*1.05/(2*360))))