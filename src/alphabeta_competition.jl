# ----------------------------------------------------------------------------------------------- #
#                                                                                                 #
#                                      SELF-CONSISTENCY                                           #
#                                                                                                 #          
# ----------------------------------------------------------------------------------------------- #

""" 
computes the order parameter σ0τzsz (in continuum basis).
Energetically favorable GS.
Competition α vs β
"""
function β_orderparameter_vsfilling(p, tol, kpoints)
    of_matsnon = []
    μ_arrnon = []
    e_arrnon = []

    of_matsa = []
    μ_arra = []
    e_arra = []

    of_matsb = []
    μ_arrb = []
    e_arrb = []
    nlist = collect(-2:0.5:2)
    count = 0
    for n in nlist
        count += 1
        println("Process: $(count/length(nlist) * 100)%")
     
        ofnon, μnon, enon = non_sc_perturbation(p, n, tol = tol, kpoints = kpoints)
        ofa, μa, ea = α_sc_perturbation(p, n, tol = tol, kpoints = kpoints)
        # ofb, μb, eb = β_sc_perturbation(p, n)

        push!(of_matsnon,ofnon)
        push!(μ_arrnon, μnon)
        push!(e_arrnon, enon)

        push!(of_matsa,ofa)
        push!(μ_arra, μa)
        push!(e_arra, ea)
    end
    return of_matsnon, of_matsa, of_matsb, μ_arrnon, μ_arra, μ_arrb, e_arrnon, e_arra, e_arrb
end


function statistics_selfconsistency(p, n, phase, its; kw...)
    ofs = []
    mus = []
    es = []
    phase_sc_pert = ifelse(phase == :sym,  non_sc_perturbation, 
        ifelse(phase == :α, α_sc_perturbation, β_sc_perturbation))
    println(phase_sc_pert)
    for i in 1:its
        o, mu, e = phase_sc_pert(p, n; kw...)
        push!(ofs, o)
        push!(mus, mu)
        push!(es, e)
    end
    return ofs, mus, es
end

function β_orderparameter_vsfilling_vsU1J(p)
    of_mats = []
    μ_arr = []
    e_arr = []
    nlist = collect(-4:0.5:4)
    ratioJoU = collect(-20:2.5:20)
    count = 0
    for rJoU in ratioJoU
        count += 1
        println("Process: $(count/length(nlist) * 100)%")
        for n in nlist
            ofα, μα, eα = α_sc_perturbation(ParamsHF(p, J = rJoU * p.U1), n)
            ofβ, μβ, eβ = β_sc_perturbation(ParamsHF(p, J = rJoU * p.U1), n)
            if eα < eβ
                push!(of_mats, ofα)
                push!(μ_arr, μα)
                push!(e_arr, eα)
            else
                push!(of_mats, ofβ)
                push!(μ_arr, μβ)
                push!(e_arr, eβ)
            end
        end
    end
    save_order_parameter(p, of_mats, μ_arr, e_arr, nlist, ratioJoU)
    return of_mats, μ_arr, e_arr, nlist, ratioJoU
end

# ----------------------------------------------------------------------------------------------- #
#                                                                                                 #
#                                    Guided filling sweeps                                        #
#                                                                                                 #          
# ----------------------------------------------------------------------------------------------- #

function guided_sweep(sym::Symbol, p::ParamsHF, config::Union{SelfConsistency_config, SelfConsistency_config_random},
    νs = collect(-4:1:4), kws...)
    of_mats = []; μ_arr = []; e_arr = []
    ground_state_sc = ifelse(sym == :sym, non_sc_perturbation, ifelse(sym == :α, α_sc_perturbation, β_sc_perturbation))
    of, μ, e = ground_state_sc(p, νs[1], config; kws...) 
    push_hf!(of_mats, μ_arr, e_arr, of, μ, e)
    if length(νs)>1
        for ν in νs[2:end]
            of, μ, e = ground_state_sc(ParamsHF(p, μ = 0), ν, of, config; kws...)
            push_hf!(of_mats, μ_arr, e_arr, of, μ, e)
        end
    end
    return of_mats, μ_arr, e_arr
end

function push_hf!(v1, v2, v3, e1, e2, e3)
    push!(v1, e1)
    push!(v2, e2)
    push!(v3, e3)
end

# ----------------------------------------------------------------------------------------------- #
#                                                                                                 #
#                                          PHASES                                                 #
#                                                                                                 #          
# ----------------------------------------------------------------------------------------------- #
"""
Groundstate that preserves the H = H0 + HI symmetries
"""
function non_sc_perturbation(p, n, config; kw...)
    of = ones(8) .* (n+4)/8
    return sc_perturbation(:sym, p, n, of, config; kw...)
end
non_sc_perturbation(p, n, of, config; kw...) = sc_perturbation(:sym, p, n, of, config; kw...)
"""
Groundstate that breaks the H = H0 + HI symmetries.
Of α type (with SU(2) of spin at each valley). 
The commonly known as valley polarized phase. Breaks TRS.
But preserves the valley.
"""
function α_sc_perturbation(p, n, config; kw...)
    a0, b0 = (n+4)/8, (n+4)/8
    a = rand(1)[1]/5 
    if a + a0 < 4
        of = [a0, a0, b0, b0, a0, a0, b0, b0] + [a, a, 0, 0, a, a, 0, 0]
    else
        of = [a0, a0, b0, b0, a0, a0, b0, b0] - [0, 0, a, a, 0, 0, a, a]
    end
    of = [0.8,0.8,0,0,0.8,0.8,0,0]
    return sc_perturbation(:α, p, n, of, config; kw...)
end
α_sc_perturbation(p, n, of, config; kw...) = sc_perturbation(:α, p, n, of, config; kw...)
"""
Groundstate that breaks the H = H0 + HI symmetries.
Of β type 
The commonly known as valley and spin polarized phase. Breaks TRS.
But preserves the valley.
"""
function β_sc_perturbation(p, n, config; kw...)  
    a, b = (n+4)/8, (n+4)/8
    of = [a, a, b, b, a, a, b, b] 
   return sc_perturbation(:β, p, n, of, config; kw...)
end
β_sc_perturbation(p, n, of, config; kw...) = sc_perturbation(:β, p, n, of, config; kw...)

function sc_perturbation(sym, p, n, of, config::SelfConsistency_config; kw...)
    method = get(kw, :method, :hybrid) # convergence method 
    func = selfconsistency_twovalleystwospins_optim_adaptive
    of, μ, _ = func(sym, p, n, of, config)
    e = total_energy(p, of, μ, evals = config.global_int_evals)
    println("Energy without penalities:", e)
    return of, μ, e
end

function sc_perturbation(sym, p, n, of, config::SelfConsistency_config_random; kw...)
    method = get(kw, :method, :hybrid) # convergence method 
    func = selfconsistency_twovalleystwospins_optim_random
    of, μ, _ = func(sym, p, n, of, config)
    e = total_energy(p, of, μ, evals = config.int_evals)
    println("Energy without penalities:", e)
    return of, μ, e
end

selfconsistency_aux(sym, p, n, of, kmesh; kw...) =
    selfconsistency_twovalleystwospins(sym, p, n, of, iterations, tol, α, kmesh, evals) # old method not used
# ----------------------------------------------------------------------------------------------- #
#                                                                                                 #
#                                     ORDER PARAMETERS                                            #
#                                       -_-__---__-_-                                             #          
# ----------------------------------------------------------------------------------------------- #
"""
order parameter is: σ0τzs0 
OUR BASIS: 
    1⬆+, 2⬆+, 1⬆-, 2⬆-, 1⬇+, 2⬇+, 1⬇-, 2⬇-"""
function α1_orderparameter(p, mu, of)
    dim = 8
    op = spzeros(ComplexF64, dim, dim)
    op[1:2,1:2] = [1 0; 0 1]
    op[3:4,3:4] = [-1 0; 0 -1]
    op[5:6,5:6] = [1 0; 0 1]
    op[7:8,7:8] = [-1 0; 0 -1]
    return sum(expected_value_op(p, mu, diagm(of), op))
end

function β1_orderparameter(p, ofmat, mus, es, nlist)
    ord_param = []
    for i in 1:length(nlist)
        println(i)
        push!(ord_param, β1_orderparameter(p, mus[i], ofmat[i]))
    end
    return nlist, ord_param
end

""" order parameter is: σ0τzsz 
# OUR BASIS:  1⬆+, 2⬆+, 1⬆-, 2⬆-, 1⬇+, 2⬇+, 1⬇-, 2⬇-
"""
function β1_orderparameter(p, mu, of)
    dim = 8
    op = spzeros(ComplexF64, dim, dim)
    op[1:2,1:2] = [1 0; 0 1]
    op[3:4,3:4] = [-1 0; 0 -1]
    op[5:6,5:6] = [-1 0; 0 -1]
    op[7:8,7:8] = [1 0; 0 1]
    return sum(expected_value_op(p, mu, diagm(of), op))
end

"""
collect the list of order parameters that are odd against PHS and
diagonal. Those of interest see description above.
σzρ0s0, σ0ρ0s0, σzρ0s0, σzρ0s0, σzρ0s0, σzρ0s0
"""
function ph_orderparameter(p, mu, of; kpoints = 101)
    op = spzeros(ComplexF64, 8, 8)

    #PHS f - electrons | 1i*(C2zTP) = -1im * σyτz (τ es valle)
    Uc = -1 .*[0 -im 0 0;
             im 0  0 0;
             0 0 0 im;
             0 0 -im 0]
    uc = [Uc 0I; 0I Uc]
    phs = expected_value_op(p, mu, diagm(of), uc, thirdMBZ(p, kpoints = kpoints))
    obs0 = expected_value_op(p, mu, diagm(of), 1I, thirdMBZ(p, kpoints = kpoints))
    # σzρ0s0
    op[1:2,1:2] = [1 0; 0 -1]
    op[3:4,3:4] = [1 0; 0 -1]
    op[5:6,5:6] = [1 0; 0 -1]
    op[7:8,7:8] = [1 0; 0 -1]
    obs1 = expected_value_op(p, mu, diagm(of), op, thirdMBZ(p, kpoints = kpoints))
    
    # σzρzs0
    op[1:2,1:2] = [1 0; 0 -1]
    op[3:4,3:4] = [-1 0; 0 1]
    op[5:6,5:6] = [1 0; 0 -1]
    op[7:8,7:8] = [-1 0; 0 1]

    obs2 = expected_value_op(p, mu, diagm(of), op, thirdMBZ(p, kpoints = kpoints))

    # σzρ0sz
    op[1:2,1:2] = [1 0; 0 -1]
    op[3:4,3:4] = [1 0; 0 -1]
    op[5:6,5:6] = [-1 0; 0 1]
    op[7:8,7:8] = [-1 0; 0 1]

    obs3 = expected_value_op(p, mu, diagm(of), op, thirdMBZ(p, kpoints = kpoints))
    
    # σzρzsz
    op[1:2,1:2] = [1 0; 0 -1]
    op[3:4,3:4] = [-1 0; 0 1]
    op[5:6,5:6] = [-1 0; 0 1]
    op[7:8,7:8] = [1 0; 0 -1]
    obs4 = expected_value_op(p, mu, diagm(of), op, thirdMBZ(p, kpoints = kpoints))
    return [phs,obs0, obs1, obs2, obs3, obs4]
end

# ----------------------------------------------------------------------------------------------- #
#                                                                                                 #
#                                     E N E R G I E S                                             #
#                                       ==========                                                #
#                                                                                                 #          
# ----------------------------------------------------------------------------------------------- #

""" Total energy of the system under the approximation that Hw and Hv as well as 
Ev and Ew are zero. Note that E_j is included inside the ϵ(k) the eigenvalues
of the diagonalized hamiltonians so they are not included in here. """
total_energy(p::ParamsHF, of, mu; kw...) = 
    occupiedbands_energy(p, diagm(of), mu; kw...) - E_U(p, diagm(of))

""" determines ∑_k ∑_n_occ (ϵ_k + μ).
There are referred with respect to μ = 0 so we can compare between different fillings"""
function occupiedbands_energy(p::ParamsHF, of, mu, ks)  # path with simple sum over the energies
    return occupiedbands_energy(real.(bands(p, ks, of)), mu)/length(ks)
end

function occupiedbands_energy(p, of, mu; evals = 100) #new parth with adaptive integral over the energies. 
    M, xmin, xmax = int_boundaries(p)
    g1, g2 = bravais_vectors(p)
    κ1 = κ(g1, g2)
    shape = hex_shape(p, norm(κ1), 0)
    integrand(k) = k_occupiedbands_energy(p, of, mu, k, shape)
    Δx = [0, xmin[2]]
    Δy = [xmax[1]/2, xmax[2]]
    val, err = Cubature.hcubature(integrand, Δx, Δy;  maxevals=evals);
    return val/( (Δy[2] - Δx[2]) * (Δy[1] - Δx[1]))
end

function k_occupiedbands_energy(p, of, mu, k,shape)
    e = eigvals!(Matrix(hf_valley_spin_hamiltonian(p, [k[1], k[2]], of)))
    sum(e[e .< mu])
end

occupiedbands_energy(es::Matrix, mu) = sum(es[es.<mu])

function refer_integral_to1MBZ(p, k)
    g1, g2 = bravais_vectors(p)
    κ1 = κ(g1, g2)
    M, xmin, xmax = int_boundaries(p)
    shape = hex_shape(p, norm(κ1), 0)
    if PolygonOps.inpolygon(k, shape; in = true, on = false, out = false)
        nothing
    elseif k[2] > 0.0
        k -= g1
    else
        k -= g2
    end
    return k
end

""" Returns E_U over N, being N the number of moire unitcells """
E_U(p::ParamsHF, of) = p.U1/2*((trace(of)-4)^2+8*(trace(of)-4)-trace(of*of)) + 
    0*3p.U2 * ((trace(of)-4)^2 + 8 * (trace(of)-4)) 
function E_J_adhoc(p, of,i,j,k) 
    - p.J*(i*(trace(of)-4)^2+j*8*(trace(of)-4)-k*trace(of*of)) 
end

#########
# Order params
########
"""
plot_four_op_vsfilling(ParamsHF(p, U1 = 10, U2 = 0, J = 6), 
    gssym, gsalpha, gsbeta, evals = 2000, kpoints = 200)
"""
function plot_four_op_vsfilling(p, gs1, gs2, gs3; kws...)
    a1, b1, c1, d1 = four_op_vsfilling(p, gs1; kws...)
    a2, b2, c2, d2 = four_op_vsfilling(p, gs2; kws...)
    a3, b3, c3, d3 = four_op_vsfilling(p, gs3; kws...)
    
    fig = Figure(resolution = (1000,1000))
    ax1 = Axis(fig[1,1], xlabel = "ν ",  ylabel = "sz")
    ax2 = Axis(fig[1,2], xlabel = "ν ",  ylabel = "τz")
    ax3 = Axis(fig[2,1], xlabel = "ν ",  ylabel = "szτz")
    ax4 = Axis(fig[2,2], xlabel = "ν ",  ylabel = "σz")
    
    scatter!(ax1, -4:4, [i for i in a1], color = :gray)
    lines!(ax1, -4:4, [i for i in a1], color = :gray)
    scatter!(ax1, -4:4, [i for i in a2], color = :red)
    lines!(ax1, -4:4, [i for i in a2], color = :red)
    scatter!(ax1, -4:4, [i for i in a3], color = :blue)
    lines!(ax1, -4:4, [i for i in a3], color = :blue)

    scatter!(ax2, -4:4, [i for i in b1], color = :gray)
    lines!(ax2, -4:4, [i for i in b1], color = :gray)
    scatter!(ax2, -4:4, [i for i in b2], color = :red)
    lines!(ax2, -4:4, [i for i in b2], color = :red)
    scatter!(ax2, -4:4, [i for i in b3], color = :blue)
    lines!(ax2, -4:4, [i for i in b3], color = :blue)

    scatter!(ax3, -4:4, [i for i in c1], color = :gray)
    lines!(ax3, -4:4, [i for i in c1], color = :gray)
    scatter!(ax3, -4:4, [i for i in c2], color = :red)
    lines!(ax3, -4:4, [i for i in c2], color = :red)
    scatter!(ax3, -4:4, [i for i in c3], color = :blue)
    lines!(ax3, -4:4, [i for i in c3], color = :blue)

    scatter!(ax4, -4:4, [i for i in d1], color = :gray)
    lines!(ax4, -4:4, [i for i in d1], color = :gray)
    scatter!(ax4, -4:4, [i for i in d2], color = :red)
    lines!(ax4, -4:4, [i for i in d2], color = :red)
    scatter!(ax4, -4:4, [i for i in d3], color = :blue)
    lines!(ax4, -4:4, [i for i in d3], color = :blue)

    ax1.xticks = collect(-4:1:4)
    ax2.xticks = collect(-4:1:4)
    ax3.xticks = collect(-4:1:4)
    ax4.xticks = collect(-4:1:4)
    fig
end


four_op_vsfilling(p, gs; kws...) = four_op_vsfilling(p, [g for g in gs[2]], [diagm(g) for g in gs[1]]; kws...)

function four_op_vsfilling(p, mus, ofs; evals = 1000, kpoints = 100)
    op_sz = []
    op_tauz = []
    op_sztauz = []
    op_sigmaz = [] 
    for i in 1:length(mus)
        a,b,c,d = four_orderparameters(p, mus[i], ofs[i]; evals = evals, kpoints = kpoints)
        push!(op_sz, a)
        push!(op_tauz, b)
        push!(op_sztauz, c)
        push!(op_sigmaz, d)
    end

    return op_sz, op_tauz, op_sztauz, op_sigmaz
end

function four_orderparameters(p, mu, of; evals = 1000, kpoints = 100)
    dim = 8
    op_sz = sparse(diagm([1, 1, 1, 1, -1, -1, -1 ,-1]))
    op_tauz = sparse(diagm([1, 1, -1, -1, 1, 1, -1 ,-1]))
    op_tauzsz = sparse(diagm([1, 1, -1, -1, -1, -1, 1 ,1]))
    op_sigmaz = sparse(diagm([1, -1, 1, -1, 1, -1, 1 ,-1]))
    
    ev_sz = sum(expected_value_op(p, mu, of, op_sz; evals = evals))
    ev_tauz = sum(expected_value_op(p, mu, of, op_tauz; evals = evals))
    ev_tauzsz = sum(expected_value_op(p, mu, of, op_tauzsz; evals = evals))
    ev_sigmaz = sum(expected_value_op(p, mu, of, op_sigmaz; evals = evals))
    return ev_sz, ev_tauz, ev_tauzsz, ev_sigmaz
end