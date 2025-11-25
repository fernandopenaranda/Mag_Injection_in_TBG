function selfconsistency(p, n; iterations = 10, tol = 0.1e-1 ,α = 0.3, kw...)
    if p.twovalleys == false && p.twovalleystwospins == false
        selfconsistency_onevalley(p, n, iterations, tol, α; kw...)
    elseif p.twovalleys == true
        selfconsistency_twovalleys(p, n, iterations, tol, α; kw...)
    elseif p.twovalleystwospins == true
        selfconsistency_twovalleystwospins(p, n, iterations, tol, α; kw...)
    end
end

function selfconsistency_onevalley(p, n, iterations, tol, α; kw...)
    nf = 1  # seed neutrality
    p = ParamsHF(p, n = n)
    hdim = ham_matrix_size(p)
    println("it: ")
    err = 1; it = 1; μ = 0.0
    while abs(err) > tol && it < 15
        print(" $(it).. ")
        (nf, μ, err) = chargedistribution(nf, p, n, hdim, α; kw...)
        it += 1
    end
    return real(nf), μ, n
end

""" REMARK important the chemical potential corresponding to a full filling of n must not see the conduction bands
in that case we would have to include the O^c and O^cf matrices. """
function selfconsistency_twovalleys(p, n, iterations, tol, α; kpoints = 101)
    of = [1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0] # seed neutrality in the orbital and valley subspace
    # the two f orbitals in the top valley at neutrality are occupied
    p = ParamsHF(p, n = n)
    hdim = ham_matrix_size(p)
    println("it: ")
    err = 1; it = 1; μ = 0.0
    kmesh = sixthMBZ(p, kpoints = kpoints) #sixth
    while abs(err) > tol && it < 30
        print(" $(it).. ")
        (of, μ, err) = chargedistribution(of, p, n, hdim, α, kmesh)
        display(of)
        it += 1
    end
    display(real(of))
    println("μ = ", μ, "n = ", n)
    return  hf_plotbands(ParamsHF(p, μ = μ), of)
end


#=
########################################################################################################

                                PATH FOR PHASE COMPETITION

########################################################################################################
=#

""" the seed can be either specified externally or use the halffillng guess
for the α phase this seed works great. This seed implictly assumes that there is 
no valley coherence. 
# note that we have imposed that Of is diagonal meaning that there 
# is no valley coherence.
# return  hf_plotbands(ParamsHF(p, μ = μ), of)"""
selfconsistency_twovalleystwospins(sym, p, n,  iterations, tol, α, kmesh, evals) = 
    selfconsistency_twovalleystwospins(p, n,[1, 1, 0, 0, 1, 1, 0, 0],
        iterations, tol, α, kmesh, evals)

function selfconsistency_twovalleystwospins(sym, p, n, ofvector::Vector, iterations, tol, α, kmesh, evals)
    of = ofvector
    p = ParamsHF(p, n = n)
    hdim = ham_matrix_size(p)
    println("ite: ")
    err = 1; it = 1; μ = 0.0
    while abs(err) > tol && it < 15
        print(" $(it).. ")
        μ = compute_mu(p, n, of, kmesh, hdim) 
        (of, err) = chargedistribution(sym, of, p, n, μ, hdim, α, evals)
        it += 1
    end
    println("μ = ", μ, "  n = ", n)
    return of, compute_mu(p, n, of, kmesh, hdim), n 
end

"""
here we use several convergence methods to reach optimal convergence in of.
We combine a global low resolution and a local high resolution method on A
constrained interval obtained from the global step.
it calls `selfconsistency_twovalleystwospins_optim`
Update global works better in all cases ignore the local path
"""
function selfconsistency_twovalleystwospins_optim_adaptive(sym, p::ParamsHF, n, ofvector::Vector, 
        config::SelfConsistency_config)
    # p = ParamsHF(p, n = n)
    println("ite: ")
    start_time = Dates.now()
    optim_of, mu, n, sampling_spacing, status = selfconsistency_twovalleystwospins_optim(
        sym, p, n, ofvector, config.global_int_evals; points = config.global_points, 
        max_evals = config.global_max_evals, 
        method = :global, penalty_type = :fixed_n, penalty_energy = config.global_penalty)
    
    upper_bound, lower_bound = upper_and_lower_bounds(sym, optim_of, sampling_spacing)
    count = 0
    mun = 0
    for i in 1:config.max_count  # local minima loop to make it closer
        if status == :MAXEVAL_REACHED
            optim_of, mun, n, sampling_spacing, status = selfconsistency_twovalleystwospins_optim(
                sym, p, n, optim_of, config.local_int_evals; points = config.local_points, 
                max_evals = config.local_max_evals, upper_bound = upper_bound, lower_bound = lower_bound, 
                method = :local, penalty_type = :of_convergence, penalty_energy = config.local_penalty)
            upper_bound, lower_bound = upper_and_lower_bounds(sym, optim_of, sampling_spacing)
            count += 1
            println("Sampling spacing of it $(count) is: $(round.(sampling_spacing, digits = 4))")
        end
    end
    mu = mun
    end_time = Dates.now()  # Record end time 
    elapsed_time = end_time - start_time
    println(" ")
    println(" TIME: $(elapsed_time)")

    println("second convergence...")
    tol = 1e-5; it = 1; μ = 0.0; err = 1
    while abs(err) > tol && it < 15
        print(" $(it).. ")
        mu = compute_mu_experimental(p, n, optim_of, ham_matrix_size(p), points = config.global_points)
        err = abs.(sum(optim_of))
        optim_of = Of(p, mu, diagm(optim_of), evals = config.global_int_evals)
        println("It: $(it) with Of: $(round.(optim_of, digits = 4))")
        err -= abs(sum(optim_of))
        it += 1
    end
    return optim_of, mu, n, sampling_spacing #type instability in sampling spacing
end

function selfconsistency_twovalleystwospins_optim_random(sym, p::ParamsHF, n, ofvector::Vector, 
    config::SelfConsistency_config_random)
    start_time = Dates.now()
    optim_of = rand(8)
    rand_seed = similar(ofvector)
    calc_of = similar(ofvector)
    store_optim_of = similar(rand_seed)
    store_mu = 1000
    e = 1e8
    for i in 1:config.rand_calls
        tt = Dates.now() - start_time
        if tt.value/1000/3600 <  5.0 # 10 hours calculatiom 
            println("Progress: $((i-1)/config.rand_calls * 100)%")
            rand_seed .= rand(8)
            println("0.. Rand seed: $(config.rand_calls)")
            calc_of, calc_mu, calc_e = naif_self_consistency_loop(p, n, rand_seed, config)
            if calc_e < e 
                store_mu = calc_mu
                store_optim_of .= calc_of
                e = calc_e
            end
        else println("time limit") end
    end
    println(" TIME: $(Dates.now() - start_time)")
    return store_optim_of, store_mu, n, 0.0
end

function naif_self_consistency_loop(p, n, optim_of, config)
    err = 1
    mu = 0.0
    e = 0.0
    start_time = Dates.now()
    for it in 1:config.conv_iterations
        tt = Dates.now() - start_time
        if tt.value/1000/3600 <  1 # 1 hours calculatiom 
            if abs(err) > config.tol
                print(" $(it).. ")
                mu, optim_of, err = naif_self_consistency(p, n, optim_of, config)
                println("It: $(it) with Of: $(round.(optim_of, digits = 4))")
            end
        else nothing end
    end
    e = total_energy(p, optim_of, mu, evals = config.int_evals)
    println("energy: $(e)")
    return optim_of, mu, e

end


# function VPvsQAH_self_consistency_loop(p, n, optim_of, config)
#     err_vp = 1
#     err_qah = 1
#     mu_vp = 0.0
#     e_vp = 0.0
#     mu_qah = 0.0
#     e_qaj = 0.0
#     vp_test = [1, 1, -1, -1, 1, 1, -1, -1]
#     start_time = Dates.now()
#     for it in 1:config.conv_iterations
#         tt = Dates.now() - start_time
#         if tt.value/1000/3600 <  1 # 1 hours calculatiom 
#             if abs(err) > config.tol
#                 print(" $(it).. ")
#                 mu, optim_of, err = naif_self_consistency(p, n, optim_of, config)
#                 if sum(optim_of .* vp) < 0.01
#                 println("It: $(it) with Of: $(round.(optim_of, digits = 4))")
#             end
#         else nothing end
#     end
#     e = total_energy(p, optim_of, mu, evals = config.int_evals)
#     println("energy: $(e)")
#     return optim_of, mu, e

# end

function naif_self_consistency(p, n, optim_of, config)
    mu = compute_mu_experimental(p, n, optim_of, ham_matrix_size(p), points = config.points)
    err = abs.(sum(optim_of))
    optim_of = Of(p, mu, diagm(optim_of), evals = config.int_evals)
    err -= abs(sum(optim_of))
    return mu, optim_of, err
end


function generate_vp_aqh_solutions(s::Self_consistent_data)
    s2 = fix_valley_convention_ofmats(s, method = "1");
    s3 = fix_spin_convention_ofmats(s2)
    s4 = fix_spin_convention2_ofmats(s3)
    s5 = permute_vp_aqh_character(s4)

    aux_vp = [zeros(Float64, 8) for i in 1:length(s4.ofmats)]
    aux_qah = [zeros(Float64, 8) for i in 1:length(s4.ofmats)]

    for i in 1:length(s.ofmats)
        # println(s4.ofmats[i][5]-s4.ofmats[i][3], "  ", s5.ofmats[i][5]-s5.ofmats[i][3]  )
        if abs(s4.ofmats[i][5]-s4.ofmats[i][3]) > 1e-2
            if s4.ofmats[i][5]-s4.ofmats[i][3] > 0
                aux_vp[i][:]  .= s4.ofmats[i]  
                aux_qah[i][:] .= s5.ofmats[i]  
            else 
                aux_vp[i][:]  .= s5.ofmats[i]  
                aux_qah[i][:] .= s4.ofmats[i]  
            end    
        else 
            aux_vp[i][:] .= s4.ofmats[i]  
            aux_qah[i][:] .= s4.ofmats[i]  
           
        end
    end
    return  Self_consistent_data(aux_vp, s.mus, s.es, s.ns),  Self_consistent_data(aux_qah, s.mus, s.es, s.ns)
end

permute_vp_aqh_character(s::Self_consistent_data) =
    Self_consistent_data([change_valley(o) for o in s.ofmats], s.mus, s.es, s.ns)

change_valley(o) = [o[1],o[2],o[5],o[6],o[3],o[4],o[7],o[8]]

# function select_degen_solution(s::Self_consistent_data)
#     aux = similar(s.ofmats)
#     vp 
#     for i in 1:length(s.ofmats)
#         if s.ofmats[i] 



# end


"""
upper and lower bounds for the local convergence procedure around
the minima found in the global minima finding procedure. Note that
η is the minimum average spacing between sampling points in the global
algorithm
"""
function upper_and_lower_bounds(sym, of, η = 0.2)
    if sym == :sym
        upper_bound = [of[1]] .+ η
        lower_bound = [of[1]] .- η
    elseif sym == :α
        upper_bound = [maximum(of), minimum(of)] .+ η
        lower_bound = [maximum(of), minimum(of)] .- η
    elseif sym == :β
        symmetrizer!(:β, of)
        nl = of
        upper_bound =  nl .+ η
        lower_bound =  nl .- η
    else nothing end
    for i in 1:length(upper_bound) # avoid boundary issues
        if lower_bound[i] < 0
            lower_bound[i] = 0
        elseif lower_bound[i] >1
            lower_bound[i] = 1
        elseif upper_bound[i] > 1
            upper_bound[i] = 1
        elseif upper_bound[i] < 0
            upper_bound[i] = 0
        else nothing end
    end
    println("   |  up_bound:  ", round.(upper_bound, digits = 2), "   |  low_bound:  ", round.(lower_bound, digits = 2))
    
    return upper_bound, lower_bound
end

function selfconsistency_twovalleystwospins_optim(sym, p, n, ofvector::Vector, evals; kw...)
    hdim = ham_matrix_size(p)
    points = get(kw, :points, 10)
    optimize_func(of) = oferr_optim(sym, of, p, n, hdim, evals; kw...)[1]
    println("starting the optimization: n = $(n) with points = $(points)")
    result = optimize_sym(sym, optimize_func, ofvector; delete!(Dict(kw),:points)...)

    # result = optimize_sym(sym, optimize_func, ofvector; process_keywords([:points]; kw...)...)
    optim_of = of_optimized(sym, result)
    μ = compute_mu_experimental(p, n, optim_of, hdim, points = points)
    #n_comp = compute_n_experimental(p, optim_of, μ, ham_matrix_size(p), points = points)
    println("optimized energy: $(result.e) || vars: $(round.(optim_of, digits = 3))")#  || n: $n_comp ")

    # println("optimized energy: $(result.e)  ||  μ = $(round.(μ, digits = 3)) || vars: $(round.(optim_of, digits = 3))")#  || n: $n_comp ")
    return optim_of, μ, n, result.average_sampling_diff, result.status
end


"""
here i set up the optimization scheme instead of the convergence strategy
passed an of of a given `sym` phase, we determine the μ at a given filling `n`.
Then we determine the total_energy as the function of `of` and `μ` that we want 
to optimize. Fixing the filling not only constrains μ but also the trace of `of`
because of this a lagrange multiplier is added constraining the variation of the 
solution. alternative:

# mu2 = compute_mu(p, n, of2 .* ones(8), kmesh, hdim)
# return e + penalty * abs(mu-mu2) 
"""
function oferr_optim(sym, of, p, n, hdim, evals; kw...)
    if sym == :sym
        oferr_opt(:sym, of .* ones(8), p, n, hdim, evals; kw...)
    elseif sym == :α
        oferr_opt(:α, [of[1], of[1], of[2], of[2], of[1], of[1], of[2], of[2]], p, n, hdim, evals; kw...)
    elseif sym == :β
        oferr_opt(:β, of, p, n, hdim, evals; kw...)
    end
end

#mu and n estimation standard and adatptive methods.

"given a density matrix of computes the position of the chemical potential 
for a given filling n,
# considers the 4 lowest eigenvalues, 
#do not cross conduction bands"
function compute_mu(p, n, of, kpoints::Vector{SVector{2, T}}, hdim) where  T <: Real 
    data = zeros(Float64, 8length(kpoints))
    for (i, k) in enumerate(kpoints)
        eigs = real.(sc_eigs(p, k, diagm(of)).values)
        data[8i-7:8i] .= eigs[hdim÷2-3:hdim÷2+4] 
    end
    ind = Int(round(n * length(kpoints)))
    auxmu = sort(data)
    return auxmu[ifelse(n != 4, 4length(kpoints) + ind + 1, 4length(kpoints) + ind)] 
end

"""
IMPORTANT I ASSUME THST THERE IS C3 Symmetry OTHERWISE CONSIDER MBZ INSTEAD OF THIRDMBZ FOR THE FIRST STEP
I Compute the chemical potential adaptively,
"""
function compute_mu_experimental(p, n, of, hdim::Number; points = 50, mu_tol = 0.05)::Float64 
    kpoints = sixthMBZ(p, kpoints = points)
    eigs = zeros(Float64, hdim)
    dataes = zeros(Float64, 8length(kpoints))
    mun = compute_mu_experimental_sub(dataes,  eigs, p, n, of, kpoints, hdim)
    return mun
end

function compute_mu_experimental_sub(dataes, eigs, p, n, of, kpoints::Vector{Vector{Float64}}, hdim)
    i = 1
    for k in kpoints
        sc_eigvals!(eigs, p, k, diagm(of))
        dataes[8i-7:8i] .= eigs[hdim÷2-3:hdim÷2+4] 
        i += 1
    end
    return dataes[sortperm(dataes)][4length(kpoints) + Int(round(n * length(kpoints)))]
end

""" returns mu for a given n and of and also returns the computed n """
function compute_mu_and_n_experimental(p, n, of, hdim::Number; points = 50, mu_tol = 0.05)::Tuple{Float64, Float64}
    kpoints = sixthMBZ(p, kpoints = points)
    eigs = zeros(Float64, hdim)
    dataes = zeros(Float64, 8length(kpoints))
    mun, n = compute_mu_and_n_experimental_sub(dataes,  eigs, p, n, of, kpoints, hdim)
    return mun, n
end

function compute_mu_and_n_experimental_sub(dataes, eigs, p, n, of, kpoints::Vector{Vector{Float64}}, hdim)
    i = 1
    for k in kpoints
        sc_eigvals!(eigs, p, k, diagm(of))
        dataes[8i-7:8i] .= eigs[hdim÷2-3:hdim÷2+4] 
        i += 1
    end
    mu = dataes[sortperm(dataes)][4length(kpoints) + Int(round(n * length(kpoints)))]
    ind = findmin(abs.(dataes) .- mu)[2]
    n = (8*length(dataes[dataes .< mu])/8length(kpoints))-4
    return mu, n

end

function compute_mu_experimental_sub(dataes, dataks, eigs, p, n, of, kpoints::Vector{Vector{Float64}}, hdim)
    i = 1
    for k in kpoints
        eigs .= real.(sc_eigs(p, k, diagm(of)).values)
        dataes[8i-7:8i] .= eigs[hdim÷2-3:hdim÷2+4] 
        dataks[8i-7:8i] .= fill(k,8)
        i += 1
    end
    ind = Int(round(n * length(kpoints)))
    indices = sortperm(dataes)
    mu = dataes[indices][4length(kpoints) + ind]
    kt = dataks[indices][4length(kpoints) + ind]
    return mu, kt
end

"""
given of and mu it determines the occupation of the flat bands (combination of f and c orbitals).
Note that this breaks if the flat bands are not inside a gap.
The remote bands for every k should be above the flat bands.
"""
function compute_n(p, of, mu, kpoints, hdim)
    data = zeros(Float64, 8length(kpoints))
    for (i, k) in enumerate(kpoints)
        eigs = real.(sc_eigs(p, k, diagm(of)).values)  # DIAGONALIZATION 1 
        data[8i-7:8i] .= eigs[hdim÷2-3:hdim÷2+4] 
    end
    auxmu = data[data .< mu]
    return (8*length(auxmu)/8length(kpoints))-4
end

function compute_n_experimental(p, of, mu, hdim::Number; points = 50, mu_tol = 0.005)::Float64
    nkpoints = func(p, kpoints = points) # init
    kpoints = MBZ(p, kpoints = points)
    nn, kt = compute_n_experimental(p, of, mu, nkpoints, hdim) # compute mu such that the filling is n
    count = 1
    n = 100
    while abs.(n - nn)/abs(nn) > mu_tol && count < 10 
        nkpoints =  [k ./ (10*count) + kt for k in kpoints] 
        n = nn
        k0 = kt
        n, kt = compute_n_experimental(p, of, n, nkpoints, hdim)
        count += 1
    end
    println("  |  n: ", nn)
    return nn
end

function compute_n_experimental(p, of, mu, kpoints::Union{Vector{SVector{2, T}},Vector{AbstractVector{T}}}, hdim) where  T <: Real 
    dataes = zeros(Float64, 8length(kpoints))
    dataks = Vector{SVector{2, Float64}}(undef, 8length(kpoints))
    eigs = zeros(Float64, hdim)
    for (i, k) in enumerate(kpoints)
        eigs .= real.(sc_eigs(p, k, diagm(of)).values)
        dataes[8i-7:8i] .= eigs[hdim÷2-3:hdim÷2+4] 
        dataks[8i-7:8i] .= [k, k, k ,k ,k ,k ,k ,k]
    end
    ind = findmin(abs.(dataes) .- mu)[2]
    auxmu = dataes[dataes .< mu]
    n = (8*length(auxmu)/8length(kpoints))-4
    k = dataks[ind]
    return n, k
end
#= 

   CONSTRAINTS (SYMMETRIZERS)

=#

function symmetrizer!(sym, of)
    if sym == :sym
        sym_symmetrizer!(of)
    elseif sym == :α
        α_symmetrizer!(of)
    elseif sym == :β
        β_symmetrizer!(of)
    end
end

"""enforce the order parameter to pick a specific (energetically degenerate) solution
within the irrep that defines it. This respects all the symmetries of the HF or continuum model"""
function sym_symmetrizer!(of::Vector)
    a = mean(of)
    if a > 1
        a = 1
    elseif a < 0 
        a = 0
    end
    of .= a
    
    return nothing
end

"""enforce the order parameter to pick a specific (energetically degenerate) solution
within the irrep that defines it. α constrains implying there is SU2g symmetry """
function α_symmetrizer!(of::Vector)
    a = mean(of[[1,5,2,6]])
    b = mean(of[[3,7,4,8]])
    if a > 1
        a = 1
    elseif b > 1
        b = 1
    elseif a < 0
        a = 0
    elseif b < 0
        b = 0
    end
    if a > b
        of[[1,2,5,6]].= a
        of[[3,4,7,8]].= b
    else
        of[[1,2,5,6]].= b
        of[[3,4,7,8]].= a
    end
    return nothing
end

"""enforce the order parameter to pick a specific (energetically degenerate) solution
within the irrep that defines it. Rhozsigmaz. 8 free parameters"""
function β_symmetrizer!(of::Vector)
    of .= of
end

function chargedistribution(sym, nf, p, n, mu, hdim, α, evals)
    if p.twovalleys == false && p.twovalleystwospins == false
        chargedistribution_onevalley(nf, p, n, mu, hdim, α )
    elseif  p.twovalleys == true
        chargedistribution_twovalleys(nf, p, n, mu, hdim, α)
    elseif  p.twovalleystwospins == true
        chargedistribution_twovalleystwospins(sym, nf, p, n, mu, hdim, α,  evals)
    end
end

""" standard convergence strategy """
chargedistribution_twovalleystwospins(sym::Symbol, of, p, n, mu, hdim, α, evals) =
     oferr(sym, of, p, n, mu, hdim,  evals)

"""
computes the density matrix corresponding to a filling n. Core of the 
self consistency.
"""
function oferr(sym::Symbol, of, p, n, mu, hdim, α, evals)
    symmetrizer! = ifelse(sym == :sym, sym_symmetrizer!, ifelse(sym == :α, 
        α_symmetrizer!, β_symmetrizer!))
    new_of = Of(p, mu, diagm(of), evals = evals) # DIAGONALIZATION
    symmetrizer!(new_of, n)
    err = maximum(abs.(of - new_of))
    of = (1-α) .* new_of + α * of
    symmetrizer!(of, n)
    return (of, err)
end

function chargedistribution_onevalley(nf, p, n, hdim, α, kpoints)
    data = zeros(Float64, 2length(kpoints)); 
    for (i, k) in enumerate(kpoints)
        eigs = real.(hartree_eigs(ParamsHF(p, nf = nf), k)[1])
        data[2i-1:2i] .= eigs[hdim÷2:hdim÷2+1]
    end
    ind = Int(round(n/4 * length(kpoints)))
    auxmu = sort(data)
    mu = auxmu[ifelse(ind == -length(kpoints), 1, length(kpoints) + ind)]

    nf_aux = 0.0
    meshdim = length(kpoints)
    for k in kpoints
        nf_aux += fdensityink(ParamsHF(p, nf = nf), k, mu, meshdim)
    end
    nf = (1-α) * nf + α * nf_aux
    println("  nf = ", real(nf), " error: ", real(nf_aux-nf))
    return (nf, mu, real(nf_aux-nf))
end


function chargedistribution_twovalleys(of, p, n, hdim, α, kpoints)
    println("points in the mesh: ", kpoints)
    data = zeros(Float64, 4length(kpoints))
    for (i, k) in enumerate(kpoints)
        eigs = real.(sc_eigs(p, k, of; kw...).values)
        data[4i-3:4i] .= eigs[hdim÷2-1:hdim÷2+2] 
    end
    ind = Int(round(n * length(kpoints)))
    auxmu = sort(data)
    mu = auxmu[ifelse(n != 2 , 2length(kpoints) + ind + 1, 2length(kpoints) + ind)] # mu for a given nf
    new_of = Of(p, mu, of, kpoints)
    of = (1-α) .* new_of + α * of
    new_of_trace = trace(new_of) # new trace
    old_of_trace = trace(of)
    println("μ: ", mu," [meV]",  " err: ", sum(abs.(of-new_of)))
    return (of, mu, sum(abs.(of-new_of)))
end


"expected value of the h_hartree operator at a given k over occupied states. 
# this 4 assumes that the two valleys contribute equally as well as the two spins 
# if there is spin/valley polarization or coherence with broken TRS it won't be simply a 4 
# you should go over the two valleys separately i thing to compute ρG
The two meshizerenorm is not to count 6 times the Gamma point
# hartree_eigs return the eigs ane corresponding eigvecs ordered by eigs
"
function chargedistribution(ρG::Number, mu, k, meshdim, p, hartree_term, hdim; kw...)
    sum = 0.0
    meshsizerenorm = (6 * meshdim - 5)/ifelse(k == [0.,0.], 6, 1)               #
    eigs, eigvecs= hartree_eigs(ParamsHF(p, ρG = real(ρG)), k; kw...)
    for nindex in 1:hdim  
        if eigs[nindex] ≤ mu
            sum += (eigvecs[:, nindex]' * hartree_term * eigvecs[:,nindex])/meshsizerenorm 
        else nothing end
    end
    return sum
end

function MBZ(p; kpoints = 42) #new
    kpoints = kpoints ÷ 6
    g1, g2 = bravais_vectors(p)
    κ1 = κ(g1, g2)
    d = norm(g1) / kpoints
    u1 = κ1 / norm(κ1)
    u2 = (g1 + g2)/2 / norm((g1 + g2)/2)
    u3 = g1 /norm(g1)

    shape = hex_shape(p, norm(κ1), d)
 
    mesh = Vector{Vector{Float64}}()
    ik = -kpoints * d * sqrt(3)
    for i in -kpoints:kpoints
        push!(mesh, u1 * ik)
        jk = -kpoints * d
        for j in -kpoints:kpoints
            push!(mesh, u1 * ik .+ jk * u2)
            push!(mesh, u1 * ik .+ jk * u3)
            jk += d
        end
        ik += d * sqrt(3)
    end
    kgrid = [mesh[i] for i in 1:size(mesh, 1)]
    inside = [PolygonOps.inpolygon(m, shape; in=true, on=false, out=false) for m in kgrid]
    pos = findall(x -> x != 0, inside)
    return uniquetol(kgrid[pos])
end

function uniquetol(A::AbstractArray{T}; kwargs...) where T
    S = Vector{T}()
    for a in A
         if !any(s -> isapprox(s, a; kwargs...), S)
             push!(S, a)
         end
    end
    return S
end

function thirdMBZ(p; kpoints = 42)
    kgrid = sixthMBZ(p; kpoints = kpoints)
    nkgrid = [[k[1],-k[2]] for k in kgrid]
    return uniquetol(vcat(kgrid, nkgrid))
end


function sixthMBZ(p; kpoints = 42) #new
    g1, g2 = bravais_vectors(p)
    κ1 = κ(g1, g2)
    d = norm(g1) / kpoints
    u1 = κ1 / norm(κ1)
    u2 = (g1 + g2)/2 / norm((g1 + g2)/2)
    u3 = g1 /norm(g1)

    shape = hex_shape(p, norm(κ1), d)
 
    mesh = Vector{Vector{Float64}}()
    ik = 0
    for i in 0:kpoints
        push!(mesh, u1 * ik)
        jk = d
        for j in 0:kpoints
            push!(mesh, u1 * ik .+ jk * u2)
            push!(mesh, u1 * ik .+ jk * u3)
            jk += d
        end
        ik += d * sqrt(3)
    end
    kgrid = [mesh[i] for i in 1:size(mesh, 1)]
    inside = [PolygonOps.inpolygon(m, shape; in=true, on=false, out=false)
        for m in kgrid]
    pos = findall(x -> x != 0, inside)
    return kgrid[pos]
end

function plotsixthMBZ(p; kpoints = 42)
    kgrid = sixthMBZ(p,kpoints = kpoints)
    fig = Figure(resolution = (600, 600));
    ax = Axis(fig[1, 1], aspect = 1.0)   
    for k in kgrid
        scatter!(ax, (k[1],k[2]), color = :blue, markersize = 3)
    end   
    fig
end

function plotMBZ(p; kpoints = 42)
    G1, G2 = bravais_vectors(p)
    κ1 = κ(G1, G2)
    d = norm(G1) / kpoints
    kgrid = MBZ(p, kpoints = kpoints)
    fig = Figure(resolution = (600, 600));
    ax = Axis(fig[1, 1], aspect = 1.0)   
    shape = hex_shape(p, norm(κ1), d)
    for i in 1:length(shape)
        scatter!(ax, (shape[i][1], shape[i][2]), markersize = 3)
    end
    for k in kgrid
        scatter!(ax, (k[1],k[2]), color = :blue, markersize = 3)
    end   
    fig
end

hex_shape(p, R, d) = rotate.([SA[0,1+d],SA[√3/2*(1+d),1/2*(1+d)],SA[√3/2*(1+d),-1/2*(1+d)],SA[0,-1*(1+d)],
    SA[-√3/2*(1+d),-1/2*(1+d)],SA[-√3/2*(1+d),1/2*(1+d)],SA[0,1*(1+d)]], θ = 0*π/2) .* (R)

rotate(v::Vector{Float64}; θ = 2π/6) = [cos(θ) -sin(θ); sin(θ) cos(θ)] * v 
rotate(v::SVector{2, Float64}; θ = 2π/6) = [cos(θ) -sin(θ); sin(θ) cos(θ)] * v 
inverty(vec) = [vec[1], -vec[2]]

function ρGvsn(p; kws...)
    nlist = collect(-4:0.5:4)
    ρlist = []
    for n in nlist
        push!(ρlist, selfconsistency(p, n; kws...)[1])
    end
    return nlist, ρlist
end

function unique_zeros(vectors)
    adjusted_vectors = [map(x -> x == -0.0 ? 0.0 : x, vec) for vec in vectors]
    return unique(adjusted_vectors)
end

function hartree_ham(p)
    hdim = ham_matrix_size(p)
    hartree_term = spzeros(ComplexF64, hdim, hdim)
    hf_hartree_mod!(hartree_term, ParamsHF(p, nf = 2)) 
    return hartree_term
end

function k_focus(p, nkpoints) # new k values inside the MBZ
    g1, g2 = bravais_vectors(p)
    κ1 = κ(g1, g2)
    inside = [PolygonOps.inpolygon(k, hex_shape(p, norm(κ1), +norm(κ1)/10); in=true, on=false, out=false) for k in nkpoints]
    pos = findall(x -> x != 0, inside)
    return nkpoints[pos]
end

function minimal_distance_to_element(vectors, index)
    ref_vector = vectors[index]  # The specific vector to compare against
    distances = [norm(v - ref_vector) for v in vectors if v != ref_vector]  # Compute distances
    return minimum(distances)  # Return the smallest distance
end

function process_keywords(sym::Symbol; kw...)
    kw_dict = Dict(kw)
    for s in sym
        delete!(kw_dict, s)
    end
    kw0 = ( ;kw_dict...)
end

""" 
Returns a `Self_consistent_data` object with a given convention for the filling of the valleys between one of the 
4 possible `methods` below. 
the two valleys are degenerate. So we choose in the unperturbed solution one fix filling so the injection current
is not oscillating between manifolds with the same dof.
Regarding the shift current there are 4 degenerate filling options:
    (1) First valley +
    (2) First valley -
    (3) First valley + for positive and - for negative
    (4) First valley - for positive and + for negative
    (5) + - + - + - + - as a function of the filling. To do alternating pattern this would break phs

# la injection es insensitive al spin por lo que no hago nada con el
GENERAL REMARKS ON DEGENERACIES OF H(of):
# [1 ,2 ,3, 4, 5, 6, 7, 8] = [3, 4, 1, 2, 7, 8, 5, 6] = [1, 2, 7, 8, 3, 4, 5, 6] = [5, 6, 3, 4, 1, 2, 7, 8] = [5, 6, 7, 8, 1, 2, 3, 4]
# [1, 2, 3, 4, 5, 6,  7, 8] ≠ [3, 4, 1, 2, 5, 6, 7, 8] Esta operación lleva la solución TRS sym a la broken-TRS sym
"""
function fix_valley_convention_ofmats(s::Self_consistent_data; method = "1")
    reshuffled_ofmats = [zeros(8) for i in 1:length(s.ofmats)]
    for i in 1:length(s.ofmats) #populate the valleys in the same fashion
        # que valle está más ocupado 
        if method == "1"
            if sum(s.ofmats[i][[1,2,5,6]]) > sum(s.ofmats[i][[3,4,7,8]])
                reshuffled_ofmats[i] .= s.ofmats[i]
            else
                reshuffled_ofmats[i] .= s.ofmats[i][[3,4,1,2,7,8,5,6]]
            end
        elseif method == "2"
            if sum(s.ofmats[i][[1,2,5,6]]) < sum(s.ofmats[i][[3,4,7,8]])
                reshuffled_ofmats[i] .= s.ofmats[i]
            else
                reshuffled_ofmats[i] .= s.ofmats[i][[3,4,1,2,7,8,5,6]]
            end
        elseif method == "3"
            if sign(s.ns[i])*sum(s.ofmats[i][[1,2,5,6]]) > sign(s.ns[i])*sum(s.ofmats[i][[3,4,7,8]])
                reshuffled_ofmats[i] .= s.ofmats[i]
            else
                reshuffled_ofmats[i] .= s.ofmats[i][[3,4,1,2,7,8,5,6]]
            end
        elseif method == "4"
            if sign(s.ns[i])*sum(s.ofmats[i][[1,2,5,6]]) < sign(s.ns[i])*sum(s.ofmats[i][[3,4,7,8]])
                reshuffled_ofmats[i] .= s.ofmats[i]
            else
                reshuffled_ofmats[i] .= s.ofmats[i][[3,4,1,2,7,8,5,6]]
            end
        end
    end
        Self_consistent_data(reshuffled_ofmats, s.mus, s.es, s.ns)
end

function fix_spin_convention_ofmats(s::Self_consistent_data; method = "1")
    reshuffled_ofmats = [zeros(8) for i in 1:length(s.ofmats)]
    for i in 1:length(s.ofmats) #populate the spinsin the same fashion
            if sum(s.ofmats[i][[1,2,3,4]]) > sum(s.ofmats[i][[5,6,7,8]])
                reshuffled_ofmats[i] .= s.ofmats[i]
            else
                reshuffled_ofmats[i] .= s.ofmats[i][[5,6,7,8,1,2,3,4]]
            end
    end
        Self_consistent_data(reshuffled_ofmats, s.mus, s.es, s.ns)
end

function fix_spin_convention2_ofmats(s::Self_consistent_data; method = "1")
    reshuffled_ofmats = [zeros(8) for i in 1:length(s.ofmats)]
    for i in 1:length(s.ofmats) #populate the spinsin the same fashion
            if sum(s.ofmats[i][[1,2]]) > sum(s.ofmats[i][[5,6]])
                reshuffled_ofmats[i] .= s.ofmats[i]
            else
                reshuffled_ofmats[i] .= s.ofmats[i][[5,6,3,4,1,2,7,8]]
            end
    end
        Self_consistent_data(reshuffled_ofmats, s.mus, s.es, s.ns)
end

""" 
Returns a `Self_consistent_data` object with a given convention for the filling of the valleys between one of the 
4 possible `methods` below. 
the two valleys are degenerate. So we choose in the unperturbed solution one fix filling so the injection current
is not oscillating between manifolds with the same dof.
Regarding the shift current there are 4 degenerate filling options:
    (1) First valley +
    (2) First valley -
    (3) First valley + for positive and - for negative
    (4) First valley - for positive and + for negative
    (5) + - + - + - + - as a function of the filling. To do alternating pattern this would break phs

# la injection es insensitive al spin por lo que no hago nada con el
GENERAL REMARKS ON DEGENERACIES OF H(of):
# [1 ,2 ,3, 4, 5, 6, 7, 8] = [3, 4, 1, 2, 7, 8, 5, 6] = [1, 2, 7, 8, 3, 4, 5, 6] = [5, 6, 3, 4, 1, 2, 7, 8] = [5, 6, 7, 8, 1, 2, 3, 4]
# [1, 2, 3, 4, 5, 6,  7, 8] ≠ [3, 4, 1, 2, 5, 6, 7, 8] Esta operación lleva la solución TRS sym a la broken-TRS sym
"""
function fix_valley_convention_ofmats(ofmats, mus, es, ns; method = "1")
    reshuffled_ofmats = [zeros(8) for i in 1:length(ofmats)]
    for i in 1:length(ofmats) #populate the valleys in the same fashion
        # que valle está más ocupado 
        if method == "1"
            if sum(ofmats[i][[1,2,5,6]]) > sum(ofmats[i][[3,4,7,8]])
                reshuffled_ofmats[i] .= ofmats[i]
            else
                reshuffled_ofmats[i] .= ofmats[i][[3,4,1,2,7,8,5,6]]
            end
        elseif method == "2"
            if sum(ofmats[i][[1,2,5,6]]) < sum(ofmats[i][[3,4,7,8]])
                reshuffled_ofmats[i] .= ofmats[i]
            else
                reshuffled_ofmats[i] .= ofmats[i][[3,4,1,2,7,8,5,6]]
            end
        elseif method == "3"
            if sign(ns[i])*sum(ofmats[i][[1,2,5,6]]) > sign(ns[i])*sum(ofmats[i][[3,4,7,8]])
                reshuffled_ofmats[i] .= ofmats[i]
            else
                reshuffled_ofmats[i] .= ofmats[i][[3,4,1,2,7,8,5,6]]
            end
        elseif method == "4"
            if sign(ns[i])*sum(ofmats[i][[1,2,5,6]]) < sign(ns[i])*sum(ofmats[i][[3,4,7,8]])
                reshuffled_ofmats[i] .= ofmats[i]
            else
                reshuffled_ofmats[i] .= ofmats[i][[3,4,1,2,7,8,5,6]]
            end
        end
    end
        return reshuffled_ofmats, mus, es, ns
end

"""
I will start filling the +↑ first, so then I have to decide what to do with the rest.
All combinations are allowed with the symmetry so I'll define different methods for that:
the order is specified by permlist, the filling order of the spin and valleys is given by a 
4 dim array where each element in a 1: list sets which of is the filling order among
1 = +↑,
2 = -↓,
3 = +↑
4 = -↓, e.g.
[1,2,3,4] -> +↑, -↓, +↑, -↓ starting from 0 to 4
[2,1,3,4] -> -↓,+↑, +↑, -↓
"""
function set_filling_order(s::Self_consistent_data; permlist = [1,2,3,4])
    reshuffled_ofmats = [zeros(8) for i in 1:length(s.ofmats)]
    for i in 1:length(s.ofmats)
        positions = sortperm(s.ofmats[i][[1,3,5,7]])
        reshuffled_ofmats[i][[1,3,5,7]] .= s.ofmats[i][[1,3,5,7]][positions][permlist]
    end
    Self_consistent_data(reshuffled_ofmats, s.mus, s.es, s.ns)
end

function tauz(s::Self_consistent_data)
    taus = zeros(Float64,length(s.ofmats))
    for i in 1:length(s.ofmats)
        taus[i] = abs(sum(s.ofmats[i] .* [1, 1, -1, -1, 1, 1, -1, -1]))
    end
    taus
end
        

# fix discontinuities by enforcing maximal TRS breakdown and minimal TRS breakdown for the states labeled as QAH and SVH/VH
function discontinuities_fix(s::Self_consistent_data, gs)
    aux_vp = [zeros(Float64, 8) for i in 1:length(s.ofmats)]
    aux_qah = [zeros(Float64, 8) for i in 1:length(s.ofmats)]
    for i in 1:length(s.ofmats)
        occ_Kup = sum(s.ofmats[i][1:2])
        occ_Kdown = sum(s.ofmats[i][3:4])
        occ_Kpup = sum(s.ofmats[i][5:6])
        occ_Kpdown = sum(s.ofmats[i][7:8])

        # indices of the orbitals blocks with larger occupancies
        r_ind = sortperm([occ_Kup, occ_Kdown, occ_Kpup, occ_Kpdown])
        # now I promote the most unequal occupation of the valleys
        aux_qah[i][:] =vcat([s.ofmats[i][2r_ind[k]-1:2r_ind[k]] for k in 1:4]...)
        # # now I promote the most TRS preserving combination
        aux_vp[i][:] =vcat([s.ofmats[i][2r_ind[k]-1:2r_ind[k]] for k in [1,4,2,3]]...)

    end
    if gs == "vp"
         return  Self_consistent_data(aux_vp, s.mus, s.es, s.ns)
    else
        return Self_consistent_data(aux_qah, s.mus, s.es, s.ns)
     end
end