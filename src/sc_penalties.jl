# path minimizing of-of2
function oferr_opt(sym, of, p, n, hdim,  evals; penalty_type = :of_convergence, penalty_energy = 500, points = 10, kw...)
    symmetrizer!(sym, of)
    mu, new_n = compute_mu_and_n_experimental(p, n, of, hdim, points = points) # first
        of2 = Of(p, mu, diagm(of), evals = evals)
        symmetrizer!(sym, of2)
        e = maximum(abs.(abs.(of)-abs.(of2)))
        e += penalty_energy * abs(sum(of-of2))
        e += penalty_energy * abs(n-new_n)    
    println("  |  Energy: ", e, "  || vars: $(round.(of, digits = 3))  ")

    return e
end

####
# Lagrange multipliers
####

""" 
finding a mimimum in the energy for a specific phase we force other solutions not to 
be energetically favoured (deprecated)
"""
function lagrange_multipliers_phases(sym, of)
    if sym == :sym
        return 0
    elseif sym == :α
        return lagrange_cost_alpha(of)
    elseif sym == :β
        return 0 
    end

end
function lagrange_cost_alpha(of)
    # Here I prevent the symmetric solution σ0s0τ0
    a = minimum(of)
    b = maximum(of)
    return 1e-6/abs.(sum(a-b))
end