using NLopt

function optimize_sym(sym, f, init_vals; kw...)::My_optim
    if sym == :sym   
        sym_ln_newuoa_bound(f, init_vals; kw...)
    elseif sym == :α
        α_ln_newuoa_bound(f, init_vals; kw...)
    elseif sym == :β
        β_ln_newuoa_bound(f, init_vals; kw...)
    end
end

function sym_alg(f, init_vals)  #method 1 UNUSED
    result = optimize(f, 0.0, 1.0, show_trace=true, rel_tol = 0.00001)
    return My_optim(Optim.minimizer(result), Optim.minimum(result))
end

"""
two derivative free (N)  methods  :LN_COBYLA and :GN_ISRES the first is local and 
the second global. The global method is passed before and its followed by the local one.
"""
function sym_ln_newuoa_bound(f, init_vals; lower_bound = [0.0], upper_bound = [1.0], 
                                max_evals = 20, method = :LN_COBYLA, kw...)  #method 2
    if method == :local 
        which_algorithm = :LN_COBYLA
    else
        which_algorithm = :GN_ISRES
    end
    evaluations = Float64[] # track evaluations 
    objective_function(x::Vector, grad::Vector) = begin 
        push!(evaluations, x[1])
        f(x)
    end
    opt = Opt(which_algorithm, 1)
    opt.lower_bounds = lower_bound
    opt.upper_bounds = upper_bound
    opt.min_objective = objective_function
    # opt.xtol_rel = 1e-6
    opt.ftol_rel = 1e-6
    opt.maxeval = max_evals
    initial_guess = [init_vals[1]]
    (minf, minx, ret) = NLopt.optimize(opt, initial_guess)
    println(" ")
    println("CONVERGENCE: ", ret)
    return My_optim(minx, minf, mean(diff(sort(evaluations))), ret)
end

"""
local (l) derivative free (n) algorithm for the α phase
"""
function α_ln_newuoa_bound(f, init_vals; lower_bound = [0.0,0.0], upper_bound = [1.0, 1.0], 
                            max_evals = 20, method = :local, kw...)
    
    evaluations = Float64[]
    objective_function(x::Vector, grad::Vector) = begin 
        push!(evaluations, x[1])
        f(x)
    end
    if method == :local 
        which_algorithm = :LN_NEWUOA_BOUND
    else
        which_algorithm = :GN_ISRES
    end
    opt = Opt(which_algorithm, 2)
    opt.lower_bounds = lower_bound
    opt.upper_bounds = upper_bound
    opt.min_objective = objective_function
    opt.ftol_rel = 1e-6  # Relative tolerance on function value
    opt.maxeval = max_evals   # Maximum number of evaluations
    initial_guess = [maximum(init_vals), minimum(init_vals)]
    (minf, minx, ret) = NLopt.optimize(opt, initial_guess)
    println(" ")
    println("CONVERGENCE: ", ret)
    return My_optim(minx, minf, mean(diff(sort(evaluations))), ret)
end

function β_ln_newuoa_bound(f::Function, init_vals; lower_bound = zeros(Float64, 8), upper_bound = ones(Float64, 8), 
                            max_evals = 20, method = :local, kw...)
    evaluations = Float64[]
    objective_function(x::Vector, grad::Vector) = begin 
        # allocations = @allocations f(x)
        # println("allocations: ", allocations)
        push!(evaluations, x[1])
        f(x)
    end
    if method == :local 
        which_algorithm = :LN_NEWUOA_BOUND
    else
        which_algorithm = :GN_ISRES
    end
    opt = Opt(which_algorithm, 8)

    opt.lower_bounds =  lower_bound 
    opt.upper_bounds = upper_bound 
    opt.min_objective = objective_function

    opt.ftol_rel = 0.5e-5   # Relative tolerance on function value, Means an error of 0.04 meV
    opt.maxeval = max_evals    # Maximum number of evaluations
    
    initial_guess = init_vals
    (minf, minx, ret) = NLopt.optimize(opt, initial_guess)
    println(" ")
    println("CONVERGENCE - β: ", ret)
    return My_optim(minx, minf, mean(diff(sort(evaluations))), ret)
end


function of_optimized(sym, result)
    if sym == :sym
        a = result.of[1]
        optim_of = a .* ones(8)
    elseif sym == :α
        a = result.of[1]
        b = result.of[2]
        optim_of = [a, a, b, b, a, a, b, b]
        symmetrizer!(:α, optim_of)
    elseif sym == :β
        a = result.of[1]
        b = result.of[2]
        c = result.of[3]
        d = result.of[4]
        e = result.of[5]
        f = result.of[6]
        g = result.of[7]
        h = result.of[8]

        optim_of = [a, b, c, d, e, f, g, h]
        symmetrizer!(:β, optim_of)
    end
    return optim_of
end