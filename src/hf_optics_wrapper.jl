
using ColorSchemes
#=
                                                JDOS
=#
function jdos(ωlist; ν = 1, μ = 0.0, v = -4.303*1e3, vp = 1.622*1e3,
    M = 3.697, γ = -24.75, θ = 1.5, λ::Union{Missing, Number} = missing, hartree = 0.0, color = :black, ylims = [-10,10], kws...)
    p = paramsHF(θ, 1, 1)
    p = ParamsHF(p, μ = μ, ν = ν, v = v, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ),  hartree = hartree)
    ωlist, jdos = hf_jdos(p, ωlist; kws...)
    plot_jdos(ωlist, jdos, color, ylims)
end

function jdos(ωlist; ν = 1, μ = 0.0, v = -4.303*1e3, vp = 1.622*1e3,
    M = 3.697, γ = -24.75, θ = 1.5, λ::Union{Missing, Number} = missing, hartree = 0.0, color = :black, kws...)
    p = paramsHF(θ, 1, 1)
    p = ParamsHF(p, μ = μ, ν = ν, v = v, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ), hartree = hartree)
    ωlist, jdos = hf_jdos(p, ωlist; kws...)
    plot_jdos(ωlist, jdos ./ maximum(jdos), color, [0,1])
end
#=
                            LINEAR CONDUCTIVITY
=#
function sweep_linear_conductivity(a, b, ωlist, which_param::Symbol, param_sweep::Array; ν = 1, μ = 0.0, v = -4.303*1e3, vp = 1.622*1e3,
        M = 3.697, γ = -24.75, θ = 1.5, λ::Union{Missing, Number} = missing, hartree = 0.0, kws...)
    cmap = ColorSchemes.viridis
    p = paramsHF(θ, 1, 1)
    p = ParamsHF(p, μ = μ, ν = ν, v = v, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ), hartree = hartree)
    p = ParamsHF(p, Dict(which_param => param_sweep[1]))

    count = 0.0
    strobs = string("Re[" ,"σ_ ", string(a),string(b), "] ")
    fig = linear_conductivity(a, b, ωlist, p, strobs, :black; kws...)

    ax = Axis(fig[1,1])
    for param in param_sweep[2:end]
        count += 1/length(param_sweep)
        p = paramsHF(θ, 1, 1)
        p = ParamsHF(p, μ = μ, ν = ν, v = v, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ), hartree = hartree)
        p = ParamsHF(p, Dict(which_param => param))
        linear_conductivity!(ax, a, b, ωlist, p, strobs, cmap[count]; kws...) 
    end
    fig
end

function linear_conductivity(a, b, ωlist, p, strobs, color;  kws...)
    ωlist, conds = σab_inter_linear(a, b, p, ωlist; kws...)
    return plot_linear_conductivity(ωlist, conds, strobs, color)
end
function linear_conductivity!(ax, a, b, ωlist, p, strobs, color;  kws...)
    ωlist, conds = σab_inter_linear(a, b, p, ωlist; kws...)
    return plot_linear_conductivity!(ax, ωlist, conds, color)
end

function linear_conductivity(a, b, ωlist; ν = 1, μ = 0.0, v = -4.303*1e3, vp = 1.622*1e3,
        M = 3.697, γ = -24.75, θ = 1.5, λ::Union{Missing, Number} = missing,  color = :black, kws...)
    p = paramsHF(θ, 1, 1)
    p = ParamsHF(p, μ = μ, ν = ν, v = v, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ))

   linear_conductivity(a, b, ωlist, p; kws...)
end

function self_consistent_linear_conductivity(a, b, p, n, ωlist, color = :black; kws...)
    nf, μ, n = selfconsistency(p, n; tol = 1e-3, α = 0.4)
    ωlist, conds = σab_inter_linear(a, b, ParamsHF(p, nf = nf, μ = μ, n= n), ωlist; kws...)
    strobs = string("Re[" ,"σ_ ", string(a),string(b), "] ")
    return plot_linear_conductivity(ωlist, conds, strobs, color)
end

function self_consistent_linear_conductivity!(ax, a, b, p, n, ωlist, color = color; kws...)
    nf, μ, n = selfconsistency(p, n; tol = 1e-3, α = 0.4)
    ωlist, conds = σab_inter_linear(a, b, ParamsHF(p, nf = nf, μ = μ, n= n), ωlist; kws...)
    strobs = string("Re[" ,"σ_ ", string(a),string(b), "] ")
    return plot_linear_conductivity!(ax, ωlist, conds, color)
end

function sweep_self_consistent_linear_conductivity(nsweep, a, b, ωlist, p; kws...)
    cmap = ColorSchemes.balance
    count = 0.0
    strobs = string("Re[" ,"σ_ ", string(a),string(b), "] ")
    fig = self_consistent_linear_conductivity(a, b, p, nsweep[1], ωlist, cmap[1]; kws...)
    ax = Axis(fig[1,1])
    for n in nsweep[2:end]
        count += 1/length(nsweep)
        self_consistent_linear_conductivity!(ax, a, b, p, n, ωlist, cmap[count]; kws...) 
    end
    fig
end

function linear_conductivity(a, b, ωlist, p; color = :black, kws...)
    strobs = string("Re[" ,"σ_ ", string(a),string(b), "] ")
    ωlist, conds = σab_inter_linear(:x, :x, p, ωlist; kws...)
    return plot_linear_conductivity(ωlist, conds, strobs, color)
end


function linear_conductivity!(fig, a, b, ωlist; ν = 1, μ = 0.0, v = -4.303*1e3, vp = 1.622*1e3,
    M = 3.697, γ = -24.75, θ = 1.5, λ::Union{Missing, Number} = missing, hartree = 0.0, color = :black, kws...)
    p = paramsHF(θ, 1, 1)
    p = ParamsHF(p, μ = μ, ν = ν, v = v, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ), hartree = hartree)
    ωlist, conds = σab_inter_linear(:x, :x, p, ωlist; kws...)
    plot_linear_conductivity!(fig, ωlist, conds, :x,:x, part = :real) 
end
#=
                                         SHIFT CURRENT
=#

function shift_current(part, a, b, c, ωlist, p; kws...)
    ωlist, vals1, vals2 = shift_current(a, b, c, p, ωlist, part; kws...)
    plot_shift_current(part, a, b, c, ωlist, vals1, vals2)
end

function shift_current(part, a, b, c, ωlist; ν = 1, μ = 0.0, v = -4.303*1e3, vp = 1.622*1e3,
    M = 3.697, γ = -24.75, θ = 1.5, λ::Union{Missing, Number} = 0.0, hartree = 0.0, sigmaz = 0.0, color = :black, kws...)
    p = paramsHF(θ, 1, 1) #1 shell
    p = ParamsHF(p, μ = μ, ν = ν, v = v, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ), sigmaz = sigmaz)
    ωlist, vals1, vals2 = shift_current(a, b, c, p, ωlist, part; kws...)
    plot_shift_current(part, a, b, c, ωlist, vals1, vals2)
end

function self_consistent_shift_current(part, a, b, c, ωlist, p, n; kws...)
    nf, μ, n = selfconsistency(p, n; tol = 1e-3, α = 0.4)
    shift_current(part, a, b, c, ωlist, ParamsHF(p, nf = nf, μ = μ, n= n); kws...)
end
function self_consistent_shift_current_2spin2valleys(part, a, b, c, ωlist, p, n; kws...)
    of = spdiagm([1, 1, .0, .0, .8, .8, .2, .2])
    display(of) 
    p = ParamsHF(p, twovalleystwospins = true, U1 = 10, J = 0)
    display(hf_plotbands(p, of))
    shift_current(a, b, c, p, of, ωlist, part; kws...)
end
    

#=
                                        INJECTION CURRENT
=#
function injection_current(part, a, b, c, ωlist; ν = 1, μ = 0.0, v = -4.303*1e3, vp = 1.622*1e3, U1 = 0,
    M = 10, γ = -24.75, θ = 1.5, λ::Union{Missing, Number} = missing, hartree = 0.0, color = :black, kws...)
    p = paramsHF(θ, 1, 1)
    p = ParamsHF(p, μ = μ, ν = ν, v = v, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ), U1 = U1)
    display(hf_plotbands(p))
    ωlist, vals1, vals2 = injection_current(a, b, c, p, ωlist, part; kws...)
    plot_injection_current(part, a, b, c, ωlist, vals1, vals2)
end


function self_consistent_injection_current_2spin2valleys(part, a, b, c, ωlist, p, n; kws...)
    #            oup, oup, uup up, od, od, ud, ud
    # of = spdiagm([1, 1, .0, .0, .8, .8, .2, .2]) # tests
    # display(of) 
    if p.twovalleystwospins == false || p.VP == false
        throw(ArgumentError("p.twovalleystwospins must be == true"))
    end
    of, μ, n_new = selfconsistency(p, n)
    # μ = ifelse(n == 0, 0, μ) # to avoid intravalley optical transitions at charge neutrality.
    display(hf_plotbands(ParamsHF(p, μ = μ, n = n_new), of))
    ωlist, vals1, _= injection_current(a, b, c, ParamsHF(p, μ = μ, n = n_new), of, ωlist, part; kws...)
    plot_injection_current(part, a, b, c, ωlist, vals1)
    save_magneticinjection_cuts(ParamsHF(p, μ = μ, n = n_new), ωlist, vals1, get(kws, :η, 0.5), get(kws, :evals, 1), a, b, c, part)
    return ωlist, vals1
end

function linecuts_self_consistent_injection_current_2spin2valleys(part, a, b, c, ωlist, p, nlist; kws...)
    for n in nlist
        println("_________________________________________________________________________________________")
        println(" ")
        println("filling: ", n)
        println(" ")
        println("_________________________________________________________________________________________")
        self_consistent_injection_current_2spin2valleys(part, a, b, c, ωlist, p, n; kws...)
    end
end

function injection_map(part, a, b, c, ωlist, μlist; v = -4.303*1e3, vp = 1.622*1e3,
    M = 3.697, γ = -24.75, θ = 1.5, λ::Union{Missing, Number} = missing, hartree = 0.0, color = :black, kws...)
    p = paramsHF(θ, 1, 1)
    p = ParamsHF(p, v = v, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ), hartree = hartree)
    # fig = hf_plotbands(μ = μ, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ),  ν = 1,  color = :gray, ylims = [-40,40])
    # hf_plotbands!(fig,μ = μ, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ) , ν = -1,  color = :orange, ylims = [-40,40])
    # show(fig)
    injection_map(part, a, b, c, p, ωlist, μlist; kws...)
end
#= 
    COMPOSITE PLOTS
=#

function hatree_linecuts(part, a, b, c, ωval, μlist, hartree_list; v = -4.303*1e3, vp = 1.622*1e3,
    M = 3.697, γ = -24.75, θ = 1.5, λ::Union{Missing, Number} = missing, color = :black, kws...)
    p = paramsHF(θ, 1, 1)
    # FIGURE
    ylab = ylabel_injection(a,b,c, part)
    scale = 1e-15
    fig = Figure(resolution = (1500, 600)); 
    ax1 = Axis(fig[1, 1]; 
        title = "+ valley, ω = $(ωval)",
        xlabel = "μ [meV]", 
        ylabel = ylab,
        xlabelsize = 24, 
        ylabelsize = 24)
    ax2 = Axis(fig[1, 2]; 
        title = "- valley, ω = $(ωval)",
        xlabel = "μ [meV]", 
        ylabel = ylab,
        xlabelsize = 24, 
        ylabelsize = 24)
    ax3 = Axis(fig[1, 3]; 
        title = "both valleys, ω = $(ωval)",
        xlabel = "μ [meV]", 
        ylabel = ylab,
        xlabelsize = 24, 
        ylabelsize = 24)
    # CALCULATION
    for (i, hartree) in enumerate(hartree_list)
        auxval1 = []
        auxval2 = []
        colors = reverse(ColorSchemes.deep[LinRange(0, 1, length(hartree_list))])
        for (j, μ) in enumerate(μlist)
            p = ParamsHF(p,  μ = μ, v = v, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ), hartree = hartree)
            _ , vals1, vals2 = injection_current(a, b, c, p, ωval, part; kws...)
            push!(auxval1, vals1)
            push!(auxval2, vals2)
        end
        println(hartree)
        plot_injection_current_vs_μ!(ax1, ax2, ax3, μlist, flatten(auxval1), flatten(auxval2), color = colors[i], ylims = [-20,20])
    end
    fig
end