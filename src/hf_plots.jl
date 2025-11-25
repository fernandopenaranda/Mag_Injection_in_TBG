function hf_plotbands(p::ParamsHF; kpoints = 100, ylims = [-30, 30], res =  (500*1, 700))
    fig = Figure(size =res)
    if p.twovalleys == false
        plotbands!(fig, real.(bands(ParamsHF(p, ν = 1), kpoints)), ylimits = ylims, color = :orange)
        plotbands!(fig, real.(bands(ParamsHF(p, ν = -1), kpoints)), ylimits = ylims, color = :gray)
    else
        plotbands!(fig, real.(bands(ParamsHF(p, ν = 1), kpoints)), ylimits = ylims, color = :black)
    end
        fig
end

function hf_plotbands(p::ParamsHF, of; kpoints = 100, ylims = [-30, 30], res =  (500*1, 700))
    fig = Figure(size = res)
    plotbands!(fig, real.(bands(p, of, kpoints)), ylimits = ylims, color = :black)
    fig
end


function hf_plotbands!(fig, p::ParamsHF, of; kpoints = 100, ylims = [-30, 30])
    plotbands!(fig, real.(bands(p, of, kpoints)), ylimits = ylims, color = :black)
    fig
end

function hf_plotbands(;nmax = 1, ν = 1, μ = 0.0, v = -4.303*1e3, vp = 1.622*1e3, M = 3.697, γ = -24.75, θ = 1.5, λ::Union{Missing, Number} = missing, hartree= 0.0, sigmaz = sigmaz,
     kpoints = 1000, color = :lightgray, ylims = [-60,60])
    p = paramsHF(θ, nmax, 1); 
    p = ParamsHF(p, μ = μ, ν = ν, v = v, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ), hartree = hartree, sigmaz = sigmaz)
    plotbands(bands(p, kpoints), ylimits = ylims, color = color)
end


function hf_plotbands!(fig ;nmax = 1, ν = 1, μ = 0.0, v = -4.303*1e3, vp = 1.622*1e3, M = 3.697, γ = -24.75, θ = 1.5, λ::Union{Missing, Number} = missing, 
    hartree = 0.0, sigmaz = 0.0, kpoints = 1000, color = :orange, ylims = [-60,60])
    p = paramsHF(θ, nmax, 1); 
    p = ParamsHF(p, μ = μ, ν = ν, v = v, vp = vp, M = M, γ = γ, λ = (isa(λ, Missing) ? lambda(θ) : λ), hartree = hartree, sigmaz = sigmaz)
    plotbands!(fig, real.(bands(p, kpoints)), ylimits = ylims, color = color)
    fig
end

plotbands(mat; kw...) = plotbands!(Figure(resolution = (800, 1100)), real.(mat); kw...)

function plotmybands!(ax::Axis, mat; dash = missing, dots = false, color = missing, ylimits = missing, xlimits = missing)
    xarr = collect(1:size(mat,2))
    num_points = size(mat,2)
        for i in 1:size(mat, 1)
            lines!(ax, xarr , (mat[i,:]), color = ifelse(isa(color, Missing), :lightgray, color), linestyle = ifelse(isa(dash,Missing), :solid, dash), linewidth = 0.9)
        end
    l = num_points /  (2√3/2+ 2/2)
    ax.xticks = [1, 1/2*l, (1/2+ √3/2) * l, (1/2+ 2√3/2) * l, (2√3/2+ 2/2) * l], ["K","-M", "Γ", "M", "K'"]
    if isa(ylimits,Missing)
        nothing
    else
        CairoMakie.ylims!(ax, ylimits[1], ylimits[2])
    end
    if isa(xlimits,Missing)
        nothing
    else
        xlims!(ax, ylimits[1], ylimits[2])
    end
end

function plotbands!(f::Figure, mat; dots = false, color = missing, ylimits = missing, xlimits = missing)
    ax = Axis(f[1, 1]; xlabel = "k", ylabel = "E [meV]", xlabelsize= 27, ylabelsize= 27, xticklabelsize = 27, yticklabelsize = 27)
    xarr = collect(1:size(mat,2))
    num_points = size(mat,2)
        for i in 1:size(mat, 1)
            lines!(ax, xarr , (mat[i,:]), color = ifelse(isa(color,Missing), :lightgray, color))
        end
    l = num_points /  (2√3/2+ 2/2)
    ax.xticks = [1, 1/2*l, (1/2+ √3/2) * l, (1/2+ 2√3/2) * l, (2√3/2+ 2/2) * l], ["K","-M", "Γ", "M", "K'"]
    if isa(ylimits,Missing)
        nothing
    else
        CairoMakie.ylims!(ax, ylimits[1], ylimits[2])
    end
    if isa(xlimits,Missing)
        nothing
    else
        xlims!(ax, ylimits[1], ylimits[2])
    end
    ylims!(ax, -40,40)
    return f
end

function plotbands!(ax::Axis, mat; dots = false, color = missing, ylimits = missing, xlimits = missing)
    xarr = collect(1:size(mat,2))
    pointsk = 2/5 * (length(xarr)-1)
    if dots == false
        for i in 1:size(mat, 1)
            lines!(ax, xarr , mat[i,:], color = ifelse(isa(color,Missing), :lightgray, color))
        end
    else
        for i in 1:size(mat, 1)
            scatter!(ax, collect(1:size(mat,2)) , mat[i,:], markersize = 5, markeralpha = 0.8)
        end
    end   
    # ax.xticks = ([1, pointsk+1, 2*pointsk+1, 5*pointsk/2 + 1], ["K1", "Γ", "M", "K2"])
    ax.xticks = ([1, pointsk+1, 2*pointsk+1, 5*pointsk/2 + 1], ["-M", "Γ", "M", "K2"])

    if isa(ylimits,Missing)
        nothing
    else
        ylims!(ax, ylimits[1], ylimits[2])
    end
    if isa(xlimits,Missing)
        nothing
    else
        xlims!(ax, ylimits[1], ylimits[2])
    end

    return f
end



function plothartreebands(μ, hart; M = 4, λ = 90)
    fig = hf_plotbands(M = M, λ = λ, μ = μ, ν = 1, nmax =1,  color = :gray, ylims = [-50,50], hartree = hart, sigmaz = 0);
    hf_plotbands!(fig, M = M, λ = λ, μ = μ, ν = -1, nmax =1,  color = :orange, ylims = [-50,50], hartree = hart, sigmaz = 0)
end

########### EXECUTABLES
    # fig = hf_plotbands()
    # hf_plotbands!(fig, λ = 0)