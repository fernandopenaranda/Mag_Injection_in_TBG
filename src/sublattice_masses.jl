# Plots for bandstructures with different mass terms. Comparison between the continuum model and the THFM
# THFM
function plotmybands!(ax::Axis, mat; dash = missing, dots = false, color = missing, ylimits = missing, 
        xlimits = missing)
    xarr = collect(1:size(mat,2))
    num_points = size(mat,2)
        for i in 1:size(mat, 1)
            lines!(ax, xarr , (mat[i,:]), color = ifelse(isa(color, Missing), :lightgray, color), 
                linestyle = ifelse(isa(dash,Missing), :solid, dash), linewidth = 1.3)
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

# Continuum
function compute_bands_continuum!(pHF::ParamsHF, pcont, ax1; mass_term = :none,  ylimits = missing)
    num_points = 277
    energies(valley) = bands_bistritzer(pHF, ParamsBM(pcont, ν = valley, pointsk = num_points), mass_term = mass_term, eigvecs = false)[1]
    plotbandws!(ax1, real.(energies(1)), color = :blue)
    plotbandws!(ax1, real.(energies(-1)), color = :red)
    l = num_points /  (2√3/2+ 2/2)
    ax1.xticks = [1, 1/2*l, (1/2+ √3/2) * l, (1/2+ 2√3/2) * l, (2√3/2+ 2/2) * l], ["K","-M", "Γ", "M", "K'"]
    if isa(ylimits, Missing)
        nothing
    else
        CairoMakie.ylims!(ax1, ylimits[1], ylimits[2])
    end
    xlims!(ax1, 1,  (2√3/2+ 2/2) * l)
end

# Continuum

function bands_bistritzer(pHF::ParamsHF, p::ParamsBM; eigvecs = false, mass_term = :none, kws...)
    km = k_path(pHF, 101) 
    kmesh = hcat([x for (x, y) in  km], [y for (x, y) in  km]) .* a0
    
    sys_dim = 4 * Int(1 + 3p.nmax * (1 + p.nmax))
    es = zeros(Float64, sys_dim, size(kmesh,1))
    for i in 1:size(kmesh,1)
        es[:,i] = real.(bistritzer_eigs(p, kmesh[i,:], mass_term = mass_term; kws...)[1])
    end
    return (es, )
end


function plotbandws!(ax, mat; dots = false, color = missing, ylimits = missing, xlimits = missing)
    xarr = collect(1:size(mat,2))
    for i in 1:size(mat, 1)
        lines!(ax,  xarr  ,1000* mat[i,:], color = (color, 0.8), linewidth = 1.5)
    end
    return f
end
