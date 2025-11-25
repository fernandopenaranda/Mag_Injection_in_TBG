
include("figures_correlated.jl")
""" 
run: `correlateddensityplot(folder, :REAL, :y, :y, :y, "shift", colrange = 0.1e5, rs = "m1", gs = "vp")`
"""
function correlateddensityplot(folder, part, a, b, c, which_current = "shift"; rs = "", gs = "", colrange = 1e5)
    fig = Figure(resolution = (700, 600), fontsize = 24)
    # with_theme(theme_latexfonts()) do
    if which_current == "shift"
        scale = 1
        func = (x -> x)
    else
        scale = 1e-15
        func = (x -> x)
        # func = (x -> abs(x))
    end
    if part == :REAL
        which_current *= "_real_"
    else 
        which_current *= "_imag_"
    end
    which_current *= string(a) * string(b) * string(c) * rs
    str = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Projects/HeavyFermion_Optics/src/data_atlas_selfconsistency/" * string(folder) * "/" * which_current
    # ns = collect(-4:0.2:0)
    ns =  flatten(parse_str_to_arr.(import_array(str*"cond_ns"*gs*".csv"))[:])
    inds = sortperm(ns)
    ωs = flatten(parse_str_to_arr.(import_array(str*"cond_omegas"*gs*".csv"))[:])
    l = parse_str_to_arr.(import_array(str*"cond_mat"*gs*".csv"))
    mat = zeros(length(l), length(l[1]))
    for i in 1:length(l)
        mat[i,:] .= l[i]
    end
    ax1 = Axis(fig[1, 1],
        # ylabel = "μ [meV]",
        ylabel = L"$\nu",
        xlabel = L"$\omega [\text{meV}]",
        xlabelsize = 17,
        ylabelsize = 17)
    hm = heatmap!(ax1 , ωs[1:length(l[1])], ns[inds],(func.(mat'[:,inds])./scale),colormap = :balance, colorrange = (-colrange,colrange), rasterize = true)#1 .*(-1,1) .* maximum(real.(mat./scale)))#, colorrange = colrange)#, colorrange = colrange, )
    ax1.yticks = -4:1:4
    abc = string(a) * string(b) * string(c)
    Colorbar(fig[1, 2], hm, width = 30, labelsize = 30) #  label = L"$\sigma_{$(xxx)}^\text{inj} [\text{nm} \mu A/V^2 \text{fs}])$",
    return fig
    # end
end

function correlateddensityplot!(ax, folder, part, a, b, c, which_current = "shift"; rs = "", gs = "", colrange = 1e5)
    # with_theme(theme_latexfonts()) do
    j = 1
    if which_current == "shift"
        scale2 = 1
        scale = colrange
        j *= -1
        s_select = 1 # global sign selector
        func = (x -> x)
    else
        scale2 = colrange
        scale = 1e-15 
        func = (x -> x)
        s_select = -1 # global sign selector
        # j = -1
        # func = (x -> abs(x))
    end
    if part == :REAL
        which_current *= "_real_"
    else 
        which_current *= "_imag_"
    end
    which_current *= string(a) * string(b) * string(c) * rs
    str = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Projects/HeavyFermion_Optics/src/data_atlas_selfconsistency/" * string(folder) * "/" * which_current
    # ns = collect(-4:0.2:0)
    ns =  flatten(parse_str_to_arr.(import_array(str*"cond_ns"*gs*".csv"))[:])
    inds = sortperm(ns)
    ωs = flatten(parse_str_to_arr.(import_array(str*"cond_omegas"*gs*".csv"))[:])
    l = parse_str_to_arr.(import_array(str*"cond_mat"*gs*".csv"))
    mat = zeros(length(l), length(l[1]))
    for i in 1:length(l)
        mat[i,:] .= s_select *  l[i]
    end
    sign_conv = ifelse(a == :x, -1, 1)
    matt = hcat(j .* func.(mat'[:,inds[end:-1:1]])./scale,func.(mat'[:,inds])./scale)
    hm = heatmap!(ax , ωs[1:length(l[1])], (vcat(-ns[inds[end:-1:1]], ns[inds])), sign_conv  .* matt,colormap = :balance, colorrange = (-scale2,scale2), rasterize = 2)
    ax.yticks = -4:1:4
    return hm
end


function bands_animation(p, s)
    fig =  hf_plotbands(ParamsHF(p, μ = s.mus[1]), s.ofmats[1], kpoints = 101)
    framerate = 30
    timestamps = range(1, length(s.mus), step=1)

    record(fig, "time_animation.mp4", timestamps;
            framerate = framerate) do t
        println(t)
        empty!(fig)
        hf_plotbands!(fig, ParamsHF(p, μ = s.mus[t]), s.ofmats[t], kpoints = 101)
        Label(fig[0, :], text="ν: $(s.ns[t])", fontsize=30, tellwidth=false, halign=:center)
    end
end


function jdos(p, s2, evals = 100, η = 0.5)
res2 = [shift_current(:x, :x, :x, ParamsHF(p, μ = s2.mus[i]), diagm(s2.ofmats[i]) ,collect(0:1:50), :REAL; η = η, evals = evals) for i in 1:length(s2.mus)];
end

function plot_myjdos(s2, res2)#, colrange =1)
    num_omegas = length(res2[1][1])
    num_ns = length(res2)
    mat = zeros(num_ns ,num_omegas)
    for i in 1:num_ns
        mat[i,:] .= res2[i][2]
    end
    fig = Figure(resolution = (700, 600), fontsize = 24)
    ax1 = Axis(fig[1, 1],
    # ylabel = "μ [meV]",
    ylabel = L"$\nu",
    xlabel = L"$\omega [\text{meV}]",
    xlabelsize = 24,
    ylabelsize = 24)

    hm = heatmap!(ax1 , res2[1][1], s2.ns, mat',colormap = Reverse(:grays))#, colorrange = (-0*colrange,colrange))#1 .*(-1,1) .* maximum(real.(mat./scale)))#, colorrange = colrange)#, colorrange = colrange, )
    return fig
end

function plot_dos_part(s, res)#, colrange =1)
    num_omegas = length(res[1][1])
    num_ns = length(res)
    mat = zeros(num_ns ,num_omegas)
    for i in 1:num_ns
        mat[i,:] .= res[i][2]
    end
    fig = Figure(size = (400, 600), fontsize = 24)
    ax1 = Axis(fig[1, 1],
    # ylabel = "μ [meV]",
    ylabel = L"$\nu",
    xlabel = L"$\omega [\text{meV}]",
    xlabelsize = 24,
    ylabelsize = 24)
    hm = heatmap!(ax1 , res[1][1], s.ns, mat',colormap = Reverse(:grays))#, colorrange = (-0*colrange,colrange))#1 .*(-1,1) .* maximum(real.(mat./scale)))#, colorrange = colrange)#, colorrange = colrange, )
    mat2 = mat'
    println(size(mat2))
    for i in 1:num_ns
        mat2[:,i] .= reverse(res[i][2])
    end 
    heatmap!(ax1 , res[1][1], -s.ns, mat2, colormap = Reverse(:grays))
    ax1.yticks = -4:4
    return fig
end



#bands_animation(p,s)


function orderparameters_fig(s, l)
    fig = Figure(resolution = (1000,1000));
    ax1 = Axis(fig[1,1], xlabel = "ν ",  ylabel = "sz")
    ax2 = Axis(fig[1,2], xlabel = "ν ",  ylabel = "τz")
    ax3 = Axis(fig[2,1], xlabel = "ν ",  ylabel = "szτz")
    ax4  = Axis(fig[2,2], xlabel = "ν ",  ylabel = "σz")

    scatterlines!(ax1, s.ns, [abs(i) for i in l[1]], color = :black);
    scatterlines!(ax2, s.ns, [abs(i) for i in l[2]], color = :black);
    scatterlines!(ax3, s.ns, [abs(i) for i in l[3]], color = :black);
    scatterlines!(ax4, s.ns, [abs(i) for i in l[4]], color = :black);
    ax1.xticks = collect(0:4)
    ax2.xticks = collect(0:4)
    ax3.xticks = collect(0:4)
    ax4.xticks = collect(0:4)
    fig
end

band_occups(folder::Int64) = band_occups(provide_folder_get_observables(folder))

function band_occups(s::Self_consistent_data)
    fig = Figure(resolution = (700, 600), fontsize = 24)
    # s2 = fix_valley_convention_ofmats(s, method = "1");
    # s3 = fix_spin_convention_ofmats(s2)
    # s4 = fix_spin_convention2_ofmats(s3)
    println("hey")
    
    labels = ["+↑", "-↑", "+↓", "-↓"]
    ax1 = Axis(fig[1,1], xlabel = "ν ",  ylabel = "⟨n_f⟩", title = labels[1], titlesize = 20)
    ax2 = Axis(fig[1,2], xlabel = "ν ",  ylabel = "⟨n_f⟩", title = labels[2], titlesize = 20)
    ax3 = Axis(fig[2,1], xlabel = "ν ",  ylabel = "⟨n_f⟩", title = labels[3], titlesize = 20)
    ax4 = Axis(fig[2,2], xlabel = "ν ",  ylabel = "⟨n_f⟩", title = labels[4], titlesize = 20)
    flat_ofmats= flatten(s.ofmats)
    indices(j) = [8i-7+(j-1) for i in 1:length(flat_ofmats) ÷ 8]
    count = 1
    

    hidexdecorations!(ax1)
    hidexdecorations!(ax2)
    hideydecorations!(ax2)
    hideydecorations!(ax4)
    # ax1.xticks = collect(0:4)
    # ax2.xticks = collect(0:4)
    ylims!(ax1, (0,1))
    ylims!(ax2, (0,1))
    ylims!(ax3, (0,1))
    ylims!(ax4, (0,1))
    hlines!(ax2, [0.5], color = :lightgray, linewidth=0.5)
    hlines!(ax4, [0.5], color = :lightgray, linewidth=0.5)
    vlines!(ax1, collect(0:4), color = :lightgray, linewidth=0.5)
    vlines!(ax2, collect(0:4), color = :lightgray, linewidth=0.5)
    lines!(ax1, s.ns, flat_ofmats[indices(1)])
    lines!(ax2, s.ns, flat_ofmats[indices(3)])
    lines!(ax3, s.ns, flat_ofmats[indices(5)])
    lines!(ax4, s.ns, flat_ofmats[indices(7)])
    # scatter!(ax1, s.ns, flat_ofmats[indices(1)])
    # scatter!(ax2, s.ns, flat_ofmats[indices(3)])
    # scatter!(ax3, s.ns, flat_ofmats[indices(5)])
    # scatter!(ax4, s.ns, flat_ofmats[indices(7)])
    # Legend(fig[1, 2], ax)
    fig
end


function band_occups(svp::Self_consistent_data, sqah)
    fig = Figure(resolution = (700, 600), fontsize = 24)
    # s2 = fix_valley_convention_ofmats(s, method = "1");
    # s3 = fix_spin_convention_ofmats(s2)
    # s4 = fix_spin_convention2_ofmats(s3)
    println("hey")
    
    labels = ["+↑", "-↑", "+↓", "-↓"]
    ax1 = Axis(fig[1,1], xlabel = "ν ",  ylabel = "⟨n_f⟩", title = labels[1],     titlesize = 20)
    ax2 = Axis(fig[1,2], xlabel = "ν ",  ylabel = "⟨n_f⟩", title = labels[2],     titlesize = 20)
    ax3 = Axis(fig[2,1], xlabel = "ν ",  ylabel = "⟨n_f⟩", title = labels[3],     titlesize = 20)
    ax4 = Axis(fig[2,2], xlabel = "ν ",  ylabel = "⟨n_f⟩", title = labels[4],     titlesize = 20)
    flat_ofmats_vp= flatten(svp.ofmats)
    flat_ofmats_qah= flatten(sqah.ofmats)
    indices(j) = [8i-7+(j-1) for i in 1:length(flat_ofmats_vp) ÷ 8]
    count = 1
    

    hidexdecorations!(ax1)
    hidexdecorations!(ax2)
    hideydecorations!(ax2)
    hideydecorations!(ax4)
    # ax1.xticks = collect(0:4)
    # ax2.xticks = collect(0:4)
    ylims!(ax1, (0,1))
    ylims!(ax2, (0,1))
    ylims!(ax3, (0,1))
    ylims!(ax4, (0,1))
    hlines!(ax2, [0.5], color = :lightgray, linewidth=0.5)
    hlines!(ax4, [0.5], color = :lightgray, linewidth=0.5)
    vlines!(ax1, collect(0:4), color = :lightgray, linewidth=0.5)
    vlines!(ax2, collect(0:4), color = :lightgray, linewidth=0.5)
    lines!(ax1, svp.ns, flat_ofmats_vp[indices(1)])
    lines!(ax2, svp.ns, flat_ofmats_vp[indices(3)])
    lines!(ax3, svp.ns, flat_ofmats_vp[indices(5)])
    lines!(ax4, svp.ns, flat_ofmats_vp[indices(7)])

    lines!(ax1, sqah.ns, flat_ofmats_qah[indices(1)])
    lines!(ax2, sqah.ns, flat_ofmats_qah[indices(3)])
    lines!(ax3, sqah.ns, flat_ofmats_qah[indices(5)])
    lines!(ax4, sqah.ns, flat_ofmats_qah[indices(7)])
    # Legend(fig[1, 2], ax)
    fig
end



# hf_plotbands(ParamsHF(p, μ = svp.mus[60], U1 = 10, U2 = 0, J = 5, sigmaz = 2.5, sigmazlayerz = 0, λ = 140, vp = 2e3), svp.ofmats[60], kpoints = 101)
# hf_plotbands(ParamsHF(p, μ = sqah.mus[60], U1 = 10, U2 = 0, J = 5, sigmaz = 2.5, sigmazlayerz = 0, λ = 140, vp = 2e3), sqah.ofmats[60], kpoints = 101)

# for this comment in shift_current_ω the indicated lines:
# joint_dos: res2 = [shift_current(:x, :x, :x, ParamsHF(p, μ = s2.mus[i]), diagm(s2.ofmats[i]) ,collect(0:0.3:30), :REAL; η = 0.5, evals = 400) for i in 1:length(s2.mus)];
# dos: res2 = [shift_current(:x, :x, :x, p, diagm(s2.ofmats[i]) ,collect(0:1:40), :REAL; η = 0.5, evals = 10) for i in 1:length(s2.ofmats)];
# o # dos: res2 = [shift_current(:x, :x, :x, ParamsHF(p, μ = s2.mus[i]), diagm(s2.ofmats[i]) ,collect(-40:1:40), :REAL; η = 0.5, evals = 50) for i in 1:length(s2.ofmats)];
# plot_myjdos(s2, res2)

# dos(folder, evals = 1000, omax = 50, steps = 0.25)
function dos(folder, p; evals = 10, omax = 50, steps = 0.25) # RUN DOS
    s = provide_folder_get_observables(folder)
    # s = set_filling_order(provide_folder_get_observables(folder), permlist = [1,4,3,2]) # I select one which has to be degenerate with the rest.
    res = [shift_current(:x, :x, :x, ParamsHF(p, μ = s.mus[i]), diagm(s.ofmats[i]) ,collect(-omax:steps:omax), :REAL; η = 0.25, evals = evals) for i in 1:length(s.ofmats)]
    save_dos(s, res, folder)
    res
end

function save_dos(s, res, folder)
    omegalist = res[1][1] 
    nlist = s.ns
    str = pwd()
    ofvec = []
    for i in 1:length(res)
        push!(ofvec, res[i][2])
    end
    save_arrays(string(str, "/dos/"*string(folder)*"nlist_2.csv"), nlist)
    save_arrays(string(str, "/dos/"*string(folder)*"omegas_2.csv"), omegalist)
    save_arrays(string(str, "/dos/"*string(folder)*"ofmats_2.csv"), ofvec)
end

function import_dos(str, folder)
    str = pwd()
    ns = import_array(str * "/dos/"*string(folder)*"nlist_2.csv")[:]
    omegas = import_array(str * "/dos/"*string(folder)*"omegas_2.csv")[:]
    mat = parse_str_to_arr.(import_array(str * "/dos/"*string(folder)*"ofmats_2.csv"))
    return ns, omegas, mat
end

function plot_dos_part!(ax, folder)
    s = provide_folder_get_observables(folder)
    str = pwd() 
    res = import_dos(str, folder)
    plot_dos_part!(ax, s, res)
end

function plot_dos_part!(ax, s, res)#, colrange =1)
    num_omegas = length(res[2])
    num_ns = length(res[1])
    mat = zeros(2num_ns-1 ,num_omegas)
    for i in 1:2num_ns-1
        if i<num_ns+1
            mat[i,:] .= reverse(res[3][num_ns+1-i] )  
        else
            mat[i,:] .= res[3][i-num_ns+1] 
        end
    end
    # fig = Figure(size = (400, 600), fontsize = 24)
    # ax1 = Axis(fig[1, 1],
    # # ylabel = "μ [meV]",
    # ylabel = L"$\nu",
    # xlabel = L"$\omega [\text{meV}]",
    # xlabelsize = 24,
    # ylabelsize = 24)
    # hm = heatmap!(ax , res[2], s.ns, mat',colormap = Reverse(:grays), rasterize = false)#, colorrange = (-0*colrange,colrange))#1 .*(-1,1) .* maximum(real.(mat./scale)))#, colorrange = colrange)#, colorrange = colrange, )
    
    hm = heatmap!(ax ,res[2], vcat(-s.ns[end:-1:2], s.ns), mat', colormap = Reverse(:grays), rasterize = 10)

    # heatmap!(ax , res[2], -s.ns[2:end], mat2[:,2:end], colormap = Reverse(:grays), rasterize = true)
    ax.yticks = -4:4
end

