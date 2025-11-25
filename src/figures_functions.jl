#____________________________________________________________________________________________________________________
# Figure 1
#____________________________________________________________________________________________________________________

function lat_fig1()
    with_theme(merge(theme_latexfonts(), Theme(fontsize = 15))) do #
        figure1()     
    end
end

function figure1(;ylims = (-40,40))
    fig = Figure(size = (350, 400))
    ax1 = Axis(fig[1,1], ylabel = L"$\text{E [meV]}$") #bands
    p = ParamsHF(paramsHF(1.05, 1,  1),
                  μ = 0,
                  nf = 1,
                  λ = 140,
                  M = 5,
                  v = -7e3,
                  vp = 2e3,
                  γ =-30,
                  sigmaz = 0,
                  sigmazlayerz = 0,
                  U1 = 0,
                  U2 = 0,   
                  J = 1,
                  VP =true,
                  twovalleystwospins= true)
    plotmybands!(ax1, real.(bands(ParamsHF(p, μ = 0), ones(8),  101)), ylimits = ylims, color = :lightgray)
    plotmybands!(ax1, real.(bands(ParamsHF(p, μ = 0, sigmaz = 2.5), ones(8), 101)), ylimits = ylims, color = :black)
    fig
end
#____________________________________________________________________________________________________________________
# Figure 2
#____________________________________________________________________________________________________________________


function lat_fig2()
    with_theme(merge(theme_latexfonts(), Theme(fontsize = 15))) do #
        fig2()
    end
end

function fig2()
    # Presets BM and THFM
    pBM = ParamsBM(paramsBM(1.05, 3,  1), mass = 0);
    pHF = ParamsHF(paramsHF(1.05, 1,  1),
        μ = 0,
        nf = 1,
        λ = 92,
        M = 3.65,
        v = -4.303e3*a0,
        vp = 1.622e3*a0,
        vpp = -0.0332*a0,
        γ = -24.75,
        sigmaz = 0.,
        sigmazlayerz = 0,
        layerz = 0,
        U1 = 0,
        U2 = 0,   
        J = 0,
        VP =true,
        twovalleystwospins= true);
    mass0 = 0.25
    ξ = (2sin(-π/sqrt(3))+sin(2π/sqrt(3)))
    mass = 4π*ξ* mass0/1000
    pBM_delta3 = ParamsBM(pBM, mass = mass)
    pHF_delta3= ParamsHF(pHF, layerz = mass0)

    #heavyfermion bands 
    bhf = bands(ParamsHF(pHF, sigmaz = 0, sigmazlayerz= 0, layerz = 0), ones(8), 101);
    bhfdelta1 = bands(ParamsHF(pHF, sigmaz = mass0, sigmazlayerz= 0, layerz = 0), ones(8), 101);
    bhfdelta2 = bands(ParamsHF(pHF, sigmaz = 0, sigmazlayerz= mass0, layerz = 0), ones(8), 101);
    bhfdelta3 = bands(ParamsHF(pHF, sigmaz = 0, sigmazlayerz= 0, layerz = mass0), ones(8), 101);
    #continuum bands
    bbm_pv = bands_bistritzer(pHF, ParamsBM(pBM, ν = 1, mass = 0mass), pointsk = 101, mass_term = :tauz, eigvecs = false)[1];
    bbm_nv = bands_bistritzer(pHF, ParamsBM(pBM, ν = -1, mass = 0mass), pointsk = 101, mass_term = :tauz, eigvecs = false)[1];
    bbm_pv_delta1 = bands_bistritzer(pHF_delta3, ParamsBM(pBM_delta3, ν = 1, mass = mass/2ξ/pi), pointsk = 101, mass_term = :normal, eigvecs = false)[1];
    bbm_nv_delta1 = bands_bistritzer(pHF_delta3, ParamsBM(pBM_delta3, ν = -1, mass = mass/2ξ/pi), pointsk = 101, mass_term = :normal, eigvecs = false)[1];
    bbm_pv_delta2 = bands_bistritzer(pHF_delta3, ParamsBM(pBM_delta3, ν = 1, mass = mass), pointsk = 101, mass_term = :sigmaztauz, eigvecs = false)[1];
    bbm_nv_delta2 = bands_bistritzer(pHF_delta3, ParamsBM(pBM_delta3, ν = -1, mass = mass), pointsk = 101, mass_term = :sigmaztauz, eigvecs = false)[1];
    bbm_pv_delta3 = bands_bistritzer(pHF_delta3, ParamsBM(pBM_delta3, ν = 1, mass = mass), pointsk = 101, mass_term = :layered, eigvecs = false)[1];
    bbm_nv_delta3 = bands_bistritzer(pHF_delta3, ParamsBM(pBM_delta3, ν = -1, mass = mass), pointsk = 101, mass_term = :layered, eigvecs = false)[1];

    fig_sublattices = Figure(size = (600, 600))
    gqah = GridLayout(fig_sublattices[1,1])
    ax1 = Axis(gqah[1,1], ylabel = L"$\text{E [meV]}$") ;
    plotmybands!(ax1, real.(bhf), ylimits = [-5,5], color = :black); 
    plotbandws!(ax1, real.(bbm_nv), color = :darkorange, ylimits =[-5,5])
    plotbandws!(ax1, real.(bbm_pv), color = :darkorange, ylimits = [-5,5])

    ax2 = Axis(gqah[1,2], ylabel = L"$\text{E [meV]}$", xlabel = L"$\text{k}$") ;
    plotmybands!(ax2, real.(bhfdelta1), ylimits = [-5,5], color = :black); 
    plotbandws!(ax2, real.(bbm_nv_delta1), color = :darkorange, ylimits =[-5,5])
    plotbandws!(ax2, real.(bbm_pv_delta1), color = :darkorange, ylimits = [-5,5])

    ax3 = Axis(gqah[2,1], ylabel = L"$\text{E [meV]}$", xlabel = L"$\text{k}$") ;
    plotmybands!(ax3, real.(bhfdelta2), ylimits = [-5,5], color = :black); 
    plotbandws!(ax3, real.(bbm_nv_delta2), color = :darkorange, ylimits =[-5,5])
    plotbandws!(ax3, real.(bbm_pv_delta2), color = :darkorange, ylimits = [-5,5])

    ax4 = Axis(gqah[2,2], ylabel = L"$\text{E [meV]}$", xlabel = L"$\text{k}$") ;
    plotmybands!(ax4, real.(bhfdelta3), ylimits = [-5,5], color = :black); 
    plotbandws!(ax4, real.(bbm_nv_delta3), color = :darkorange, ylimits =[-4,4])
    plotbandws!(ax4, real.(bbm_pv_delta3), color = :darkorange, ylimits = [-4,4])

    hidexdecorations!(ax1, grid = false)
    hidexdecorations!(ax2, grid = false)
    hideydecorations!(ax2, grid = false)
    hideydecorations!(ax4, grid = false)

    add_label!(ax, label) = text!(ax, label, position = (2, 3.4), fontsize = 15) 
    add_label!(ax1, "(a)")
    add_label!(ax2, "(b)")
    add_label!(ax3, "(c)")
    add_label!(ax4, "(d)")
    lines!(ax1, [1], [1000], label = "THFM", color = :black)
    lines!(ax1, [1], [1000], label = "Continuum", color = :darkorange)  
    lg = axislegend(ax1,  labelsize= 13, boxed = false)
    lg.framevisible = false  
    lg.halign = :right
    rowgap!(gqah, 15)
    colgap!(gqah, 15)
    xlims!(ax1, 1,277)
    xlims!(ax2, 1,277)
    xlims!(ax3, 1,277)
    xlims!(ax4, 1,277)
    ylims!(ax1, -4,4)
    ylims!(ax2, -4,4)
    ylims!(ax3, -4,4)
    ylims!(ax4, -4,4)
    fig_sublattices
end
#____________________________________________________________________________________________________________________
# Figure 3
#____________________________________________________________________________________________________________________

function lat_fig3(folder)
    with_theme(merge(theme_latexfonts(), Theme(fontsize = 17))) do #
        figure3(folder)     
    end
end

function figure3(folder)
    fig = Figure()
    lsize = 15
    ax_dos = Axis(fig[1:3,1], xlabel = L"\text{E [meV]}", ylabel = L"$\nu$", xticklabelsize = lsize, yticklabelsize = lsize)
    ga = GridLayout(fig[1:3,3])
    ax_vh = Axis(ga[2,1], xlabel = L"$\nu$", ylabel = L"$\langle n_i \rangle$", xticklabelsize = lsize, yticklabelsize = lsize)
    ax_svh = Axis(ga[3,1], xlabel = L"$\nu$", ylabel = L"$\langle n_i \rangle$", xticklabelsize = lsize, yticklabelsize = lsize)
    ax_qah = Axis(ga[4,1], xlabel = L"$\nu$", ylabel = L"$\langle n_i \rangle$", xticklabelsize = lsize, yticklabelsize = lsize)
    ax_tau = Axis(fig[1:3,2], xlabel = L"$| \langle \tau_z \rangle |$", xticklabelsize = lsize, yticklabelsize = lsize)
    text!(ax_vh, L"VH$$", position = (-0,0.8))
    text!(ax_svh, L"SVH$$", position = (-0,0.8))
    text!(ax_qah, L"QAH$$", position = (-0,0.8))
    plot_dos_part!(ax_dos, folder)

    s_vh = set_filling_order(provide_folder_get_observables(folder), permlist = [1,4,3,2])
    s_svh = set_filling_order(provide_folder_get_observables(folder), permlist = [1,2,3,4])
    s_qah = set_filling_order(provide_folder_get_observables(folder), permlist = [1,3,2,4])
    plot_occup!(ax_vh, s_vh)
    plot_occup!(ax_svh, s_svh)
    plot_occup!(ax_qah, s_qah)
    hidexdecorations!(ax_vh, grid = false)
    hidexdecorations!(ax_svh, grid = false)

    plot_tau!(ax_tau, s_vh, s_svh, s_qah)
    colsize!(fig.layout, 2, Auto(0.5))

    Label(fig.layout[1, 1, TopLeft()], "(a)",
    fontsize = 15,
    padding = (0, -90, 5, 0),
    halign = :right)

    Label(fig.layout[1, 2, Top()], "(b)",
    fontsize = 15,
    padding = (0, 28, 5, 0),
    halign = :right)

    Label(fig.layout[1, 3, TopLeft()], "(c)",
    fontsize = 15,
    padding = (0, 35, -70, 0),
    halign = :right)

    Label(fig.layout[2, 3, TopLeft()], "(d)",
    fontsize = 15,
    padding = (0, 35, -5, 0),
    halign = :right)

    Label(fig.layout[3, 3, TopLeft()], "(e)",
    fontsize = 15,
    padding = (0, 35, -5, 0),
    halign = :right)
    elem1 = MarkerElement(color = :red, marker=:rect, markersize = 9, strokecolor = :black)
    elem2 = MarkerElement(color = :blue, marker=:rect, markersize = 9, strokecolor = :black)
    elem3 = LineElement(color = :black,  linewidth = 1.5)
    elem4 = LineElement(color = :black, linestyle = :dash, linewidth = 1.5)
    legenditems = [elem1, elem2, elem3, elem4]

    terms = [L"$s$", L"$\bar{s}$", L"$\eta$", L"$\bar{\eta}$"]
    Legend(ga[1, 1], legenditems, terms, orientation = :horizontal, tellwidth = true, height =20, framewidth = .5,  textsize = 8)#, padding = (0,0,0,0))
    fig
end

function plot_tau!(ax, s1, s2, s3)
    t1 = tauz(s1)
    t2 = tauz(s2)
    t3 = tauz(s3)
    lines!(ax,  t3, s1.ns,  color = :black, label = "QAH")
    lines!(ax,  t1, s1.ns, color = :gray, linestyle = :dash, label = "VH/SVH")
    lines!(ax,  t1, -s1.ns, color = :gray, linestyle = :dash)
    lines!(ax,  t2, s1.ns,  color = :gray, linestyle = :dash)
    lines!(ax,  t2, -s1.ns, color = :gray, linestyle = :dash)
    lines!(ax,  t3, -s1.ns, color = :black)
    ax.yticks = -4:4
    ax.xticks = 0:2
    xlims!(-0.01,2)
    hideydecorations!(ax, grid = false)
    ylims!(-4,4)
    axislegend(ax; position= (2, 0.5), align = :right, patchsize=(12,10), labelsize = 11, backgroundcolor = (:white, 0.6), height = 40, width = 80 )
end

function plot_occup!(ax, s)
    flat_ofmats= flatten(s.ofmats)
    indices(j) = [8i-7+(j-1) for i in 1:length(flat_ofmats) ÷ 8]
    ylims!(ax, (0,1))
    lines!(ax, s.ns, flat_ofmats[indices(1)], color = :red, linestyle = :dash)
    lines!(ax, s.ns, flat_ofmats[indices(3)], color = :red)
    lines!(ax, s.ns, flat_ofmats[indices(5)], color = :blue, linestyle = :dash)
    lines!(ax, s.ns, flat_ofmats[indices(7)], color = :blue)
end
#____________________________________________________________________________________________________________________
# Figure 4
#____________________________________________________________________________________________________________________
function lat_fig4(folder)
    with_theme(merge(theme_latexfonts(), Theme(fontsize = 13))) do #
        figure4(folder)     
    end
end

function figure4(folder)
    fig = Figure(size = (1800/2, 700/2))
    valshift = 1e4
    valinj = 4

    g1 = GridLayout(fig[1,1])
    g2 = GridLayout(fig[1,2])
    g3 = GridLayout(fig[1,3])

    ax_1 = Axis(g1[2,1], xlabel = L"$\omega$ [meV]", ylabel = L"$\nu$")
    ax_2 = Axis(g2[2,1], xlabel = L"$\omega$ [meV]", ylabel = L"$\nu$")
    ax_3 = Axis(g3[2,1], xlabel = L"$\omega$ [meV]", ylabel = L"$\nu$")
    
    hm1 = correlateddensityplot!(ax_1, folder, :REAL, :y, :y, :y, "shift", rs = "m1", gs = "vp", colrange = valshift) #da igual vp o vh
    hm2 = correlateddensityplot!(ax_2, folder, :REAL, :x, :x, :x ,  "inj", rs = "m1", gs = "vp", colrange = valinj) #da igual vp o vh
    hm3 = correlateddensityplot!(ax_3, folder, :REAL,  :x, :x, :x, "inj", rs = "m1", gs = "qah", colrange = valinj) #da igual vp o vh

    cb = Colorbar(g1[1,1], hm1,  vertical = false, flipaxis = false, ticklabelpad= -2)
    Label(g1[1,1, Top()], halign = :center,  L" $\sigma_{yyy}^\text{sh} \text{ }  [\text{nm} \mu A/V^2]$")
    Label(g1[2,1, Top()], halign = :center,  L"\text{VH/VSH/QAH}")
    cb.height = 10
    cb.ticks = [-1,1]
    cb.tellwidth = true
    rowgap!(g1, -10)

    cb = Colorbar(g2[1,1], hm2,  vertical = false, flipaxis = false, ticklabelpad= -2)
    Label(g2[1,1, Top()], halign = :center,  L" $\sigma_{xxx}^\text{inj}/\tau \text{ }  [\text{nm} \mu A/V^2 \text{fs}]$")
    Label(g2[2,1, Top()], halign = :center,  L"\text{QAH}")
    cb.height = 10
    cb.ticks = [-valinj, valinj]
    cb.tellwidth = true
    rowgap!(g2, -10)
    
    cb = Colorbar(g3[1,1], hm3,  vertical = false, flipaxis = false, ticklabelpad= -2)
    Label(g3[1,1, Top()], halign = :center,  L" $\sigma_{xxx}^\text{inj}/\tau \text{ } [\text{nm} \mu A/V^2 \text{fs}]$")
    Label(g3[2,1, Top()], halign = :center,  L"\text{VH/VSH}")
    cb.height = 10
    cb.ticks = [-valinj, valinj]
    cb.tellwidth = true
    rowgap!(g3, -10)
    
    valstr = Int(log10(valshift));
    Label(g1[1, 1, Top()], halign = :right, latexstring("\$\\times 10^{$(valstr)}\$"),  padding = (0, 0, 0, 0))
    hideydecorations!(ax_2, grid = false)
    hideydecorations!(ax_3,grid = false)

    Label(fig.layout[1, 1, Top()], "(a)",
    fontsize = 12,
    padding = (0, 0, 5, 0),
    halign = :left)

    Label(fig.layout[1, 2, Top()], "(b)",
    fontsize = 12,
    padding = (0, 0, 5, 0),
    halign = :left)

    Label(fig.layout[1, 3, TopLeft()], "(c)",
    fontsize = 12,
    padding = (0, -15, 5, 0),
    halign = :left)
    resize_to_layout!(fig)
    fig
end

#____________________________________________________________________________________________________________________
# Figure 5
#____________________________________________________________________________________________________________________

function lat_fig5(folder)
    with_theme(merge(theme_latexfonts(), Theme(fontsize = 15))) do #
        figure5(folder)     
    end
end

function figure5(folder) 
    fig = Figure(resolution = (1800/1.5, 1000/1.5))
    gqah = GridLayout(fig[1,1])
    gvh = GridLayout(fig[1,2])

    Label(gqah[1, 1:2, Top()], "QAH", valign = :bottom,
        padding = (0, 20, 5, 0), fontsize =22)

    Label(gvh[1, 1:2, Top()], "VH/SVH", valign = :bottom,
        padding = (0, 20, 5, 0), fontsize =22)

    ax_1 = Axis(gqah[1,1], xlabel = L"$\omega$ [meV]", ylabel = L"$\nu$")
    ax_2 = Axis(gqah[1,2], xlabel = L"$\omega$ [meV]", ylabel = L"$\nu$")
    ax_3 = Axis(gqah[2,1], xlabel = L"$\omega$ [meV]", ylabel = L"$\nu$")
    ax_4 = Axis(gqah[2,2], xlabel = L"$\omega$ [meV]", ylabel = L"$\nu$")
    
    ax_5 = Axis(gvh[1,1], xlabel = L"$\omega$ [meV]", ylabel = L"$\nu$")
    ax_6 = Axis(gvh[1,2], xlabel = L"$\omega$ [meV]", ylabel = L"$\nu$")
    ax_7 = Axis(gvh[2,1], xlabel = L"$\omega$ [meV]", ylabel = L"$\nu$")
    ax_8 = Axis(gvh[2,2], xlabel = L"$\omega$ [meV]", ylabel = L"$\nu$")

    valshift = 0.1e5

    hm1 = correlateddensityplot!(ax_1, folder, :REAL, :x, :x, :x, "shift", rs = "m1", gs = "vp", colrange = valshift)
    hm2 = correlateddensityplot!(ax_2, folder, :REAL, :y, :y, :y, "shift", rs = "m1", gs = "vp", colrange = valshift)
    hm3 = correlateddensityplot!(ax_3, folder, :REAL, :x, :x, :x, "inj", rs = "m1", gs = "vp", colrange = 4)
    hm4 = correlateddensityplot!(ax_4, folder, :REAL, :y, :y, :y, "inj", rs = "m1", gs = "vp", colrange = 4)

    hm5 = correlateddensityplot!(ax_5, folder, :REAL, :x, :x, :x, "shift", rs = "m1", gs = "qah", colrange = valshift)
    hm6 = correlateddensityplot!(ax_6, folder, :REAL, :y, :y, :y, "shift", rs = "m1", gs = "qah", colrange = valshift)
    hm7 = correlateddensityplot!(ax_7, folder, :REAL, :x, :x, :x, "inj", rs = "m1", gs = "qah", colrange = 4)
    hm8 = correlateddensityplot!(ax_8, folder, :REAL, :y, :y, :y, "inj", rs = "m1", gs = "qah", colrange = 4)

    Label(gqah[1, 1, Top()], L"$\sigma_{xxx}^\text{sh}$", valign = :bottom, padding = (0, 0, -30, 0), fontsize =25)
    Label(gqah[1, 2, Top()], L"$\sigma_{yyy}^\text{sh}$", valign = :bottom, padding = (0, 0, -30, 0), fontsize =25)
    Label(gvh[1, 1, Top()], L"$\sigma_{xxx}^\text{sh}$", valign = :bottom, padding = (0, 0, -30, 0), fontsize =25)
    Label(gvh[1, 2, Top()], L"$\sigma_{yyy}^\text{sh}$", valign = :bottom, padding = (0, 0, -30, 0), fontsize =25)

    Label(gqah[2, 1, Top()], L"$\sigma_{xxx}^\text{inj}$", valign = :bottom, padding = (0, 0, -30, 0), fontsize =25)
    Label(gqah[2, 2, Top()], L"$\sigma_{yyy}^\text{inj}$", valign = :bottom, padding = (0, 0, -30, 0), fontsize =25)
    Label(gvh[2, 1, Top()], L"$\sigma_{xxx}^\text{inj}$", valign = :bottom, padding = (0, 0, -30, 0), fontsize =25)
    Label(gvh[2, 2, Top()], L"$\sigma_{yyy}^\text{inj}$", valign = :bottom, padding = (0, 0, -30, 0), fontsize =25)

    colgap!(gvh, 10)
    colgap!(gqah, 10)

    rowgap!(gvh, 10)
    rowgap!(gqah, 10)

    hideydecorations!(ax_2, grid = false)
    hideydecorations!(ax_4, grid = false)
    hideydecorations!(ax_6, grid = false)
    hideydecorations!(ax_8, grid = false)
    hidexdecorations!(ax_1, grid = false)
    hidexdecorations!(ax_2, grid = false)
    hidexdecorations!(ax_5, grid = false)
    hidexdecorations!(ax_6, grid = false)

    valstr = Int(log10(valshift));

    cb = Colorbar(gvh[1,3], hm5, label = L" $\sigma_{iii}^\text{sh} \text{ }  [10^{4} \times  \text{nm} \mu A/V^2]$", vertical = true, flipaxis = true, ticklabelpad= -1)
    cb = Colorbar(gvh[2,3], hm3, label = L" $\sigma_{iii}^\text{inj}/\tau \text{ } [\text{nm} \mu A/V^2 \text{fs}]$", vertical = true, flipaxis = true, ticklabelpad= -1)
  
    Label(gqah[1, 1, Top()], "(a)",
    fontsize = 16,
    padding = (3, 0, -50, 0),
    halign = :left)

    Label(gqah[1, 2, Top()], "(b)",
    fontsize = 16,
    padding = (3, 0, -50, 0),
    halign = :left)

    Label(gqah[2, 1, Top()], "(c)",
    fontsize = 16,
    padding = (3, 0, -30, 0),
    halign = :left)

    Label(gqah[2, 2, Top()], "(d)",
    fontsize = 16,
    padding = (3, 0, -30, 0),
    halign = :left)

    Label(gvh[1, 1, Top()], "(e)",
    fontsize = 16,
    padding = (3, 0, -50, 0),
    halign = :left)

    Label(gvh[1, 2, Top()], "(f)",
    fontsize = 16,
    padding = (3, 0, -50, 0),
    halign = :left)


    Label(gvh[2, 1, Top()], "(g)",
    fontsize = 16,
    padding = (3, 0, -30, 0),
    halign = :left)

    Label(gvh[2, 2, Top()], "(h)",
    fontsize = 16,
    padding = (3, 0, -30, 0),
    halign = :left)
    fig
end
#____________________________________________________________________________________________________________________
# Figure 6
#____________________________________________________________________________________________________________________

function lat_fig6(folder,  ν= 2.5,ω = 24)
    with_theme(merge(theme_latexfonts(), Theme(fontsize = 18))) do #
        figure6(folder, ν ,ω)     
    end
end

function figure6(folder, ν ,ω)
    fig = Figure()
    gqah = GridLayout(fig[1,1])

    Label(gqah[1, 1, Top()], "QAH", valign = :bottom,
        padding = (0, 20, 5, 0), fontsize =17)

    Label(gqah[1, 2, Top()], "VH/SVH", valign = :bottom,
        padding = (0, 20, 5, 0), fontsize =17)
    val = ν
    latex_str = L"$\nu = %$(val)$"
    latex_str2 = L"$\omega = %$(ω) \text{ } [meV]$"

    ax_1 = Axis(gqah[1,1], xlabel = L"$\omega$ [meV]", ylabel = L"$\sigma_{xxx} \text{ }  [10^{4} \times  \text{nm} \mu A/V^2]")
    ax_2 = Axis(gqah[1,2], xlabel = L"$\omega$ [meV]", ylabel = L"$\sigma_{xxx} \text{ }  [10^{4} \times  \text{nm} \mu A/V^2]")
    
    valshift = 0.1e5
    filling = ν
    ωs, shxxx = correlatedlcplot!(folder, filling, :REAL, :x, :x, :x, "shift", rs = "m1", gs = "vp", colrange = valshift)
    _, injxxx = correlatedlcplot!(folder, filling, :REAL, :x, :x, :x, "inj", rs = "m1", gs = "vp")

    lines!(ax_1, ωs, -shxxx, color = :gray, label = L"$\sigma_{xxx}^\text{shift}")
    lines!(ax_1, ωs, -injxxx .* (200/1e4), color = :orange, label = L"$\sigma_{xxx}^\text{inj}")
    lines!(ax_1, ωs, -shxxx + -injxxx .* (200/1e4), color = :black, label = L"$\sigma_{xxx}")
    ωs, shxxx = correlatedlcplot!(folder, filling, :REAL, :x, :x, :x, "shift", rs = "m1", gs = "qah", colrange = valshift)
    _, injxxx = correlatedlcplot!(folder, filling, :REAL, :x, :x, :x, "inj", rs = "m1", gs = "qah")
    lines!(ax_2, ωs, -shxxx, color = :gray)
    lines!(ax_2, ωs, -injxxx .* (200/1e4), color = :orange)
    lines!(ax_2, ωs, -shxxx + -injxxx .* (200/1e4), color = :black)
    freq = ω
    ωs, shxxx = correlatedlcplotycut!(folder, freq, :REAL, :x, :x, :x, "shift", rs = "m1", gs = "vp", colrange = valshift)
    _, injxxx = correlatedlcplotycut!(folder, freq, :REAL, :x, :x, :x, "inj", rs = "m1", gs = "vp")
  
    Label(gqah[1, 1, Top()], "(a)",
    fontsize = 14,
    padding = (3, 0, -45, 0),
    halign = :left)

    Label(gqah[1, 2, Top()], "(b)",
    fontsize = 14,
    padding = (3, 0, -45, 0),
    halign = :left)

    axislegend(ax_1, orientation = :horizontal,   
    labelsize = 18,                               
    framevisible = true,                          
    patchsize = (12, 6), colgap = 5,
    position = (1.08,1.02),
    framecolor = (:black, 0.5),
    backgroundcolor = (:white, 0.7),
    padding = 1.5 .* (2, 2, 2, 2))

    hideydecorations!(ax_2, grid = false)
    xlims!(ax_1, (0,40))
    xlims!(ax_2, (0,40))
    ylims = (-1,1)
    ylims!(ax_1, ylims)
    ylims!(ax_2, ylims)
    fig
end

function correlatedlcplotycut!(folder, freq, part, a, b, c, which_current = "shift"; rs = "", gs = "", colrange = 1e5)
    j = 1
    if which_current == "shift"
        scale2 = 1
        scale = colrange
        j *= -1
        func = (x -> x)
    else
        scale2 = colrange
        scale = 1e-15
        func = (x -> x)
    end
    if part == :REAL
        which_current *= "_real_"
    else 
        which_current *= "_imag_"
    end
    which_current *= string(a) * string(b) * string(c) * rs
    str = string(folder) * "/" * which_current
    ns =  flatten(parse_str_to_arr.(import_array(str*"cond_ns"*gs*".csv"))[:])
    inds = sortperm(ns)
    ωs = flatten(parse_str_to_arr.(import_array(str*"cond_omegas"*gs*".csv"))[:])
    l = parse_str_to_arr.(import_array(str*"cond_mat"*gs*".csv"))
    mat = zeros(length(l), length(l[1]))
    for i in 1:length(l)
        mat[i,:] .= l[i]
    end

    matt = hcat(j .* func.(mat'[:,inds[end:-1:1]])./scale,func.(mat'[:,inds])./scale)
    nss = vcat(-ns[inds[end:-1:1]], ns[inds])
    idx = argmin(abs.(ωs .- freq))
    return nss, matt[idx,:]
end


