function plot_jdos(ωlist, jdos, color, ylims)
    println("color ", color)
    fig = Figure(); 
    ax = Axis(fig[1, 1]; 
        xlabel = "ω [meV]", 
        ylabel = "JDOS (arb. units.)",
        xlabelsize = 24, 
        ylabelsize = 24)
    lines!(ax, ωlist, jdos ./ maximum(jdos), color = color)
    fig
end

function plot_jdos!(ax, ωlist, jdos, color, ylims)
    lines!(ax, ωlist, jdos ./ maximum(jdos) , color = color)
    hidedecorations!(ax)
end

plot_linear_conductivity(ωlist, conds, a, b; part = :real, kw...) = 
    plot_linear_conductivity(ωlist, conds,
        string((part == :real ? "Re[" : "Im["),"σ_ ", string(a),string(b), "] "); kw...)

function plot_linear_conductivity(ωlist, conds, strobs, color)
    println("color ", color)
    fig = Figure(); 
    ax = Axis(fig[1, 1]; 
        xlabel = "ω [meV]", 
        ylabel = string(strobs, "/σ_mono", "(ω)" ),
        xlabelsize = 24, 
        ylabelsize = 24)
    lines!(ax, ωlist, conds, color = color)
    fig
end

plot_linear_conductivity!(fig, ωlist, conds, a, b; part = :real, kw...) = 
    plot_linear_conductivity!(fig, ωlist, conds,
        string((part == :real ? "Re[" : "Im["),"σ_ ", string(a),string(b), "] "); kw...)

function plot_linear_conductivity!(ax, ωlist, conds, color)
    lines!(ax, ωlist, conds, color = color)
    hidedecorations!(ax)
end


function plot_shift_current(part, a,b,c, ωlist, vals1, vals2)
    ylab = ylabel_shift(a,b,c, part)
    println(maximum(vals1))
    scale = ()
    fig = Figure(); 
    ax = Axis(fig[1, 1]; 
        xlabel = "ω [meV]", 
        ylabel = ylab,
        xlabelsize = 24, 
        ylabelsize = 24)
    lines!(ax, ωlist, vals1 , color = :orange) # ./1000 meV -> eV
    lines!(ax, ωlist, vals2, color = :gray)
    lines!(ax, ωlist, (vals1+vals2) , color = :black)
    # ylims!(ax, [-0.01,0.3])
    fig
end

function ylabel_shift(a,b,c, part) 
    if part == :REAL
        if a == b == c && a == :x
            ylab = L"$Re[\sigma_{xxx}]$ [nm $\mu$A/$V^2$]"
        elseif a == b == c && a == :y
            ylab = L"$Re[\sigma_{yyy}]$ [nm $\mu$A/$V^2$]"
        elseif a == :x && b == :x && c == :y
            ylab = L"$Re[\sigma_{xxy}]$ [nm $\mu$A/$V^2$]"
        elseif a == :x && b == :y && c == :y
            ylab = L"$Re[\sigma_{xyy}]$ [nm $\mu$A/$V^2$]"
        elseif a == :x && b == :y && c == :x
            ylab = L"$Re[\sigma_{xyx}]$ [nm $\mu$A/$V^2$]"
        elseif a == :y && b == :x && c == :y
            ylab = L"$Re[\sigma_{yxx}]$ [nm $\mu$A/$V^2$]"
        elseif a == :y && b == :x && c == :y
            ylab = L"$Re[\sigma_{yxy}]$ [nm $\mu$A/$V^2$]"
        elseif a == :y && b == :x && c == :x
            ylab = L"$Re[\sigma_{yxx}]$ [nm $\mu$A/$V^2$]"
        elseif a == :y && b == :y && c == :x
            ylab = L"$Re[\sigma_{yyx}]$ [nm $\mu$A/$V^2$]"
        else nothing end
    else
        if a == b == c && a == :x
            ylab = L"$Im[\sigma_{xxx}]$ [nm $\mu$A/$V^2$]"
        elseif a == b == c && a == :y
            ylab = L"$Im[\sigma_{yyy}]$ [nm $\mu$A/$V^2$]"
        elseif a == :x && b == :x && c == :y
            ylab = L"$Im[\sigma_{xxy}]$ [nm $\mu$A/$V^2$]"
        elseif a == :x && b == :y && c == :y
            ylab = L"$Im[\sigma_{xyy}]$ [nm $\mu$A/$V^2$]"
        elseif a == :x && b == :y && c == :x
            ylab = L"$Im[\sigma_{xyx}]$ [nm $\mu$A/$V^2$]"
        elseif a == :y && b == :x && c == :y
            ylab = L"$Im[\sigma_{yxx}]$ [nm $\mu$A/$V^2$]"
        elseif a == :y && b == :x && c == :y
            ylab = L"$Im[\sigma_{yxy}]$ [nm $\mu$A/$V^2$]"
        elseif a == :y && b == :y && c == :x
            ylab = L"$Im[\sigma_{yyx}]$ [nm $\mu$A/$V^2$]"
        else nothing end
    end
end

function plot_injection_current(part, a,b,c, ωlist, vals1)
    ylab = ylabel_injection(a,b,c, part)
    println(maximum(vals1))
    scale = 1e-15
    fig = Figure(); 
    ax = Axis(fig[1, 1]; 
        xlabel = "ω [meV]", 
        ylabel = ylab,
        xlabelsize = 24, 
        ylabelsize = 24)
    lines!(ax, ωlist, vals1 ./ (scale), color = :black) 
    ylims!(ax, [-20,20])
    fig
end

function plot_injection_current!(ax, part, a,b,c, ωlist, vals1; color = :black, label = false)
    scale = 1e-15
    lines!(ax, ωlist, vals1 ./ (scale), color = color, label = label)
    # scatter!(ax, ωlist, vals1 ./ (scale), color = color) 
    ylims!(ax, [-1200,1200])
end

function plot_injection_current(part, a,b,c, ωlist, vals1, vals2)
    ylab = ylabel_injection(a,b,c, part)
    println(maximum(vals1))
    scale = 1e-15
    fig = Figure(); 
    ax = Axis(fig[1, 1]; 
        xlabel = "ω [meV]", 
        ylabel = ylab,
        xlabelsize = 24, 
        ylabelsize = 24)
    lines!(ax, ωlist, vals1 ./ (scale), color = :orange) 
    lines!(ax, ωlist, vals2 ./ (scale), color = :gray)
    lines!(ax, ωlist, (vals1+vals2) ./(scale), color = :black)
    ylims!(ax, [-20,20])
    fig
end

function plot_injection_current_vs_μ!(ax1, ax2, ax3, μlist, vals1, vals2; color = :gray, ylims = [-20,20], scale = 1e-15)
    lines!(ax1, μlist, vals1 ./ (scale), color = color)
    lines!(ax2, μlist, vals2 ./ (scale), color = color)
    lines!(ax3, μlist, (vals1+vals2) ./(scale), color = color)
    ylims!(ax1, ylims)
    ylims!(ax2, ylims)
    ylims!(ax3, ylims)
end


function ylabel_injection(a,b,c, part) 
    if part == :REAL
        if a == b == c && a == :x
            ylab = L"$[\sigma^{inj}_{xxx}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == b == c && a == :y
            ylab = L"$Re[\sigma^{inj}_{yyy}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == :x && b == :x && c == :y
            ylab = L"$Re[\sigma^{inj}_{xxy}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == :x && b == :y && c == :y
            ylab = L"$Re[\sigma^{inj}_{xyy}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == :x && b == :y && c == :x
            ylab = L"$Re[\sigma^{inj}_{xyx}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == :y && b == :x && c == :y
            ylab = L"$Re[\sigma^{inj}_{yxx}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == :y && b == :x && c == :y
            ylab = L"$Re[\sigma^{inj}_{yxy}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == :y && b == :y && c == :x
            ylab = L"$Re[\sigma^{inj}_{yyx}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        else nothing end
    else
        if a == b == c && a == :x
            ylab = L"$Im[\eta^{inj_{xxx}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == b == c && a == :y
            ylab = L"$Im[\eta^{inj_{yyy}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == :x && b == :x && c == :y
            ylab = L"$Im[\eta^{inj_{xxy}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == :x && b == :y && c == :y
            ylab = L"$Im[\eta^{inj_{xyy}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == :x && b == :y && c == :x
            ylab = L"$Im[\eta^{inj_{xyx}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == :y && b == :x && c == :y
            ylab = L"$Im[\eta^{inj_{yxx}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == :y && b == :x && c == :y
            ylab = L"$Im[\eta^{inj_{yxy}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        elseif a == :y && b == :y && c == :x
            ylab = L"$Im[\eta^{inj_{yyx}/\tau]$ [nm $\mu$A/$(V^2 fs)$]"
        else nothing end
    end
end

#= 

    DENSITY PLOTS

=#
function densityplot(str1, str2, str3, str4; colrange = (-20,20), scale = 1e-15)
    fig = Figure(resolution = (2000, 600), fontsize = 24)
    with_theme(theme_latexfonts()) do
    μs = Matrix(CSV.read(str1, DataFrame))
    ωs = Matrix(CSV.read(str2, DataFrame)) 
    mat1 = Matrix(CSV.read(str3, DataFrame))
    mat2 = Matrix(CSV.read(str4, DataFrame))
    ax1 = Axis(fig[1, 1],
        title = "+ valley",
        # ylabel = "μ [meV]",
        ylabel = "ν",
        xlabel = "ω [meV]",
        xlabelsize = 24,
        ylabelsize = 24
    )
    ax2 = Axis(fig[1, 2],
    title = "- valley",
    # ylabel = "μ [meV]",
    ylabel = "ν",
    xlabel = "ω [meV]",
    xlabelsize = 24,
    ylabelsize = 24
)
    ax3 = Axis(fig[1, 3],
    title = "both valleys",
    # ylabel = "μ [meV]",
    ylabel = "ν",
    xlabel = "ω [meV]",
    xlabelsize = 24,
    ylabelsize = 24
    )
    hideydecorations!(ax2)
    hideydecorations!(ax3)
    colrange = colrange
    hm = heatmap!(ax1 , ωs, -μs , (mat1./scale), colorrange = colrange, colormap = :balance)
    hm = heatmap!(ax2 , ωs, -μs , (mat2./scale), colorrange = colrange, colormap = :balance)
    hm = heatmap!(ax3 , ωs, -μs , ((mat1+mat2)./scale), colorrange = colrange, colormap = :balance)
    Colorbar(fig[1, 4], hm, width = 30, label = L"$\sigma_{yyy}^\text{shift} [\text{nm} \mu A/V^2])$", labelsize = 30)
    return fig
    end
end

function densityplot(str1, str2, str3; colrange = (-20,20), scale = 1e-15)
    fig = Figure(resolution = (700, 600), fontsize = 24)
    with_theme(theme_latexfonts()) do
    μs = Matrix(CSV.read(str1, DataFrame))
    ωs = Matrix(CSV.read(str2, DataFrame)) 
    mat1 = Matrix(CSV.read(str3, DataFrame))
    ax1 = Axis(fig[1, 1],
        title = "+ valley",
        ylabel = "ν",
        xlabel = "ω [meV]",
        xlabelsize = 24,
        ylabelsize = 24)
    colrange = colrange
    hm = heatmap!(ax1 , ωs, -μs , (mat1./scale), colorrange = colrange, colormap = :grays)
    Colorbar(fig[1, 2], hm, width = 30, label = L"$\sigma_{xxx}^\text{shift} [\text{nm} \mu A/V^2])$", labelsize = 30)
    return fig
    end
end