#using Implicit3DPlotting #https://github.com/matthiashimmelmann/Implicit3DPlotting.jl
# using GLMakie

# function ϵ_k(p, of)
#     eigs = real.(sc_eigs(p, k, diagm(of))[1])
#     return eigs[hdim÷2-3:hdim÷2+4] 
# end

function klims(p)
 return (-0.026746521118468435, 0.026746521118468435), 2/sqrt(3) .* (-0.026746521118468435, 0.026746521118468435)
end

""" 
plots the flat bands 
"""
function plot_ϵk(p, of; points = 20)
    xs = LinRange(klims(p)[1][1], klims(p)[1][2], points)
    ys = LinRange(klims(p)[2][1], klims(p)[2][2], points)
    zs = zeros(Float64,length(xs), length(ys))
    fig = Figure(size = (2000,800))
    ϵ_k(k, i) = real.(sc_eigs(p, k, diagm(of)).values)[120÷2-3:120÷2+4][i] 
    axlisttop = [Axis3(fig[1, i], xlabel = "k_x", ylabel = "k_y", zlabel = "E [meV]") for i in 1:4]
    axlistbot = [Axis3(fig[2, i-4],xlabel = "k_x", ylabel = "k_y", zlabel = "E [meV]") for i in 5:8]
    for i in 1:4
        zs .= [ϵ_k([x,y],i) for x in xs, y in ys]
        surface!(axlisttop[i],xs, ys, zs, colormap = :deep)
        hidedecorations!(axlisttop[i])

        zs .= [ϵ_k([x,y],i+4) for x in xs, y in ys] 
        surface!(axlistbot[i],xs, ys, zs, colormap = :deep)
        hidedecorations!(axlistbot[i])

    end
    
    return fig
end


function plot_ϵk(p, of, it, points = 20)
    xs = LinRange(klims(p)[1][1], klims(p)[1][2], points)
    ys = LinRange(klims(p)[2][1], klims(p)[2][2], points)
    zs = zeros(Float64,length(xs), length(ys))
    fig = Figure()
    ϵ_k(k, i) = real.(sc_eigs(p, k, diagm(of)).values)[120÷2-3:120÷2+4][i] 
    ax = Axis3(fig[1, 1], title = "Band $(it)", xlabel = "X", ylabel = "Y", zlabel = "Z")
    zs .+= [ϵ_k([x,y],it) for x in xs, y in ys]
    surface!(ax,xs, ys, zs, colormap = :deep)
    return fig
end


"at which μ, n is n(μ)"
function  n(p, of, ;ϵ = 4e-1, evals = 4000)
    Δn(mu) = ν - n(p, of, mu, ϵ = ϵ, evals = evals)[1][1]
    result = optimize(Δn, -60, 60; show_trace = true, abs_tol = 1e-3, iterations = 10)
    return Optim.minimizer(result)
end

"""
occupation for a given mu integrating the flat bands
"""
function n(p, of, μ; ϵ = 4e-1, evals = 4000)
    M, xmin, xmax = int_boundaries(p)
    hdim = ham_matrix_size(p)
    ϵ_ki(k) = real.(sc_eigs(p, k, diagm(of))[1])[hdim÷2-3:hdim÷2+4][1] # just one band
    dirac_delta(k, E) =1/(√(4π*ϵ)) * exp(-0.5*(E-ϵ_ki(k))^2/ϵ^2)
    
    Δx = [0, xmin[2]]
    minE = ϵ_ki([0.,0.]) # assumption that the minimum in energy is at the gamma point
    push!(Δx,minE + minE/5 )
    Δy = [xmax[1]/2, xmax[2]]
    push!(Δy, μ)

    integrand(vec) = dirac_delta([vec[1], vec[2]], vec[3])
    val, err = Cubature.hcubature(1, (x,v) -> v[:] .= integrand(x), Δx, Δy; maxevals = evals)
    return ϵ/abs(Δx[3]-Δy[3])*1/abs.(xmin[2]*(xmax[2]-xmax[1]/2)) * val, err
end

