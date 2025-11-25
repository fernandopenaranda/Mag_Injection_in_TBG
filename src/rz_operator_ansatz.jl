# https://www.youtube.com/watch?v=EjZR7z9VH-w real space k-space
# implementation of the non-local rz position operator. 

# Since the WF at K (K') are polarized in the top (bottom) layer 
# rz(K) (rz(K')) must return 1 (-1) provided that at K and K' the lower bands
# are purely composed of f-electrons.
"""
it only affects the f'electrons. This is possibly not entirely true.
Also there is one 1/ξ possibly missing, be aware of that.
"""
function hf_rz(p, k)
    ham_dim = ham_matrix_size(gs_indices(p))
    rz_mat = spzeros(ComplexF64, ham_dim, ham_dim)
    rz_mat[1:2,1:2] .= σ0 .* non_local(p, k)
    return rz_mat
end
rz(ψs, p, k) = ψs' * hf_rz(p,k) * ψs

function hf_rz(p, k, of)
    ham_dim = ham_matrix_size(gs_indices(p))
    rz_mat = spzeros(ComplexF64, ham_dim, ham_dim)
    rz_mat[1:2,1:2] .= σ0 .* non_local(p, k)
    return [rz_mat 0I 0I 0I; 0I -rz_mat 0I 0I; 0I 0I rz_mat 0I; 0I 0I 0I -rz_mat] 
end
rz(ψs, p, k, of) = ψs' * hf_rz(p, k, of) * ψs

function non_local(p, k)
    g1, g2 = bravais_vectors(p)
    K1 =  κ(g1, g2)
    mod_a_n = 0.542*π*a0/(2*1.703*sin(1.05*π/360))
    val = 0.0im
    for n in [0,1,2]
        a_n = mod_a_n .* [real(exp(1im*(n*2π/3+π/2))), imag(exp(1im*(n*2π/3+π/2)))] 
        # three real space vectors perpendicular to the moir'e momentum transfer
        val += sin(k' * a_n)
    end
    return  val * non_local_pref(p)
end

function d_non_local(p, k, a)
    mod_a_n = 0.542*π*a0/(2*1.703*sin(1.05*π/360))
    val = 0.0im
    if a == :x 
        ind = 1
    elseif a == :y
        ind = 2
    else ind = 0 end
    for n in [0,1,2] # first neighbours
        a_n = mod_a_n .* [real(exp(1im*(n*2π/3+π/2))), imag(exp(1im*(n*2π/3+π/2)))] 
        # three real space vectors perpendicular
        # to the moir'e momentum transfer
        val += cos(k' * a_n) * a_n[ind]
    end
    return  val * non_local_pref(p)
end


non_local_pref(p) = p.ν /(2sin(-π/sqrt(3))+sin(2π/sqrt(3)))/1.0787840086239031

function hf_rz_matrix_elements_in_K_Kp(p)
    g1, g2 = bravais_vectors(p)
    K1 =  κ(g1, g2)
    K2 =  κ(g2, g1)
    vecsK1 = f_eigenstates(p, K1) 
    vecsK2 = f_eigenstates(p, K2) 
    matK1 = hf_rz(p, K1)
    matK2 = hf_rz(p, K2)
    matelements_K = [vecsK1[1]' * matK1 * vecsK1[1] vecsK1[1]' * matK1 * vecsK1[2]; 
        vecsK1[2]' * matK1 * vecsK1[1]  vecsK1[2]' * matK1 * vecsK1[2]]
    matelements_Kp = [vecsK2[1]' * matK2 * vecsK2[1] vecsK2[1]' * matK2 * vecsK2[2]; 
        vecsK2[2]' * matK2 * vecsK2[1]  vecsK2[2]' * matK2 * vecsK2[2]]
    return matelements_K, matelements_Kp
end

function f_eigenstates(p, k)
    es, evecs = eigen(Matrix(hf_hamiltonian(p, [k[1], k[2]])))
    half_dim = size(evecs,1)÷2
    return [evecs[:,i] for i in half_dim:half_dim+1] # lowest_evecs
end

##########
# Function to define a regular hexagon centered at the origin with distance a0 from center to vertex
 hexagon_vertices(K1) = [norm(K1) * Complex(cos(π/3 * i+π/2), sin(π/3 * i+π/2)) for i in 0:5]

# Function to check if a point is inside the hexagon
function inside_hexagon(x, y, K)
    x_abs, y_abs = abs(x), abs(y)
    return x_abs <= K[1]  && y_abs <= norm(K) && (sqrt(3) * y_abs + x_abs <= sqrt(3) * norm(K))
end

function hexagon_plot(p)
    f = Figure(resolution = (700, 700))
    ax = Axis(f[1,1], xlabel = "kx", ylabel = "ky")
    g1, g2 = bravais_vectors(p)
    K1 =  κ(g1, g2)
    K2 =  κ(g2, g1)
    mod_a_n = 2π/abs(K1[1])
    re_vertices = real.(hexagon_vertices(K1))
    imag_vertices = imag.(hexagon_vertices(K1))
    lines!(ax, push!(re_vertices, re_vertices[1]), push!(imag_vertices, imag_vertices[1]),
        color = :black, linewidth = 0.5)
    x_min, x_max = -norm(K1), norm(K1)
    y_min, y_max = -norm(K1),norm(K1)

    x_grid = range(x_min, x_max; length=200)
    y_grid = range(y_min, y_max; length=200)
    
    points = [(x, y) for x in x_grid, y in y_grid if  inside_hexagon(x, y, K1)]
    x_values = [p[1] for p in points]
    y_values = [p[2] for p in points]
    z_values = [non_local(p, [kx, ky]) for (kx, ky) in points]
    scatter!(ax, x_values, y_values, markersize=7, color=real.(z_values), colormap=:redblue)
    f
end

function z_covariant(zmat, rmat)
    dim = size(zmat, 1)
    mat = zeros(ComplexF64, size(rmat,1),size(rmat,1))
    @inbounds begin
        for n in 1:dim
            for m in 1:dim
                if m != n && n > m 
                    for p in 1:dim
                        if  p == n 
                            val = zmat[n,p] * rmat[p,m] 
                        elseif p == m 
                            val = -rmat[n,p] * zmat[p, m]
                        else
                            val = zmat[n,p] * rmat[p,m] - rmat[n,p] * zmat[p, m]
                        end
                        mat[n,m] += val    
                    end
                else nothing end
            end
        end
    end
    return -1im .* mat
end