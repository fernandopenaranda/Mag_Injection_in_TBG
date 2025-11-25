"""
Kramers intervalley coherent groundstate at nu =0
We take the one-shot approach. 
cc -> J σyτy
ff -> U1 σyτy
"""
function kivc!(h_mat, p)
    dim = size(h_mat, 1) ÷ 2
    h_mat[1:2, dim+1:dim+2] .+= -1im * p.U1/2 .* σy #\sigmay \tauy
    h_mat[dim+1:dim+2, 1:2] .+= 1im * p.U1/2 .* σy
    g = zeros(Float64, 2)
    ham_dim = size(h_mat, 1)
    g1, g2 = bravais_vectors(p)
    cc_mat = 1im * p.J .* σy
    for i in 3:4:dim-3
        h_mat[i+2:i+3, dim+i+2:dim+i+3] .+= cc_mat
        h_mat[dim+i+2:dim+i+3,i+2:i+3] .+= -cc_mat
    end
end