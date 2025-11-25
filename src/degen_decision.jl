# after running: run_random.jl
# for the sake of the computation of the injection current we must choose among the degenerate solutions
# at a given filling. The Hamiltonian has a symmetry with respect to valley exchange but the injection current
# does not. It changes sign upon this operation, so we choose how to fix the valley after compoting the f-density
# matrices.
""" 
Returns a `Self_consistent_data` object with a given convention for the filling of the valleys between one of the 
4 possible `methods` below. 
the two valleys are degenerate. So we choose in the unperturbed solution one fix filling so the injection current
is not oscillating between manifolds with the same dof.
Regarding the shift current there are 4 degenerate filling options:
    (1) First valley +
    (2) First valley -
    (3) First valley + for positive and - for negative
    (4) First valley - for positive and + for negative
    (5) + - + - + - + - as a function of the filling. To do alternating pattern this would break phs

# la injection es insensitive al spin por lo que no hago nada con el
GENERAL REMARKS ON DEGENERACIES OF H(of):
# [1 ,2 ,3, 4, 5, 6, 7, 8] = [3, 4, 1, 2, 7, 8, 5, 6] = [1, 2, 7, 8, 3, 4, 5, 6] = [5, 6, 3, 4, 1, 2, 7, 8] = [5, 6, 7, 8, 1, 2, 3, 4]
# [1, 2, 3, 4, 5, 6,  7, 8] ≠ [3, 4, 1, 2, 5, 6, 7, 8] Esta operación lleva la solución TRS sym a la broken-TRS sym
"""
function fix_valley_convention_ofmats(s::Self_consistent_data; method = "1")
    reshuffled_ofmats = [zeros(8) for i in 1:length(s.ofmats)]
    for i in 1:length(s.ofmats) #populate the valleys in the same fashion
        if method == "1"
            if sum(s.ofmats[i][[1,2,5,6]]) > sum(s.ofmats[i][[3,4,7,8]])
                reshuffled_ofmats[i] .= s.ofmats[i]
            else
                reshuffled_ofmats[i] .= s.ofmats[i][[3,4,1,2,7,8,5,6]]
            end
        elseif method == "2"
            if sum(s.ofmats[i][[1,2,5,6]]) < sum(s.ofmats[i][[3,4,7,8]])
                reshuffled_ofmats[i] .= s.ofmats[i]
            else
                reshuffled_ofmats[i] .= s.ofmats[i][[3,4,1,2,7,8,5,6]]
            end
        elseif method == "3"
            if sign(s.ns[i])*sum(s.ofmats[i][[1,2,5,6]]) > sign(s.ns[i])*sum(s.ofmats[i][[3,4,7,8]])
                reshuffled_ofmats[i] .= s.ofmats[i]
            else
                reshuffled_ofmats[i] .= s.ofmats[i][[3,4,1,2,7,8,5,6]]
            end
        elseif method == "4"
            if sign(s.ns[i])*sum(s.ofmats[i][[1,2,5,6]]) < sign(s.ns[i])*sum(s.ofmats[i][[3,4,7,8]])
                reshuffled_ofmats[i] .= s.ofmats[i]
            else
                reshuffled_ofmats[i] .= s.ofmats[i][[3,4,1,2,7,8,5,6]]
            end
        end
    end
        Self_consistent_data(reshuffled_ofmats, s.mus, s.es, s.ns)
end


function provide_folder_get_observables_atlas(folder)
    str = pwd() * "/output" * "/Data" * "/" * string(folder)
    es_path = str * "/joined_es.csv"
    mus_path = str * "/joined_mus.csv"
    ns_path = str * "/joined_ns.csv"
    ofmats_path = str * "/joined_ofmats.csv"
    s = return_variables(es_path,mus_path, ns_path, ofmats_path )
    Self_consistent_data(s[1], s[2], s[3], s[4])
end

function save_sbp(PID, s, rs = "")
    str = "/scratch/ferpe/HFCorrelated/cluster_temp/src/output/Data"
    str2 = string(str,"/",  PID)
    if !isdir(str2)
        mkdir(str2)
    end
     str3 = str2
    ofvec = []
    for mat in s.ofmats
        push!(ofvec, mat)
    end
    save_arrays(string(str3, "/joined_ns_"*rs*".csv"), s.ns)
    save_arrays(string(str3, "/joined_mus_"*rs*".csv"), s.mus)
    save_arrays(string(str3, "/joined_es_"*rs*".csv"), s.es)
    save_arrays(string(str3, "/joined_ofmats_"*rs*".csv"), ofvec)
end

function store_reshuffled(folder)
	s = provide_folder_get_observables_atlas(folder)
	s2 = fix_valley_convention_ofmats(s; method = "1")
	save_sbp(folder, s2)
end

function store_vpqah(folder)
        s = provide_folder_get_observables_atlas(folder)
        svp, sqah = generate_vp_aqh_solutions(s)
        save_sbp(folder, svp, "vp")
	save_sbp(folder, sqah, "qah")
end

function generate_vp_aqh_solutions(s::Self_consistent_data)
    s2 = fix_valley_convention_ofmats(s, method = "1");
    s3 = fix_spin_convention_ofmats(s2)
    s4 = fix_spin_convention2_ofmats(s3)
    s5 = permute_vp_aqh_character(s4)

    aux_vp = [zeros(Float64, 8) for i in 1:length(s4.ofmats)]
    aux_qah = [zeros(Float64, 8) for i in 1:length(s4.ofmats)]

    for i in 1:length(s.ofmats)
        if abs(s4.ofmats[i][5]-s4.ofmats[i][3]) > 1e-2
            if s4.ofmats[i][5]-s4.ofmats[i][3] > 0
                aux_vp[i][:]  .= s4.ofmats[i]  
                aux_qah[i][:] .= s5.ofmats[i]  
            else 
                aux_vp[i][:]  .= s5.ofmats[i]  
                aux_qah[i][:] .= s4.ofmats[i]  
            end    
        else 
            aux_vp[i][:] .= s4.ofmats[i]  
            aux_qah[i][:] .= s4.ofmats[i]  
           
        end
    end
    return  Self_consistent_data(aux_vp, s.mus, s.es, s.ns),  Self_consistent_data(aux_qah, s.mus, s.es, s.ns)
end

permute_vp_aqh_character(s::Self_consistent_data) =
    Self_consistent_data([change_valley(o) for o in s.ofmats], s.mus, s.es, s.ns)

change_valley(o) = [o[1],o[2],o[5],o[6],o[3],o[4],o[7],o[8]]

function fix_spin_convention_ofmats(s::Self_consistent_data; method = "1")
    reshuffled_ofmats = [zeros(8) for i in 1:length(s.ofmats)]
    for i in 1:length(s.ofmats) #populate the spinsin the same fashion
            if sum(s.ofmats[i][[1,2,3,4]]) > sum(s.ofmats[i][[5,6,7,8]])
                reshuffled_ofmats[i] .= s.ofmats[i]
            else
                reshuffled_ofmats[i] .= s.ofmats[i][[5,6,7,8,1,2,3,4]]
            end
    end
        Self_consistent_data(reshuffled_ofmats, s.mus, s.es, s.ns)
end

function fix_spin_convention2_ofmats(s::Self_consistent_data; method = "1")
    reshuffled_ofmats = [zeros(8) for i in 1:length(s.ofmats)]
    for i in 1:length(s.ofmats) #populate the spinsin the same fashion
            if sum(s.ofmats[i][[1,2]]) > sum(s.ofmats[i][[5,6]])
                reshuffled_ofmats[i] .= s.ofmats[i]
            else
                reshuffled_ofmats[i] .= s.ofmats[i][[5,6,3,4,1,2,7,8]]
            end
    end
        Self_consistent_data(reshuffled_ofmats, s.mus, s.es, s.ns)
end

