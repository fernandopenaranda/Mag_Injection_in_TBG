function self_consistent_densitymap(p, ωlist, nlist, η = 0.5, evals = 10000)
    start_time = time_ns()
    mat1 = zeros(Float64, length(ωlist), length(nlist))
    count = 1;
    for n in nlist
        nf, μ, n = selfconsistency(p, n; tol = 1e-3, α = 0.4)
        _, vals = hf_jdos(ParamsHF(p, n = n, nf = nf, μ = μ), ωlist, η = η, evals = 10000) 
        mat1[:,count] .= vals
        count +=1
    end
    save_data(nlist, ωlist, mat1)
    elapsed_time = (time_ns() - start_time) / 10^9
    println("Elapsed time: ", elapsed_time, " seconds")
    return nlist, ωlist, mat1
end


function self_consistent_linear_map(a, b, p, ωlist, nlist; kws...)
    start_time = time_ns()
    mat1 = zeros(Float64, length(ωlist), length(nlist))
    mat2 = zeros(Float64, length(ωlist), length(nlist))
    count = 1;
    for n in nlist
        nf, μ, n = selfconsistency(p, n; tol = 1e-3, α = 0.4)
        _ , posval= σab_inter_linear(a, b, ParamsHF(p, ν = 1, nf = nf, μ = μ, n= n), ωlist; kws...)
        _ , negval= σab_inter_linear(a, b, ParamsHF(p, ν = -1, nf = nf, μ = μ, n= n), ωlist; kws...)
        mat1[:,count] .= posval
        mat2[:,count] .= negval
        count +=1
    end
    save_data(nlist, ωlist, mat1, mat2)
    elapsed_time = (time_ns() - start_time) / 10^9
    println("Elapsed time: ", elapsed_time, " seconds")
    return nlist, ωlist, mat1, mat2
end

function self_consistent_shift_map(part, a, b, c, p, ωlist, nlist; kws...)
    start_time = time_ns()
    mat1 = zeros(Float64, length(ωlist), length(nlist))
    mat2 = zeros(Float64, length(ωlist), length(nlist))
    count = 1;
    for n in nlist
        # print(" ", i/length(nlist))
        nf, μ, n = selfconsistency(p, n; tol = 1e-3, α = 0.4)
        _ , posval, negval = shift_current(a, b, c, ParamsHF(p, ν = 1, nf = nf, μ = μ, n= n), ωlist, part; kws...)
        mat1[:,count] .= posval
        mat2[:,count] .= negval
        count +=1
    end
    save_data(nlist, ωlist, mat1, mat2)
    elapsed_time = (time_ns() - start_time) / 10^9
    println("Elapsed time: ", elapsed_time, " seconds")
    return nlist, ωlist, mat1, mat2
end


function self_consistent_shift_map_onevalley(part, a, b, c, p, ωlist, nlist; kws...)
    start_time = time_ns()
    mat1 = zeros(Float64, length(ωlist), length(nlist))
    count = 1;
    for n in nlist
        nf, μ, n = selfconsistency(p, n; tol = 1e-3, α = 0.4)
        _ , posval, negval = shift_current(a, b, c, ParamsHF(p, ν = 1, nf = nf, μ = μ, n= n), ωlist, part; kws...)
        mat1[:,count] .= posval
        mat2[:,count] .= negval
        count +=1
    end
    save_data(nlist, ωlist, mat1, mat2)
    elapsed_time = (time_ns() - start_time) / 10^9
    println("Elapsed time: ", elapsed_time, " seconds")
    return nlist, ωlist, mat1, mat2
end


injection_map(part, a, b, c, p, ωlist, μlist; kws...) = 
    response_map(part, a, b, c, p, ωlist, μlist, injection_current; kws...)

function response_map(part, a, b, c, p, ωs, μlist, func; kws...) 
    start_time = time_ns()
    mat1 = zeros(Float64, length(ωs), length(μlist))
    mat2 = zeros(Float64, length(ωs), length(μlist))
    for (i,μ) in enumerate(μlist)
        print(" ", i/length(μlist))
        _ , posval, negval = injection_current(a, b, c, ParamsHF(p, μ = μ) , ωs, part; kws...)
        mat1[:,i] .= posval
        mat2[:,i] .= negval
    end
    save_data(μlist, ωs, mat1, mat2)
    elapsed_time = (time_ns() - start_time) / 10^9
    println("Elapsed time: ", elapsed_time, " seconds")
    return μlist, ωs, mat1, mat2
end

#= 

                SAVE DATA

=#
function save_data(array1, array2, matrix1, matrix2)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
    df_array1 = DataFrame(omegas = vec(array1))
    df_array2 = DataFrame(omegas = vec(array2))

    # Define file names
    str = "~/"
    file_array1 = str * "array1_$timestamp.csv"
    file_array2 = str * "array2_$timestamp.csv"
    file_matrix1 = str * "posvalley_$timestamp.csv"
    file_matrix2 = str * "negvalley_$timestamp.csv"
    CSV.write(file_array1, df_array1)
    CSV.write(file_array2, df_array2)
    CSV.write(file_matrix1,  Tables.table(matrix1), writeheader=false)
    CSV.write(file_matrix2,  Tables.table(matrix2), writeheader=false)

    # Write to CSV files
    println("Files saved:")
    println("  - ", file_array1)
    println("  - ", file_array2)
    println("  - ", file_matrix1)
    println("  - ", file_matrix2)
end

function save_data(array1, array2, matrix1)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
    df_array1 = DataFrame(omegas = vec(array1))
    df_array2 = DataFrame(omegas = vec(array2))

    # Define file names
    str = "~/"
    file_array1 = str * "array1_$timestamp.csv"
    file_array2 = str * "array2_$timestamp.csv"
    file_matrix1 = str * "densities_$timestamp.csv"
    CSV.write(file_array1, df_array1)
    CSV.write(file_array2, df_array2)
    CSV.write(file_matrix1,  Tables.table(matrix1), writeheader=false)

    # Write to CSV files
    println("Files saved:")
    println("  - ", file_array1)
    println("  - ", file_array2)
    println("  - ", file_matrix1)
end