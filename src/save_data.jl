
function save_arrays(filename::AbstractString, array1)
    data = DataFrame(array1=array1)
    CSV.write(filename, data)
end

function save_arrays(filename::AbstractString, array1, array2)
    data = DataFrame(array1=array1, array2=array2)
    CSV.write(filename, data)
end


function store_struct_data(filename::AbstractString, data)
    field_names = fieldnames(typeof(data))
    field_values = [getfield(data, field) for field in field_names] 
    data_dict = Dict([(field_names[i], field_values[i]) for i in 1:length(field_names)])
    df = DataFrame(data_dict)
    CSV.write(filename, df)
end

function read_struct_data(filename)
    df = CSV.File(filename) |> DataFrame
    p = paramsHF(1.05, 1,  1)
    col_names = names(df)
    for i in 1:length(col_names)
        dc = Dict(Symbol(col_names[i]) => df[1, i])
        p = ParamsHF(p, dc)
    end
    return p
end

function import_data(file_array1, file_array2, file_matrix1, file_matrix2)

    df_array1 = CSV.read(file_array1, DataFrame)
    df_array2 = CSV.read(file_array2, DataFrame)
    df_matrix1 = CSV.read(file_matrix1, DataFrame)
    df_matrix2 = CSV.read(file_matrix2, DataFrame)

    # Convert DataFrames to arrays
    array1 = Matrix(df_array1)
    array2 = Matrix(df_array2)

    # Convert DataFrames to matrices
    matrix1 = Matrix(df_matrix1)
    matrix2 = Matrix(df_matrix2)

    return array1, array2, matrix1, matrix2
end

import_array(file_arrays) =  Matrix(CSV.read(file_arrays, DataFrame))
    
    
    