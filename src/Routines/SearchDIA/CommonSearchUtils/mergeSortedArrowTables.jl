#
#=
Have N of these tables. Need to combine into one sorted Arrow table without loading all tables
into memory at once. 
julia> DataFrame(Arrow.Table(readdir(second_quant_folder, join = true)[1]))
280488×11 DataFrame
    Row │ precursor_idx  prob      weight         target  irt_obs    missed_cleavage  isotopes_captured  scan_idx  ms_file_idx  peak_area   new_best_scan 
        │ UInt32         Float32   Float32        Bool    Float32    UInt8            Tuple{Int8, Int8}  UInt32    Int64        Float32     UInt32        
────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
      1 │       1468734  0.808784    201.467        true   0.659221                1  (0, 3)                49270            1     6.03541          50180
      2 │        262434  0.989585   2696.17         true   0.659221                0  (0, 3)                76753            1   121.201            76753
=#

getColNames(at::Arrow.Table) = keys(at)
getColTypes(at::Arrow.Table) = [eltype(at[col]) for col in getColNames(at)]
function getEmptyDF(at::Arrow.Table, N::Int)
    df = DataFrame()
    [df[!,Symbol(col)] = Vector{coltype}(undef, N) for (col, coltype) in zip(getColNames(at), getColTypes(at))]
    return df
end
    
function addPrecursorToHeap!(
    precursor_heap::BinaryMinHeap{Tuple{UInt32, UInt32, Int64}},
    first_sort_key::AbstractVector{UInt32},
    second_sort_key::AbstractVector{UInt32},
    table_idx::Int64,
    row_idx::Int64)
    push!(
        precursor_heap,
        (
        first_sort_key[row_idx],
        second_sort_key[row_idx],
        table_idx
        )
    )
end

function addPrecursorToHeap!(
    precursor_heap::BinaryMinHeap{Tuple{S, UInt32, Int64}},
    first_sort_key::AbstractVector{S},
    second_sort_key::AbstractVector{UInt32},
    table_idx::Int64,
    row_idx::Int64) where {S<:AbstractString}
    push!(
        precursor_heap,
        (
        first_sort_key[row_idx],
        second_sort_key[row_idx],
        table_idx
        )
    )
end

# Handle AbstractVector with Union{Missing, String} elements
function addPrecursorToHeap!(
    precursor_heap::BinaryMinHeap{Tuple{Union{Missing, String}, UInt32, Int64}},
    first_sort_key::AbstractVector{Union{Missing, String}},
    second_sort_key::AbstractVector{UInt32},
    table_idx::Int64,
    row_idx::Int64)
    push!(
        precursor_heap,
        (
        first_sort_key[row_idx],
        second_sort_key[row_idx],
        table_idx
        )
    )
end

function mergeSortedArrowTables(
    input_dir::String, 
    output_path::String,
    sort_keys::Tuple{Symbol, Symbol};
    N = 1000000
)

    function fillColumn!(
        peptide_batch_col::Vector{R},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    ) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::R
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{Union{Missing, R}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    ) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::Union{Missing, R}
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{String},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    )
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::String
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{Union{Missing, String}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    )
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::Union{String,Missing}
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{Tuple{R, R}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    ) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::Tuple{R, R}
        end
    end

    #Get all .arrow files in the input 
    input_paths = [path for path in readdir(input_dir, join=true) if endswith(path, ".arrow")]
    #Keep track of which tables have 
    tables = [Arrow.Table(path) for path in input_paths]
    table_idxs = ones(Int64, length(tables))

    peptide_batch = getEmptyDF(first(tables), N)
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, N)
    peptide_batch_names = Symbol.(names(peptide_batch))
    precursor_heap = BinaryMinHeap{Tuple{eltype(first(tables)[first(sort_keys)]), UInt32, Int64}}()
    for (i, table) in enumerate(tables)
        addPrecursorToHeap!(
            precursor_heap,
            table[first(sort_keys)],
            table[last(sort_keys)],
            i,
            1
        )
    end
    i = 1
    n_writes = 0
    while length(precursor_heap) > 0
        _, _, table_idx = pop!( precursor_heap)
        table = tables[table_idx]
        idx = table_idxs[table_idx]
        sorted_tuples[i] = (table_idx, idx)
        table_idxs[table_idx] += 1
        idx = table_idxs[table_idx]
        if (idx > length(table[1]))
            continue
        end
        addPrecursorToHeap!(
            precursor_heap,
            table[first(sort_keys)],
            table[last(sort_keys)],
            table_idx,
            idx
        )
        i += 1
        if i > N
            for col in peptide_batch_names
                fillColumn!(
                    peptide_batch[!,col],
                    col,
                    sorted_tuples,
                    tables,
                    i
                )
            end
            if iszero(n_writes)
                if isfile(output_path)
                    rm(output_path)
                end
                open(output_path, "w") do io
                    Arrow.write(io, peptide_batch; file=false)  # file=false creates stream format
                end
            else
                Arrow.append(
                    output_path,
                    peptide_batch
                )
            end
            n_writes += 1
            i = 1
        end
    end
    for col in peptide_batch_names
        fillColumn!(
            peptide_batch[!,col],
            col,
            sorted_tuples,
            tables,
            i
        )
    end
    Arrow.append(
        output_path,
        peptide_batch[range(1, max(1, i - 1)),:]
    )
    return nothing
end

# 4-key sorting support for protein group analysis
function addPrecursorToHeap!(
    precursor_heap::BinaryMinHeap{Tuple{S, Bool, UInt8, UInt32, Int64}},
    first_sort_key::AbstractVector{S},
    second_sort_key::AbstractVector{Bool},
    third_sort_key::AbstractVector{UInt8},
    fourth_sort_key::AbstractVector{UInt32},
    table_idx::Int64,
    row_idx::Int64) where {S<:AbstractString}
    push!(
        precursor_heap,
        (
        first_sort_key[row_idx],
        second_sort_key[row_idx],
        third_sort_key[row_idx],
        fourth_sort_key[row_idx],
        table_idx
        )
    )
end

# Handle Union{Missing, String} for 4-key sorting
function addPrecursorToHeap!(
    precursor_heap::BinaryMinHeap{Tuple{Union{Missing, String}, Bool, UInt8, UInt32, Int64}},
    first_sort_key::AbstractVector{Union{Missing, String}},
    second_sort_key::AbstractVector{Bool},
    third_sort_key::AbstractVector{UInt8},
    fourth_sort_key::AbstractVector{UInt32},
    table_idx::Int64,
    row_idx::Int64)
    push!(
        precursor_heap,
        (
        first_sort_key[row_idx],
        second_sort_key[row_idx],
        third_sort_key[row_idx],
        fourth_sort_key[row_idx],
        table_idx
        )
    )
end

function addPrecursorToHeap!(
    precursor_heap::BinaryMinHeap{Tuple{Float32, Int64}},
    first_sort_key::AbstractVector{Float32},
    table_idx::Int64,
    row_idx::Int64)
    push!(
        precursor_heap,
        (
        first_sort_key[row_idx],
        table_idx
        )
    )
end

function addPrecursorToHeap!(
    precursor_heap::BinaryMinHeap{Tuple{Float32, Int64}},
    first_sort_key::AbstractVector{Union{Missing, Float32}},
    table_idx::Int64,
    row_idx::Int64)
    push!(
        precursor_heap,
        (
        coalesce(first_sort_key[row_idx], 0.0f0),
        table_idx
        )
    )
end

function mergeSortedArrowTables(
    input_dir::String, 
    output_path::String,
    sort_key::Symbol;
    N = 1000000)

    function fillColumn!(
        peptide_batch_col::Vector{R},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    ) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::R
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{Union{Missing, R}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    ) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::Union{Missing, R}
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{String},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    )
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::String
        end
    end


    function fillColumn!(
        peptide_batch_col::Vector{Union{Missing, String}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    )
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::String
        end
    end


    function fillColumn!(
        peptide_batch_col::Vector{Tuple{R, R}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::Tuple{R, R}
        end
    end

    #Get all .arrow files in the input 
    input_paths = [path for path in readdir(input_dir, join=true) if endswith(path, ".arrow")]
    #Keep track of which tables have 
    tables = [Arrow.Table(path) for path in input_paths]
    table_idxs = ones(Int64, length(tables))

    peptide_batch = getEmptyDF(first(tables), N)
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, N)
    peptide_batch_names = Symbol.(names(peptide_batch))
    precursor_heap = BinaryMinHeap{Tuple{Float32, Int64}}()
    for (i, table) in enumerate(tables)
        addPrecursorToHeap!(
            precursor_heap,
            table[sort_key],
            i,
            1
        )
    end
    i = 1
    n_writes = 0
    while length(precursor_heap) > 0
        _, table_idx = pop!(precursor_heap)
        table = tables[table_idx]
        idx = table_idxs[table_idx]
        sorted_tuples[i] = (table_idx, idx)
        table_idxs[table_idx] += 1
        idx = table_idxs[table_idx]
        if (idx > length(table[1]))
            continue
        end
        addPrecursorToHeap!(
            precursor_heap,
            table[sort_key],
            table_idx,
            idx
        )
        i += 1
        if i > N
            for col in peptide_batch_names
                fillColumn!(
                    peptide_batch[!,col],
                    col,
                    sorted_tuples,
                    tables,
                    i
                )
            end

            if iszero(n_writes)
                if isfile(output_path)
                    rm(output_path)
                end
                open(output_path, "w") do io
                    Arrow.write(io, peptide_batch; file=false)  # file=false creates stream format
                end
            else
                Arrow.append(
                    output_path,
                    peptide_batch
                )
            end
            n_writes += 1
            i = 1
        end
    end
    for col in peptide_batch_names
        fillColumn!(
            peptide_batch[!,col],
            col,
            sorted_tuples,
            tables,
            i
        )
    end
    Arrow.append(
        output_path,
        peptide_batch[range(1, max(1, i - 1)),:]
    )
    return nothing
end


# 4-key sorting version for protein group analysis
function mergeSortedArrowTables(
    input_dir::String, 
    output_path::String,
    sort_keys::NTuple{4, Symbol};
    N = 1000000
)

    function fillColumn!(
        peptide_batch_col::Vector{R},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    ) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::R
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{Union{Missing, R}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    ) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::Union{Missing, R}
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{String},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    )
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::String
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{Union{Missing, String}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    )
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::Union{String,Missing}
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{Tuple{R, R}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    ) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::Tuple{R, R}
        end
    end

    #Get all .arrow files in the input 
    input_paths = [path for path in readdir(input_dir, join=true) if endswith(path, ".arrow")]
    #Keep track of which tables have 
    tables = [Arrow.Table(path) for path in input_paths]
    table_idxs = ones(Int64, length(tables))

    peptide_batch = getEmptyDF(first(tables), N)
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, N)
    peptide_batch_names = Symbol.(names(peptide_batch))
    
    # Create heap with appropriate type based on the first table's column types
    first_table = first(tables)
    T1 = eltype(first_table[sort_keys[1]])  # String or Union{Missing, String}
    T2 = eltype(first_table[sort_keys[2]])  # Bool
    T3 = eltype(first_table[sort_keys[3]])  # UInt8
    T4 = eltype(first_table[sort_keys[4]])  # UInt32
    
    precursor_heap = BinaryMinHeap{Tuple{T1, T2, T3, T4, Int64}}()
    
    for (i, table) in enumerate(tables)
        addPrecursorToHeap!(
            precursor_heap,
            table[sort_keys[1]],
            table[sort_keys[2]],
            table[sort_keys[3]],
            table[sort_keys[4]],
            i,
            1
        )
    end
    
    i = 1
    n_writes = 0
    while length(precursor_heap) > 0
        _, _, _, _, table_idx = pop!(precursor_heap)
        table = tables[table_idx]
        idx = table_idxs[table_idx]
        sorted_tuples[i] = (table_idx, idx)
        table_idxs[table_idx] += 1
        idx = table_idxs[table_idx]
        if (idx > length(table[1]))
            continue
        end
        addPrecursorToHeap!(
            precursor_heap,
            table[sort_keys[1]],
            table[sort_keys[2]],
            table[sort_keys[3]],
            table[sort_keys[4]],
            table_idx,
            idx
        )
        i += 1
        if i > N
            for col in peptide_batch_names
                fillColumn!(
                    peptide_batch[!,col],
                    col,
                    sorted_tuples,
                    tables,
                    i
                )
            end
            if iszero(n_writes)
                if isfile(output_path)
                    rm(output_path)
                end
                open(output_path, "w") do io
                    Arrow.write(io, peptide_batch; file=false)  # file=false creates stream format
                end
            else
                Arrow.append(
                    output_path,
                    peptide_batch
                )
            end
            n_writes += 1
            i = 1
        end
    end
    for col in peptide_batch_names
        fillColumn!(
            peptide_batch[!,col],
            col,
            sorted_tuples,
            tables,
            i
        )
    end
    Arrow.append(
        output_path,
        peptide_batch[range(1, max(1, i - 1)),:]
    )
    return nothing
end
