function addRowToHeap!(
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

function mergeSortedProteinGroups(
    input_dir::String, 
    output_path::String,
    sort_key::Symbol;
    N = 1000000 ) #N -> batch size 

    function fillColumn!(
        pg_batch_col::Vector{R},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n) where {R<:Real}
        for i in range(1, max(min(length(pg_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            pg_batch_col[i] = tables[table_idx][col][idx]::R
        end
    end

    #Get all .arrow files in the input 
    input_paths = [path for path in readdir(input_dir, join=true) if endswith(path, ".arrow")]
    #Keep track of which tables have 
    tables = [Arrow.Table(path) for path in input_paths]
    table_idxs = ones(Int64, length(tables))

    #pre-allocate section of merged table 
    pg_batch = getEmptyDF(first(tables), N) #See /src/utils/mergePsmTables.jl for definition
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, N)
    #Names of columns for merge table 
    pg_batch_names = Symbol.(names(pg_batch))
    #Keeps track of the talbe with the highest ranked row
    pg_heap = BinaryMinHeap{Tuple{Float32, Int64}}()
    for (i, table) in enumerate(tables)
        addRowToHeap!(
            pg_heap,
            table[sort_key],
            i,
            1
        )
    end
    i = 1
    while length(pg_heap) > 0
        _, table_idx = pop!(pg_heap)
        table = tables[table_idx]
        idx = table_idxs[table_idx]
        sorted_tuples[i] = (table_idx, idx)
        table_idxs[table_idx] += 1
        idx = table_idxs[table_idx]
        if (idx > length(table[1]))
            continue
        end
        addRowToHeap!(
            pg_heap,
            table[sort_key],
            table_idx,
            idx
        )
        i += 1
        if i > N
            for col in pg_batch_names
                fillColumn!(
                    pg_batch[!,col],
                    col,
                    sorted_tuples,
                    tables,
                    i
                )
            end
            Arrow.append(
                output_path,
                pg_batch
            )
            i = 1
        end
    end
    for col in pg_batch_names
        fillColumn!(
            pg_batch[!,col],
            col,
            sorted_tuples,
            tables,
            i
        )
    end
    Arrow.append(
        output_path,
        pg_batch[range(1, max(1, i - 1)),:]
    )
    return nothing
end

function getPgScoreThreshold(
    scores::AbstractVector{Float32},
    targets::AbstractVector{Bool},
    q_value_threshold::Float32
)
    # Second pass: Find the probability threshold
    targets_above = 0
    decoys_above = 0

    for (score, target) in zip(scores, targets)

        targets_above += target
        decoys_above += (1 - target)
        
        current_q_value = decoys_above / (targets_above + decoys_above)
        
        if current_q_value > q_value_threshold
           return score#return prob  # This is the probability threshold we're looking for
        end
    end

    return scores[end]

end

function getProteinGroups(
    passing_psms_folder::String,
    protein_groups_folder::String,
    temp_folder::String,
    accession_numbers::AbstractVector{String},
    accession_number_to_id::Dict{String, UInt32},
    precursor_sequence::AbstractVector{String};
    min_peptides = 2,
    protein_q_val_threshold::Float32 = 0.01f0)

    function getProteinGroupsDict(
        psms_path::String,
        accession_numbers::AbstractVector{String},
        accession_number_to_id::Dict{String, UInt32},
        precursor_sequence::AbstractVector{String};
        min_peptides::Int64 = 2)

        psms_table = Arrow.Table(psms_path)

        protein_groups = Dictionary{@NamedTuple{protein_idx::UInt32, target::Bool},
        @NamedTuple{
            max_pg_score::Float32, 
            peptides::Set{String}}
        }()

        for i in range(1, length(psms_table[1]))
            precursor_idx = psms_table[:precursor_idx][i]
            sequence = precursor_sequence[precursor_idx]
            score = psms_table[:prob][i]
            protein_idx = accession_number_to_id[accession_numbers[precursor_idx]]
            keyname = (protein_idx = protein_idx, target = psms_table[:target][i])
            if haskey(protein_groups, keyname)
                max_pg_score, peptides = protein_groups[keyname]
                if score > max_pg_score
                    max_pg_score = score
                end
                push!(peptides, sequence)
                protein_groups[keyname] = (max_pg_score = max_pg_score, peptides = peptides)
            else
                sequences = Set{String}((sequence,))
                insert!(protein_groups,
                keyname,
                        (max_pg_score = score,
                        peptides = sequences)
                )
            end
        end
        filter!(x->length(x[:peptides])>=min_peptides, protein_groups)
        #modify the table
        psms_table = DataFrame(Tables.columntable(psms_table))
        psms_table[!,:max_pg_score] = Vector{Union{Missing, Float32}}(undef, size(psms_table, 1))
        for i in range(1, size(psms_table, 1))
            precursor_idx = psms_table[i,:precursor_idx]
            protein_idx = accession_number_to_id[accession_numbers[precursor_idx]]
            key = (protein_idx = protein_idx, target = psms_table[i,:target])
            if haskey(protein_groups, key)
                psms_table[i,:max_pg_score] = protein_groups[key][:max_pg_score]
            else
                psms_table[i,:max_pg_score] = missing
            end
        end
        Arrow.write(
            psms_path,
            psms_table
        )
        return protein_groups
    end

    function writeProteinGroups(
                                    protein_groups::Dictionary{
                                    @NamedTuple{protein_idx::UInt32, target::Bool},
                                    @NamedTuple{max_pg_score::Float32,  peptides::Set{String}}
                                    },
                                    protein_groups_path::String)
        # Extract keys and values
        keys_array = keys(protein_groups)
        values_array = values(protein_groups)

        # Create vectors for each column
        protein_idx = [k[:protein_idx] for k in keys_array]
        target = [k[:target] for k in keys_array]
        max_pg_score = [v[:max_pg_score] for v in values_array]
        #peptides = [join(v[:peptides], ";") for v in values_array]  # Convert Set to String

        # Create DataFrame
        df = DataFrame((
            protein_idx = protein_idx,
            target = target,
            max_pg_score = max_pg_score,
        )
        )
        sort!(df, :max_pg_score, rev = true)
        # Convert DataFrame to Arrow.Table
        Arrow.write(protein_groups_path, df)
        return size(df, 1)
    end


    pg_count = 0
    for file_path in readdir(passing_psms_folder, join=true)
        _, extention = splitext(file_path)
        if extention != ".arrow"
            continue
        end
        protein_groups_path = joinpath(protein_groups_folder, basename(file_path))
        protein_groups =  getProteinGroupsDict(
            file_path,
            accession_numbers,
            accession_number_to_id,
            precursor_sequence;
            min_peptides = min_peptides
        )
        pg_count += writeProteinGroups(
            protein_groups,
            protein_groups_path
        )
    end
    #=
    max_pg_scores = Vector{@NamedTuple{score::Float32, target::Bool}}(undef, pg_count)
    n = 1
    protein_groups_paths = readdir(protein_groups_folder, join=true)
    for file_path in protein_groups_paths
        pg_table = Arrow.Table(file_path)
        for i in range(1, length(pg_table[:max_pg_score]))
            max_pg_scores[n] = (score = pg_table[:max_pg_score][i], target = pg_table[:target][i])
            n += 1
        end
    end
    sort!(max_pg_scores, alg=QuickSort, by=x->x[:score],rev = true)
    # Second pass: Find the probability threshold
    pg_score_threshold = find_score_threshold(max_pg_scores, protein_q_val_threshold)
    =#
    sorted_pg_scores_path = joinpath(temp_folder, "sorted_pg_scores.arrow")
    mergeSortedProteinGroups(
        protein_groups_folder,
        sorted_pg_scores_path,
        :max_pg_score,
        N = 1000000
    )
    sorted_pg_scores = Arrow.Table(sorted_pg_scores_path)
    return getPgScoreThreshold(
            sorted_pg_scores[:max_pg_score],
            sorted_pg_scores[:target],
            protein_q_val_threshold
        )
end
#test_table = Arrow.Table("/Users/n.t.wamsley/Desktop/testresults/RESULTS/Search/temp/passing_psms/E10H50Y40_02.arrow")

