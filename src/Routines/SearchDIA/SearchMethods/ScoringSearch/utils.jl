#==========================================================
Trace Selection  
==========================================================#

"""
    get_best_traces(second_pass_psms_paths::Vector{String}, min_prob::Float32=0.75f0) 
    -> Set{@NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}}

Identify best scoring isotope traces for each precursor.

# Process
1. Accumulates scores across files
2. Selects highest scoring trace per precursor
3. Returns set of best precursor-isotope combinations
"""
function get_best_traces(
    second_pass_psms_paths::Vector{String},
    min_prob::Float32 = 0.75f0
)
    psms_trace_scores = Dictionary{
            @NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}, Float32}()

    for file_path in second_pass_psms_paths
        if splitext(file_path)[end] != ".arrow"
            continue
        end
        row_score = zero(Float32)
        psms_table = Arrow.Table(file_path)
        for i in range(1, length(psms_table[1]))
            psms_key = (precursor_idx = psms_table[:precursor_idx][i],  isotopes_captured = psms_table[:isotopes_captured][i])

            if (psms_table[:prob][i]>min_prob)#&(psms_table[:fraction_transmitted]>0.25)
                row_score = psms_table[:weight][i]
            else
                row_score = zero(Float32)
            end

            row_score = log2(psms_table[:prob][i])
            if haskey(psms_trace_scores, psms_key)
                psms_trace_scores[psms_key] = psms_trace_scores[psms_key] + row_score
            else
                insert!(
                    psms_trace_scores,
                    psms_key,
                    row_score
                )
            end
        end
    end

    psms_trace_df = DataFrame(
    (precursor_idx = [key[:precursor_idx] for key in keys(psms_trace_scores)],
    isotopes_captured = [key[:isotopes_captured] for key in keys(psms_trace_scores)],
    score = [val for val in values(psms_trace_scores)])
    );
    psms_trace_df[!,:best_trace] .= false;
    gpsms = groupby(psms_trace_df,:precursor_idx)
    for (precursor_idx, psms) in pairs(gpsms)
        psms[argmax(psms[!,:score]),:best_trace] = true
    end
    filter!(x->x.best_trace, psms_trace_df);
    traces_passing = Set([(precursor_idx = x.precursor_idx, isotopes_captured = x.isotopes_captured) for x in eachrow(psms_trace_df)]);
    return traces_passing
end

"""
    sort_and_filter_quant_tables(second_pass_psms_paths::Vector{String},
                                merged_quant_path::String,
                                best_traces::Set{@NamedTuple})

Filter PSM tables to retain only best traces and sort by probability.

Modifies files in place to optimize storage.
"""
function sort_and_filter_quant_tables(
    second_pass_psms_paths::Vector{String},
    merged_quant_path::String,
    best_traces::Set{@NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}}
)

    #Remove if present 
    if isfile(merged_quant_path)
        rm(merged_quant_path)
    end
    #file_paths = [fpath for fpath in readdir(quant_psms_folder,join=true) if endswith(fpath,".arrow")]
    #Sort and filter each psm table 
    for fpath in second_pass_psms_paths
        psms_table = DataFrame(Tables.columntable(Arrow.Table(fpath)))

        #Indicator variable of whether each psm is from the best trace 
        psms_table[!,:best_trace] = zeros(Bool, size(psms_table, 1))
        for i in range(1, size(psms_table, 1))
            key = (precursor_idx = psms_table[i, :precursor_idx], isotopes_captured = psms_table[i, :isotopes_captured])
            if key ∈ best_traces
                psms_table[i,:best_trace]=true
            end
        end
        #Filter out unused traces 
        filter!(x->x.best_trace,psms_table)
        #Sort in descending order of probability
        sort!(psms_table, :prob, rev = true, alg=QuickSort)
        #write back
        writeArrow(fpath, psms_table)
    end
    return nothing
end

#==========================================================
PSM score merging and processing 
==========================================================#
"""
    merge_sorted_psms_scores(input_paths::Vector{String}, output_path::String; N=10000000)

Merge sorted PSM scores from multiple files.

Uses heap-based merging for memory efficiency with batched processing.
"""
function merge_sorted_psms_scores(
                    input_paths::Vector{String}, 
                    output_path::String,
                    ; N = 10000000
)
    
    function fillColumn!(
        peptide_batch_col::Vector{R},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::R
        end
    end

    function addPrecursorToHeap!(
        psms_heap::BinaryMaxHeap{Tuple{Float32, Int64}},
        sort_key::AbstractVector{Float32},
        table_idx::Int64,
        row_idx::Int64)
        push!(
            psms_heap,
            (
            sort_key[row_idx],
            table_idx
            )
        )
    end

    ##println("psms input_paths $input_paths")
    #input_paths = [path for path in readdir(input_dir,join=true) if endswith(path,".arrow")]
    tables = [Arrow.Table(path) for path in input_paths]
    table_idxs = ones(Int64, length(tables))
    psms_batch = DataFrame()
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, N)
    psms_heap = BinaryMaxHeap{Tuple{Float32, Int64}}()
    psms_batch[!,:prob] = zeros(Float32, N)
    psms_batch[!,:target] = zeros(Bool, N)

    
    #Load psms_heap
    for (i, table) in enumerate(tables)
        ##println("i $i")
        ##println("size(table) ", size(DataFrame(table)))
        addPrecursorToHeap!(
            psms_heap,
            table[:prob],
            i,
            1
        )
    end
    i = 1
    n_writes = 0
    while length(psms_heap)>0
        _, table_idx = pop!(psms_heap)
        table = tables[table_idx]
        idx = table_idxs[table_idx]
        sorted_tuples[i] = (table_idx, idx)
        table_idxs[table_idx] += 1
        idx = table_idxs[table_idx]
        if (idx > length(table[1]))
            continue
        end
        addPrecursorToHeap!(
            psms_heap,
            table[:prob],
            table_idx,
            idx
        )
        i += 1
        if i > N
            for col in [:prob,:target]
                fillColumn!(
                    psms_batch[!,col],
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
                    Arrow.write(io, psms_batch; file=false)  # file=false creates stream format
                end
            else
                Arrow.append(
                    output_path,
                    psms_batch
                )
            end
            n_writes += 1
            i = 1
        end
    end
    for col in [:prob,:target]
        fillColumn!(
            psms_batch[!,col],
            col,
            sorted_tuples,
            tables,
            i
        )
    end
    Arrow.append(
        output_path,
        psms_batch[range(1, max(1, i - 1)),:]
    )
    #println("output_path $output_path")
    return nothing
end


"""
    get_pep_spline(merged_psms_path::String, score_col::Symbol;
                   min_pep_points_per_bin=5000, n_spline_bins=20) -> UniformSpline

Create posterior error probability spline from merged scores.

Returns spline for PEP calculation based on target/decoy distributions.
"""
function get_pep_spline(
                            merged_psms_path::String,
                            score_col::Symbol;
                            min_pep_points_per_bin = 5000,
                            n_spline_bins = 20
)

    psms_scores = Arrow.Table(merged_psms_path)
    Q = length(psms_scores[1])
    M = ceil(Int, Q / min_pep_points_per_bin)
    bin_target_fraction, bin_mean_prob = Vector{Float32}(undef, M), Vector{Float32}(undef, M)
    bin_size = 0
    bin_idx = 0
    mean_prob, targets = 0.0f0, 0
    for i in range(1, Q)
        bin_size += 1
        targets += psms_scores[:target][i]
        mean_prob += psms_scores[score_col][i]
        if bin_size == min_pep_points_per_bin
            bin_idx += 1
            bin_target_fraction[bin_idx] = targets/bin_size
            bin_mean_prob[bin_idx] = mean_prob/bin_size
            bin_size, targets, mean_prob = zero(Int64), zero(Int64), zero(Float32)
        end
    end
    bin_target_fraction[end] = targets/max(bin_size, 1)
    bin_mean_prob[end] = mean_prob/max(bin_size, 1)
    try 
        if length(bin_target_fraction)<20
            @warn "Less than 20 bins to estimate PEP. PEP results suspect..."
        end
        return UniformSpline(bin_target_fraction, bin_mean_prob, 3, 3)
    catch
        @warn "Failed to estimate PEP spline"
        return UniformSpline(SVector{4, Float32}([0, 0, 0, 0]), 3, 0.0f0, 1.0f0, 100.0f0)
    end
end

"""
    get_qvalue_spline(merged_psms_path::String, score_col::Symbol;
                      min_pep_points_per_bin=1000) -> Interpolation

Create q-value interpolation function from merged scores.

Returns interpolation for mapping scores to q-values.
"""
function get_qvalue_spline(
                            merged_psms_path::String,
                            score_col::Symbol;
                            min_pep_points_per_bin = 1000
)

    psms_scores = Arrow.Table(merged_psms_path)
    Q = length(psms_scores[1])
    M = ceil(Int, Q / min_pep_points_per_bin)
    bin_qval, bin_mean_prob = Vector{Float32}(undef, M), Vector{Float32}(undef, M)
    bin_size = 0
    bin_idx = 0
    mean_prob, targets, decoys = 0.0f0, 0, 0
    targets, decoys = 0, 0
    for i in range(1, Q)
        targets += psms_scores[:target][i]
        decoys += (1 - psms_scores[:target][i])
    end

    min_q_val = typemax(Float32)
    for i in reverse(range(1, Q))
        bin_size += 1
        targets -= psms_scores[:target][i]
        decoys -= (1 - psms_scores[:target][i])
        mean_prob += psms_scores[score_col][i]
        if bin_size == min_pep_points_per_bin
            bin_idx += 1
            qval = decoys/(targets) #not decoys/(targets + decoys)
            if qval > min_q_val
                bin_qval[bin_idx] = min_q_val
            else
                min_q_val = qval
                bin_qval[bin_idx] = qval
            end
            bin_mean_prob[bin_idx] = mean_prob/bin_size
            bin_size, mean_prob = zero(Int64), zero(Float32)
        end
    end
    bin_qval[end] = targets/bin_size
    bin_mean_prob[end] = mean_prob/bin_size
    prepend!(bin_qval, 1.0f0)
    prepend!(bin_mean_prob, 0.0f0)
    bin_qval = bin_qval[isnan.(bin_mean_prob).==false]
    bin_mean_prob = bin_mean_prob[isnan.(bin_mean_prob).==false]
    #return bin_qval, bin_mean_prob
    return linear_interpolation(
        Interpolations.deduplicate_knots!(bin_mean_prob, move_knots=true),
        bin_qval, extrapolation_bc=Line())
end

"""
    get_psms_passing_qval(passing_psms_paths::Vector{String}, passing_psms_folder::String,
                         second_pass_psms_paths::Vector{String}, pep_spline::UniformSpline,
                         qval_interp::Interpolations.Extrapolation, q_val_threshold::Float32)

Filter PSMs that pass a given q-value threshold and calculate error probabilities.

# Arguments
- `passing_psms_paths`: Vector to store paths of filtered PSM files
- `passing_psms_folder`: Folder to store passing PSM files
- `second_pass_psms_paths`: Paths to second pass PSM files
- `pep_spline`: Spline for posterior error probability calculation
- `qval_interp`: Interpolation function for q-value calculation
- `q_val_threshold`: Q-value threshold for filtering

# Process
1. For each file:
   - Calculates q-values and PEP for each PSM
   - Filters PSMs below q-value threshold
   - Selects relevant columns for output
   - Saves filtered PSMs to new file
   - Updates passing_psms_paths with new file locations

# Selected Columns
- precursor_idx, prob, qval, pep
- weight, target, irt_obs
- missed_cleavage, isotopes_captured
- scan_idx, ms_file_idx
"""
function get_psms_passing_qval(
                            passing_psms_paths::Vector{String},
                            passing_psms_folder::String,
                            second_pass_psms_paths::Vector{String}, 
                            pep_spline::UniformSpline,
                            qval_interp::Interpolations.Extrapolation,
                            q_val_threshold::Float32,
)
    
    #file_paths = [fpath for fpath in readdir(quant_psms_folder, join=true) if endswith(fpath,".arrow")]
    for (ms_file_idx, file_path) in enumerate(second_pass_psms_paths)
        # Read the Arrow table
        passing_psms = DataFrame(Tables.columntable(Arrow.Table(file_path)))
        passing_psms[!,:qval] = qval_interp.(passing_psms[!,:prob])
        passing_psms[!,:pep] = pep_spline.(passing_psms[!,:prob])
        # Sample the rows and convert to DataFrame
        select!(passing_psms,
        [
            :precursor_idx,
            :prob,
            :qval,
            :pep,
            :weight,
            :target,
            :irt_obs,
            :missed_cleavage,
            :isotopes_captured,
            :scan_idx,
            :ms_file_idx])
        ###println("size(passing_psms) ", size(passing_psms))
        ###println("q_val_threshold $q_val_threshold")
        filter!(x->x.qval<=q_val_threshold, passing_psms)
        ##println("size(passing_psms) ", size(passing_psms))
        # Append to the result DataFrame
        ##println("joinpath(passing_psms_folder, basename(file_path)) ", joinpath(passing_psms_folder, basename(file_path)))
        Arrow.write(
            joinpath(passing_psms_folder, basename(file_path)),
            passing_psms
        )
        passing_psms_paths[ms_file_idx] = joinpath(passing_psms_folder, basename(file_path))
    end

    return
end


#==========================================================
Protein group analysis
==========================================================#
"""
    get_protein_groups(passing_psms_paths::Vector{String}, passing_pg_paths::Vector{String},
                      protein_groups_folder::String, temp_folder::String,
                      precursors::BasicLibraryPrecursors; min_peptides=2,
                      protein_q_val_threshold::Float32=0.01f0) -> String

Create and score protein groups from passing PSMs.

Returns path to sorted protein group scores.
"""
function get_protein_groups(
    passing_psms_paths::Vector{String},
    passing_pg_paths::Vector{String},
    protein_groups_folder::String,
    temp_folder::String,
    precursors::BasicLibraryPrecursors;
    min_peptides = 2,
    protein_q_val_threshold::Float32 = 0.01f0
)

    function getProteinGroupsDict(
        protein_inference_dict::Dictionary{String, Tuple{String, Bool}},
        psm_precursor_idx::AbstractVector{UInt32},
        psm_score::AbstractVector{Float32},
        psm_is_target::AbstractVector{Bool},
        precursors::BasicLibraryPrecursors;
        min_peptides::Int64 = 2)

        accession_numbers = getAccessionNumbers(precursors)
        precursor_sequence = getSequence(precursors)
        protein_groups = Dictionary{@NamedTuple{protein_name::String, target::Bool},
        @NamedTuple{
            max_pg_score::Float32, 
            peptides::Set{String}}
        }()

        for i in range(1, length(psm_precursor_idx))
            precursor_idx = psm_precursor_idx[i]
            sequence = precursor_sequence[precursor_idx]
            #Exclude peptide 
            if last(protein_inference_dict[sequence])==false
                continue
            end
            score = psm_score[i]
            protein_name = first(protein_inference_dict[sequence])
            keyname = (protein_name = protein_name, target = psm_is_target[i])
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
        #psms_table = DataFrame(Tables.columntable(psms_table))
        max_pg_score = Vector{Union{Missing, Float32}}(undef, length(psm_precursor_idx))
        for i in range(1, length(psm_precursor_idx))
            precursor_idx = psm_precursor_idx[i]
            sequence = precursor_sequence[precursor_idx]
            protein_name = first(protein_inference_dict[sequence])
            #protein_idx = accession_number_to_id[accession_numbers[precursor_idx]]
            key = (protein_name = protein_name, target = psm_is_target[i])
            if haskey(protein_groups, key)
                max_pg_score[i] = protein_groups[key][:max_pg_score]
            else
                max_pg_score[i] = missing
            end
        end
        #Arrow.write(
        #    psms_path,
        #    psms_table
        #)
        return max_pg_score, protein_groups
    end

    function writeProteinGroups(
                                    protein_groups::Dictionary{
                                    @NamedTuple{protein_name::String, target::Bool},
                                    @NamedTuple{max_pg_score::Float32,  peptides::Set{String}}
                                    },
                                    protein_groups_path::String)
        # Extract keys and values
        keys_array = keys(protein_groups)
        values_array = values(protein_groups)

        # Create vectors for each column
        protein_name = [k[:protein_name] for k in keys_array]
        target = [k[:target] for k in keys_array]
        max_pg_score = [v[:max_pg_score] for v in values_array]
        #peptides = [join(v[:peptides], ";") for v in values_array]  # Convert Set to String

        # Create DataFrame
        df = DataFrame((
            protein_name = protein_name,
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
    #Concatenate psms 
    ##########
    #Protein inference
    ##########
    passing_psms = Arrow.Table(passing_psms_paths)
    protein_peptide_rows = Set{Tuple{String, String}}()
    passing_precursor_idx = passing_psms[:precursor_idx]
    accession_numbers = getAccessionNumbers(precursors)
    sequences = getSequence(precursors)
    for pid in passing_precursor_idx
        push!(
            protein_peptide_rows, 
            (
                sequences[pid],
                accession_numbers[pid]
            )
        )
    end
    protein_peptide_rows = collect(protein_peptide_rows)
    peptides = [first(row) for row in protein_peptide_rows]
    proteins = [last(row) for row in protein_peptide_rows]
    @time protein_inference_dict = infer_proteins(
        proteins,
        peptides
        )
    #Returns Dictionary{String, Tuple{String, Bool}}
    #Key is a peptide base sequence (no mods or charge state)
    #Value is 1) protein group name 2) Whether to exclude the peptide from protien quant and exclude
    #for purposes of protein scoring as well .

    for (ms_file_idx, file_path) in enumerate(passing_psms_paths)#readdir(passing_psms_folder, join=true)
        _, extention = splitext(file_path)
        if extention != ".arrow"
            continue
        end
        protein_groups_path = joinpath(protein_groups_folder, basename(file_path))
        passing_pg_paths[ms_file_idx] = protein_groups_path
        psms_table = Arrow.Table(file_path)
        max_pg_score, protein_groups = getProteinGroupsDict(
            protein_inference_dict,
            psms_table[:precursor_idx],
            psms_table[:prob],
            psms_table[:target],
            precursors;
            min_peptides = min_peptides
        )
        psms_table = DataFrame(Tables.columntable(psms_table))
        psms_table[!,:max_pg_score] = max_pg_score
        writeArrow(file_path, psms_table)
        pg_count += writeProteinGroups(
            protein_groups,
            protein_groups_path
        )
    end

    sorted_pg_scores_path = joinpath(temp_folder, "sorted_pg_scores.arrow")
    if isfile(sorted_pg_scores_path)
        rm(sorted_pg_scores_path)
    end
    merge_sorted_protein_groups(
        protein_groups_folder,
        sorted_pg_scores_path,
        :max_pg_score,
        N = 1000000
    )
    return sorted_pg_scores_path, protein_inference_dict
end

function add_protein_inferrence_col(
    passing_psms_paths::Vector{String},
    protein_inference_dict::Dictionary{String, Tuple{String, Bool}},
    precursor_sequences::AbstractVector{S}
) where {S<:AbstractString}

    for (ms_file_idx, file_path) in enumerate(passing_psms_paths)
        passing_psms = DataFrame(Tables.columntable(Arrow.Table(file_path)))
        precursor_idx = passing_psms[!,:precursor_idx]
        inferred_protein_group = Vector{String}(undef, length(precursor_idx))
        use_for_protein_quant = zeros(Bool, length(precursor_idx))
        for (i, pid) in enumerate(precursor_idx)
            protein_group, use_for_inference = protein_inference_dict[precursor_sequences[pid]]
            inferred_protein_group[i] = protein_group
            use_for_protein_quant[i] = use_for_inference
        end
        passing_psms[!,:use_for_protein_quant] = use_for_protein_quant
        passing_psms[!,:inferred_protein_group] = inferred_protein_group
        writeArrow(
            file_path,
            passing_psms            
        )
    end
end

"""
    merge_sorted_protein_groups(input_dir::String, output_path::String,
                              sort_key::Symbol; N=1000000)

Merge sorted protein group scores from multiple files.

Uses heap-based merging for memory efficiency.
"""
function merge_sorted_protein_groups(
    input_dir::String, 
    output_path::String,
    sort_key::Symbol;
    N = 1000000 
) #N -> batch size 

    function addRowToHeap!(
        precursor_heap::BinaryMaxHeap{Tuple{Float32, Int64}},
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

    function fillColumn!(
        pg_batch_col::Vector{S},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n) where {S<:AbstractString}
        for i in range(1, max(min(length(pg_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            pg_batch_col[i] = tables[table_idx][col][idx]::S
        end
    end

    #Get all .arrow files in the input 
    input_paths = [path for path in readdir(input_dir, join=true) if endswith(path, ".arrow")]
    #println("input_paths protein ", input_paths)
    #Keep track of which tables have 
    tables = [Arrow.Table(path) for path in input_paths]
    table_idxs = ones(Int64, length(tables))

    #pre-allocate section of merged table 
    pg_batch = getEmptyDF(first(tables), N) #See /src/utils/mergePsmTables.jl for definition
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, N)
    #Names of columns for merge table 
    pg_batch_names = Symbol.(names(pg_batch))
    #Keeps track of the talbe with the highest ranked row
    pg_heap = BinaryMaxHeap{Tuple{Float32, Int64}}()
    for (i, table) in enumerate(tables)
        addRowToHeap!(
            pg_heap,
            table[sort_key],
            i,
            1
        )
    end
    i = 1
    n_writes = 0
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
            if iszero(n_writes)
                if isfile(output_path)
                    rm(output_path)
                end
                open(output_path, "w") do io
                    Arrow.write(io, pg_batch; file=false)  # file=false creates stream format
                end
            else
                Arrow.append(
                    output_path,
                    pg_batch
                )
            end
            n_writes += 1
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