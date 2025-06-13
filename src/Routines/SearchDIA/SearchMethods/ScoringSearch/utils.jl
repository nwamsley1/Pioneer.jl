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
    @info "[PERF] get_best_traces: Starting" n_files=length(second_pass_psms_paths)
    start_time = time()
    initial_memory = Base.gc_live_bytes() / 1024^2  # MB
    
    psms_trace_scores = Dictionary{
            @NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}, Float32}()

    # Track aggregated stats
    total_rows_processed = 0
    files_processed = 0
    
    for (file_idx, file_path) in enumerate(second_pass_psms_paths)
        if splitext(file_path)[end] != ".arrow"
            continue
        end
        row_score = zero(Float32)
        psms_table = Arrow.Table(file_path)
        n_rows = length(psms_table[1])
        total_rows_processed += n_rows
        files_processed += 1
        
        for i in range(1, n_rows)
            psms_key = (precursor_idx = psms_table[:precursor_idx][i],  isotopes_captured = psms_table[:isotopes_captured][i])

            row_score = psms_table[:prob][i]
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
    
    elapsed = time() - start_time
    final_memory = Base.gc_live_bytes() / 1024^2
    memory_used = final_memory - initial_memory
    @info "[PERF] get_best_traces: Completed" elapsed=round(elapsed, digits=3) memory_used_MB=round(memory_used, digits=1) n_traces=length(traces_passing) files_processed=files_processed total_rows=total_rows_processed avg_rows_per_file=round(total_rows_processed/max(files_processed,1), digits=1)
    
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
    isotope_trace_type::IsotopeTraceType,
    prob_col::Symbol,
    best_traces::Set{@NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}}
)
    @info "[PERF] sort_and_filter_quant_tables: Starting" n_files=length(second_pass_psms_paths)
    start_time = time()
    total_rows_processed = 0
    total_rows_kept = 0
    files_processed = 0

    #Remove if present 
    if isfile(merged_quant_path)
        safeRm(merged_quant_path)
    end
    #file_paths = [fpath for fpath in readdir(quant_psms_folder,join=true) if endswith(fpath,".arrow")]
    #Sort and filter each psm table 
    for (file_idx, fpath) in enumerate(second_pass_psms_paths)
        # Read file
        psms_table = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        initial_rows = size(psms_table, 1)
        total_rows_processed += initial_rows
        files_processed += 1

        if seperateTraces(isotope_trace_type)
            #Indicator variable of whether each psm is from the best trace 
            psms_table[!,:best_trace] = zeros(Bool, size(psms_table, 1))
            for i in range(1, size(psms_table, 1))
               key = (precursor_idx = psms_table[i, :precursor_idx], isotopes_captured = psms_table[i, :isotopes_captured])
               if key ∈ best_traces
                   psms_table[i,:best_trace]=true
               end
            end
        else
            transform!(
                groupby(psms_table, :precursor_idx),
                :prob => (p -> p .== maximum(p)) => :best_trace
            )
        end

        #Filter out unused traces
        filter!(x->x.best_trace,psms_table)
        filtered_rows = size(psms_table, 1)
        total_rows_kept += filtered_rows
        
        #Sort in descending order of probability
        sort!(psms_table, [prob_col, :target], rev = [true, true], alg=QuickSort)
        
        #write back
        writeArrow(fpath, psms_table)
    end
    
    elapsed = time() - start_time
    @info "[PERF] sort_and_filter_quant_tables: Completed" elapsed=round(elapsed, digits=3) files_processed=files_processed total_rows=total_rows_processed kept_rows=total_rows_kept retention_rate=round(total_rows_kept/max(total_rows_processed,1), digits=3) avg_rows_per_file=round(total_rows_processed/max(files_processed,1), digits=1)
    
    return nothing
end

function sort_quant_tables(
    second_pass_psms_paths::Vector{String},
    merged_quant_path::String,
    prob_col::Symbol,
)

    #Remove if present 
    if isfile(merged_quant_path)
        safeRm(merged_quant_path)
    end
    #file_paths = [fpath for fpath in readdir(quant_psms_folder,join=true) if endswith(fpath,".arrow")]
    #Sort and filter each psm table 
    for fpath in second_pass_psms_paths
        psms_table = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        #Sort in descending order of probability
        sort!(psms_table, [prob_col, :target], rev = [true, true], alg=QuickSort)
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
                    prob_col::Symbol
                    ; N = 10000000
)
    @info "[PERF] merge_sorted_psms_scores: Starting" n_files=length(input_paths) batch_size=N
    start_time = time()
    initial_memory = Base.gc_live_bytes() / 1024^2
    total_batches = 0
    
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

    function addPrecursorToHeap!(
        psms_heap::BinaryMaxHeap{Tuple{Float32, Int64}},
        sort_key::AbstractVector{S},
        table_idx::Int64,
        row_idx::Int64) where {S<:AbstractString}
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
    
    # Load all tables
    load_start = time()
    tables = [Arrow.Table(path) for path in input_paths]
    table_sizes = [length(table[1]) for table in tables]
    total_rows = sum(table_sizes)
    load_time = time() - load_start
    # Tables loaded
    
    table_idxs = ones(Int64, length(tables))
    psms_batch = DataFrame()
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, N)
    psms_heap = BinaryMaxHeap{Tuple{Float32, Int64}}()
    psms_batch[!,prob_col] = zeros(Float32, N)
    psms_batch[!,:target] = zeros(Bool, N)
    psms_batch[!,:precursor_idx] = zeros(UInt32, N)

    
    #Load psms_heap
    for (i, table) in enumerate(tables)
        ##println("i $i")
        ##println("size(table) ", size(DataFrame(table)))
        addPrecursorToHeap!(
            psms_heap,
            table[prob_col],
            i,
            1
        )
    end
    i = 1
    n_writes = 0
    rows_written = 0
    if Sys.iswindows()
        writeArrow(output_path, DataFrame())
    else
        safeRm(output_path, force=true)
    end
    
    merge_start = time()
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
            table[prob_col],
            table_idx,
            idx
        )
        i += 1
        if i > N
            for col in [prob_col,:target,:precursor_idx]
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
                    safeRm(output_path, force=true)
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
            rows_written += (i - 1)
            total_batches += 1
            
            # Check memory and GC if needed
            current_memory = Base.gc_live_bytes() / 1024^2
            if current_memory > initial_memory * 1.5
                GC.gc()
            end
            
            i = 1
        end
    end
    for col in [prob_col,:target,:precursor_idx]
        fillColumn!(
            psms_batch[!,col],
            col,
            sorted_tuples,
            tables,
            i
        )
    end
    if n_writes > 0
        Arrow.append(
            output_path,
            psms_batch[range(1, max(1, i - 1)),:]
        )
    else
        writeArrow(
            output_path,
            psms_batch[range(1, max(1, i - 1)),:]
        )
    end
    rows_written += max(1, i - 1)
    elapsed = time() - start_time
    final_memory = Base.gc_live_bytes() / 1024^2
    @info "[PERF] merge_sorted_psms_scores: Completed" elapsed=round(elapsed, digits=3) n_files=length(input_paths) total_batches=total_batches total_rows=rows_written rows_per_sec=round(rows_written/elapsed, digits=1) memory_used_MB=round(final_memory-initial_memory, digits=1) avg_rows_per_file=round(total_rows/length(input_paths), digits=1)
    
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
                      min_pep_points_per_bin=1000, fdr_scale_factor=1.0f0) -> Interpolation

Create q-value interpolation function from merged scores.

Returns interpolation for mapping scores to q-values.
"""
function get_qvalue_spline(
                            merged_psms_path::String,
                            score_col::Symbol,
                            use_unique::Bool;
                            min_pep_points_per_bin = 1000,
                            fdr_scale_factor::Float32 = 1.0f0
)


    psms_scores = DataFrame(Arrow.Table(merged_psms_path))
    
    if use_unique
        psms_scores = unique(psms_scores)
    end

    Q = size(psms_scores, 1)
    M = ceil(Int, Q / min_pep_points_per_bin)
    bin_qval, bin_mean_prob = Vector{Float32}(undef, M), Vector{Float32}(undef, M)
    bin_size = 0
    bin_idx = 0
    mean_prob, targets, decoys = 0.0f0, 0, 0
    targets, decoys = 0, 0
    for i in range(1, Q)
        targets += psms_scores[!, :target][i]
        decoys += (1 - psms_scores[!, :target][i])
    end

    min_q_val = typemax(Float32)
    for i in reverse(range(1, Q))
        bin_size += 1
        targets -= psms_scores[!, :target][i]
        decoys -= (1 - psms_scores[!, :target][i])
        mean_prob += psms_scores[!, score_col][i]
        if bin_size == min_pep_points_per_bin
            bin_idx += 1
            # Apply FDR scale factor to correct for library target/decoy ratio
            qval = (decoys * fdr_scale_factor) / max(targets, 1)
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
    # Apply FDR scale factor to final bin calculation
    bin_qval[end] = (decoys * fdr_scale_factor) / max(targets, 1)
    bin_mean_prob[end] = mean_prob/bin_size
    prepend!(bin_qval, 1.0f0)
    prepend!(bin_mean_prob, 0.0f0)
    bin_qval = bin_qval[isnan.(bin_mean_prob).==false]
    bin_mean_prob = bin_mean_prob[isnan.(bin_mean_prob).==false]
    
    # Sort and deduplicate knot vectors
    # First, create paired array for sorting
    paired = collect(zip(bin_mean_prob, bin_qval))
    # Sort by probability (x-values)
    sort!(paired, by = x -> x[1])
    
    # Manual deduplication keeping the first occurrence of each x value
    xs = Float32[]
    ys = Float32[]
    prev_x = NaN32
    for (x, y) in paired
        if x != prev_x && !isnan(x)
            push!(xs, x)
            push!(ys, y)
            prev_x = x
        end
    end
    
    # Ensure we have at least 2 points for interpolation
    if length(xs) < 2
        @warn "Insufficient unique points for q-value interpolation, using default"
        xs = Float32[0.0, 1.0]
        ys = Float32[1.0, 0.0]
    end
    
    # Final check for sorted and unique
    if length(xs) > 1
        for i in 2:length(xs)
            if xs[i] <= xs[i-1]
                error("Failed to properly deduplicate knot vectors. Debug info: xs=$xs")
            end
        end
    end
    #=
    xs = copy(bin_mean_prob)
    ys = copy(bin_qval)
    # remove duplicates and keep the **first** occurrence
    mask = [i == 1 || xs[i] != xs[i-1] for i in eachindex(xs)]
    xs   = xs[mask]
    ys   = ys[mask]
    =#
    return linear_interpolation(xs, ys; extrapolation_bc = Line())
end

"""
    get_psms_passing_qval(precursors::LibraryPrecursors, passing_psms_paths::Vector{String}, 
                         passing_psms_folder::String, second_pass_psms_paths::Vector{String}, 
                         pep_spline::UniformSpline, qval_interp::Interpolations.Extrapolation, 
                         q_val_threshold::Float32)

Filter PSMs that pass a given q-value threshold and calculate error probabilities.

# Arguments
- `precursors::LibraryPrecursors`: Library precursor information
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
                            precursors::LibraryPrecursors,
                            passing_psms_paths::Vector{String},
                            passing_psms_folder::String,
                            second_pass_psms_paths::Vector{String}, 
                            #pep_spline::UniformSpline,
                            qval_interp_global::Interpolations.Extrapolation,
                            qval_interp_experiment_wide::Interpolations.Extrapolation,
                            prob_col_global::Symbol,
                            prob_col_experiment_wide::Symbol,
                            qval_col_global::Symbol,
                            qval_col_experiment_wide::Symbol,
                            q_val_threshold::Float32,)
    for (ms_file_idx, file_path) in enumerate(second_pass_psms_paths)
        # Read the Arrow table
        passing_psms = DataFrame(Tables.columntable(Arrow.Table(file_path)))
        passing_psms[!,qval_col_global] = qval_interp_global.(passing_psms[!,prob_col_global])
        passing_psms[!,qval_col_experiment_wide] = qval_interp_experiment_wide.(passing_psms[!,prob_col_experiment_wide])
        #passing_psms[!,:pep] = pep_spline.(passing_psms[!,prob_col])

        cols = [
            :precursor_idx,
            :global_prob,
            :prec_prob,
            :prob,
            :global_qval,
            :run_specific_qval,
            :prec_mz,
            #:pep,
            :weight,
            :target,
            :rt,
            :irt_obs,
            :missed_cleavage,
            :isotopes_captured,
            :scan_idx,
            :entrapment_group_id,
            :ms_file_idx
        ]
        available_cols = intersect(cols, Symbol.(names(passing_psms)))

        # Sample the rows and convert to DataFrame
        select!(passing_psms,available_cols)
        filter!(x->(x[qval_col_global]<=q_val_threshold)&(x[qval_col_experiment_wide]<=q_val_threshold), passing_psms)
        # Append to the result DataFrame
        #Arrow.write(
        #    joinpath(passing_psms_folder, basename(file_path)),
        #    passing_psms
        #)
        writeArrow(joinpath(passing_psms_folder, basename(file_path)),
        passing_psms)
        passing_psms_paths[ms_file_idx] = joinpath(passing_psms_folder, basename(file_path))
    end
    return
end

"""
    get_psms_passing_qval(precursors::PlexedLibraryPrecursors, passing_psms_paths::Vector{String}, 
                         passing_psms_folder::String, second_pass_psms_paths::Vector{String}, 
                         pep_spline::UniformSpline, qval_interp::Interpolations.Extrapolation, 
                         q_val_threshold::Float32)

Filter PSMs that pass a given group-level q-value threshold and calculate error probabilities for plexed library precursors.

# Arguments
- `precursors::PlexedLibraryPrecursors`: Plexed library precursor information
- `passing_psms_paths`: Vector to store paths of filtered PSM files
- `passing_psms_folder`: Folder to store passing PSM files
- `second_pass_psms_paths`: Paths to second pass PSM files
- `pep_spline`: Spline for posterior error probability calculation
- `qval_interp`: Interpolation function for q-value calculation
- `q_val_threshold`: Q-value threshold for filtering

# Process
1. Reads all PSMs from second pass files
2. Adds sequence and structural modification information from precursors
3. Calculates group-level q-values by:
   - Grouping by ms_file_idx, sequence, and structural_mods
   - Finding maximum probability for each group
   - Computing FDR and q-values at the group level
4. Calculates PSM-level q-values and PEP
5. Filters PSMs below the group q-value threshold
6. Saves filtered PSMs by file
7. Updates passing_psms_paths with new file locations

# Selected Columns
- precursor_idx, prob, qval, group_qvalue
- weight, target, irt_obs
- missed_cleavage, isotopes_captured
- scan_idx, ms_file_idx
"""
function get_psms_passing_qval(
                            precursors::PlexedLibraryPrecursors,
                            passing_psms_paths::Vector{String},
                            passing_psms_folder::String,
                            second_pass_psms_paths::Vector{String}, 
                            pep_spline::UniformSpline,
                            qval_interp::Interpolations.Extrapolation,
                            q_val_threshold::Float32)

    function groupLevelQval(df::DataFrame)
        # Step 1: Group by the specified columns and get max probability for each group
        grouped_df = combine(groupby(df, [:ms_file_idx, :sequence, :structural_mods, :charge]),
                            :prob => maximum => :group_max_prob,
                            :decoy => first => :group_decoy)
        
        # Step 2: Calculate FDR and q-values based on the max probability
        # Sort by descending probability
        sort!(grouped_df, :group_max_prob, rev=true)
        
        # Calculate FDR for each row
        grouped_df.group_fdr = zeros(Float32, size(grouped_df, 1))
        cumulative_decoys = 0
        cumulative_targets = 0
        
        for i in 1:nrow(grouped_df)
            if grouped_df.group_decoy[i]
                cumulative_decoys += 1
            else
                cumulative_targets += 1
            end
            
            # FDR = (decoys * scale_factor) / targets
            if cumulative_targets == 0
                grouped_df.group_fdr[i] = 1.0
            else
                grouped_df.group_fdr[i] = min(1.0, (cumulative_decoys) / cumulative_targets)
            end
        end
        
        # Step 3: Convert FDR to q-value (monotonically non-decreasing values)
        grouped_df.group_qvalue = copy(grouped_df.group_fdr)
        
        # Ensure q-values are monotonically non-decreasing from bottom to top
        for i in (nrow(grouped_df)-1):-1:1
            grouped_df.group_qvalue[i] = min(grouped_df.group_qvalue[i], 
                                              grouped_df.group_qvalue[i+1])
        end
        
        # Step 4: Join the group-level information back to the original dataframe
        # Create a temporary dataframe with just the keys and new columns
        temp_df = select(grouped_df, 
                         [:ms_file_idx, :sequence, :structural_mods, 
                          :group_max_prob, :group_qvalue])
        
        # Join this information back to the original dataframe
        result_df = leftjoin(df, temp_df, 
                             on=[:ms_file_idx, :sequence, :structural_mods, :charge])
        
        return result_df
    end

    #If at least one plex passed 
    function minQval(df::DataFrame)
        # Group by sequence and structural_mods (across all plexes/files)
        # Find the minimum q-value for each peptide group
        min_qval_df = combine(
            groupby(df, [:ms_file_idx, :sequence, :structural_mods, :charge]),
            :run_specific_qval => minimum => :min_qval
        )
        
        temp_df = select(min_qval_df, 
        [:ms_file_idx, :sequence, :structural_mods, 
         :charge, :min_qval])

        # Join this information back to the original dataframe
        result_df = leftjoin(
            df, 
            temp_df, 
            on = [:ms_file_idx, :sequence, :structural_mods, :charge]
        )
        
        return result_df
    end
    psms = DataFrame(Tables.columntable(Arrow.Table(second_pass_psms_paths)))
    psms[!,:decoy] = psms[!,:target].==false
    psms[!,:sequence] = [getSequence(precursors)[pid] for pid in psms[!,:precursor_idx]]
    psms[!,:structural_mods] = [getStructuralMods(precursors)[pid] for pid in psms[!,:precursor_idx]]
    #psms = groupLevelQval(psms)
    psms[!,:run_specific_qval] = qval_interp.(psms[!,:prob])
    psms = minQval(psms)
    psms[!,:pep] = pep_spline.(psms[!,:prob])
    select!(psms,
    [
        :precursor_idx,
        :prob,
        :run_specific_qval,
        :min_qval,
        :prec_mz,
        #:pep,
        :weight,
        :target,
        :rt,
        :irt_obs,
        :missed_cleavage,
        :isotopes_captured,
        :scan_idx,
        :ms_file_idx])
    #filter!(x->x.group_qvalue<=q_val_threshold, psms)
    filter!(x->x.min_qval<=q_val_threshold, psms)
    psms_by_file = groupby(psms, :ms_file_idx)
    for (ms_file_idx, file_path) in enumerate(second_pass_psms_paths)
        file_psms = psms_by_file[(ms_file_idx = ms_file_idx,)]
        Arrow.write(
            joinpath(passing_psms_folder, basename(file_path)),
            file_psms
        )
        passing_psms_paths[ms_file_idx] = joinpath(passing_psms_folder, basename(file_path))
    end

    return
end

#==========================================================
Protein group analysis
==========================================================#
"""
    getProteinGroupsDict(protein_inference_dict, psm_precursor_idx, psm_score, 
                        psm_is_target, psm_entrapment_id, precursors; min_peptides=2)

Create protein groups from PSMs and calculate group scores.

# Arguments
- `protein_inference_dict`: Maps peptides to inferred protein groups
- `psm_precursor_idx`: Precursor indices from PSMs
- `psm_score`: PSM probability scores
- `psm_is_target`: Boolean array indicating targets
- `psm_entrapment_id`: Entrapment group IDs
- `precursors`: Library precursor information
- `min_peptides`: Minimum peptides required per group

# Returns
- `pg_score`: Protein group scores for each PSM
- `inferred_protein_group_names`: Protein names for each PSM
- `protein_groups`: Dictionary of protein groups with scores and peptide sets
"""
function getProteinGroupsDict(
    protein_inference_dict::Dictionary{NamedTuple{(:peptide, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}, NamedTuple{(:protein_name, :decoy, :entrap_id, :retain), Tuple{String, Bool, UInt8, Bool}}},
    psm_precursor_idx::AbstractVector{UInt32},
    psm_score::AbstractVector{Float32},
    psm_is_target::AbstractVector{Bool},
    psm_entrapment_id::AbstractVector{UInt8},
    precursors::LibraryPrecursors;
    min_peptides::Int64 = 1)

    #accession_numbers = getAccessionNumbers(precursors)
    precursor_sequence = getSequence(precursors)
    protein_groups = Dictionary{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},
    @NamedTuple{
        pg_score::Float32, 
        peptides::Set{String}}
    }()

    #inferred_protein_group_names = Vector{Union{Missing, String}}(undef, length(psm_precursor_idx))
    for i in range(1, length(psm_precursor_idx))
        precursor_idx = psm_precursor_idx[i]
        sequence = precursor_sequence[precursor_idx]
        
        # Create key for protein_inference_dict lookup
        peptide_key = (peptide = sequence, decoy = !psm_is_target[i], entrap_id = psm_entrapment_id[i])
        
        # Check if this peptide exists in our protein inference dictionary
        if !haskey(protein_inference_dict, peptide_key)
            throw("Peptide key not found error!")
            continue
        end
        #inferred_protein_group_names[i] = protein_inference_dict[peptide_key][:protein_name]
        # Exclude peptide 
        if protein_inference_dict[peptide_key][:retain] == false
            continue
        end
        
        score = psm_score[i]
        protein_name = protein_inference_dict[peptide_key][:protein_name]
        keyname = (protein_name = protein_name, target = psm_is_target[i], entrap_id = psm_entrapment_id[i])
        
        if haskey(protein_groups, keyname)
            pg_score, peptides = protein_groups[keyname]
            pg_score += log1p(-score)
            push!(peptides, sequence)
            protein_groups[keyname] = (pg_score = pg_score, peptides = peptides)
        else
            sequences = Set{String}((sequence,))
            insert!(protein_groups,
                keyname,
                (pg_score = log1p(-score),
                peptides = sequences)
            )
        end
    end
    
    filter!(x->length(x[:peptides])>=min_peptides, protein_groups)

    for key in keys(protein_groups)
        pg_score, peptides = protein_groups[key]
        pg_score = -pg_score
        protein_groups[key] = (pg_score = pg_score, peptides = peptides)
    end
    
    return protein_groups#inferred_protein_group_names, protein_groups#pg_score, inferred_protein_group_names, protein_groups
end


"""
    writeProteinGroups(acc_to_max_pg_score, protein_groups, 
                        protein_to_possible_peptides, protein_groups_path)

Write protein groups with features to Arrow file.

# Arguments
- `acc_to_max_pg_score`: Maximum scores across runs for each protein
- `protein_groups`: Dictionary of protein groups with scores and peptides
- `protein_to_possible_peptides`: All possible peptides for each protein
- `protein_groups_path`: Output file path

# Returns
- Number of protein groups written

# Output columns
- Basic: protein_name, target, entrap_id, pg_score, global_pg_score
- Features: n_peptides, total_peptide_length, n_possible_peptides, peptide_coverage
"""
function writeProteinGroups(
                                protein_groups::Dictionary{
                                    @NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},
                                    @NamedTuple{pg_score::Float32,  peptides::Set{String}}
                                },
                                protein_accession_to_possible_peptides::Dict{
                                    @NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},
                                    Set{String}
                                },
                                protein_groups_path::String)
    # Extract keys and values
    keys_array = keys(protein_groups)
    values_array = values(protein_groups)

    # Create vectors for each column
    protein_name = [k[:protein_name] for k in keys_array]
    target = [k[:target] for k in keys_array]
    entrap_id = [k[:entrap_id] for k in keys_array]
    pg_score = [v[:pg_score] for v in values_array]
    #peptides = [join(v[:peptides], ";") for v in values_array]  # Convert Set to String
    
    # New feature columns
    n_peptides = [UInt16(length(unique(v[:peptides]))) for v in values_array]  # Number of unique peptides
    total_peptide_length = [UInt16(sum(length(pep) for pep in v[:peptides])) for v in values_array]  # Total length of all peptides
    
    # Calculate possible peptides and peptide coverage
    # Handle protein groups with multiple proteins separated by semicolons
    n_possible_peptides = zeros(UInt16, length(keys_array))
    for (i, k) in enumerate(keys_array)
        # Split the protein group name by semicolons.
        # Key distinction from protein group and accession number is_decoy
        # a protein group can be a list of accession numbers concatenated with ';' delimeters
        protein_names_in_group = split(k[:protein_name], ';')
        
        # Union of all peptide sets from proteins in the group
        all_possible_peptides = Set{String}()
        for individual_protein in protein_names_in_group
            # Create key for each individual protein
            individual_key = (protein_name = String(individual_protein), 
                            target = k[:target], 
                            entrap_id = k[:entrap_id])
            # Get the set of peptides for this protein and union with existing
            if haskey(protein_accession_to_possible_peptides, individual_key)
                union!(all_possible_peptides, protein_accession_to_possible_peptides[individual_key])
            end
        end
        
        # Count unique peptides across all proteins in the group
        n_possible_peptides[i] = UInt16(max(length(all_possible_peptides), 1))
    end
    
    peptide_coverage = Float32[n_pep / n_poss for (n_pep, n_poss) in zip(n_peptides, n_possible_peptides)]
    
    # Calculate derived features for probit regression
    log_n_possible_peptides = log.(n_possible_peptides .+ 1)
    log_binom_coeff = [
        if n_obs <= n_poss && n_obs >= 0
            lgamma(n_poss + 1) - lgamma(n_obs + 1) - lgamma(n_poss - n_obs + 1)
        else
            0.0
        end
        for (n_poss, n_obs) in zip(n_possible_peptides, n_peptides)
    ]
    
    # No global scores yet
    df = DataFrame((
        protein_name = protein_name,
        target = target,
        entrap_id = entrap_id,
        pg_score = pg_score,
        n_peptides = n_peptides,
        total_peptide_length = total_peptide_length,
        n_possible_peptides = n_possible_peptides,
        peptide_coverage = peptide_coverage,
        log_n_possible_peptides = log_n_possible_peptides,
        log_binom_coeff = log_binom_coeff
    ))
    sort!(df, [:pg_score,:target], rev = [true,true])
    # Convert DataFrame to Arrow.Table
    Arrow.write(protein_groups_path, df)
    return
end

"""
    count_protein_peptides(precursors::LibraryPrecursors)

Count all possible peptides for each protein in the library.

# Arguments
- `precursors`: Library precursors

# Returns
- Dictionary mapping protein keys to sets of peptide sequences
"""
function count_protein_peptides(precursors::LibraryPrecursors)
    peptide_count_start = time()
    protein_to_possible_peptides = Dict{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8}, Set{String}}()
    
    # Count all peptides in the library for each protein
    all_accession_numbers = getAccessionNumbers(precursors)
    all_sequences = getSequence(precursors)
    all_decoys = getIsDecoy(precursors)
    all_entrap_ids = getEntrapmentGroupId(precursors)
    n_precursors = length(all_accession_numbers)
    @info "[PERF] count_protein_peptides: Starting" n_precursors=n_precursors
    
    for i in 1:n_precursors
        protein_names = split(all_accession_numbers[i], ';')  # Handle shared peptides
        is_decoy = all_decoys[i]
        entrap_id = all_entrap_ids[i]
        
        for protein_name in protein_names
            key = (protein_name = String(protein_name), target = !is_decoy, entrap_id = entrap_id)
            if !haskey(protein_to_possible_peptides, key)
                protein_to_possible_peptides[key] = Set{String}()
            end
            push!(protein_to_possible_peptides[key], all_sequences[i])
        end
    end
    
    peptide_count_time = time() - peptide_count_start
    @info "[PERF] count_protein_peptides: Completed" elapsed=round(peptide_count_time, digits=3) n_proteins=length(protein_to_possible_peptides)
    
    return protein_to_possible_peptides
end

"""
    perform_protein_inference(passing_psms_paths::Vector{String},
                                passing_pg_paths::Vector{String},
                                protein_groups_folder::String,
                                precursors::LibraryPrecursors,
                                protein_to_possible_peptides::Dict;
                                min_peptides::Int64 = 2)

Perform protein inference and initial scoring (first and second pass).

# Arguments
- `passing_psms_paths`: Paths to PSM files
- `passing_pg_paths`: Paths to protein group files (will be populated)
- `protein_groups_folder`: Folder to write protein groups
- `precursors`: Library precursors
- `protein_to_possible_peptides`: Dictionary of protein to peptide mappings
- `min_peptides`: Minimum peptides per protein group

# Returns
- `psm_to_pg_path`: Dictionary mapping PSM file paths to protein group file paths
"""
function perform_protein_inference(
    passing_psms_paths::Vector{String},
    passing_pg_paths::Vector{String},
    protein_groups_folder::String,
    precursors::LibraryPrecursors,
    protein_to_possible_peptides::Dict{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8}, Set{String}};
    min_peptides::Int64 = 1
)
    start_time = time()
    
    # Get sequences and other data from precursors
    all_sequences = getSequence(precursors)
    all_decoys = getIsDecoy(precursors)
    all_entrap_ids = getEntrapmentGroupId(precursors)
    
    # Create bidirectional mappings between PSM and protein group files
    psm_to_pg_path = Dict{String, String}()

    # Single pass: compute protein groups and write immediately
    @info "[PERF] perform_protein_inference: Starting combined inference and write pass"
    total_psms_processed = 0
    total_peptides_inferred = 0
    total_protein_groups_created = 0
    files_processed = 0
    
    for (ms_file_idx, file_path) in enumerate(passing_psms_paths)
        _, extention = splitext(file_path)
        if extention != ".arrow"
            continue
        end
        protein_groups_path = joinpath(protein_groups_folder, basename(file_path))
        passing_pg_paths[ms_file_idx] = protein_groups_path
        psms_table = Arrow.Table(file_path)
        n_psms = length(psms_table[:precursor_idx])
        total_psms_processed += n_psms
        files_processed += 1
        
        # Prepare for protein inference. Could combine build_protein_peptide_rows_for_file and infer_proteins
        #in one infer_proteins call w/ a wrapper function 
        # Get list of protein/peptide pairs in the psms file. 
        protein_peptide_rows = build_protein_peptide_rows_for_file(psms_table, precursors)
        peptides = [row.sequence for row in protein_peptide_rows]
        proteins = [(protein_name = row.protein_name, decoy = row.decoy, entrap_id = row.entrap_id) for row in protein_peptide_rows]
        total_peptides_inferred += length(peptides)
        #=
        returns: ::Dictionary{
                NamedTuple{(:peptide, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}, 
                NamedTuple{(:protein_name, :decoy, :entrap_id, :retain), Tuple{String, Bool, UInt8, Bool}}
        }
        Maps peptide names to their infered protein groups. The 'protein_name' of the protein group can be 
        multiple accession numbers/proein names seperated by semi-colons and sorted alpha-numerically.  
        'retain' refers to whether or not the peptide key should be used to score and quantify the protein group. 
        =#
        protein_inference_dict = infer_proteins(proteins, peptides)

        # Convert and write PSMs with protein inference info
        psms_table = DataFrame(Tables.columntable(psms_table))

        # Add use_for_protein_quant column based on protein inference results
        precursor_idx = psms_table[!,:precursor_idx]
        use_for_protein_quant = zeros(Bool, length(precursor_idx))
        inferred_protein_names = Vector{String}(undef, length(precursor_idx))
        for (i, pid) in enumerate(precursor_idx)
            inferred_prot = protein_inference_dict[(peptide = all_sequences[pid], decoy = all_decoys[pid], entrap_id = all_entrap_ids[pid])]
            use_for_protein_quant[i] = inferred_prot.retain
            inferred_protein_names[i] = inferred_prot.protein_name
        end
        psms_table[!,:use_for_protein_quant] = use_for_protein_quant
        psms_table[!,:inferred_protein_group] = inferred_protein_names
        filter!(x->x.retain, protein_inference_dict)
        #protein groups maps protein names to scores and peptide set 
        #=
        like so 
        protein_groups = Dictionary{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},
        @NamedTuple{
            pg_score::Float32, 
            peptides::Set{String}}
        }()
        =#
        #inferred_protein_group_names, 
        pg_to_peptides_dict = getProteinGroupsDict(
            protein_inference_dict,
            psms_table[!,:precursor_idx],
            psms_table[!,:prob],
            psms_table[!,:target],
            psms_table[!,:entrapment_group_id],
            precursors;
            min_peptides = min_peptides
        )
        writeArrow(file_path, psms_table)
        
        # Write protein groups immediately (without holding in memory)
        writeProteinGroups(
            pg_to_peptides_dict,
            protein_to_possible_peptides,
            protein_groups_path
        )
        total_protein_groups_created += length(protein_groups)
        
        # Store bidirectional mapping
        psm_to_pg_path[file_path] = protein_groups_path
        
        # Clear protein_groups from memory immediately
        protein_groups = nothing
    end
    
    total_elapsed = time() - start_time
    @info "[PERF] perform_protein_inference: Completed" total_elapsed=round(total_elapsed, digits=3) files_processed=files_processed total_psms=total_psms_processed total_peptides_inferred=total_peptides_inferred total_protein_groups=total_protein_groups_created avg_psms_per_file=round(total_psms_processed/max(files_processed,1), digits=1)
    
    return psm_to_pg_path
end

"""
    perform_protein_probit_regression(passing_pg_paths::Vector{String},
                                    max_psms_in_memory::Int64,
                                    qc_folder::String)

Perform probit regression on protein groups.

# Arguments
- `passing_pg_paths`: Paths to protein group files
- `max_psms_in_memory`: Memory limit for in-memory vs OOM processing
- `qc_folder`: Folder for QC plots
"""
function perform_protein_probit_regression(
    passing_pg_paths::Vector{String},
    max_psms_in_memory::Int64,
    qc_folder::String
)
    # First, count the total number of protein groups across all files
    total_protein_groups = 0
    for pg_path in passing_pg_paths
        if isfile(pg_path) && endswith(pg_path, ".arrow")
            table = Arrow.Table(pg_path)
            total_protein_groups += length(table[:protein_name])  # Get row count from any column
        end
    end
    
    # Set protein group limit to 5x the precursor limit
    max_protein_groups_in_memory_limit = 5 * max_psms_in_memory
    
    if total_protein_groups > max_protein_groups_in_memory_limit
        @info "Using out-of-memory probit regression (exceeds limit of $max_protein_groups_in_memory_limit)"
        perform_probit_analysis_oom(passing_pg_paths, total_protein_groups, max_protein_groups_in_memory_limit, qc_folder)
    else
        @info "Using in-memory probit regression"
        # Load all protein group tables into a single DataFrame
        all_protein_groups = DataFrame()
        for pg_path in passing_pg_paths
            if isfile(pg_path) && endswith(pg_path, ".arrow")
                append!(all_protein_groups, DataFrame(Tables.columntable(Arrow.Table(pg_path))))
            end
        end
        
        # Only perform analysis if we have both targets and decoys
        n_targets = sum(all_protein_groups.target)
        n_decoys = sum(.!all_protein_groups.target)
        
        if n_targets > 0 && n_decoys > 0 && nrow(all_protein_groups) > 10
            # Perform probit regression analysis
            @info "Performing Probit Analysis (targets: $n_targets, decoys: $n_decoys)"
            perform_probit_analysis(all_protein_groups, qc_folder, passing_pg_paths)
        else
            @info "Skipping Probit analysis: insufficient data (targets: $n_targets, decoys: $n_decoys)"
        end
    end
end

"""
    calculate_global_protein_scores(passing_pg_paths::Vector{String})

Calculate global protein scores from current pg_scores.

# Arguments
- `passing_pg_paths`: Paths to protein group files

# Returns
- Dictionary mapping protein keys to max pg_score across files
"""
function calculate_global_protein_scores(passing_pg_paths::Vector{String})
    third_pass_start = time()
    @info "[PERF] calculate_global_protein_scores: Starting"
    
    # Build acc_to_max_pg_score from current pg_scores (either original or probit-overwritten)
    acc_to_max_pg_score = Dict{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},Float32}()
    
    for pg_path in passing_pg_paths
        if isfile(pg_path) && endswith(pg_path, ".arrow")
            table = Arrow.Table(pg_path)
            n_rows = length(table[:protein_name])
            
            for i in 1:n_rows
                key = (protein_name = table[:protein_name][i], 
                        target = table[:target][i], 
                        entrap_id = table[:entrap_id][i])
                pg_score = table[:pg_score][i]  # Either original or probit score
                
                old = get(acc_to_max_pg_score, key, -Inf32)
                acc_to_max_pg_score[key] = max(pg_score, old)
            end
        end
    end
    
    # Add global_pg_score to all protein group files
    fourth_pass_start = time()
    @info "[PERF] calculate_global_protein_scores: Adding global_pg_score to protein group files"
    
    for pg_path in passing_pg_paths
        if isfile(pg_path) && endswith(pg_path, ".arrow")
            df = DataFrame(Tables.columntable(Arrow.Table(pg_path)))
            
            # Add global_pg_score column
            global_pg_scores = Vector{Float32}(undef, nrow(df))
            for i in 1:nrow(df)
                key = (protein_name = df[i, :protein_name],
                        target = df[i, :target],
                        entrap_id = df[i, :entrap_id])
                global_pg_scores[i] = get(acc_to_max_pg_score, key, df[i, :pg_score])
            end
            df[!, :global_pg_score] = global_pg_scores
            
            # Re-sort by global_pg_score
            sort!(df, [:global_pg_score, :target], rev = [true, true])
            
            # Write back
            writeArrow(pg_path, df)
        end
    end
    
    fourth_pass_time = time() - fourth_pass_start
    @info "[PERF] calculate_global_protein_scores: Completed" elapsed=round(fourth_pass_time, digits=3) n_protein_groups=length(acc_to_max_pg_score)
    
    return acc_to_max_pg_score
end

"""
    add_global_scores_to_psms(passing_psms_paths::Vector{String},
                                acc_to_max_pg_score::Dict)

Add global_pg_score to PSM files.

# Arguments
- `passing_psms_paths`: Paths to PSM files
- `acc_to_max_pg_score`: Dictionary mapping protein keys to global scores
"""
function add_global_scores_to_psms(
    passing_psms_paths::Vector{String},
    acc_to_max_pg_score::Dict{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},Float32}
)
    fifth_pass_start = time()
    @info "[PERF] add_global_scores_to_psms: Starting"
    
    for (ms_file_idx, file_path) in enumerate(passing_psms_paths)
        if isfile(file_path) && endswith(file_path, ".arrow")
            psms_table = DataFrame(Tables.columntable(Arrow.Table(file_path)))
            
            # Add global_pg_score column by looking up each protein group
            global_pg_scores = Vector{Float32}(undef, nrow(psms_table))
            for i in 1:nrow(psms_table)
                key = (protein_name = psms_table[i, :inferred_protein_group],
                        target = psms_table[i, :target],
                        entrap_id = psms_table[i, :entrapment_group_id])
                global_pg_scores[i] = get(acc_to_max_pg_score, key, 0.0f0)
            end
            psms_table[!, :global_pg_score] = global_pg_scores
            
            # Write back
            writeArrow(file_path, psms_table)
        end
    end
    
    fifth_pass_time = time() - fifth_pass_start
    @info "[PERF] add_global_scores_to_psms: Completed" elapsed=round(fifth_pass_time, digits=3)
end

"""
    update_psms_with_probit_scores(psm_to_pg_path::Dict{String,String},
                                  acc_to_max_pg_score::Dict)

Update PSMs with probit-scored pg_score values after regression.

# Arguments
- `psm_to_pg_path`: Dictionary mapping PSM file paths to protein group file paths
- `acc_to_max_pg_score`: Dictionary mapping protein keys to global scores

# Process
1. For each PSM file:
   - Load corresponding protein group file (has probit-scored pg_score)
   - Create lookup dictionary for probit scores
   - Update PSM pg_score values with probit scores
   - Add global_pg_score column
   - Write back updated PSM file
"""
function update_psms_with_probit_scores(
    ptable::Any,
    psm_to_pg_path::Dict{String, String},
    acc_to_max_pg_score::Dict{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},Float32},
    pg_score_to_qval::Interpolations.Extrapolation,
    global_pg_score_to_qval::Interpolations.Extrapolation
)
    start_time = time()
    @info "[PERF] update_psms_with_probit_scores: Starting"
    
    total_psms_updated = 0
    files_processed = 0
    
    for (psm_path, pg_path) in psm_to_pg_path
        # Verify both files exist
        if !isfile(psm_path) || !isfile(pg_path)
            @warn "Missing file" psm_path=psm_path pg_path=pg_path
            continue
        end
        
        file_start = time()
        
        # Load protein groups to get probit-scored pg_score values
        pg_table = Arrow.Table(pg_path)
        
        # Create lookup dictionary: (protein_name, target, entrap_id) -> probit pg_score
        pg_score_lookup = Dict{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8}, Float32}()
        n_pg_rows = length(pg_table[:protein_name])
        
        for i in 1:n_pg_rows
            key = (protein_name = pg_table[:protein_name][i],
                   target = pg_table[:target][i],
                   entrap_id = pg_table[:entrap_id][i])
            pg_score_lookup[key] = pg_table[:pg_score][i]  # This is now the probit score
        end
        
        # Load PSMs
        psms_df = DataFrame(Tables.columntable(Arrow.Table(psm_path)))
        n_psms = nrow(psms_df)
        
        # Update pg_score column with probit scores
        probit_pg_scores = Vector{Float32}(undef, n_psms)
        global_pg_scores = Vector{Float32}(undef, n_psms)
        
        for i in 1:n_psms
            # Skip if missing inferred protein group
            if ismissing(psms_df[i, :inferred_protein_group])
                throw("Missing Inferred Protein Group!!!")
                probit_pg_scores[i] = 0.0f0
                global_pg_scores[i] = 0.0f0
                continue
            end
            
            # Create key for lookup
            key = (protein_name = psms_df[i, :inferred_protein_group],
                   target = psms_df[i, :target],
                   entrap_id = psms_df[i, :entrapment_group_id])
            
            # Get probit pg_score
            if !haskey(pg_score_lookup, key)
                println("key $key")
                pid = psms_df[i,:precursor_idx]
                println("getAccessionNumbers(ptable)[pid] ", getAccessionNumbers(ptable)[pid])
                throw("Missing pg score lookup key!!!")
            end
            probit_pg_scores[i] = get(pg_score_lookup, key, 0.0f0)
            if !haskey(acc_to_max_pg_score, key)
                throw("Missing global pg score lookup key!!!")
            end
            # Get global pg_score
            global_pg_scores[i] = get(acc_to_max_pg_score, key, 0.0f0)
        end
        
        # Update columns
        psms_df[!, :pg_score] = probit_pg_scores
        psms_df[!, :global_pg_score] = global_pg_scores
        
        # Calculate q-values using interpolation functions
        pg_qvals = Vector{Float32}(undef, n_psms)
        global_pg_qvals = Vector{Float32}(undef, n_psms)
        
        for i in 1:n_psms
            if ismissing(probit_pg_scores[i])
                pg_qvals[i] = 1.0f0
                global_pg_qvals[i] = 1.0f0
            else
                pg_qvals[i] = pg_score_to_qval(probit_pg_scores[i])
                global_pg_qvals[i] = global_pg_score_to_qval(global_pg_scores[i])
            end
        end
        
        # Add q-value columns
        psms_df[!, :pg_qval] = pg_qvals
        psms_df[!, :global_qval_pg] = global_pg_qvals
        
        # Write back
        writeArrow(psm_path, psms_df)
        
        total_psms_updated += n_psms
        files_processed += 1
        
        file_time = time() - file_start
        @debug "Processed file" file=basename(psm_path) n_psms=n_psms elapsed=round(file_time, digits=3)
    end
    
    elapsed = time() - start_time
    @info "[PERF] update_psms_with_probit_scores: Completed" elapsed=round(elapsed, digits=3) files_processed=files_processed total_psms_updated=total_psms_updated avg_psms_per_file=round(total_psms_updated/max(files_processed,1), digits=1)
end

"""
    build_protein_peptide_rows_for_file(psms_table, precursors)

Build protein-peptide mappings from a single file's PSMs.

# Arguments
- `psms_table`: Arrow table containing PSMs from one file
- `precursors`: Library precursor information

# Returns
- `protein_peptide_rows`: Collection of unique peptide-protein mappings
"""
function build_protein_peptide_rows_for_file(psms_table, precursors)
    protein_peptide_rows = Set{NamedTuple{(:sequence, :protein_name, :decoy, :entrap_id), Tuple{String, String, Bool, UInt8}}}()
    
    # Get data from PSMs
    passing_precursor_idx = psms_table[:precursor_idx]
    
    # Get other data from precursors
    accession_numbers = getAccessionNumbers(precursors)
    decoys = getIsDecoy(precursors)
    entrap_ids = getEntrapmentGroupId(precursors)
    sequences = getSequence(precursors)
    
    for pid in passing_precursor_idx
        push!(
            protein_peptide_rows, 
            (
                sequence = sequences[pid],
                protein_name = accession_numbers[pid],
                decoy = decoys[pid],
                entrap_id = entrap_ids[pid]
            )
        )
    end
    
    return collect(protein_peptide_rows)
end
"""
    get_protein_groups(passing_psms_paths::Vector{String}, passing_pg_paths::Vector{String},
                      protein_groups_folder::String, temp_folder::String,
                      precursors::LibraryPrecursors; min_peptides=2,
                      protein_q_val_threshold::Float32=0.01f0)

Create and score protein groups from passing PSMs with per-file protein inference.

# Arguments
- `passing_psms_paths`: Paths to PSM files that passed FDR threshold
- `passing_pg_paths`: Output paths for protein group files
- `protein_groups_folder`: Folder to store protein group results
- `temp_folder`: Temporary folder for intermediate files
- `precursors`: Library precursor information
- `min_peptides`: Minimum peptides required for a protein group (default: 2)
- `protein_q_val_threshold`: Q-value threshold for protein groups (default: 0.01)

# Returns
- Nothing (protein inference is performed per-file)

# Process
1. Counts all possible peptides for each protein in library
2. For each file:
   a. Performs protein inference for that file's peptides
   b. Creates protein groups from inferred proteins
   c. Calculates protein group scores and features
3. Performs probit regression analysis if sufficient data
4. Generates QC plots with decision boundaries

# Note
Protein inference is now performed independently for each file, which reduces
memory usage but may result in different protein groupings across files.
"""
function get_protein_groups(
    passing_psms_paths::Vector{String},
    passing_pg_paths::Vector{String},
    protein_groups_folder::String,
    temp_folder::String,
    precursors::LibraryPrecursors;
    min_peptides = 2,
    protein_q_val_threshold::Float32 = 0.01f0,
    max_psms_in_memory::Int64 = 10000000  # Default value if not provided
)
    @info "[PERF] get_protein_groups: Starting" n_files=length(passing_psms_paths) min_peptides=min_peptides
    start_time = time()
    initial_memory = Base.gc_live_bytes() / 1024^2

    pg_count = 0
    
    # Count all possible peptides for each protein in the library
    protein_to_possible_peptides = count_protein_peptides(precursors)
    
    # Perform protein inference and initial scoring
    pg_count = perform_protein_inference(
        passing_psms_paths,
        passing_pg_paths,
        protein_groups_folder,
        precursors,
        protein_to_possible_peptides;
        min_peptides = min_peptides
    )

    # Perform Probit Regression Analysis on protein groups
    qc_folder = joinpath(dirname(temp_folder), "qc_plots")
    !isdir(qc_folder) && mkdir(qc_folder)
    perform_protein_probit_regression(passing_pg_paths, max_psms_in_memory, qc_folder)

    # Calculate global protein scores
    acc_to_max_pg_score = calculate_global_protein_scores(passing_pg_paths)
    
    # Add global scores to PSM files
    add_global_scores_to_psms(passing_psms_paths, acc_to_max_pg_score)
    
    total_elapsed = time() - start_time
    final_memory = Base.gc_live_bytes() / 1024^2
    @info "[PERF] get_protein_groups: Completed" total_elapsed=round(total_elapsed, digits=3) memory_used_MB=round(final_memory-initial_memory, digits=1) total_protein_groups=pg_count
    
    # Return nothing since protein inference is now per-file
    return nothing
end

"""
    perform_probit_analysis_oom(pg_paths::Vector{String}, total_protein_groups::Int, 
                               max_protein_groups_in_memory::Int, qc_folder::String)

Perform out-of-memory probit regression analysis on protein groups.

# Arguments
- `pg_paths`: Vector of paths to protein group Arrow files
- `total_protein_groups`: Total number of protein groups across all files
- `max_protein_groups_in_memory`: Maximum number of protein groups to hold in memory
- `qc_folder`: Folder for QC plots

# Process
1. Calculate sampling ratio for each file
2. Sample protein groups from each file proportionally
3. Fit probit model on sampled data
4. Apply model to all protein groups file by file
5. Calculate and report performance metrics
"""
function perform_probit_analysis_oom(pg_paths::Vector{String}, total_protein_groups::Int, 
                                    max_protein_groups_in_memory::Int, qc_folder::String)
    
    # Calculate sampling ratio
    sampling_ratio = max_protein_groups_in_memory / total_protein_groups
    @info "Sampling ratio: $(round(sampling_ratio * 100, digits=2))%"
    
    # Sample protein groups from each file
    sampled_protein_groups = DataFrame()
    for pg_path in pg_paths
        if isfile(pg_path) && endswith(pg_path, ".arrow")
            table = Arrow.Table(pg_path)
            n_rows = length(table[:protein_name])  # Get row count from any column
            n_sample = ceil(Int, n_rows * sampling_ratio)
            
            if n_sample > 0
                # Sample indices
                sample_indices = sort(sample(1:n_rows, n_sample, replace=false))
                
                # Load sampled rows
                df = DataFrame(Tables.columntable(table))
                append!(sampled_protein_groups, df[sample_indices, :])
            end
        end
    end
    
    @info "Sampled $(nrow(sampled_protein_groups)) protein groups for training"
    
    # Define features to use
    feature_names = [:pg_score, :peptide_coverage, :n_possible_peptides, :log_binom_coeff]
    X = Matrix{Float64}(sampled_protein_groups[:, feature_names])
    y = sampled_protein_groups.target
    
    # Fit probit model on sampled data
    β_fitted, X_mean, X_std = fit_probit_model(X, y)
    
    @info "Fitted probit model coefficients: $β_fitted"
    
    # Now apply the model to all protein groups file by file
    # and collect statistics
    total_targets = 0
    total_decoys = 0
    
    # Process each file
    for pg_path in pg_paths
        if isfile(pg_path) && endswith(pg_path, ".arrow")
            # Load file
            df = DataFrame(Tables.columntable(Arrow.Table(pg_path)))
            
            # Calculate probit scores
            X_file = Matrix{Float64}(df[:, feature_names])
            prob_scores = calculate_probit_scores(X_file, β_fitted, X_mean, X_std)
            
            # Overwrite pg_score with probit scores
            df[!, :pg_score] = Float32.(prob_scores)
            
            # Sort by pg_score and target in descending order
            sort!(df, [:pg_score, :target], rev = [true, true])
            
            # Count statistics
            total_targets += sum(df.target)
            total_decoys += sum(.!df.target)
            
            # Write back the file with probit scores
            writeArrow(pg_path, df)
        end
    end
    
    # Report results
    @info "Out-of-Memory Probit Regression Results:"
    @info "  Total protein groups: $(total_targets + total_decoys) (targets: $total_targets, decoys: $total_decoys)"
    @info "  Training sample size: $(nrow(sampled_protein_groups))"
    @info "  Model features: $(join(feature_names, ", "))"
    @info "  Model coefficients: $(round.(β_fitted, digits=3))"
end

"""
    perform_probit_analysis(all_protein_groups::DataFrame, qc_folder::String, 
                          pg_paths::Union{Vector{String}, Nothing} = nothing)

Perform probit regression analysis on protein groups with comparison to baseline.

# Arguments
- `all_protein_groups::DataFrame`: Protein group data with features
- `qc_folder::String`: Folder for QC plots
- `pg_paths::Union{Vector{String}, Nothing}`: Paths to protein group files for re-processing

# Process
1. Fits probit model with multiple features
2. Calculates performance metrics
3. Compares to pg_score-only baseline
4. Generates decision boundary plots
5. Re-processes individual files with probit scores if pg_paths provided
"""
function perform_probit_analysis(all_protein_groups::DataFrame, qc_folder::String,
                               pg_paths::Union{Vector{String}, Nothing} = nothing)
    n_targets = sum(all_protein_groups.target)
    n_decoys = sum(.!all_protein_groups.target)
    
    # Define features to use
    feature_names = [:pg_score, :peptide_coverage, :n_possible_peptides, :log_binom_coeff]
    X = Matrix{Float64}(all_protein_groups[:, feature_names])
    y = all_protein_groups.target
    
    # Fit probit model
    β_fitted, X_mean, X_std = fit_probit_model(X, y)

    # Report basic model statistics
    @info "Probit Regression completed:"
    @info "  Total protein groups: $(n_targets + n_decoys) (targets: $n_targets, decoys: $n_decoys)"
    @info "  Model features: $(join(feature_names, ", "))"
    @info "  Model coefficients: $(round.(β_fitted, digits=3))"
    
    # Create decision boundary plots
    # TODO: Fix plotting type error
    # plot_probit_decision_boundary(all_protein_groups, β_fitted, X_mean, X_std, feature_names, qc_folder)
    
    # Re-process individual files if paths are provided
    if pg_paths !== nothing
        @info "Re-processing individual protein group files with probit scores"
        for pg_path in pg_paths
            if isfile(pg_path) && endswith(pg_path, ".arrow")
                # Load file
                df = DataFrame(Tables.columntable(Arrow.Table(pg_path)))
                
                # Calculate probit scores
                X_file = Matrix{Float64}(df[:, feature_names])
                prob_scores = calculate_probit_scores(X_file, β_fitted, X_mean, X_std)
                
                # Overwrite pg_score with probit scores
                df[!, :pg_score] = Float32.(prob_scores)
                
                # Sort by pg_score and target in descending order
                sort!(df, [:pg_score, :target], rev = [true, true])
                
                # Write back the file with probit scores
                writeArrow(pg_path, df)
            end
        end
    end
end


"""
    fit_probit_model(X::Matrix{Float64}, y::Vector{Bool})

Fit a probit regression model for protein group classification.

# Arguments
- `X::Matrix{Float64}`: Feature matrix
- `y::Vector{Bool}`: Target labels (true for targets, false for decoys)

# Returns
- `β_fitted`: Fitted coefficients
- `X_mean`: Feature means for standardization
- `X_std`: Feature standard deviations for standardization
"""
function fit_probit_model(X::Matrix{Float64}, y::Vector{Bool})
    # Standardize features
    X_mean = mean(X, dims=1)
    X_std = std(X, dims=1)
    X_std[X_std .== 0] .= 1.0  # Avoid division by zero
    X_standardized = (X .- X_mean) ./ X_std
    
    # Add intercept column
    X_with_intercept = hcat(ones(size(X_standardized, 1)), X_standardized)
    X_df = DataFrame(X_with_intercept, [:intercept; Symbol.("feature_", 1:size(X, 2))])
    
    # Initialize coefficients
    β = zeros(Float64, size(X_with_intercept, 2))
    
    # Create data chunks for parallel processing
    n_chunks = max(1, Threads.nthreads())
    chunk_size = max(1, ceil(Int, length(y) / n_chunks))
    data_chunks = Iterators.partition(1:length(y), chunk_size)
    
    # Fit probit model
    β_fitted = Pioneer.ProbitRegression(β, X_df, y, data_chunks, max_iter=30)
    
    return β_fitted, vec(X_mean), vec(X_std)
end

"""
    calculate_probit_scores(X::Matrix{Float64}, β::Vector{Float64}, X_mean::Vector{Float64}, X_std::Vector{Float64})

Calculate probit probability scores for new data.

# Arguments
- `X::Matrix{Float64}`: Feature matrix
- `β::Vector{Float64}`: Fitted coefficients
- `X_mean::Vector{Float64}`: Feature means from training
- `X_std::Vector{Float64}`: Feature standard deviations from training

# Returns
- `Vector{Float64}`: Probability scores
"""
function calculate_probit_scores(X::Matrix{Float64}, β::Vector{Float64}, X_mean::Vector{Float64}, X_std::Vector{Float64})
    # Standardize using training statistics
    X_standardized = (X .- X_mean') ./ X_std'
    
    # Add intercept
    X_with_intercept = hcat(ones(size(X_standardized, 1)), X_standardized)
    X_df = DataFrame(X_with_intercept, [:intercept; Symbol.("feature_", 1:size(X, 2))])
    
    # Create data chunks
    n_chunks = max(1, Threads.nthreads())
    chunk_size = max(1, ceil(Int, size(X, 1) / n_chunks))
    data_chunks = Iterators.partition(1:size(X, 1), chunk_size)
    
    # Calculate probabilities
    prob_scores = zeros(Float64, size(X, 1))
    Pioneer.ModelPredictProbs!(prob_scores, X_df, β, data_chunks)
    
    return prob_scores
end

"""
    calculate_qvalues_from_scores(scores::Vector{<:Real}, labels::Vector{Bool})

Calculate q-values from scores and labels.

# Arguments
- `scores::Vector{<:Real}`: Scores (higher = more likely to be target)
- `labels::Vector{Bool}`: True labels (true = target, false = decoy)

# Returns
- `Vector{Float64}`: Q-values for each score
"""
function calculate_qvalues_from_scores(scores::Vector{<:Real}, labels::Vector{Bool})
    # Sort by score in descending order
    sort_idx = sortperm(scores, rev=true)
    sorted_labels = labels[sort_idx]
    
    # Calculate q-values
    qvalues = zeros(Float64, length(scores))
    cumulative_targets = 0
    cumulative_decoys = 0
    
    for i in 1:length(sorted_labels)
        if sorted_labels[i]
            cumulative_targets += 1
        else
            cumulative_decoys += 1
        end
        
        # FDR = decoys / targets
        if cumulative_targets == 0
            qvalues[i] = 1.0
        else
            qvalues[i] = cumulative_decoys / cumulative_targets
        end
    end
    
    # Convert FDR to q-value (ensure monotonicity)
    for i in (length(qvalues)-1):-1:1
        qvalues[i] = min(qvalues[i], qvalues[i+1])
    end
    
    # Map q-values back to original order
    unsorted_qvalues = zeros(Float64, length(scores))
    unsorted_qvalues[sort_idx] = qvalues
    
    return unsorted_qvalues
end

"""
    plot_probit_decision_boundary(all_protein_groups::DataFrame, β::Vector{Float64}, 
                                  X_mean::Vector{Float64}, X_std::Vector{Float64},
                                  feature_names::Vector{Symbol}, qc_folder::String)

Create scatter plots showing the probit decision boundary.

# Arguments
- `all_protein_groups::DataFrame`: Protein group data
- `β::Vector{Float64}`: Fitted probit coefficients
- `X_mean::Vector{Float64}`: Feature means from training
- `X_std::Vector{Float64}`: Feature standard deviations from training
- `feature_names::Vector{Symbol}`: Names of features used
- `qc_folder::String`: Output folder for plots
"""
function plot_probit_decision_boundary(all_protein_groups::DataFrame, β::Vector{Float64}, 
                                     X_mean::Vector{Float64}, X_std::Vector{Float64},
                                     feature_names::Vector{Symbol}, qc_folder::String)
    
    # Create output folder if it doesn't exist
    protein_ml_folder = joinpath(qc_folder, "protein_ml_plots")
    !isdir(protein_ml_folder) && mkdir(protein_ml_folder)
    
    # Calculate probit scores for all data
    X = Matrix{Float64}(all_protein_groups[:, feature_names])
    prob_scores = calculate_probit_scores(X, β, X_mean, X_std)
    
    # Separate targets and decoys
    target_mask = all_protein_groups.target
    decoy_mask = .!all_protein_groups.target
    
    # Plot pg_score vs peptide_coverage with decision boundary coloring
    p = scatter(all_protein_groups[target_mask, :peptide_coverage], 
               all_protein_groups[target_mask, :pg_score],
               label="Targets",
               color=:blue,
               alpha=0.3,
               markersize=3,
               xlabel="Peptide Coverage",
               ylabel="PG Score",
               title="Protein Group Classification with Probit Decision Boundary",
               legend=:bottomright,
               yscale=:log10,
               ylims=(1, maximum(all_protein_groups.pg_score) * 1.5))
    
    scatter!(p, all_protein_groups[decoy_mask, :peptide_coverage],
            all_protein_groups[decoy_mask, :pg_score],
            label="Decoys",
            color=:red,
            alpha=0.3,
            markersize=3)
    
    # Add decision boundary contour
    # Create a grid for the contour
    x_range = range(minimum(all_protein_groups.peptide_coverage), 
                    maximum(all_protein_groups.peptide_coverage), length=100)
    y_range = range(log10(1), log10(maximum(all_protein_groups.pg_score) * 1.5), length=100)
    
    # For each point in the grid, calculate the probit score
    z_grid = zeros(length(y_range), length(x_range))
    for (i, x) in enumerate(x_range)
        for (j, y_log) in enumerate(y_range)
            y = 10^y_log
            # Create feature vector for this point
            # Assuming features are [pg_score, peptide_coverage, n_possible_peptides, log_binom_coeff]
            # We need to estimate n_possible_peptides and log_binom_coeff for the grid point
            # For visualization, we'll use median values
            median_n_possible = median(all_protein_groups.n_possible_peptides)
            n_observed = x * median_n_possible  # peptide_coverage * n_possible
            log_binom = if n_observed <= median_n_possible && n_observed >= 0
                lgamma(median_n_possible + 1) - lgamma(n_observed + 1) - lgamma(median_n_possible - n_observed + 1)
            else
                0.0
            end
            
            X_point = reshape([y, x, median_n_possible, log_binom], 1, :)  # Create 1x4 matrix
            prob_score = calculate_probit_scores(X_point, β, X_mean, X_std)[1]
            z_grid[j, i] = prob_score
        end
    end
    
    # Add contour line at 0.5 probability
    contour!(p, x_range, 10 .^ y_range, z_grid', 
            levels=[0.5], 
            color=:black, 
            linewidth=2, 
            linestyle=:dash,
            label="Decision Boundary (p=0.5)")
    
    savefig(p, joinpath(protein_ml_folder, "pg_score_vs_peptide_coverage_with_boundary.png"))
    
    # Create a simpler 2D plot with just pg_score colored by probit probability
    p2 = scatter(all_protein_groups[target_mask, :pg_score],
                prob_scores[target_mask],
                label="Targets",
                color=:blue,
                alpha=0.5,
                markersize=4,
                xlabel="PG Score",
                ylabel="Probit Probability",
                title="Probit Model Predictions",
                legend=:bottomright,
                xscale=:log10,
                xlims=(1, maximum(all_protein_groups.pg_score) * 1.5),
                ylims=(0, 1))
    
    scatter!(p2, all_protein_groups[decoy_mask, :pg_score],
            prob_scores[decoy_mask],
            label="Decoys",
            color=:red,
            alpha=0.5,
            markersize=4)
    
    # Add horizontal line at p=0.5
    hline!(p2, [0.5], color=:black, linestyle=:dash, linewidth=2, label="Decision Boundary")
    
    savefig(p2, joinpath(protein_ml_folder, "probit_probability_vs_pg_score.png"))
    
    @info "Protein ML plots saved to: $protein_ml_folder"
end

"""
    add_protein_inference_col(passing_psms_paths, protein_inference_dict, 
                             precursor_sequences, precursor_is_decoy, precursors_entrap_id)

Add protein inference columns to passing PSM files.

# Arguments
- `passing_psms_paths`: Paths to PSM files
- `protein_inference_dict`: Protein inference results
- `precursor_sequences`: Peptide sequences for each precursor
- `precursor_is_decoy`: Decoy status for each precursor
- `precursors_entrap_id`: Entrapment group IDs

# Added columns
- `use_for_protein_quant`: Whether peptide should be used for quantification
- `inferred_protein_group`: Inferred protein group name
"""
function add_protein_inference_col(
    passing_psms_paths::Vector{String},
    protein_inference_dict::Dictionary{@NamedTuple{peptide::String, decoy::Bool, entrap_id::UInt8}, @NamedTuple{protein_name::String, decoy::Bool, entrap_id::UInt8, retain::Bool}},
    precursor_sequences::AbstractVector{S},
    precursor_is_decoy::AbstractVector{Bool},
    precursors_entrap_id::AbstractVector{UInt8}
) where {S<:AbstractString}

    for (ms_file_idx, file_path) in enumerate(passing_psms_paths)
        passing_psms = DataFrame(Tables.columntable(Arrow.Table(file_path)))
        precursor_idx = passing_psms[!,:precursor_idx]
        inferred_protein_group = Vector{String}(undef, length(precursor_idx))
        use_for_protein_quant = zeros(Bool, length(precursor_idx))
        for (i, pid) in enumerate(precursor_idx)
            #protein_group, decoy, entrap_id, use_for_inference = protein_inference_dict[(peptide = precursor_sequences[pid], decoy = precursor_is_decoy[pid], entrap_id = precursors_entrap_id[pid])]
            inferred_prot = protein_inference_dict[(peptide = precursor_sequences[pid], decoy = precursor_is_decoy[pid], entrap_id = precursors_entrap_id[pid])]
            inferred_protein_group[i] = inferred_prot.protein_name::String
            use_for_protein_quant[i] = inferred_prot.retain::Bool
        end
        passing_psms[!,:use_for_protein_quant] = use_for_protein_quant
        passing_psms[!,:inferred_protein_group] = inferred_protein_group
        writeArrow(
            file_path,
            passing_psms            
        )
    end
end

function sort_protein_tables(
    protein_groups_paths::Vector{String},
    merged_pgs_path::String,
    prob_col::Symbol,
)

    #Remove if present 
    if isfile(merged_pgs_path)
        safeRm(merged_pgs_path)
    end
    #file_paths = [fpath for fpath in readdir(quant_psms_folder,join=true) if endswith(fpath,".arrow")]
    #Sort and filter each psm table 
    for fpath in protein_groups_paths
        pgs_table = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        #Sort in descending order of probability
        sort!(pgs_table, prob_col, rev = true, alg=QuickSort)
        #write back
        writeArrow(fpath,  pgs_table)
    end
    return nothing
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
    @info "[PERF] merge_sorted_protein_groups: Starting" input_dir=input_dir batch_size=N
    start_time = time()
    initial_memory = Base.gc_live_bytes() / 1024^2
    total_batches = 0 

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
    # Found files
    
    #Keep track of which tables have 
    load_start = time()
    tables = [Arrow.Table(path) for path in input_paths]
    table_sizes = [length(table[1]) for table in tables]
    total_rows = sum(table_sizes)
    load_time = time() - load_start
    # Tables loaded
    
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
    rows_written = 0
    merge_start = time()
    
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
                    safeRm(output_path)
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
            rows_written += (i - 1)
            total_batches += 1
            
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
    if n_writes > 0
        Arrow.append(
            output_path,
            pg_batch[range(1, max(1, i - 1)),:]
        )
    else
        writeArrow(
            output_path,
            pg_batch[range(1, max(1, i - 1)),:]
        )
    end
    
    rows_written += max(1, i - 1)
    elapsed = time() - start_time
    final_memory = Base.gc_live_bytes() / 1024^2
    @info "[PERF] merge_sorted_protein_groups: Completed" elapsed=round(elapsed, digits=3) n_files=length(input_paths) total_batches=total_batches total_rows=rows_written rows_per_sec=round(rows_written/elapsed, digits=1) memory_used_MB=round(final_memory-initial_memory, digits=1) avg_rows_per_file=round(total_rows/max(length(input_paths),1), digits=1)
    
    return nothing
end

function get_proteins_passing_qval(
    passing_proteins_folder::String,
    global_qval_interp::Interpolations.Extrapolation,
    experiment_wide_qval_interp::Interpolations.Extrapolation,
    global_prob_col::Symbol,
    experiment_wide_prob_col::Symbol,
    global_qval_col::Symbol,
    experiment_wide_qval_col::Symbol,
    q_val_threshold::Float32)

    #Get all .arrow files in the input 
    input_paths = [path for path in readdir(passing_proteins_folder, join=true) if endswith(path, ".arrow")]

    for file_path in input_paths
        # Read the Arrow table
        passing_proteins = DataFrame(Tables.columntable(Arrow.Table(file_path)))
        passing_proteins[!,global_qval_col] = global_qval_interp.(passing_proteins[!,global_prob_col])
        passing_proteins[!,experiment_wide_qval_col] = experiment_wide_qval_interp.(passing_proteins[!,experiment_wide_prob_col])

        # Sample the rows and convert to DataFrame
        filter!(x->(x[global_qval_col]<=q_val_threshold)&(x[experiment_wide_qval_col]<=q_val_threshold), passing_proteins)
        # Append to the result DataFrame
        writeArrow(            joinpath(passing_proteins_folder, basename(file_path)),
        passing_proteins)
    end
    
    return
end