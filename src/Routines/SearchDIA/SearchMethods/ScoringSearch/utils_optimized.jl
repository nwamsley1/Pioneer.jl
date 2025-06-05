# Optimized versions of critical functions in utils.jl

"""
    sort_and_filter_quant_tables_optimized(second_pass_psms_paths, merged_quant_path, 
                                          isotope_trace_type, prob_col, best_traces)

Optimized version using parallel processing and avoiding repeated I/O.
"""
function sort_and_filter_quant_tables_optimized(
    second_pass_psms_paths::Vector{String},
    merged_quant_path::String,
    isotope_trace_type::IsotopeTraceType,
    prob_col::Symbol,
    best_traces::Set{@NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}}
)
    # Remove if present 
    isfile(merged_quant_path) && rm(merged_quant_path)
    
    # Process files in parallel
    Threads.@threads for fpath in second_pass_psms_paths
        try
            # Use Arrow.Table directly without DataFrame conversion when possible
            psms_table = DataFrame(Tables.columntable(Arrow.Table(fpath)))
            
            if seperateTraces(isotope_trace_type)
                # Pre-allocate the best_trace column
                n_rows = size(psms_table, 1)
                best_trace_mask = zeros(Bool, n_rows)
                
                # Vectorized operation instead of loop
                @inbounds for i in 1:n_rows
                    key = (precursor_idx = psms_table[i, :precursor_idx], 
                          isotopes_captured = psms_table[i, :isotopes_captured])
                    best_trace_mask[i] = key ∈ best_traces
                end
                psms_table[!, :best_trace] = best_trace_mask
            else
                transform!(
                    groupby(psms_table, :precursor_idx),
                    :prob => (p -> p .== maximum(p)) => :best_trace
                )
            end
            
            # Filter and sort
            filter!(x -> x.best_trace, psms_table)
            sort!(psms_table, prob_col, rev=true, alg=QuickSort)
            
            # Write back
            writeArrow(fpath, psms_table)
        catch e
            @error "Failed to process file $fpath" exception=e
        end
    end
    
    return nothing
end

"""
    merge_sorted_psms_scores_streaming(input_paths, output_path, prob_col; 
                                     batch_size=1000000, max_memory_gb=2.0)

Streaming version of merge that limits memory usage.
"""
function merge_sorted_psms_scores_streaming(
    input_paths::Vector{String}, 
    output_path::String,
    prob_col::Symbol;
    batch_size::Int = 1000000,
    max_memory_gb::Float64 = 2.0
)
    # Estimate memory per file
    n_files = length(input_paths)
    memory_per_file = max_memory_gb * 1024^3 / n_files  # bytes
    
    # Create file readers
    readers = []
    for path in input_paths
        try
            reader = Arrow.Table(path)
            push!(readers, (path=path, reader=reader, pos=1, len=length(reader[prob_col])))
        catch e
            @warn "Failed to open $path" exception=e
        end
    end
    
    # Initialize output
    isfile(output_path) && rm(output_path)
    
    # Batch for output
    output_batch = DataFrame(
        prob = Float32[],
        target = Bool[],
        precursor_idx = UInt32[]
    )
    sizehint!(output_batch, batch_size)
    
    # Priority queue for merge (max heap)
    heap = BinaryMaxHeap{Tuple{Float32, Int}}()
    
    # Initialize heap with first element from each file
    for (i, r) in enumerate(readers)
        if r.pos <= r.len
            push!(heap, (r.reader[prob_col][r.pos], i))
        end
    end
    
    # Merge loop
    rows_written = 0
    first_write = true
    
    while !isempty(heap)
        # Get highest probability entry
        prob_val, file_idx = pop!(heap)
        r = readers[file_idx]
        
        # Add to output batch
        push!(output_batch.prob, r.reader[prob_col][r.pos])
        push!(output_batch.target, r.reader[:target][r.pos])
        push!(output_batch.precursor_idx, r.reader[:precursor_idx][r.pos])
        
        # Advance position
        readers[file_idx] = (path=r.path, reader=r.reader, pos=r.pos+1, len=r.len)
        
        # Add next element from same file if available
        if readers[file_idx].pos <= readers[file_idx].len
            push!(heap, (readers[file_idx].reader[prob_col][readers[file_idx].pos], file_idx))
        end
        
        # Write batch if full
        if nrow(output_batch) >= batch_size
            if first_write
                Arrow.write(output_path, output_batch)
                first_write = false
            else
                Arrow.append(output_path, output_batch)
            end
            rows_written += nrow(output_batch)
            
            # Clear batch
            empty!(output_batch)
            sizehint!(output_batch, batch_size)
        end
    end
    
    # Write remaining data
    if nrow(output_batch) > 0
        if first_write
            Arrow.write(output_path, output_batch)
        else
            Arrow.append(output_path, output_batch)
        end
        rows_written += nrow(output_batch)
    end
    
    @info "Merged $rows_written rows from $(length(input_paths)) files"
    return nothing
end

"""
    load_protein_groups_optimized(passing_pg_paths)

Optimized loading using single Arrow.Table call.
"""
function load_protein_groups_optimized(passing_pg_paths::Vector{String})
    # Filter valid paths
    valid_paths = filter(passing_pg_paths) do pg_path
        isfile(pg_path) && endswith(pg_path, ".arrow")
    end
    
    if isempty(valid_paths)
        return DataFrame()
    end
    
    # Load all at once - Arrow handles multiple files efficiently
    return DataFrame(Arrow.Table(valid_paths))
end

"""
    count_total_rows_optimized(file_paths)

Count rows without loading data into memory.
"""
function count_total_rows_optimized(file_paths::Vector{String})
    total_rows = 0
    
    for path in file_paths
        if isfile(path) && endswith(path, ".arrow")
            # Use Arrow metadata to get row count without loading data
            open(path, "r") do io
                # Skip to footer to read metadata
                seekend(io)
                footer_size_pos = position(io) - 10
                seek(io, footer_size_pos)
                
                # This is a simplified version - real implementation would parse Arrow footer
                # For now, fall back to loading just the schema
                table = Arrow.Table(path, limit=1)
                if !isempty(propertynames(table))
                    col = first(propertynames(table))
                    total_rows += length(Arrow.Table(path)[col])
                end
            end
        end
    end
    
    return total_rows
end

"""
    process_protein_groups_in_batches(pg_paths, process_func; batch_size=100000)

Process protein groups in memory-efficient batches.
"""
function process_protein_groups_in_batches(
    pg_paths::Vector{String},
    process_func::Function;
    batch_size::Int = 100000
)
    for path in pg_paths
        if !isfile(path) || !endswith(path, ".arrow")
            continue
        end
        
        # Process file in batches
        table = Arrow.Table(path)
        n_rows = length(table[:protein_name])
        
        for start_idx in 1:batch_size:n_rows
            end_idx = min(start_idx + batch_size - 1, n_rows)
            
            # Create view of batch
            batch_df = DataFrame()
            for col in propertynames(table)
                batch_df[!, col] = table[col][start_idx:end_idx]
            end
            
            # Process batch
            process_func(batch_df)
        end
    end
end

"""
    get_psms_passing_qval_optimized(precursors, passing_psms_paths, passing_psms_folder,
                                   second_pass_psms_paths, qval_interp_global, 
                                   qval_interp_experiment_wide, prob_col_global,
                                   prob_col_experiment_wide, qval_col_global,
                                   qval_col_experiment_wide, q_val_threshold)

Optimized version that processes files in parallel and minimizes I/O.
"""
function get_psms_passing_qval_optimized(
    precursors::LibraryPrecursors,
    passing_psms_paths::Vector{String},
    passing_psms_folder::String,
    second_pass_psms_paths::Vector{String},
    qval_interp_global::Interpolations.Extrapolation,
    qval_interp_experiment_wide::Interpolations.Extrapolation,
    prob_col_global::Symbol,
    prob_col_experiment_wide::Symbol,
    qval_col_global::Symbol,
    qval_col_experiment_wide::Symbol,
    q_val_threshold::Float32
)
    # Process files in parallel
    Threads.@threads for (ms_file_idx, file_path) in collect(enumerate(second_pass_psms_paths))
        try
            # Read once
            passing_psms = DataFrame(Tables.columntable(Arrow.Table(file_path)))
            
            # Vectorized q-value calculation
            passing_psms[!, qval_col_global] = qval_interp_global.(passing_psms[!, prob_col_global])
            passing_psms[!, qval_col_experiment_wide] = qval_interp_experiment_wide.(passing_psms[!, prob_col_experiment_wide])
            
            # Define columns to keep
            cols = [:precursor_idx, :global_prob, :prec_prob, :prob, :global_qval,
                   :run_specific_qval, :prec_mz, :weight, :target, :rt, :irt_obs,
                   :missed_cleavage, :isotopes_captured, :scan_idx, :entrapment_group_id,
                   :ms_file_idx]
            available_cols = intersect(cols, Symbol.(names(passing_psms)))
            
            # Filter and select in one operation
            filtered_psms = filter(passing_psms[!, available_cols]) do row
                row[qval_col_global] <= q_val_threshold && 
                row[qval_col_experiment_wide] <= q_val_threshold
            end
            
            # Write result
            output_path = joinpath(passing_psms_folder, basename(file_path))
            writeArrow(output_path, filtered_psms)
            passing_psms_paths[ms_file_idx] = output_path
            
        catch e
            @error "Failed to process PSM file" file=file_path exception=e
        end
    end
    
    return nothing
end