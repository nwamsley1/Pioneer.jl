"""
    get_best_precursors_accross_runs_optimized(psms_paths::Vector{String}, 
                                              prec_mzs::AbstractVector{Float32},
                                              rt_irt::Dict{Int64, RtConversionModel};
                                              max_q_val::Float32=0.01f0,
                                              max_concurrent_files::Int=4) 
    -> Dictionary{UInt32, NamedTuple}

Optimized version with controlled parallelism to avoid memory bandwidth saturation.
Processes files in small batches with thread-local accumulation before merging.

# Key Optimizations
- Controlled concurrent Arrow file reading (max_concurrent_files parameter)
- Thread-local dictionary accumulation to avoid synchronization
- Efficient dictionary merging at the end
- Memory bandwidth aware - limits concurrent I/O operations

# Arguments
- `psms_paths`: Paths to PSM files from first pass search
- `prec_mzs`: Vector of precursor m/z values
- `rt_irt`: Dictionary mapping file indices to RT-iRT conversion models
- `max_q_val`: Maximum q-value threshold for considering PSMs
- `max_concurrent_files`: Maximum concurrent Arrow file reads (default: 4)

# Returns
Dictionary mapping precursor indices to NamedTuple containing best precursor statistics
"""
function get_best_precursors_accross_runs_optimized(
                         psms_paths::Vector{String},
                         prec_mzs::AbstractVector{Float32},
                         rt_irt::Dict{Int64, RtConversionModel};
                         max_q_val::Float32 = 0.01f0,
                         max_concurrent_files::Int = 4
                         )

    # Thread-safe version of readPSMs! that works with local dictionaries
    function readPSMs_local!(
        prec_to_best_prob::Dictionary{UInt32, @NamedTuple{ best_prob::Float32, 
                                                    best_ms_file_idx::UInt32,
                                                    best_scan_idx::UInt32,
                                                    best_irt::Float32, 
                                                    mean_irt::Union{Missing, Float32}, 
                                                    var_irt::Union{Missing, Float32}, 
                                                    n::Union{Missing, UInt16}, 
                                                    mz::Float32}},
        precursor_idxs::AbstractVector{UInt32},
        q_values::AbstractVector{Float16},
        probs::AbstractVector{Float32},
        rts::AbstractVector{Float32},
        scan_idxs::AbstractVector{UInt32},
        ms_file_idxs::AbstractVector{UInt32},
        rt_irt::SplineRtConversionModel,
        max_q_val::Float32)

        for row in eachindex(precursor_idxs)
            # Extract current PSM information
            q_value = q_values[row]
            precursor_idx = precursor_idxs[row]
            prob = probs[row]
            irt = rt_irt(rts[row])
            scan_idx = UInt32(scan_idxs[row])
            ms_file_idx = UInt32(ms_file_idxs[row])

            # Initialize running statistics
            passed_q_val = (q_value <= max_q_val)
            n = passed_q_val ? one(UInt16) : zero(UInt16)
            mean_irt = passed_q_val ? irt : zero(Float32)
            var_irt = zero(Float32)
            mz = prec_mzs[precursor_idx]

            # Has the precursor been encountered in a previous raw file?
            if haskey(prec_to_best_prob, precursor_idx)
                # Update existing precursor entry
                best_prob, best_ms_file_idx, best_scan_idx, best_irt, old_mean_irt, var_irt, old_n, mz = prec_to_best_prob[precursor_idx] 
                
                # Update best match if current is better
                if (best_prob < prob)
                    best_prob = prob
                    best_irt = irt
                    best_scan_idx = scan_idx
                    best_ms_file_idx = ms_file_idx
                end

                # Update running statistics
                mean_irt += old_mean_irt
                n += old_n 
                prec_to_best_prob[precursor_idx] = (
                                                best_prob = prob, 
                                                best_ms_file_idx = best_ms_file_idx,
                                                best_scan_idx = best_scan_idx,
                                                best_irt = irt,
                                                mean_irt = mean_irt, 
                                                var_irt = var_irt,
                                                n = n,
                                                mz = mz)
            else
                # Create new precursor entry
                val = (best_prob = prob,
                        best_ms_file_idx = ms_file_idx,
                        best_scan_idx = scan_idx, 
                        best_irt = irt,
                        mean_irt = mean_irt, 
                        var_irt = var_irt,
                        n = n,
                        mz = mz)
                insert!(prec_to_best_prob, precursor_idx, val)
            end
        end
    end

    # Merge two dictionaries, keeping the best values for each precursor
    function merge_dictionaries!(
        target::Dictionary{UInt32, @NamedTuple{ best_prob::Float32, 
                                                best_ms_file_idx::UInt32,
                                                best_scan_idx::UInt32,
                                                best_irt::Float32, 
                                                mean_irt::Union{Missing, Float32}, 
                                                var_irt::Union{Missing, Float32}, 
                                                n::Union{Missing, UInt16}, 
                                                mz::Float32}},
        source::Dictionary{UInt32, @NamedTuple{ best_prob::Float32, 
                                                best_ms_file_idx::UInt32,
                                                best_scan_idx::UInt32,
                                                best_irt::Float32, 
                                                mean_irt::Union{Missing, Float32}, 
                                                var_irt::Union{Missing, Float32}, 
                                                n::Union{Missing, UInt16}, 
                                                mz::Float32}})

        for (precursor_idx, source_val) in pairs(source)
            if haskey(target, precursor_idx)
                # Merge the values, keeping best probability
                target_val = target[precursor_idx]
                
                best_prob = max(target_val.best_prob, source_val.best_prob)
                best_ms_file_idx = target_val.best_prob >= source_val.best_prob ? target_val.best_ms_file_idx : source_val.best_ms_file_idx
                best_scan_idx = target_val.best_prob >= source_val.best_prob ? target_val.best_scan_idx : source_val.best_scan_idx
                best_irt = target_val.best_prob >= source_val.best_prob ? target_val.best_irt : source_val.best_irt
                
                # Combine statistics
                mean_irt = target_val.mean_irt + source_val.mean_irt
                n = target_val.n + source_val.n
                var_irt = target_val.var_irt + source_val.var_irt
                mz = target_val.mz  # Should be the same
                
                target[precursor_idx] = (
                    best_prob = best_prob,
                    best_ms_file_idx = best_ms_file_idx,
                    best_scan_idx = best_scan_idx,
                    best_irt = best_irt,
                    mean_irt = mean_irt,
                    var_irt = var_irt,
                    n = n,
                    mz = mz
                )
            else
                # Just copy the value
                target[precursor_idx] = source_val
            end
        end
    end

    # Process files in batches to control memory bandwidth
    n_files = length(psms_paths)
    batch_size = min(max_concurrent_files, n_files)
    n_precursors_vec = Vector{UInt64}(undef, n_files)
    
    # Initialize global dictionary
    prec_to_best_prob = Dictionary{UInt32, @NamedTuple{ best_prob::Float32, 
                                                        best_ms_file_idx::UInt32,
                                                        best_scan_idx::UInt32,
                                                        best_irt::Float32, 
                                                        mean_irt::Union{Missing, Float32}, 
                                                        var_irt::Union{Missing, Float32}, 
                                                        n::Union{Missing, UInt16}, 
                                                        mz::Float32}}()

    # First pass: collect best matches and mean iRT in batches
    for batch_start in 1:batch_size:n_files
        batch_end = min(batch_start + batch_size - 1, n_files)
        batch_indices = batch_start:batch_end
        
        # Process batch in parallel (limited concurrency)
        tasks = map(batch_indices) do key
            Threads.@spawn begin
                # Thread-local dictionary
                local_dict = Dictionary{UInt32, @NamedTuple{ best_prob::Float32, 
                                                            best_ms_file_idx::UInt32,
                                                            best_scan_idx::UInt32,
                                                            best_irt::Float32, 
                                                            mean_irt::Union{Missing, Float32}, 
                                                            var_irt::Union{Missing, Float32}, 
                                                            n::Union{Missing, UInt16}, 
                                                            mz::Float32}}()
                
                psms_path = psms_paths[key]
                psms = Arrow.Table(psms_path)
                
                # Store n_precursors
                n_precursors = length(psms[:precursor_idx])
                
                # Process PSMs into local dictionary
                readPSMs_local!(
                    local_dict,
                    psms[:precursor_idx],
                    psms[:q_value],
                    psms[:prob],
                    psms[:rt],
                    psms[:scan_idx],
                    psms[:ms_file_idx],
                    rt_irt[key],
                    max_q_val
                )
                
                return (key, n_precursors, local_dict)
            end
        end
        
        # Collect results and merge into global dictionary
        for task in tasks
            key, n_precursors, local_dict = fetch(task)
            n_precursors_vec[key] = n_precursors
            merge_dictionaries!(prec_to_best_prob, local_dict)
        end
    end

    max_precursors = maximum(n_precursors_vec)
    # Filter to top N precursors by probability
    sort!(prec_to_best_prob, by = x->x[:best_prob], alg=PartialQuickSort(1:max_precursors), rev = true);
    N = 0
    for key in collect(keys(prec_to_best_prob))
        N += 1
        if N > max_precursors
            delete!(prec_to_best_prob, key)
        end
    end

    # For the variance calculation, we need a different approach since we can't 
    # easily parallelize variance calculations. We'll use the original sequential approach.
    for (key, psms_path) in enumerate(psms_paths)
        psms = Arrow.Table(psms_path)
        
        for row in eachindex(psms[:precursor_idx])
            q_value = psms[:q_value][row]
            if q_value > max_q_val
                continue
            end
            
            precursor_idx = psms[:precursor_idx][row]
            irt = rt_irt[key](psms[:rt][row])
            
            if haskey(prec_to_best_prob, precursor_idx)
                # Update variance calculation
                best_prob, best_ms_file_idx, best_scan_idx, best_irt, mean_irt, var_irt, n, mz = prec_to_best_prob[precursor_idx] 
                var_irt += (irt - mean_irt/n)^2
                prec_to_best_prob[precursor_idx] = (
                    best_prob = best_prob, 
                    best_ms_file_idx = best_ms_file_idx,
                    best_scan_idx = best_scan_idx,
                    best_irt = best_irt,
                    mean_irt = mean_irt, 
                    var_irt = var_irt,
                    n = n,
                    mz = mz)
            end
        end
    end
    
    return prec_to_best_prob
end