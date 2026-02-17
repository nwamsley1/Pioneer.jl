# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
    get_best_precursors_accross_runs(psms_paths, prec_mzs, rt_irt, prec_is_decoy;
                                    max_q_val=0.01f0, min_pep=Float16(0.9),
                                    fdr_scale_factor=1.0f0)
    -> Dictionary{UInt32, NamedTuple}

Identify and collect best precursor matches across multiple runs for retention time calibration.
Pools PSMs from all runs and computes a **global PEP** via isotonic regression (`get_PEP!`),
giving a more stable score→PEP mapping than per-run estimates. A precursor passes if its
minimum global PEP across all runs is ≤ the threshold.

# Arguments
- `psms_paths`: Paths to PSM files from first pass search
- `prec_mzs`: Vector of precursor m/z values
- `rt_irt`: Dictionary mapping file indices to RT-iRT conversion models
- `prec_is_decoy`: Boolean vector indicating decoy status per precursor index
- `max_q_val`: Maximum q-value threshold for iRT statistics
- `min_pep`: Global PEP threshold — precursors with min global PEP ≤ this pass to second pass
- `fdr_scale_factor`: Scale factor for library target/decoy ratio correction

# Returns
Dictionary mapping precursor indices to NamedTuple containing:
- `best_prob`: Highest probability score
- `best_ms_file_idx`: File index with best match
- `best_scan_idx`: Scan index of best match
- `best_irt`: iRT value of best match
- `mean_irt`: Mean iRT across qualifying matches
- `var_irt`: Variance in iRT across qualifying matches
- `n`: Number of qualifying matches
- `mz`: Precursor m/z value

# Process
1. First pass: Collects best matches, mean iRT, and pools (precursor_idx, prob) for global PEP
2. Computes global PEP via isotonic regression on pooled data
3. Filters to precursors whose minimum global PEP ≤ min_pep
4. Second pass: Calculates iRT variance for remaining precursors
"""
function get_best_precursors_accross_runs(
                         psms_paths::Vector{String},
                         prec_mzs::AbstractVector{Float32},
                         rt_irt::Dict{Int64, RtConversionModel},
                         prec_is_decoy::AbstractVector{Bool};
                         max_q_val::Float32 = 0.01f0,
                         min_pep::Float16 = Float16(0.9),
                         fdr_scale_factor::Float32 = 1.0f0
                         )

    function readPSMs!(
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
        rt_irt::RtConversionModel,
        max_q_val::Float32)

        for row in eachindex(precursor_idxs)
            # Extract current PSM information
            q_value = q_values[row]
            precursor_idx = precursor_idxs[row]
            prob = probs[row]
            irt =  rt_irt(rts[row])
            scan_idx = UInt32(scan_idxs[row])
            ms_file_idx = UInt32(ms_file_idxs[row])

            # Initialize running statistics
            passed_q_val = (q_value <= max_q_val)
            n = passed_q_val ? one(UInt16) : zero(UInt16)
            mean_irt = passed_q_val ? irt : zero(Float32)
            var_irt = zero(Float32)
            mz = prec_mzs[precursor_idx]

            #Has the precursor been encountered in a previous raw file?
            #Keep a running mean irt for instances below q-val threshold
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
                                                best_prob = best_prob, 
                                                best_ms_file_idx = best_ms_file_idx,
                                                best_scan_idx = best_scan_idx,
                                                best_irt = best_irt,
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
    function getVariance!(
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
        rts::AbstractVector{Float32},
        rt_irt::RtConversionModel,
        max_q_val::Float32)
        for row in eachindex(precursor_idxs)
            # Skip PSMs that don't pass q-value threshold
            q_value = q_values[row]

            # Get precursor info
            precursor_idx = precursor_idxs[row]
            irt =  rt_irt(rts[row])

            if q_value > max_q_val
                continue
            end
            if haskey(prec_to_best_prob, precursor_idx)
                # Update variance calculation
                best_prob, best_ms_file_idx, best_scan_idx, best_irt, mean_irt, var_irt, n, mz = prec_to_best_prob[precursor_idx] 
                var_irt += (irt - mean_irt/n)^2
                prec_to_best_prob[precursor_idx] = (
                    best_prob = best_prob, 
                    best_ms_file_idx= best_ms_file_idx,
                    best_scan_idx = best_scan_idx,
                    best_irt = best_irt,
                    mean_irt = mean_irt, 
                    var_irt = var_irt,
                    n = n,
                    mz = mz)

            end
        end
    end
    # Initialize dictionary to store best precursor matches
    prec_to_best_prob = Dictionary{UInt32, @NamedTuple{ best_prob::Float32, 
                                                        best_ms_file_idx::UInt32,
                                                        best_scan_idx::UInt32,
                                                        best_irt::Float32, 
                                                        mean_irt::Union{Missing, Float32}, 
                                                        var_irt::Union{Missing, Float32}, 
                                                        n::Union{Missing, UInt16}, 
                                                        mz::Float32}}()

    # First pass: collect best matches, mean iRT, and pool (precursor_idx, prob) for global PEP

    pooled_precursor_idxs = Vector{UInt32}()
    pooled_probs = Vector{Float32}()
    n_valid_files = 0
    for psms_path in psms_paths #For each data frame
        psms = Arrow.Table(psms_path)

        # Get the original file index from the PSM data
        if isempty(psms[:ms_file_idx])
            continue  # Skip empty files
        end
        file_idx = first(psms[:ms_file_idx])  # All PSMs in a file should have the same ms_file_idx

        # Check if RT model exists for this file
        if !haskey(rt_irt, file_idx)
            @warn "No RT model found for file index $file_idx, skipping"
            continue
        end

        n_valid_files += 1

        # Per-file PSM count diagnostics
        pep_col = psms[:PEP]
        target_col = [!prec_is_decoy[pid] for pid in psms[:precursor_idx]]
        n_total = length(pep_col)
        fname = basename(psms_path)
        parts = String[]
        for thresh in (Float16(0.1), Float16(0.5), Float16(0.75), Float16(0.9))
            nt = count(i -> pep_col[i] <= thresh && target_col[i], eachindex(pep_col))
            nd = count(i -> pep_col[i] <= thresh && !target_col[i], eachindex(pep_col))
            push!(parts, "≤$thresh: T=$nt D=$nd")
        end
        @info "File $fname: $n_total total PSMs | if filtered — " * join(parts, " | ") * "\n"

        #One row for each precursor
        readPSMs!(
            prec_to_best_prob,
            psms[:precursor_idx],
            psms[:q_value],
            psms[:prob],
            psms[:rt],
            psms[:scan_idx],
            psms[:ms_file_idx],
            rt_irt[file_idx],
            max_q_val
        )

        # Pool (precursor_idx, prob) for global PEP computation
        prec_col = psms[:precursor_idx]
        prob_col = psms[:prob]
        for i in eachindex(prec_col)
            push!(pooled_precursor_idxs, prec_col[i])
            push!(pooled_probs, prob_col[i])
        end
    end

    # Handle case where no valid files were processed
    if n_valid_files == 0
        @warn "No valid PSM files found for cross-run analysis"
        return prec_to_best_prob
    end

    # Compute global PEP via isotonic regression on pooled data
    n_pooled = length(pooled_probs)

    # Derive target labels from prec_is_decoy
    pooled_targets = Vector{Bool}(undef, n_pooled)
    for i in eachindex(pooled_precursor_idxs)
        pooled_targets[i] = !prec_is_decoy[pooled_precursor_idxs[i]]
    end

    global_peps = Vector{Float32}(undef, n_pooled)
    get_PEP!(pooled_probs, pooled_targets, global_peps;
             doSort=true, fdr_scale_factor=fdr_scale_factor)

    # Global q-values for diagnostics
    global_qvals = Vector{Float32}(undef, n_pooled)
    get_qvalues!(pooled_probs, pooled_targets, global_qvals;
                 doSort=true, fdr_scale_factor=fdr_scale_factor)

    n_targets = count(pooled_targets)
    n_decoys = n_pooled - n_targets
    n_1pct = count(<=(0.01f0), global_qvals)
    n_5pct = count(<=(0.05f0), global_qvals)
    @info "Global PEP: pooled $n_pooled PSMs ($n_targets targets, $n_decoys decoys) from $n_valid_files files\n"
    @info "  PSMs at 1% global FDR: $n_1pct, at 5% global FDR: $n_5pct\n"

    # For each precursor, keep its minimum global PEP across all runs
    prec_min_global_pep = Dictionary{UInt32, Float32}()
    for i in eachindex(pooled_precursor_idxs)
        pid = pooled_precursor_idxs[i]
        gpep = global_peps[i]
        if haskey(prec_min_global_pep, pid)
            if gpep < prec_min_global_pep[pid]
                prec_min_global_pep[pid] = gpep
            end
        else
            insert!(prec_min_global_pep, pid, gpep)
        end
    end

    # Diagnostics: pre-filter counts
    n_pre_filter = length(prec_to_best_prob)
    n_pre_targets = count(pid -> !prec_is_decoy[pid], keys(prec_to_best_prob))
    n_pre_decoys = n_pre_filter - n_pre_targets
    @info "Pre-filter pool: $n_pre_filter precursors ($n_pre_targets targets, $n_pre_decoys decoys)\n"

    # Global unique precursors by min global PEP threshold
    parts = String[]
    for thresh in (0.1f0, 0.5f0, 0.75f0, 0.9f0)
        nt = count(pid -> prec_min_global_pep[pid] <= thresh && !prec_is_decoy[pid], keys(prec_min_global_pep))
        nd = count(pid -> prec_min_global_pep[pid] <= thresh && prec_is_decoy[pid], keys(prec_min_global_pep))
        push!(parts, "≤$thresh: T=$nt D=$nd")
    end
    @info "Global unique precursors by min global PEP — " * join(parts, " | ") * "\n"

    # Filter: keep only precursors with min global PEP ≤ threshold
    for key in collect(keys(prec_to_best_prob))
        if !haskey(prec_min_global_pep, key) || prec_min_global_pep[key] > Float32(min_pep)
            delete!(prec_to_best_prob, key)
        end
    end

    n_post_filter = length(prec_to_best_prob)
    n_post_targets = count(pid -> !prec_is_decoy[pid], keys(prec_to_best_prob))
    n_post_decoys = n_post_filter - n_post_targets
    n_removed = n_pre_filter - n_post_filter
    n_removed_targets = n_pre_targets - n_post_targets
    n_removed_decoys = n_pre_decoys - n_post_decoys
    @info "Global PEP ≤ $(min_pep) filter: kept $n_post_filter ($n_post_targets targets, $n_post_decoys decoys), " *
          "removed $n_removed ($n_removed_targets targets, $n_removed_decoys decoys)\n"

    # Second pass: calculate iRT variance for remaining precursors
    for psms_path in psms_paths #For each data frame 
        psms = Arrow.Table(psms_path)
        
        # Get the original file index from the PSM data
        if isempty(psms[:ms_file_idx])
            continue  # Skip empty files
        end
        file_idx = first(psms[:ms_file_idx])  # All PSMs in a file should have the same ms_file_idx
        
        # Check if RT model exists for this file
        if !haskey(rt_irt, file_idx)
            continue  # Skip files without RT models (already warned in first pass)
        end
        
        #One row for each precursor 
        getVariance!(
            prec_to_best_prob,
            psms[:precursor_idx],
            psms[:q_value],
            psms[:rt],
            rt_irt[file_idx],
            max_q_val
        )
    end 

    
    return prec_to_best_prob #[(prob, idx) for (idx, prob) in sort(collect(top_probs), rev=true)]
end
