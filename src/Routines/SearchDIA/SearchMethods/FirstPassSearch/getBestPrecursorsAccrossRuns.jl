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
                                    max_q_val=0.01f0, min_pep=Float16(0.5))
    -> Dictionary{UInt32, NamedTuple}

Identify and collect best precursor matches across multiple runs for retention time calibration.
Uses PEP-based filtering: any precursor with PEP ≤ min_pep in at least one run is kept.

# Arguments
- `psms_paths`: Paths to PSM files from first pass search
- `prec_mzs`: Vector of precursor m/z values
- `rt_irt`: Dictionary mapping file indices to RT-iRT conversion models
- `prec_is_decoy`: Boolean vector indicating decoy status per precursor index
- `max_q_val`: Maximum q-value threshold for iRT statistics
- `min_pep`: PEP threshold — precursors with PEP ≤ this in any run pass to second pass

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
1. First pass: Collects best matches and mean iRT; builds set of PEP-passing precursors
2. Filters to precursors that achieved PEP ≤ min_pep in at least one run
3. Second pass: Calculates iRT variance for remaining precursors
"""
function get_best_precursors_accross_runs(
                         psms_paths::Vector{String},
                         prec_mzs::AbstractVector{Float32},
                         rt_irt::Dict{Int64, RtConversionModel},
                         prec_is_decoy::AbstractVector{Bool};
                         max_q_val::Float32 = 0.01f0,
                         min_pep::Float16 = Float16(0.5)
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

    # First pass: collect best matches, mean iRT, and PEP-passing precursors

    pep_passing_precs = Set{UInt32}()
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

        # Collect precursors with PEP ≤ threshold in this run
        pep_col = psms[:PEP]
        prec_col = psms[:precursor_idx]
        for i in eachindex(prec_col)
            if pep_col[i] <= min_pep
                push!(pep_passing_precs, prec_col[i])
            end
        end
    end

    # Handle case where no valid files were processed
    if n_valid_files == 0
        @warn "No valid PSM files found for cross-run analysis"
        return prec_to_best_prob
    end

    # Diagnostics: pre-filter counts
    n_pre_filter = length(prec_to_best_prob)
    n_pre_targets = count(pid -> !prec_is_decoy[pid], keys(prec_to_best_prob))
    n_pre_decoys = n_pre_filter - n_pre_targets
    @info "Pre-filter pool: $n_pre_filter precursors ($n_pre_targets targets, $n_pre_decoys decoys)"

    # Filter to precursors that passed PEP ≤ min_pep in at least one run
    for key in collect(keys(prec_to_best_prob))
        if !(key in pep_passing_precs)
            delete!(prec_to_best_prob, key)
        end
    end

    n_post_filter = length(prec_to_best_prob)
    n_post_targets = count(pid -> !prec_is_decoy[pid], keys(prec_to_best_prob))
    n_post_decoys = n_post_filter - n_post_targets
    n_removed = n_pre_filter - n_post_filter
    n_removed_targets = n_pre_targets - n_post_targets
    n_removed_decoys = n_pre_decoys - n_post_decoys
    @info "PEP ≤ $(min_pep) filter: kept $n_post_filter ($n_post_targets targets, $n_post_decoys decoys), " *
          "removed $n_removed ($n_removed_targets targets, $n_removed_decoys decoys)"

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
