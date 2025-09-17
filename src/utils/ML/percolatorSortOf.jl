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

#############################################################################
# Target-Decoy Pairing Implementation
#############################################################################

const PAIRING_RANDOM_SEED = 1844  # Fixed seed for reproducible pairing
const IRT_BIN_SIZE = 1000

#############################################################################
# Running Statistics Helper Functions
#############################################################################

"""
    update_pair_statistics(current_stats, new_prob::Float32)

Updates running statistics for MBR pair probabilities in a memory-efficient manner.
Maintains best 2, worst 2, running mean, and count without storing all values.
"""
function update_pair_statistics(current_stats, new_prob::Float32)
    # Update count and running mean
    new_count = current_stats.count_pairs + 1
    new_mean = (current_stats.mean_prob * current_stats.count_pairs + new_prob) / new_count

    # Update best probabilities (existing logic enhanced)
    new_best_1, new_best_2 = if new_prob > current_stats.best_prob_1
        (new_prob, current_stats.best_prob_1)
    elseif new_prob > current_stats.best_prob_2
        (current_stats.best_prob_1, new_prob)
    else
        (current_stats.best_prob_1, current_stats.best_prob_2)
    end

    # Update worst probabilities (new logic)
    new_worst_1, new_worst_2 = if new_count == 1
        (new_prob, zero(Float32))  # First probability
    elseif new_count == 2
        (min(current_stats.worst_prob_1, new_prob), max(current_stats.worst_prob_1, new_prob))
    elseif new_prob < current_stats.worst_prob_1
        (new_prob, current_stats.worst_prob_1)  # New minimum
    elseif new_prob < current_stats.worst_prob_2
        (current_stats.worst_prob_1, new_prob)  # New second minimum
    else
        (current_stats.worst_prob_1, current_stats.worst_prob_2)  # No change
    end

    return merge(current_stats, (
        best_prob_1 = new_best_1,
        best_prob_2 = new_best_2,
        worst_prob_1 = new_worst_1,
        worst_prob_2 = new_worst_2,
        mean_prob = new_mean,
        count_pairs = new_count
    ))
end

#=
"""
    assign_random_target_decoy_pairs!(psms::DataFrame)

Assigns random target-decoy pairings within iRT bins to PSMs, overriding library-based pair_id values.

This function implements a minority-to-majority pairing strategy where precursors from the smaller 
set (targets or decoys) are randomly paired with precursors from the larger set within iRT bins.
Excess decoys are paired with unpaired targets from nearby bins using cross-bin overflow handling.

# Arguments
- `psms::DataFrame`: PSM data with columns: precursor_idx, target, irts, isotopes_captured, pair_id

# Algorithm
1. Creates iRT bins with ~1000 targets per bin for stratified pairing
2. Randomly pairs targets/decoys within bins (minority gets 100% pairing)
3. Handles overflow by pairing excess decoys with nearby unpaired targets
4. Updates pair_id column with new assignments (missing for unpaired)
5. Validates 1:1 precursor relationships with comprehensive diagnostics

# Returns
- Modifies psms DataFrame in-place by updating the pair_id column
- Uses fixed random seed (1844) for reproducible results
"""
function assign_random_target_decoy_pairs!(psms::DataFrame)
    @info "Starting iRT-stratified target-decoy pairing..."
    
    # Set random seed for reproducibility
    Random.seed!(PAIRING_RANDOM_SEED)
    
    # 1. Create iRT bins
    irt_bins = create_irt_bins(psms)
    @info "Created $(length(irt_bins)) iRT bins"
    
    # 2. Extract unique precursors with their iRT values
    target_precursors = get_unique_precursors_with_irt(psms, true)
    decoy_precursors = get_unique_precursors_with_irt(psms, false)
    
    @info "Found $(length(target_precursors)) unique targets, $(length(decoy_precursors)) unique decoys"
    
    # 3. Create stratified random pairings within iRT bins
    pairings = create_stratified_pairings(target_precursors, decoy_precursors, irt_bins)
    
    @info "Created $(length(pairings)) target-decoy pairs across iRT bins"
    
    # 4. Add pair_idx column and assign values
    assign_pair_indices!(psms, pairings)
    
    # 5. Report pairing statistics
    report_pairing_statistics(psms)
    
    # 6. Assign unique pair_id values to unpaired PSMs to prevent downstream errors
    unpaired_mask = psms.pair_id .== 0
    n_unpaired = sum(unpaired_mask)
    
    if n_unpaired > 0
        @info "Assigning unique pair_id values to $n_unpaired unpaired PSMs"
        max_pair_id = length(pairings) > 0 ? maximum(psms.pair_id[psms.pair_id .> 0]) : UInt32(0)
        
        # Group unpaired PSMs by precursor_idx and assign unique pair_ids
        unpaired_precursors = unique(psms.precursor_idx[unpaired_mask])
        for (i, precursor_idx) in enumerate(unpaired_precursors)
            precursor_mask = (psms.precursor_idx .== precursor_idx) .& unpaired_mask
            psms.pair_id[precursor_mask] .= max_pair_id + UInt32(i)
        end
        
        @info "Assigned unique pair_ids $(max_pair_id + 1) to $(max_pair_id + length(unpaired_precursors)) for unpaired precursors"
    end
    
    return psms
end

function get_unique_precursors_with_irt(psms::DataFrame, is_target::Bool)
    mask = psms.target .== is_target
    
    # Group by precursor_idx and take first iRT value (they should all be the same for a given precursor)
    precursor_groups = groupby(psms[mask, [:precursor_idx, :irts]], :precursor_idx)
    
    precursor_irt_map = Tuple{UInt32, Float64}[]
    for group in precursor_groups
        precursor_idx = group.precursor_idx[1]
        irt_value = group.irts[1]  # Take first iRT value for this precursor
        push!(precursor_irt_map, (precursor_idx, irt_value))
    end
    
    return precursor_irt_map
end

function create_irt_bins(psms::DataFrame)
    # Get all target precursors with their iRT values
    target_precursors = get_unique_precursors_with_irt(psms, true)
    n_targets = length(target_precursors)
    
    if n_targets < TARGET_PRECURSORS_PER_BIN
        @info "Only $n_targets targets available, using single iRT bin"
        return [(-Inf, Inf)]  # Single bin containing all
    end
    
    # Calculate number of bins needed
    n_bins = max(1, div(n_targets, TARGET_PRECURSORS_PER_BIN))
    
    # Sort targets by iRT and determine bin boundaries
    sorted_irts = sort([irt for (_, irt) in target_precursors])
    
    # Create bin boundaries at quantiles
    bin_edges = Vector{Float64}(undef, n_bins + 1)
    bin_edges[1] = -Inf
    bin_edges[end] = Inf
    
    for i in 2:n_bins
        quantile_pos = (i - 1) / n_bins
        idx = max(1, min(length(sorted_irts), round(Int, quantile_pos * length(sorted_irts))))
        bin_edges[i] = sorted_irts[idx]
    end
    
    # Convert to (min, max) tuples
    bins = [(bin_edges[i], bin_edges[i+1]) for i in 1:n_bins]
    
    @info "Created $n_bins iRT bins with ~$(div(n_targets, n_bins)) targets per bin"
    return bins
end

function create_stratified_pairings(target_precursors::Vector{Tuple{UInt32, Float64}}, 
                                   decoy_precursors::Vector{Tuple{UInt32, Float64}}, 
                                   irt_bins::Vector{Tuple{Float64, Float64}})
    n_targets = length(target_precursors)
    n_decoys = length(decoy_precursors)
    
    if n_targets == 0 || n_decoys == 0
        @warn "No target-decoy pairs possible (one set is empty)"
        return Tuple{UInt32, UInt32}[]
    end
    
    # Pre-allocate result array
    all_pairings = Vector{Tuple{UInt32, UInt32}}()
    sizehint!(all_pairings, min(n_targets, n_decoys))
    
    # Assign precursors to bins
    target_bins = assign_precursors_to_bins(target_precursors, irt_bins)
    decoy_bins = assign_precursors_to_bins(decoy_precursors, irt_bins)
    
    # Determine global pairing strategy
    minority_is_target = n_targets <= n_decoys
    @info "Pairing strategy: Each $(minority_is_target ? "target" : "decoy") paired with random $(minority_is_target ? "decoy" : "target") within iRT bins"
    
    # Track unpaired precursors per bin for smart overflow handling
    unpaired_targets_by_bin = Vector{Vector{UInt32}}(undef, length(irt_bins))
    unpaired_decoys_by_bin = Vector{Vector{UInt32}}(undef, length(irt_bins))
    
    # Diagnostic counters
    bins_with_excess_targets = 0
    bins_with_excess_decoys = 0
    
    # Pair within each bin
    for (bin_idx, (min_irt, max_irt)) in enumerate(irt_bins)
        bin_targets = get(target_bins, bin_idx, UInt32[])
        bin_decoys = get(decoy_bins, bin_idx, UInt32[])
        
        @info "Bin $bin_idx analysis: $(length(bin_targets)) targets, $(length(bin_decoys)) decoys (iRT: $(round(min_irt, digits=2)) - $(round(max_irt, digits=2)))"
        
        # Diagnose bin imbalance
        if length(bin_targets) > length(bin_decoys)
            bins_with_excess_targets += 1
            @info "  → Excess targets in bin $bin_idx: $(length(bin_targets) - length(bin_decoys)) targets will remain unpaired"
        elseif length(bin_decoys) > length(bin_targets)
            bins_with_excess_decoys += 1
            @info "  → Excess decoys in bin $bin_idx: $(length(bin_decoys) - length(bin_targets)) decoys need pairing with targets from other bins"
        end
        
        bin_pairings, unpaired_t, unpaired_d = pair_within_bin(bin_targets, bin_decoys)
        append!(all_pairings, bin_pairings)
        
        unpaired_targets_by_bin[bin_idx] = unpaired_t
        unpaired_decoys_by_bin[bin_idx] = unpaired_d
        
        @info "  → Bin $bin_idx result: $(length(bin_pairings)) pairs, $(length(unpaired_t)) unpaired targets, $(length(unpaired_d)) unpaired decoys"
    end
    
    @info "Bin summary: $bins_with_excess_targets bins with excess targets, $bins_with_excess_decoys bins with excess decoys"
    
    # Handle overflow pairing: prioritize pairing excess decoys with nearby unpaired targets
    overflow_pairings = handle_cross_bin_overflow(unpaired_targets_by_bin, unpaired_decoys_by_bin, irt_bins)
    append!(all_pairings, overflow_pairings)
    
    if length(overflow_pairings) > 0
        @info "Overflow pairing across bins: $(length(overflow_pairings)) additional pairs"
    end
    
    total_unpaired = (minority_is_target ? length(unpaired_decoys_by_bin) : length(unpaired_targets_by_bin)) - length(overflow_pairings)
    @info "Final result: $(length(all_pairings)) total pairs, $total_unpaired unpaired precursors"
    
    return all_pairings
end

function assign_precursors_to_bins(precursors::Vector{Tuple{UInt32, Float64}}, 
                                 irt_bins::Vector{Tuple{Float64, Float64}})
    bin_assignments = Dict{Int, Vector{UInt32}}()
    
    # Pre-allocate bin vectors
    for i in 1:length(irt_bins)
        bin_assignments[i] = Vector{UInt32}()
    end
    
    for (precursor_idx, irt_val) in precursors
        # Find which bin this precursor belongs to
        bin_idx = findfirst(((min_irt, max_irt),) -> min_irt <= irt_val < max_irt, irt_bins)
        
        if bin_idx === nothing
            # Handle edge case - assign to last bin if >= max value
            bin_idx = length(irt_bins)
        end
        
        push!(bin_assignments[bin_idx], precursor_idx)
    end
    
    return bin_assignments
end

function pair_within_bin(targets::Vector{UInt32}, decoys::Vector{UInt32})
    n_targets = length(targets)
    n_decoys = length(decoys)
    
    if n_targets == 0 || n_decoys == 0
        return Tuple{UInt32, UInt32}[], targets, decoys
    end
    
    # Always try to pair as many as possible within bin
    n_pairs = min(n_targets, n_decoys)
    
    # Use randperm for efficient random selection without copying full arrays
    target_indices = randperm(n_targets)
    decoy_indices = randperm(n_decoys)
    
    # Create pairings
    pairings = [(targets[target_indices[i]], decoys[decoy_indices[i]]) for i in 1:n_pairs]
    
    # Determine unpaired precursors
    unpaired_targets = targets[target_indices[(n_pairs+1):end]]
    unpaired_decoys = decoys[decoy_indices[(n_pairs+1):end]]
    
    return pairings, unpaired_targets, unpaired_decoys
end

function handle_cross_bin_overflow(unpaired_targets_by_bin::Vector{Vector{UInt32}}, 
                                  unpaired_decoys_by_bin::Vector{Vector{UInt32}},
                                  irt_bins::Vector{Tuple{Float64, Float64}})
    cross_bin_pairings = Vector{Tuple{UInt32, UInt32}}()
    n_bins = length(irt_bins)
    
    # Priority: pair excess decoys with unpaired targets in nearby bins
    for bin_idx in 1:n_bins
        excess_decoys = unpaired_decoys_by_bin[bin_idx]
        
        if isempty(excess_decoys)
            continue
        end
        
        @info "Processing $(length(excess_decoys)) excess decoys from bin $bin_idx"
        
        # Find unpaired targets in nearby bins (search outward from current bin)
        for distance in 0:(n_bins-1)
            if isempty(excess_decoys)
                break
            end
            
            # Check bins at this distance
            nearby_bins = Int[]
            if distance == 0
                push!(nearby_bins, bin_idx)
            else
                # Add bins at +/- distance
                if bin_idx - distance >= 1
                    push!(nearby_bins, bin_idx - distance)
                end
                if bin_idx + distance <= n_bins
                    push!(nearby_bins, bin_idx + distance)
                end
            end
            
            for target_bin_idx in nearby_bins
                if isempty(excess_decoys)
                    break
                end
                
                available_targets = unpaired_targets_by_bin[target_bin_idx]
                if isempty(available_targets)
                    continue
                end
                
                # Pair as many as possible
                n_pairs = min(length(excess_decoys), length(available_targets))
                
                if n_pairs > 0
                    @info "  → Cross-bin pairing: $(n_pairs) decoys from bin $bin_idx with targets from bin $target_bin_idx (distance: $distance)"
                    
                    # Use randperm for random selection
                    decoy_indices = randperm(length(excess_decoys))[1:n_pairs]
                    target_indices = randperm(length(available_targets))[1:n_pairs]
                    
                    for i in 1:n_pairs
                        push!(cross_bin_pairings, (available_targets[target_indices[i]], excess_decoys[decoy_indices[i]]))
                    end
                    
                    # Remove paired precursors
                    excess_decoys = excess_decoys[setdiff(1:length(excess_decoys), decoy_indices)]
                    unpaired_targets_by_bin[target_bin_idx] = available_targets[setdiff(1:length(available_targets), target_indices)]
                end
            end
        end
        
        if !isempty(excess_decoys)
            @info "  → Warning: $(length(excess_decoys)) decoys from bin $bin_idx could not be paired (no available targets)"
        end
    end
    
    return cross_bin_pairings
end

function assign_pair_indices!(psms::DataFrame, pairings::Vector{Tuple{UInt32, UInt32}})
    # Initialize pair_id column with zeros (will be updated later for unpaired)
    psms[!, :pair_id] = zeros(UInt32, nrow(psms))
    
    # Assign pair_id for each target-decoy pair
    for (pair_id, (target_precursor, decoy_precursor)) in enumerate(pairings)
        # All instances of this target precursor get this pair_id
        target_mask = (psms.precursor_idx .== target_precursor) .& (psms.target .== true)
        psms.pair_id[target_mask] .= UInt32(pair_id)
        
        # All instances of this decoy precursor get this pair_id  
        decoy_mask = (psms.precursor_idx .== decoy_precursor) .& (psms.target .== false)
        psms.pair_id[decoy_mask] .= UInt32(pair_id)
    end
end

function report_pairing_statistics(psms::DataFrame)
    # Count paired vs unpaired rows (before final assignment)
    paired_mask = psms.pair_id .> 0
    n_paired = sum(paired_mask)
    n_unpaired = sum(.!paired_mask)
    
    @info "Initial pairing: $n_paired PSMs paired, $n_unpaired PSMs unpaired"
    
    # Count by target/decoy
    paired_targets = sum(paired_mask .& psms.target)
    paired_decoys = sum(paired_mask .& .!psms.target)
    unpaired_targets = sum(.!paired_mask .& psms.target)
    unpaired_decoys = sum(.!paired_mask .& .!psms.target)
    
    @info "Targets: $paired_targets paired, $unpaired_targets unpaired"
    @info "Decoys: $paired_decoys paired, $unpaired_decoys unpaired"
    
    # Validate pairing balance (only for paired PSMs)
    validate_pairing_balance(psms)
end

function validate_pairing_balance(psms::DataFrame)
    paired_psms = psms[psms.pair_id .> 0, :]
    
    if nrow(paired_psms) == 0
        @warn "No paired PSMs found"
        return
    end
    
    # Count targets and decoys for each pair_id
    pair_counts = combine(groupby(paired_psms, :pair_id)) do group
        (
            n_targets = sum(group.target),
            n_decoys = sum(.!group.target),
            target_precursors = length(unique(group.precursor_idx[group.target])),
            decoy_precursors = length(unique(group.precursor_idx[.!group.target]))
        )
    end
    
    @info "Validation: $(nrow(pair_counts)) pairs created"
    
    # Check for 1:1 precursor pairing
    invalid_pairs = pair_counts[(pair_counts.target_precursors .!= 1) .| (pair_counts.decoy_precursors .!= 1), :]
    if nrow(invalid_pairs) > 0
        @error "Invalid pairing detected: some pairs don't have exactly 1 target and 1 decoy precursor"
        println(invalid_pairs)
    else
        @info "✓ All pairs have exactly 1 target and 1 decoy precursor"
    end
end
=#

function getIrtBins(irts::AbstractVector{R}) where {R <:Real}
    sort_idx = sortperm(irts)
    bin_idx, bin_count = zero(UInt32), zero(UInt32)
    bin_idxs = similar(irts, UInt32, length(irts))
    for idx in sort_idx
        bin_count += one(UInt32)
        bin_idxs[idx] = bin_idx
        if bin_count >= IRT_BIN_SIZE
            bin_idx += one(UInt32)
            bin_count = zero(UInt32)
        end
    end
    return bin_idxs 
end


function getIrtBins!(psms::AbstractDataFrame)
    psms[!, :irt_bin_idx] = getIrtBins(psms.irt_pred)
    return psms
end


function assign_random_target_decoy_pairs!(psms::DataFrame)
    last_pair_id = zero(UInt32)
    psms[!,:pair_id] = zeros(UInt32, nrow(psms))  # Initialize pair_id column
    psms[!,:irt_bin_idx] = getIrtBins(psms.irt_pred)  # Ensure irt_bin_idx column exists
    for (irt_bin_idx, sub_psms) in pairs(groupby(psms, :irt_bin_idx))
        last_pair_id = assignPairIds!(sub_psms, last_pair_id)
    end
end


function assignPairIds!(psms::AbstractDataFrame, last_pair_id::UInt32)
    psms[!,:pair_id], last_pair_id = assign_pair_ids(
        psms.target, psms.decoy, psms.precursor_idx, psms.irt_bin_idx, last_pair_id
    )
    return last_pair_id
end

function assign_pair_ids(
    target::AbstractVector{Bool}, decoy::AbstractVector{Bool},
    precursor_idx::AbstractVector{UInt32}, irt_bin_idx::AbstractVector{UInt32},
    last_pair_id::UInt32
)
    targets = unique(precursor_idx[target])
    decoys = unique(precursor_idx[decoy])
    target_perm = randperm(MersenneTwister(PAIRING_RANDOM_SEED), length(targets))
    precursor_idx_to_pair_id = Dict{UInt32,UInt32}()  # Map from precursor_idx to pair_id
    pair_ids = similar(precursor_idx, UInt32)
    for i in range(1, min(length(targets), length(decoys)))
        last_pair_id += one(UInt32)
        precursor_idx_to_pair_id[targets[target_perm[i]]] = last_pair_id
        precursor_idx_to_pair_id[decoys[i]] = last_pair_id
    end
    if length(decoys) < length(targets)
        #@user_warn "Fewer decoy precursors ($(length(decoys))) than target precursors ($(length(targets))) in iRT bin $(first(irt_bin_idx)). Some targets will remain unpaired."
        for i in range(length(decoys)+1, length(targets))
            last_pair_id += one(UInt32)
            precursor_idx_to_pair_id[targets[target_perm[i]]] = last_pair_id
        end
    elseif length(targets) < length(decoys)
        @user_warn "Fewer target precursors ($(length(targets))) than decoy precursors ($(length(decoys))) in iRT bin $(first(irt_bin_idx)). Some decoys will remain unpaired."
        for i in range(length(targets)+1, length(decoys))
            last_pair_id += one(UInt32)
            precursor_idx_to_pair_id[decoys[i]] = last_pair_id
        end
    end
    for (row_idx, precursor_idx) in enumerate(precursor_idx)
        pair_ids[row_idx] = precursor_idx_to_pair_id[precursor_idx]
    end
    return pair_ids, last_pair_id
end

function sort_of_percolator_in_memory!(psms::DataFrame, 
                  features::Vector{Symbol},
                  match_between_runs::Bool = true;
                  max_q_value_xgboost_rescore::Float32 = 0.01f0,
                  max_q_value_xgboost_mbr_rescore::Float32 = 0.20f0,
                  min_PEP_neg_threshold_xgboost_rescore = 0.90f0,
                  colsample_bytree::Float64 = 0.5,
                  eta::Float64 = 0.15,
                  min_child_weight::Int = 1,
                  subsample::Float64 = 0.5,
                  gamma::Float64 = 0.0,
                  max_depth::Int = 10,
                  iter_scheme::Vector{Int} = [100, 200, 200],
                  print_importance::Bool = false,
                  show_progress::Bool = true,
                  verbose_logging::Bool = false)

    # Apply random target-decoy pairing before ML training
    assign_random_target_decoy_pairs!(psms)
    
    #Faster if sorted first (handle missing pair_id values)
    sort!(psms, [:pair_id, :isotopes_captured])
    # Display target/decoy/entrapment counts for training dataset
    if verbose_logging
        n_targets = sum(psms.target)
        n_decoys = sum(.!psms.target)
        n_entrapments = hasproperty(psms, :entrapment) ? sum(psms.entrapment) : 0
        n_total = nrow(psms)

        if n_entrapments > 0
            @user_info "ML Training Dataset: $n_targets targets, $n_decoys decoys, $n_entrapments entrapments (total: $n_total PSMs)"
        else
            @user_info "ML Training Dataset: $n_targets targets, $n_decoys decoys (total: $n_total PSMs)"
        end
    end

    prob_test   = zeros(Float32, nrow(psms))  # final CV predictions
    prob_train  = zeros(Float32, nrow(psms))  # temporary, used during training
    MBR_estimates = zeros(Float32, nrow(psms)) # optional MBR layer
    nonMBR_estimates  = zeros(Float32, nrow(psms)) # keep track of last nonMBR test scores

    unique_cv_folds = unique(psms[!, :cv_fold])
    models = Dict{UInt8, Vector{EvoTrees.EvoTree}}()
    mbr_start_iter = length(iter_scheme)

    cv_fold_col = psms[!, :cv_fold]
    fold_indices = Dict(fold => findall(==(fold), cv_fold_col) for fold in unique_cv_folds)
    train_indices = Dict(fold => findall(!=(fold), cv_fold_col) for fold in unique_cv_folds)

    Random.seed!(1776)
    non_mbr_features = [f for f in features if !startswith(String(f), "MBR_")]

    pbar = show_progress ? ProgressBar(total=length(unique_cv_folds) * length(iter_scheme)) : nothing

    for test_fold_idx in unique_cv_folds

        initialize_prob_group_features!(psms, match_between_runs)

        train_idx = train_indices[test_fold_idx]
        test_idx  = fold_indices[test_fold_idx]
        
        psms_train = @view psms[train_idx, :]
        psms_test  = @view psms[test_idx, :]

        # Display counts for this CV fold
        if verbose_logging
            n_train_targets = sum(psms_train.target)
            n_train_decoys = sum(.!psms_train.target)
            n_test_targets = sum(psms_test.target)
            n_test_decoys = sum(.!psms_test.target)

            if hasproperty(psms, :entrapment)
                n_train_entrapments = sum(psms_train.entrapment)
                n_test_entrapments = sum(psms_test.entrapment)
                @user_info "Fold $test_fold_idx - Train: $n_train_targets targets, $n_train_decoys decoys, $n_train_entrapments entrapments | Test: $n_test_targets targets, $n_test_decoys decoys, $n_test_entrapments entrapments"
            else
                @user_info "Fold $test_fold_idx - Train: $n_train_targets targets, $n_train_decoys decoys | Test: $n_test_targets targets, $n_test_decoys decoys"
            end
        end

        fold_models = Vector{EvoTrees.EvoTree}(undef, length(iter_scheme))

        for (itr, num_round) in enumerate(iter_scheme)
            psms_train_itr = get_training_data_for_iteration!(psms_train,
                                                              itr,
                                                              match_between_runs,
                                                              max_q_value_xgboost_rescore,
                                                              max_q_value_xgboost_mbr_rescore,
                                                              min_PEP_neg_threshold_xgboost_rescore,
                                                              itr >= mbr_start_iter)

            train_feats = itr < mbr_start_iter ? non_mbr_features : features
            
            bst = train_booster(psms_train_itr, train_feats, num_round;
                               colsample=colsample_bytree,
                               eta=eta,
                               min_child_weight=min_child_weight,
                               subsample=subsample,
                               gamma=gamma,
                               max_depth=max_depth)
                               
            fold_models[itr] = bst

            # Print feature importances for each iteration and fold
            if print_importance
                importances = EvoTrees.importance(bst)
                @user_info "Feature Importances - Fold $(test_fold_idx), Iteration $(itr) ($(length(importances)) features):"
                for i in 1:10:length(importances)
                    chunk = importances[i:min(i+9, end)]
                    feat_strs = ["$(feat):$(round(score, digits=3))" for (feat, score) in chunk]
                    @user_info "  " * join(feat_strs, " | ")
                end
            end

            #predict_fold!(bst, psms_train, psms_test, train_feats)
            # **temporary predictions for training only**
            prob_train[train_idx] = predict(bst, psms_train)
            psms_train[!,:prob] = prob_train[train_idx]
            get_qvalues!(psms_train.prob, psms_train.target, psms_train.q_value)

            # **predict held-out fold**
            prob_test[test_idx] = predict(bst, psms_test)
            psms_test[!,:prob] = prob_test[test_idx]

            if itr == (mbr_start_iter - 1)
			    nonMBR_estimates[test_idx] = prob_test[test_idx]
            end

            if match_between_runs
                update_mbr_features!(psms_train, psms_test, prob_test,
                                     test_idx, itr, mbr_start_iter,
                                     max_q_value_xgboost_rescore)
            end

            show_progress && update(pbar)

            if (!match_between_runs) && itr == (mbr_start_iter - 1)
                break
            end
        end
        # Make predictions on hold out data.
        if match_between_runs
            MBR_estimates[test_idx] = psms_test.prob
        else
            prob_test[test_idx] = psms_test.prob
        end
        # Store models for this fold
        models[test_fold_idx] = fold_models
    end

    if match_between_runs
        # Determine which precursors failed the q-value cutoff prior to MBR
        qvals_prev = Vector{Float32}(undef, length(nonMBR_estimates))
        get_qvalues!(nonMBR_estimates, psms.target, qvals_prev)
        pass_mask = (qvals_prev .<= max_q_value_xgboost_rescore)
        prob_thresh = any(pass_mask) ? minimum(nonMBR_estimates[pass_mask]) : typemax(Float32)
        # Label as transfer candidates only those failing the q-value cutoff but
        # whose best matched pair surpassed the passing probability threshold.
        psms[!, :MBR_transfer_candidate] .= .!pass_mask .&
                                            (psms.MBR_max_pair_prob .>= prob_thresh)

        # Use the final MBR probabilities for all precursors
        psms[!, :prob] = MBR_estimates
    else
        psms[!, :prob] = prob_test
    end
    
    return models
end

function sort_of_percolator_out_of_memory!(psms::DataFrame, 
                    file_paths::Vector{String},
                    features::Vector{Symbol},
                    match_between_runs::Bool = true; 
                    max_q_value_xgboost_rescore::Float32 = 0.01f0,
                    max_q_value_xgboost_mbr_rescore::Float32 = 0.20f0,
                    min_PEP_neg_threshold_xgboost_rescore::Float32 = 0.90f0,
                    colsample_bytree::Float64 = 0.5, 
                    eta::Float64 = 0.15, 
                    min_child_weight::Int = 1, 
                    subsample::Float64 = 0.5, 
                    gamma::Float64 = 0.0, 
                    max_depth::Int = 10,
                    iter_scheme::Vector{Int} = [100, 200, 200],
                    print_importance::Bool = false)

    function getBestScorePerPrec!(
        prec_to_best_score_new::Dictionary,
        file_paths::Vector{String},
        models::Dictionary{UInt8,EvoTrees.EvoTree},
        features::Vector{Symbol},
        match_between_runs::Bool;
        is_last_iteration::Bool = false)
    
        # Reset counts for new scores
        reset_precursor_scores!(prec_to_best_score_new)
            
        for file_path in file_paths
            psms_subset = DataFrame(Arrow.Table(file_path))
            
            probs = predict_cv_models(models, psms_subset, features)
            
            if match_between_runs && !is_last_iteration
                #Update maximum probabilities for tracked precursors 
                qvals = zeros(Float32, nrow(psms_subset))
                get_qvalues!(probs, psms_subset.target, qvals)

                for (i, pair_id) in enumerate(psms_subset[!,:pair_id])
                    prob = probs[i]
                    key = (pair_id = pair_id, isotopes = psms_subset[i,:isotopes_captured])
                    if haskey(prec_to_best_score_new, key)
                        scores = prec_to_best_score_new[key]

                        # Update running statistics with new probability
                        updated_stats = update_pair_statistics(scores, prob)

                        if prob > scores.best_prob_1
                           new_scores = merge(updated_stats, (
                                # replace best_prob_2 with best_prob_1
                                best_prob_2                     = scores.best_prob_1,
                                best_log2_weights_2             = scores.best_log2_weights_1,
                                best_irts_2                     = scores.best_irts_1,
                                best_weight_2                   = scores.best_weight_1,
                                best_log2_intensity_explained_2 = scores.best_log2_intensity_explained_1,
                                best_ms_file_idx_2              = scores.best_ms_file_idx_1,
                                is_best_decoy_2                 = scores.is_best_decoy_1,
                                # overwrite best_prob_1
                                best_prob_1                     = prob,
                                best_log2_weights_1             = log2.(psms_subset.weights[i]),
                                best_irts_1                     = psms_subset.irts[i],
                                best_weight_1                   = psms_subset.weight[i],
                                best_log2_intensity_explained_1 = psms_subset.log2_intensity_explained[i],
                                best_ms_file_idx_1              = psms_subset.ms_file_idx[i],
                                is_best_decoy_1                 = psms_subset.decoy[i]
                            ))
                            prec_to_best_score_new[key] = new_scores

                        elseif prob > scores.best_prob_2
                            # overwrite best_prob_2
                            new_scores = merge(updated_stats, (
                                best_prob_2                     = prob,
                                best_log2_weights_2             = log2.(psms_subset.weights[i]),
                                best_irts_2                     = psms_subset.irts[i],
                                best_weight_2                   = psms_subset.weight[i],
                                best_log2_intensity_explained_2 = psms_subset.log2_intensity_explained[i],
                                best_ms_file_idx_2              = psms_subset.ms_file_idx[i],
                                is_best_decoy_2                 = psms_subset.decoy[i]
                            ))
                            prec_to_best_score_new[key] = new_scores
                        else
                            # No change to best/second best, but update running stats
                            prec_to_best_score_new[key] = updated_stats
                        end

                        if qvals[i] <= max_q_value_xgboost_rescore
                            push!(scores.unique_passing_runs, psms_subset.ms_file_idx[i])
                        end

                    else
                        insert!(prec_to_best_score_new, key, (
                                best_prob_1                     = prob,
                                best_prob_2                     = zero(Float32),
                                worst_prob_1                    = prob,
                                worst_prob_2                    = zero(Float32),
                                mean_prob                       = prob,
                                count_pairs                     = Int32(1),
                                best_log2_weights_1             = log2.(psms_subset.weights[i]),
                                best_log2_weights_2             = Vector{Float32}(),
                                best_irts_1                     = psms_subset.irts[i],
                                best_irts_2                     = Vector{Float32}(),
                                best_weight_1                   = psms_subset.weight[i],
                                best_weight_2                   = zero(Float32),
                                best_log2_intensity_explained_1 = psms_subset.log2_intensity_explained[i],
                                best_log2_intensity_explained_2 = zero(Float32),
                                best_ms_file_idx_1              = psms_subset.ms_file_idx[i],
                                best_ms_file_idx_2              = zero(UInt32),
                                is_best_decoy_1                 = psms_subset.decoy[i],
                                is_best_decoy_2                 = false,
                                unique_passing_runs             = ( qvals[i] <= max_q_value_xgboost_rescore ?
                                                                    Set{UInt16}([psms_subset.ms_file_idx[i]]) :
                                                                    Set{UInt16}() )
                            ))
                    end
                end
            end


            if is_last_iteration
                if match_between_runs
                    update_mbr_probs!(psms_subset, probs, max_q_value_xgboost_rescore)
                else
                    psms_subset.prob = probs
                end

            end
        end
    

        # Compute probs and features for next round
        for file_path in file_paths
            psms_subset = DataFrame(Tables.columntable(Arrow.Table(file_path)))
            probs = predict_cv_models(models, psms_subset, features)

            for (i, pair_id) in enumerate(psms_subset[!,:pair_id])
                psms_subset[i,:prob] = probs[i]
                
                    
                if match_between_runs && !is_last_iteration
                    key = (pair_id = pair_id, isotopes = psms_subset[i,:isotopes_captured])
                    if haskey(prec_to_best_score_new, key)
                        scores = prec_to_best_score_new[key]

                        psms_subset.MBR_num_runs[i] = length(scores.unique_passing_runs)

                        best_log2_weights = Float32[]
                        best_irts = Float32[]
                        best_weight = zero(Float32)
                        best_log2_ie = zero(Float32)

                        if (scores.best_ms_file_idx_1 != psms_subset.ms_file_idx[i]) &&
                           (!isempty(scores.best_log2_weights_1))
                            best_log2_weights                   = scores.best_log2_weights_1
                            best_irts                           = scores.best_irts_1
                            best_weight                         = scores.best_weight_1
                            best_log2_ie                        = scores.best_log2_intensity_explained_1
                            psms_subset.MBR_max_pair_prob[i]    = scores.best_prob_1
                            psms_subset.MBR_mean_pair_prob[i]   = scores.mean_prob
                            psms_subset.MBR_min_pair_prob[i]    = scores.worst_prob_1
                            psms_subset.MBR_num_pairs[i]        = scores.count_pairs
                            MBR_is_best_decoy                   = scores.is_best_decoy_1
                        elseif (scores.best_ms_file_idx_2 != psms_subset.ms_file_idx[i]) &&
                               (!isempty(scores.best_log2_weights_2))
                            best_log2_weights                   = scores.best_log2_weights_2
                            best_irts                           = scores.best_irts_2
                            best_weight                         = scores.best_weight_2
                            best_log2_ie                        = scores.best_log2_intensity_explained_2
                            psms_subset.MBR_max_pair_prob[i]    = scores.best_prob_2
                            psms_subset.MBR_mean_pair_prob[i]   = scores.mean_prob
                            psms_subset.MBR_min_pair_prob[i]    = scores.worst_prob_1
                            psms_subset.MBR_num_pairs[i]        = scores.count_pairs
                            MBR_is_best_decoy                   = scores.is_best_decoy_2
                        else
                            psms_subset.MBR_best_irt_diff[i]        = -1.0f0
                            psms_subset.MBR_rv_coefficient[i]       = -1.0f0
                            psms_subset.MBR_is_best_decoy[i]        = true
                            psms_subset.MBR_max_pair_prob[i]        = -1.0f0
                            psms_subset.MBR_mean_pair_prob[i]       = -1.0f0
                            psms_subset.MBR_min_pair_prob[i]        = -1.0f0
                            psms_subset.MBR_num_pairs[i]            = 0
                            psms_subset.MBR_log2_weight_ratio[i]    = -1.0f0
                            psms_subset.MBR_log2_explained_ratio[i] = -1.0f0
                            psms_subset.MBR_is_missing[i]           = true
                            continue
                        end

                        best_log2_weights_padded, weights_padded = pad_equal_length(best_log2_weights, log2.(psms_subset.weights[i]))
                        best_iRTs_padded, iRTs_padded = pad_rt_equal_length(best_irts, psms_subset.irts[i])

                        best_irt_at_apex = best_irts[argmax(best_log2_weights)]
                        psms_subset.MBR_best_irt_diff[i] = abs(best_irt_at_apex - psms_subset.irts[i][argmax(psms_subset.weights[i])])
                        psms_subset.MBR_rv_coefficient[i] = MBR_rv_coefficient(best_log2_weights_padded, best_iRTs_padded, weights_padded, iRTs_padded)
                        psms_subset.MBR_log2_weight_ratio[i] = log2(psms_subset.weight[i] / best_weight)
                        psms_subset.MBR_log2_explained_ratio[i] = psms_subset.log2_intensity_explained[i] - best_log2_ie
                        psms_subset.MBR_is_best_decoy[i] = MBR_is_best_decoy
                    end
                end
            end

            write_subset(
                file_path,
                psms_subset,
                probs,
                match_between_runs,
                max_q_value_xgboost_rescore;
                dropVectors = is_last_iteration,
            )
        end
        
        return prec_to_best_score_new
    end


    unique_cv_folds = unique(psms[!, :cv_fold])
    #Train the model for 1:K-1 cross validation folds and apply to the held-out fold
    models = Dictionary{UInt8, Vector{EvoTrees.EvoTree}}()
    pbar = ProgressBar(total=length(unique_cv_folds)*length(iter_scheme))
    Random.seed!(1776);
    non_mbr_features = [f for f in features if !startswith(String(f), "MBR_")]

    for test_fold_idx in unique_cv_folds#(0, 1)#range(1, n_folds)
        #Clear prob stats 
        initialize_prob_group_features!(psms, match_between_runs)
        # Get training data
        psms_train = @view(psms[findall(x -> x != test_fold_idx, psms[!, :cv_fold]), :])

        for (itr, num_round) in enumerate(iter_scheme)

            psms_train_itr = get_training_data_for_iteration!(psms_train, 
                                                                itr,
                                                                match_between_runs, 
                                                                max_q_value_xgboost_rescore,
                                                                max_q_value_xgboost_mbr_rescore,
                                                                min_PEP_neg_threshold_xgboost_rescore,
                                                                itr >= length(iter_scheme))
            ###################
            #Train a model on the n-1 training folds.
            train_feats = itr < length(iter_scheme) ? non_mbr_features : features
            bst = train_booster(psms_train_itr, train_feats, num_round;
                               colsample=colsample_bytree,
                               eta=eta,
                               min_child_weight=min_child_weight,
                               subsample=subsample,
                               gamma=gamma,
                               max_depth=max_depth)
            if !haskey(models, test_fold_idx)
                insert!(
                    models,
                    test_fold_idx,
                    Vector{EvoTrees.EvoTree}([bst])
                )
            else
                push!(models[test_fold_idx], bst)
            end
            # Print feature importances
            if print_importance
                importances = EvoTrees.importance(bst)
                @user_info "Feature Importances ($(length(importances)) features):"
                for i in 1:10:length(importances)
                    chunk = importances[i:min(i+9, end)]
                    feat_strs = ["$(feat):$(round(score, digits=3))" for (feat, score) in chunk]
                    @user_info "  " * join(feat_strs, " | ")
                end
            end

            # Get probabilities for training sample so we can get q-values
            psms_train[!,:prob] = EvoTrees.predict(bst, psms_train)
            
            if match_between_runs
                summarize_precursors!(psms_train, q_cutoff = max_q_value_xgboost_rescore)
            end

            show_progress && update(pbar)
        end
    end
    
    pbar = ProgressBar(total=length(iter_scheme))
    prec_to_best_score = Dictionary{@NamedTuple{pair_id::UInt32,
                                                isotopes::Tuple{Int8,Int8}},
                                    @NamedTuple{best_prob_1::Float32,
                                                best_prob_2::Float32,
                                                worst_prob_1::Float32,
                                                worst_prob_2::Float32,
                                                mean_prob::Float32,
                                                count_pairs::Int32,
                                                best_log2_weights_1::Vector{Float32},
                                                best_log2_weights_2::Vector{Float32},
                                                best_irts_1::Vector{Float32},
                                                best_irts_2::Vector{Float32},
                                                best_weight_1::Float32,
                                                best_weight_2::Float32,
                                                best_log2_intensity_explained_1::Float32,
                                                best_log2_intensity_explained_2::Float32,
                                                best_ms_file_idx_1::UInt32,
                                                best_ms_file_idx_2::UInt32,
                                                is_best_decoy_1::Bool,
                                                is_best_decoy_2::Bool,
                                                unique_passing_runs::Set{UInt16}}}()

    for (train_iter, num_round) in enumerate(iter_scheme)
        models_for_iter = Dictionary{UInt8,EvoTrees.EvoTree}()
        for test_fold_idx in unique_cv_folds
            insert!(models_for_iter, test_fold_idx, models[test_fold_idx][train_iter])
        end
        prec_to_best_score = getBestScorePerPrec!(
            prec_to_best_score,
            file_paths,
            models_for_iter,
            features,
            match_between_runs;
            is_last_iteration = (train_iter == length(iter_scheme))
        )
        update(pbar)
    end

    return models
end

function train_booster(psms::AbstractDataFrame, features, num_round;
                       colsample::Float64,
                       eta::Float64,
                       min_child_weight::Int,
                       subsample::Float64,
                       gamma::Float64,
                       max_depth::Int)

    config = EvoTreeRegressor(
        loss=:logloss,
        nrounds = num_round,
        max_depth = max_depth,
        min_weight = min_child_weight,
        rowsample = subsample,
        colsample = colsample,
        eta = eta,
        gamma = gamma
    )
    model = fit(config, psms; target_name = :target, feature_names = features, verbosity = 0)
    return model
end

function predict_fold!(bst, psms_train::AbstractDataFrame,
                       psms_test::AbstractDataFrame, features)
    psms_test[!, :prob] = predict(bst, psms_test)
    psms_train[!, :prob] = predict(bst, psms_train)
    get_qvalues!(psms_train.prob, psms_train.target, psms_train.q_value)
end

function update_mbr_features!(psms_train::AbstractDataFrame,
                              psms_test::AbstractDataFrame,
                              prob_test::Vector{Float32},
                              test_fold_idxs,
                              itr::Int,
                              mbr_start_iter::Int,
                              max_q_value_xgboost_rescore::Float32)
    if itr >= mbr_start_iter - 1
        get_qvalues!(psms_test.prob, psms_test.target, psms_test.q_value)
        summarize_precursors!(psms_test, q_cutoff = max_q_value_xgboost_rescore)
        summarize_precursors!(psms_train, q_cutoff = max_q_value_xgboost_rescore)
    end
    if itr == mbr_start_iter - 1
        prob_test[test_fold_idxs] = psms_test.prob
    end
end

function summarize_precursors!(psms::AbstractDataFrame; q_cutoff::Float32 = 0.01f0)
    # Compute pair specific features that rely on decoys and chromatograms
    pair_groups = collect(pairs(groupby(psms, [:pair_id, :isotopes_captured])))
    Threads.@threads for idx in eachindex(pair_groups)
        _, sub_psms = pair_groups[idx]
        
        # Efficient way to find the top 2 precursors so we can do MBR on the 
        # best precursor match that isn't itself. It's always one of the top 2.

        # single pass: record the best PSM index & prob per run
        offset = Int(minimum(sub_psms.ms_file_idx))
        range_len = Int(maximum(sub_psms.ms_file_idx)) - offset + 1
        best_i = zeros(Int, range_len)
        best_p = fill(-Inf, range_len)
        for (i, run) in enumerate(sub_psms.ms_file_idx)
            idx = Int(run) - offset + 1
            p = sub_psms.prob[i]
            if p > best_p[idx]
                best_p[idx] = p
                best_i[idx] = i
            end
        end

        # if more than one run, find the global top-2 runs by their best-PSM prob
        run_best_indices = zeros(Int, range_len)
        runs = findall(!=(0), best_i)
        if length(runs) > 1
            # track top two runs (r1 > r2)
            r1 = 0; p1 = -Inf
            r2 = 0; p2 = -Inf
            for r in runs
                p = best_p[r]
                if p > p1
                    r2, p2 = r1, p1
                    r1, p1 = r, p
                elseif p > p2
                    r2, p2 = r, p
                end
            end

            # assign, for each run, the best index in “any other” run
            for r in runs
                run_best_indices[r] = (r == r1 ? best_i[r2] : best_i[r1])
            end
        end

        # Compute MBR features
        num_runs_passing = length(sub_psms.ms_file_idx[sub_psms.q_value .<= q_cutoff])
        for i in 1:nrow(sub_psms)
            sub_psms.MBR_num_runs[i] = num_runs_passing - (sub_psms.q_value[i] .<= q_cutoff)

            idx = Int(sub_psms.ms_file_idx[i]) - offset + 1
            best_idx = run_best_indices[idx]
            if best_idx == 0 || sub_psms.MBR_num_runs[i] == 0
                sub_psms.MBR_best_irt_diff[i]           = -1.0f0
                sub_psms.MBR_rv_coefficient[i]          = -1.0f0
                sub_psms.MBR_is_best_decoy[i]           = true
                sub_psms.MBR_log2_weight_ratio[i]       = -1.0f0
                sub_psms.MBR_log2_explained_ratio[i]    = -1.0f0
                sub_psms.MBR_max_pair_prob[i]           = -1.0f0
                sub_psms.MBR_mean_pair_prob[i]          = -1.0f0
                sub_psms.MBR_min_pair_prob[i]           = -1.0f0
                sub_psms.MBR_num_pairs[i]               = 0
                sub_psms.MBR_is_missing[i]              = true
                continue
            end

            best_log2_weights = log2.(sub_psms.weights[best_idx])
            best_iRTs = sub_psms.irts[best_idx]
            best_log2_weights_padded, weights_padded = pad_equal_length(best_log2_weights, log2.(sub_psms.weights[i]))
            best_iRTs_padded, iRTs_padded = pad_rt_equal_length(best_iRTs, sub_psms.irts[i])
            
            best_irt_at_apex = sub_psms.irts[best_idx][argmax(best_log2_weights)]
            sub_psms.MBR_max_pair_prob[i] = sub_psms.prob[best_idx]
            # Note: In this algorithm path, we only have single best match per run
            # So mean/min are the same as max, and count is 1
            sub_psms.MBR_mean_pair_prob[i] = sub_psms.prob[best_idx]
            sub_psms.MBR_min_pair_prob[i] = sub_psms.prob[best_idx]
            sub_psms.MBR_num_pairs[i] = 1
            sub_psms.MBR_best_irt_diff[i] = abs(best_irt_at_apex - sub_psms.irts[i][argmax(sub_psms.weights[i])])
            sub_psms.MBR_rv_coefficient[i] = MBR_rv_coefficient(best_log2_weights_padded, best_iRTs_padded, weights_padded, iRTs_padded)
            sub_psms.MBR_log2_weight_ratio[i] = log2(sub_psms.weight[i] / sub_psms.weight[best_idx])
            sub_psms.MBR_log2_explained_ratio[i] = sub_psms.log2_intensity_explained[i] - sub_psms.log2_intensity_explained[best_idx]
            sub_psms.MBR_is_best_decoy[i] = sub_psms.decoy[best_idx]
        end
    end
end

function initialize_prob_group_features!(
    psms::AbstractDataFrame,
    match_between_runs::Bool
)
    n = nrow(psms)
    psms[!, :prob]      = zeros(Float32, n)
    psms[!, :q_value]   = zeros(Float64, n)

    if match_between_runs
        psms[!, :MBR_max_pair_prob]             = zeros(Float32, n)
        psms[!, :MBR_mean_pair_prob]            = zeros(Float32, n)
        psms[!, :MBR_min_pair_prob]             = zeros(Float32, n)
        psms[!, :MBR_num_pairs]                 = zeros(Int32, n)
        psms[!, :MBR_best_irt_diff]             = zeros(Float32, n)
        psms[!, :MBR_log2_weight_ratio]         = zeros(Float32, n)
        psms[!, :MBR_log2_explained_ratio]      = zeros(Float32, n)
        psms[!, :MBR_rv_coefficient]            = zeros(Float32, n)
        psms[!, :MBR_is_best_decoy]             = trues(n)
        psms[!, :MBR_num_runs]                  = zeros(Int32, n)
        psms[!, :MBR_transfer_candidate]        = falses(n)
        psms[!, :MBR_is_missing]                = falses(n)
    end

    return psms
end

function get_training_data_for_iteration!(
    psms_train::AbstractDataFrame,
    itr::Int,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32,
    last_iter::Bool
)
   
    if itr == 1
        # Train on all precursors during first iteration. 
        return copy(psms_train)
    else
        # Do a shallow copy to avoid overwriting target/decoy labels
        psms_train_itr = copy(psms_train)

        # Convert the worst-scoring targets to negatives using PEP estimate
        order = sortperm(psms_train_itr.prob, rev=true)
        sorted_scores  = psms_train_itr.prob[order]
        sorted_targets = psms_train_itr.target[order]
        PEPs = Vector{Float32}(undef, length(order))
        get_PEP!(sorted_scores, sorted_targets, PEPs; doSort=false)

        idx_cutoff = findfirst(x -> x >= min_PEP_neg_threshold_xgboost_rescore, PEPs)
        if !isnothing(idx_cutoff)
            worst_idxs = order[idx_cutoff:end]
            psms_train_itr.target[worst_idxs] .= false
        end

        # Also train on top scoring MBR candidates if requested
        if match_between_runs && last_iter
            # Determine prob threshold for precursors passing the q-value threshold
            max_prob_threshold = minimum(
                psms_train_itr.prob[
                    psms_train_itr.target .& (psms_train_itr.q_value .<= max_q_value_xgboost_rescore)
                ]
            )

            # Hacky way to ensure anything passing the initial q-value threshold
            # will pass the next q-value threshold
            psms_train_itr.q_value[psms_train_itr.q_value .<= max_q_value_xgboost_rescore] .= 0.0
            psms_train_itr.q_value[psms_train_itr.q_value .> max_q_value_xgboost_rescore]  .= 1.0

            # Must have at least one precursor passing the q-value threshold,
            # and the best precursor can't be a decoy
            psms_train_mbr = subset(
                psms_train_itr,
                [:MBR_is_best_decoy, :MBR_max_pair_prob, :prob, :MBR_is_missing] => ByRow((d, mp, p, im) ->
                    (!im && !d && mp >= max_prob_threshold && p < max_prob_threshold)
                );
                view = true
            )

            # Compute MBR q-values.
            get_qvalues!(psms_train_mbr[!,:prob], psms_train_mbr[!,:target], psms_train_mbr[!,:q_value])

            # Take all decoys and targets passing q_thresh (all 0's now) or mbr_q_thresh
            psms_train_itr = subset(
                psms_train_itr,
                [:target, :q_value, :MBR_is_best_decoy, :MBR_is_missing] => ByRow((t, q, MBR_d, im) -> (!t) || (t && !im && !MBR_d && q <= max_q_value_xgboost_mbr_rescore))
            )
        else
            # Take all decoys and targets passing q_thresh
            psms_train_itr = subset(
                psms_train_itr,
                [:target, :q_value] => ByRow((t,q) -> (!t) || (t && q <= max_q_value_xgboost_rescore))
            )
        end

        return psms_train_itr
    end
end

function dropVectorColumns!(df)
    to_drop = String[]
    for col in names(df)
        if eltype(df[!, col]) <: AbstractVector
            push!(to_drop, col)
        end
    end
    # 2) Drop those columns in place
    select!(df, Not(to_drop))
end

"""
    reset_precursor_scores!(dict)

Set all values of `dict` to an empty precursor score tuple.  This helps reuse
the same dictionary between iterations without reallocating memory.
"""
function reset_precursor_scores!(dict)
    for key in keys(dict)
        dict[key] = (
            best_prob_1 = zero(Float32),
            best_prob_2 = zero(Float32),
            best_log2_weights_1 = Vector{Float32}(),
            best_log2_weights_2 = Vector{Float32}(),
            best_irts_1 = Vector{Float32}(),
            best_irts_2 = Vector{Float32}(),
            best_weight_1 = zero(Float32),
            best_weight_2 = zero(Float32),
            best_log2_intensity_explained_1 = zero(Float32),
            best_log2_intensity_explained_2 = zero(Float32),
            best_ms_file_idx_1 = zero(UInt32),
            best_ms_file_idx_2 = zero(UInt32),
            is_best_decoy_1 = false,
            is_best_decoy_2 = false,
            unique_passing_runs = Set{UInt16}(),
        )
    end
    return dict
end

"""
    predict_cv_models(models, df, features)

Return a vector of probabilities for `df` using the cross validation `models`.
"""
function predict_cv_models(models::Dictionary{UInt8,EvoTrees.EvoTree},
                           df::AbstractDataFrame,
                           features::Vector{Symbol})
    probs = zeros(Float32, nrow(df))
    for (fold_idx, bst) in pairs(models)
        fold_rows = findall(==(fold_idx), df[!, :cv_fold])
        if !isempty(fold_rows)
            probs[fold_rows] = predict(bst, df[fold_rows, :])
        end
    end
    return probs
end

"""
    update_mbr_probs!(df, probs, qval_thresh)

Store final MBR probabilities and mark transfer candidates as those
failing the pre-MBR q-value threshold but whose best matched pair passed
the corresponding probability cutoff.
"""
function update_mbr_probs!(
    df::AbstractDataFrame,
    probs::AbstractVector{Float32},
    qval_thresh::Float32,
)
    prev_qvals = similar(df.prob)
    get_qvalues!(df.prob, df.target, prev_qvals)
    pass_mask = (prev_qvals .<= qval_thresh) .& df.target
    prob_thresh = any(pass_mask) ? minimum(df.prob[pass_mask]) : typemax(Float32)
    df[!, :MBR_transfer_candidate] = (prev_qvals .> qval_thresh) .&
                                     (df.MBR_max_pair_prob .>= prob_thresh)
    df[!, :prob] = probs
    return df
end

"""
    write_subset(file_path, df, probs, match_between_runs, qval_thresh; dropVectors=false)

Write the updated subset to disk, optionally dropping vector columns.
The `qval_thresh` argument is used to mark transfer candidates when
`match_between_runs` is true and `dropVectors` is set.
"""
function write_subset(
    file_path::String,
    df::DataFrame,
    probs::AbstractVector{Float32},
    match_between_runs::Bool,
    qval_thresh::Float32;
    dropVectors::Bool=false,
)
    if dropVectors
        if match_between_runs
            update_mbr_probs!(df, probs, qval_thresh)
        else
            df[!, :prob] = probs
        end
        writeArrow(file_path, dropVectorColumns!(df))
    else
        df[!, :prob] = probs
        writeArrow(file_path, convert_subarrays(df))
    end
end

function MBR_rv_coefficient(weights_A::AbstractVector{<:Real},
    times_A::AbstractVector{<:Real},
    weights_B::AbstractVector{<:Real},
    times_B::AbstractVector{<:Real})

    # Construct two Nx2 matrices, each row is (weight, time)
    X = hcat(collect(weights_A), collect(times_A))
    Y = hcat(collect(weights_B), collect(times_B))

    # Compute cross-products (Gram matrices)
    Sx = X' * X
    Sy = Y' * Y

    # Numerator: trace(Sx * Sy)
    numerator = tr(Sx * Sy)

    # Denominator: sqrt( trace(Sx*Sx)* trace(Sy*Sy) )
    denominator = sqrt(tr(Sx * Sx) * tr(Sy * Sy))

    # Protect against zero in denominator (e.g. if X or Y is all zeros)
    if denominator == 0
        return 0.0
    end

    return numerator / denominator
end

function pad_equal_length(x::AbstractVector, y::AbstractVector)

    function pad(data, d)
        left  = div(d, 2)
        right = d - left        # put extra on the right if d is odd
        padded_data = vcat(zeros(eltype(data), left), data, zeros(eltype(data), right))
        return padded_data
    end

    d = length(y) - length(x)
    if d == 0
        # No padding needed
        return (x, y)
    elseif d > 0
        # Pad x
        return (pad(x, d), y)
    else
        # Pad y
        return (x, pad(y, abs(d)))
    end
end


function pad_rt_equal_length(x::AbstractVector, y::AbstractVector)

    function pad(data, n, d, rt_step)
        left  = div(d, 2) 
        right = d - left        # put extra on the right if d is odd
        padded_data = vcat(zeros(eltype(data), left), data, zeros(eltype(data), right))
        @inbounds @fastmath for i in range(1, left)
            padded_data[i] = padded_data[left+1] - (rt_step * ((left+1) - i))
        end 
        @inbounds @fastmath for i in range(1, right)
            padded_data[i+left+n] = padded_data[left+n] + (rt_step * i)
        end 
        return padded_data
    end

    nx = length(x)
    ny = length(y)
    d = ny - nx
    rt_step_x = (x[end] - x[1]) / nx
    rt_step_y = (y[end] - y[1]) / ny
    rt_step = nx > 1 ? rt_step_x : rt_step_y

    if d == 0
        # No padding needed
        return (x, y)
    elseif d > 0
        # Pad x
        return (pad(x, nx, d, rt_step), y)
    else
        # Pad y
        return (x, pad(y, ny, abs(d), rt_step))
    end
end

function summarize_prob(probs::AbstractVector{Float32})
    minimum(probs), maximum(probs), mean(probs)
end

function convert_subarrays(df::DataFrame)
    for col_name in names(df)
        col_data = df[!, col_name]
        # If it's a SubArray, let's convert it to a plain vector:
        if eltype(col_data) <: AbstractVector
           df[!, col_name] = [Vector(col_data[i]) for i in eachindex(col_data)]
        end
    end
    return df
end
