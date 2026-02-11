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

"""
    sort_of_percolator!(psms, features, match_between_runs; kwargs...) -> Models

Generic PSM scoring function using LightGBM with cross-validation.
Delegates to the trait-based `percolator_scoring!` function.

# Arguments
- `psms::AbstractPSMContainer`: PSMs to score (modified in-place)
- `features::Vector{Symbol}`: Feature columns to use
- `match_between_runs::Bool`: Whether to enable MBR features

# Keyword Arguments
- `max_q_value_lightgbm_rescore::Float32`: Q-value threshold for training (default: 0.01)
- `min_PEP_neg_threshold_itr::Float32`: PEP threshold for negative mining (default: 0.90)
- `feature_fraction::Float64`: LightGBM feature fraction (default: 0.5)
- `learning_rate::Float64`: LightGBM learning rate (default: 0.15)
- `min_data_in_leaf::Int`: Minimum samples per leaf (default: 1)
- `bagging_fraction::Float64`: LightGBM bagging fraction (default: 0.5)
- `min_gain_to_split::Float64`: Minimum gain to split (default: 0.0)
- `max_depth::Int`: Maximum tree depth (default: 10)
- `num_leaves::Int`: Maximum leaves per tree (default: 63)
- `iter_scheme::Vector{Int}`: Boosting rounds per iteration (default: [100, 200, 200])
- `show_progress::Bool`: Show progress bar (default: true)
- `verbose_logging::Bool`: Enable verbose output (default: false)

# Returns
- `Dict{UInt8, Vector{Any}}`: Dictionary mapping CV fold to trained models
"""
function sort_of_percolator!(psms::AbstractPSMContainer,
                  features::Vector{Symbol},
                  match_between_runs::Bool = true;
                  max_q_value_lightgbm_rescore::Float32 = 0.01f0,
                  max_q_value_mbr_itr::Float32 = 0.20f0,
                  min_PEP_neg_threshold_itr = 0.90f0,
                  feature_fraction::Float64 = 0.5,
                  learning_rate::Float64 = 0.15,
                  min_data_in_leaf::Int = 1,
                  bagging_fraction::Float64 = 0.5,
                  min_gain_to_split::Float64 = 0.0,
                  max_depth::Int = 10,
                  num_leaves::Int = 63,
                  iter_scheme::Vector{Int} = [100, 200, 200],
                  print_importance::Bool = false,
                  show_progress::Bool = true,
                  verbose_logging::Bool = false)

    # Separate base features from MBR features
    base_features = [f for f in features if !startswith(String(f), "MBR_")]
    mbr_features = [f for f in features if startswith(String(f), "MBR_")]

    # Build ScoringConfig from keyword arguments
    config = ScoringConfig(
        LightGBMScorer(Dict{Symbol,Any}(
            :feature_fraction => feature_fraction,
            :learning_rate => learning_rate,
            :min_data_in_leaf => min_data_in_leaf,
            :bagging_fraction => bagging_fraction,
            :min_gain_to_split => min_gain_to_split,
            :max_depth => max_depth,
            :num_leaves => num_leaves
        )),
        RandomPairing(),
        QValueNegativeMining(max_q_value_lightgbm_rescore, Float32(min_PEP_neg_threshold_itr)),
        IterativeFeatureSelection(base_features, mbr_features, length(iter_scheme)),
        FixedIterationScheme(iter_scheme),
        match_between_runs ? PairBasedMBR(max_q_value_lightgbm_rescore) : NoMBR()
    )

    # Delegate to trait-based implementation
    return percolator_scoring!(psms, config; show_progress=show_progress, verbose=verbose_logging)
end

function summarize_precursors!(psms::AbstractDataFrame; q_cutoff::Float32 = 0.01f0)
    # Diagnostic: Show isotope and pairing interaction
    n_unique_pairs = length(unique(psms.pair_id))
    unique_isotopes = unique(psms.isotopes_captured)
    n_unique_isotopes = length(unique_isotopes)

    # Compute pair specific features that rely on decoys and chromatograms
    pair_groups = collect(pairs(groupby(psms, [:pair_id, :isotopes_captured])))
    n_pair_isotope_groups = length(pair_groups)

    @debug_l2 "MBR Feature Computation: $n_unique_pairs unique pair_ids × $n_unique_isotopes isotope combinations = $n_pair_isotope_groups groups"
    @debug_l2 "Isotope combinations present: $unique_isotopes"

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
            p = sub_psms.trace_prob[i]
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

            # assign, for each run, the best index in "any other" run
            for r in runs
                run_best_indices[r] = (r == r1 ? best_i[r2] : best_i[r1])
            end
        end

        # Compute MBR features
        num_runs_passing = length(unique(sub_psms.ms_file_idx[sub_psms.q_value .<= q_cutoff]))
        for i in 1:nrow(sub_psms)
            current_ms_file_idx = sub_psms.ms_file_idx[i]
            any_passing_in_current_ms_file_idx = any(sub_psms.q_value[sub_psms.ms_file_idx .== current_ms_file_idx] .<= q_cutoff)
            sub_psms.MBR_num_runs[i] = num_runs_passing - any_passing_in_current_ms_file_idx

            idx = Int(sub_psms.ms_file_idx[i]) - offset + 1
            best_idx = run_best_indices[idx]
            if best_idx == 0 || sub_psms.MBR_num_runs[i] == 0
                sub_psms.MBR_best_irt_diff[i]           = -1.0f0
                sub_psms.MBR_rv_coefficient[i]          = -1.0f0
                sub_psms.MBR_is_best_decoy[i]           = true
                sub_psms.MBR_log2_weight_ratio[i]       = -1.0f0
                sub_psms.MBR_log2_explained_ratio[i]    = -1.0f0
                sub_psms.MBR_max_pair_prob[i]           = -1.0f0
                sub_psms.MBR_is_missing[i]              = true
                continue
            end

            best_log2_weights = log2.(sub_psms.weights[best_idx])
            best_iRTs = sub_psms.irts[best_idx]
            best_log2_weights_padded, weights_padded = pad_equal_length(best_log2_weights, log2.(sub_psms.weights[i]))
            best_iRTs_padded, iRTs_padded = pad_rt_equal_length(best_iRTs, sub_psms.irts[i])

            sub_psms.MBR_max_pair_prob[i] = sub_psms.trace_prob[best_idx]
            best_residual = irt_residual(sub_psms, best_idx)
            current_residual = irt_residual(sub_psms, i)
            sub_psms.MBR_best_irt_diff[i] = abs(best_residual - current_residual)
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
    psms[!, :trace_prob]      = zeros(Float32, n)
    psms[!, :q_value]   = zeros(Float64, n)

    if match_between_runs
        psms[!, :MBR_max_pair_prob]             = zeros(Float32, n)
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


# MBR utility functions (irt_residual, MBR_rv_coefficient, pad_equal_length,
# pad_rt_equal_length) are defined in mbr_update.jl (loaded earlier in include order)
