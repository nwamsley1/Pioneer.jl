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

#==========================================================
Pipeline Operations for Protein Inference
==========================================================#

"""
    add_peptide_metadata(precursors::LibraryPrecursors)

Add peptide sequence and protein information from precursor library.
Enriches PSM data with metadata needed for inference.
"""
function add_peptide_metadata(precursors::LibraryPrecursors)
    desc = "add_peptide_metadata"
    
    op = function(df)
        # Get library data
        all_sequences = getSequence(precursors)::AbstractVector{String}
        all_accessions = getAccessionNumbers(precursors)::AbstractVector{String}
        all_is_decoys = getIsDecoy(precursors)::AbstractVector{Bool}
        all_entrap_ids = getEntrapmentGroupId(precursors)::AbstractVector{UInt8}
        
        # Extract precursor_idx column with type assertion for performance
        precursor_idx = df.precursor_idx::AbstractVector{UInt32}
        n_rows = length(precursor_idx)
        
        # Add columns if they don't exist
        if !hasproperty(df, :sequence)
            sequences = Vector{String}(undef, n_rows)
            for i in 1:n_rows
                sequences[i] = all_sequences[precursor_idx[i]]
            end
            df.sequence = sequences
        end
        
        if !hasproperty(df, :accession_numbers)
            accessions = Vector{String}(undef, n_rows)
            for i in 1:n_rows
                accessions[i] = all_accessions[precursor_idx[i]]
            end
            df.accession_numbers = accessions
        end
        
        # Add target/decoy status if needed
        if !hasproperty(df, :is_decoy)
            is_decoy_vec = Vector{Bool}(undef, n_rows)
            for i in 1:n_rows
                is_decoy_vec[i] = all_is_decoys[precursor_idx[i]]
            end
            df.is_decoy = is_decoy_vec
        end
        
        # Add entrapment group if needed
        if !hasproperty(df, :entrap_id)
            entrap_vec = Vector{UInt8}(undef, n_rows)
            for i in 1:n_rows
                entrap_vec[i] = all_entrap_ids[precursor_idx[i]]
            end
            df.entrap_id = entrap_vec
        end
        
        return df
    end
    
    return desc => op
end

"""
    validate_peptide_data()

Ensure all required columns exist for protein inference.
"""
function validate_peptide_data()
    desc = "validate_peptide_data"
    
    op = function(df)
        required_cols = [:sequence, :accession_numbers, :is_decoy, 
                        :entrap_id, :precursor_idx]
        
        missing_cols = setdiff(required_cols, Symbol.(names(df)))
        if !isempty(missing_cols)
            error("Missing required columns for protein inference: $missing_cols")
        end
        
        return df
    end
    
    return desc => op
end

#==========================================================
Inference Wrapper Functions
==========================================================#

"""
    apply_inference_to_dataframe(df::DataFrame, precursors::LibraryPrecursors)

Apply the core infer_proteins algorithm to a prepared DataFrame.
Returns structured inference results using type-safe protein inference.
"""
function apply_inference_to_dataframe(df::DataFrame, precursors::LibraryPrecursors)
    if nrow(df) == 0
        # Return empty result for empty input
        return InferenceResult(
            Dictionary{PeptideKey, ProteinKey}()
        )
    end
    
    # Extract unique peptide-protein pairs
    unique_pairs = unique(df, [:sequence, :accession_numbers, :is_decoy, :entrap_id])
    
    # Convert to format expected by infer_proteins
    proteins_vec = Vector{ProteinKey}(undef, nrow(unique_pairs))
    peptides_vec = Vector{PeptideKey}(undef, nrow(unique_pairs))
    
    for (i, row) in enumerate(eachrow(unique_pairs))
        proteins_vec[i] = ProteinKey(
            row.accession_numbers,
            !row.is_decoy,  # Convert is_decoy to is_target
            row.entrap_id
        )
        
        peptides_vec[i] = PeptideKey(
            row.sequence,
            !row.is_decoy,  # Convert is_decoy to is_target
            row.entrap_id
        )
    end
    
    # Call inference algorithm
    return infer_proteins(proteins_vec, peptides_vec)
end

#==========================================================
Post-Inference Pipeline Operations
==========================================================#

"""
    add_inferred_protein_column(inference_result::InferenceResult)

Add inferred protein group assignments to PSMs based on inference results.
"""
function add_inferred_protein_column(inference_result::InferenceResult)
    desc = "add_inferred_protein_column"
    
    op = function(df)
        # Extract columns with type assertions for performance
        sequences = df.sequence::AbstractVector{String}
        is_decoy = df.is_decoy::AbstractVector{Bool}
        entrap_ids = df.entrap_id::AbstractVector{UInt8}
        accession_numbers = df.accession_numbers::AbstractVector{String}
        
        n_rows = length(sequences)
        inferred_proteins = Vector{Union{Missing, String}}(undef, n_rows)
        
        for i in 1:n_rows
            pep_key = PeptideKey(
                sequences[i],
                !is_decoy[i],
                entrap_ids[i]
            )
            
            if haskey(inference_result.peptide_to_protein, pep_key)
                protein_key = inference_result.peptide_to_protein[pep_key]
                inferred_proteins[i] = protein_key.name
            else
                # Fallback to original protein
                inferred_proteins[i] = missing#accession_numbers[i]
            end
        end
        
        df.inferred_protein_group = inferred_proteins
        return df
    end
    
    return desc => op
end

"""
    add_quantification_flag(inference_result::InferenceResult)

Add flag indicating whether peptide should be used for quantification.

Peptides present in the inference result (unique peptides assigned to minimal protein set)
are marked as usable for quantification. Peptides not in the result (shared peptides or
filtered out) are marked as NOT usable for quantification.
"""
function add_quantification_flag(inference_result::InferenceResult)
    desc = "add_quantification_flag"

    op = function(df)
        # Extract columns with type assertions for performance
        sequences = df.sequence::AbstractVector{String}
        is_decoy = df.is_decoy::AbstractVector{Bool}
        entrap_ids = df.entrap_id::AbstractVector{UInt8}

        n_rows = length(sequences)
        use_for_quant = Vector{Bool}(undef, n_rows)

        for i in 1:n_rows
            pep_key = PeptideKey(
                sequences[i],
                !is_decoy[i],
                entrap_ids[i]
            )

            # Peptide is usable for quantification if it's in the inference results
            # (i.e., it's a unique peptide assigned to a protein)
            use_for_quant[i] = haskey(inference_result.peptide_to_protein, pep_key)
        end

        df.use_for_protein_quant = use_for_quant
        return df
    end

    return desc => op
end

function _median_abs_deviation(values::Vector{Float64})::Float64
    isempty(values) && return 0.0
    med = Statistics.median(values)
    abs_devs = Vector{Float64}(undef, length(values))
    @inbounds for i in eachindex(values)
        abs_devs[i] = abs(values[i] - med)
    end
    return Statistics.median(abs_devs)
end

"""
    estimate_weight_detection_model(df::DataFrame)

Estimate per-file weight calibration for protein coverage surprise features.
Uses unique inferred peptides with valid positive weight values.
"""
function estimate_weight_detection_model(df::DataFrame)
    max_rank = 12
    min_rank_support = 8
    default_model = (
        log_threshold = 0.0f0,
        sigma_log = 1.0f0,
        n_unique_peptides = 0,
        used_fallback = true,
        rank_drop_profile = Float32[0.0f0],
        rank_scale_profile = Float32[1.0f0],
        profiled_rank_count = 1,
        rank_profile_examples = 0,
        used_rank_profile_fallback = true
    )

    n_rows = nrow(df)
    if n_rows == 0
        return default_model
    end

    best_weight_by_protein_peptide = Dict{Tuple{String, Bool, UInt8, String}, Float64}()

    for i in 1:n_rows
        if df.use_for_protein_quant[i] != true
            continue
        end

        protein_name = df.inferred_protein_group[i]
        if ismissing(protein_name)
            continue
        end

        weight_val = Float64(df.weight[i])
        if !isfinite(weight_val) || weight_val <= 0.0
            continue
        end

        target_val = hasproperty(df, :target) ? Bool(df.target[i]) : true
        entrap_val = hasproperty(df, :entrap_id) ? UInt8(df.entrap_id[i]) : UInt8(0)
        key = (String(protein_name), target_val, entrap_val, String(df.sequence[i]))
        if haskey(best_weight_by_protein_peptide, key)
            if weight_val > best_weight_by_protein_peptide[key]
                best_weight_by_protein_peptide[key] = weight_val
            end
        else
            best_weight_by_protein_peptide[key] = weight_val
        end
    end

    if isempty(best_weight_by_protein_peptide)
        return default_model
    end

    log_weights = Float64[log(weight) for weight in values(best_weight_by_protein_peptide)]
    log_threshold = Float64(Statistics.quantile(log_weights, 0.05))

    protein_to_log_weights = Dict{Tuple{String, Bool, UInt8}, Vector{Float64}}()
    for ((protein_name, target_val, entrap_val, _), weight) in best_weight_by_protein_peptide
        protein_key = (protein_name, target_val, entrap_val)
        push!(get!(protein_to_log_weights, protein_key, Float64[]), log(weight))
    end

    pooled_residuals = Float64[]
    for log_vals in values(protein_to_log_weights)
        if length(log_vals) < 2
            continue
        end
        med = Statistics.median(log_vals)
        for v in log_vals
            push!(pooled_residuals, v - med)
        end
    end

    sigma_log = 1.0
    used_fallback = true
    if !isempty(pooled_residuals)
        mad_val = _median_abs_deviation(pooled_residuals)
        if isfinite(mad_val) && mad_val > 0.0
            sigma_log = clamp(mad_val / 0.6744897501960817, 0.25, 2.5)
            used_fallback = false
        end
    end

    rank_drop_profile = fill(0.0f0, max_rank)
    rank_scale_profile = fill(Float32(sigma_log), max_rank)
    profiled_rank_count = 1
    rank_profile_examples = 0
    used_rank_profile_fallback = true
    rank_drop_samples = [Float64[] for _ in 1:max_rank]

    for ((_, target_val, _), log_vals) in protein_to_log_weights
        target_val || continue
        length(log_vals) >= 2 || continue

        sort!(log_vals; rev = true)
        top_log_weight = log_vals[1]
        local_max_rank = min(length(log_vals), max_rank)
        for rank_idx in 2:local_max_rank
            push!(rank_drop_samples[rank_idx], log_vals[rank_idx] - top_log_weight)
        end
    end

    for rank_idx in 2:max_rank
        samples = rank_drop_samples[rank_idx]
        if length(samples) < min_rank_support
            break
        end

        drop_median = Statistics.median(samples)
        residuals = Float64[]
        sizehint!(residuals, length(samples))
        for sample in samples
            push!(residuals, sample - drop_median)
        end

        drop_scale = sigma_log
        mad_val = _median_abs_deviation(residuals)
        if isfinite(mad_val) && mad_val > 0.0
            drop_scale = clamp(mad_val / 0.6744897501960817, 0.25, 2.5)
        end

        rank_drop_profile[rank_idx] = Float32(drop_median)
        rank_scale_profile[rank_idx] = Float32(drop_scale)
        profiled_rank_count = rank_idx
        rank_profile_examples += length(samples)
        used_rank_profile_fallback = false
    end

    return (
        log_threshold = Float32(log_threshold),
        sigma_log = Float32(sigma_log),
        n_unique_peptides = length(log_weights),
        used_fallback = used_fallback,
        rank_drop_profile = rank_drop_profile[1:profiled_rank_count],
        rank_scale_profile = rank_scale_profile[1:profiled_rank_count],
        profiled_rank_count = profiled_rank_count,
        rank_profile_examples = rank_profile_examples,
        used_rank_profile_fallback = used_rank_profile_fallback
    )
end

"""
    add_weight_observation_features(calibration)

Add weight-based protein coverage features using an empirical top-rank companion model.
"""
function add_weight_observation_features(calibration::NamedTuple)
    desc = "add_weight_observation_features"

    op = function(df)
        n_rows = nrow(df)
        expected_additional_from_top = Vector{Float32}(undef, n_rows)
        coverage_match_from_top = Vector{Float32}(undef, n_rows)

        log_threshold = Float64(calibration.log_threshold)
        std_normal = Distributions.Normal()
        rank_drop_profile = Float64.(calibration.rank_drop_profile)
        rank_scale_profile = Float64.(calibration.rank_scale_profile)
        profiled_rank_count = max(Int(calibration.profiled_rank_count), 1)

        for i in 1:n_rows
            top_weight = Float64(df.top_pep_weight[i])
            N_total = max(Int(df.n_possible_peptides[i]), 0)
            k_obs = max(Int(df.n_peptides[i]), 0)

            if top_weight <= 0.0
                expected_additional_from_top[i] = 0.0f0
                coverage_match_from_top[i] = 1.0f0
                continue
            end

            log_top_weight = log(top_weight)
            effective_rank_count = min(max(N_total, 1), profiled_rank_count)
            if effective_rank_count <= 1
                expected_additional_from_top[i] = 0.0f0
                coverage_match_from_top[i] = 1.0f0
                continue
            end

            expected_additional = 0.0
            for rank_idx in 2:effective_rank_count
                rank_scale = rank_scale_profile[rank_idx]
                rank_scale = (isfinite(rank_scale) && rank_scale > 0.0) ? rank_scale : 1.0
                expected_log_weight = log_top_weight + rank_drop_profile[rank_idx]
                rank_detect_prob = clamp(
                    Distributions.cdf(
                        std_normal,
                        (expected_log_weight - log_threshold) / rank_scale
                    ),
                    1e-6,
                    1.0 - 1e-6
                )
                expected_additional += rank_detect_prob
            end

            observed_additional = min(max(k_obs - 1, 0), effective_rank_count - 1)
            coverage_deficit = max(expected_additional - observed_additional, 0.0)
            match_denominator = max(expected_additional, 1.0)

            expected_additional_from_top[i] = Float32(expected_additional)
            coverage_match_from_top[i] = Float32(clamp(1.0 - (coverage_deficit / match_denominator), 0.0, 1.0))
        end

        df.expected_additional_from_top = expected_additional_from_top
        df.coverage_match_from_top = coverage_match_from_top
        return df
    end

    return desc => op
end

"""
    add_pg_score_interaction_features()

Add pg_score interaction terms for the current coverage-surprise feature family.
"""
function add_pg_score_interaction_features()
    desc = "add_pg_score_interaction_features"

    op = function(df)
        pg_score = df.pg_score
        coverage_match_from_top = df.coverage_match_from_top

        n_rows = length(pg_score)
        pg_score_x_coverage_match_from_top = Vector{Float32}(undef, n_rows)

        @inbounds for i in eachindex(
            pg_score,
            coverage_match_from_top
        )
            pg_score_val = Float32(pg_score[i])
            pg_score_x_coverage_match_from_top[i] =
                pg_score_val * Float32(clamp(Float64(coverage_match_from_top[i]), 0.0, 1.0))
        end

        df.pg_score_x_coverage_match_from_top = pg_score_x_coverage_match_from_top
        return df
    end

    return desc => op
end

"""
    add_consensus_precursor_relative_weight_interaction_feature()

Add the interaction terms between `pg_score` and the consensus precursor
relative-weight support and enrichment features.
"""
function add_consensus_precursor_relative_weight_interaction_feature()
    desc = "add_consensus_precursor_relative_weight_interaction_feature"

    op = function(df)
        pg_score = df.pg_score
        consensus_precursor_relative_weight_support = df.consensus_precursor_relative_weight_support
        consensus_precursor_relative_weight_enrichment = df.consensus_precursor_relative_weight_enrichment

        n_rows = length(pg_score)
        pg_score_x_consensus_precursor_relative_weight_support = Vector{Float32}(undef, n_rows)
        pg_score_x_consensus_precursor_relative_weight_enrichment = Vector{Float32}(undef, n_rows)

        @inbounds for i in eachindex(
            pg_score,
            consensus_precursor_relative_weight_support,
            consensus_precursor_relative_weight_enrichment
        )
            pg_score_x_consensus_precursor_relative_weight_support[i] =
                Float32(pg_score[i]) * Float32(consensus_precursor_relative_weight_support[i])
            pg_score_x_consensus_precursor_relative_weight_enrichment[i] =
                Float32(pg_score[i]) * Float32(consensus_precursor_relative_weight_enrichment[i])
        end

        df.pg_score_x_consensus_precursor_relative_weight_support = pg_score_x_consensus_precursor_relative_weight_support
        df.pg_score_x_consensus_precursor_relative_weight_enrichment = pg_score_x_consensus_precursor_relative_weight_enrichment
        return df
    end

    return desc => op
end

function add_consensus_precursor_rank_interaction_feature()
    return add_consensus_precursor_relative_weight_interaction_feature()
end

function _protein_group_probability_column(df::DataFrame)
    if hasproperty(df, :MBR_boosted_prec_prob)
        return :MBR_boosted_prec_prob
    else
        return :prec_prob
    end
end

function _quant_peptides_and_pg_score(
    gdf::AbstractDataFrame,
    quant_mask::AbstractVector{Bool},
    prob_col::Symbol
)
    quant_peptides = unique(gdf[quant_mask, :sequence])
    if isempty(quant_peptides)
        return quant_peptides, 0.0f0
    end

    unique_pep_probs = Vector{Float32}(undef, length(quant_peptides))
    for (i, pep) in pairs(quant_peptides)
        pep_mask = (gdf.sequence .== pep) .& quant_mask
        unique_pep_probs[i] = maximum(gdf[pep_mask, prob_col])
    end

    return quant_peptides, -sum(log.(1.0f0 .- unique_pep_probs))
end

function _direct_quant_mask(df::AbstractDataFrame)
    quant_mask = df.use_for_protein_quant .== true
    if hasproperty(df, :MBR_candidate)
        return quant_mask .& .!df.MBR_candidate
    end
    return quant_mask
end

const CONSENSUS_PRECURSOR_MAX_RUNS = 5
const CONSENSUS_PRECURSOR_RUN_DECAY_RATE = 1.0
const ConsensusRunVote = @NamedTuple{
    pg_score::Float32,
    run_order::Int64,
    normalized_precursors::Vector{Pair{UInt32, Float32}}
}

@inline function _consensus_run_decay(selected_rank::Int)::Float64
    return exp(-CONSENSUS_PRECURSOR_RUN_DECAY_RATE * Float64(selected_rank - 1))
end

"""
    build_precursor_relative_weight_consensus(psm_refs::Vector{PSMFileReference})

Build a dataset-level precursor relative-weight consensus within each inferred
protein group. Consensus weights are derived from direct (non-MBR) quant
precursors after normalizing each precursor by the max precursor weight in that
protein run. Each run's vote is weighted by the run's current protein
`pg_score`. Only the top `CONSENSUS_PRECURSOR_MAX_RUNS` runs per protein group
are retained, with an exponential decay applied by run rank within that top
set.
"""
function build_precursor_relative_weight_consensus(psm_refs::Vector{PSMFileReference})
    consensus_weight_sums = Dict{Tuple{String, Bool, UInt8, UInt32}, Float64}()
    protein_total_vote = Dict{Tuple{String, Bool, UInt8}, Float64}()
    protein_run_votes = Dict{Tuple{String, Bool, UInt8}, Vector{ConsensusRunVote}}()

    for (run_order, psm_ref) in enumerate(psm_refs)
        if !exists(psm_ref)
            continue
        end

        df = load_dataframe(psm_ref)
        prob_col = _protein_group_probability_column(df)

        for gdf in groupby(df, [:inferred_protein_group, :target, :entrap_id])
            quant_mask = gdf.use_for_protein_quant .== true
            _, pg_score = _quant_peptides_and_pg_score(gdf, quant_mask, prob_col)

            direct_mask = _direct_quant_mask(gdf)
            best_weight_by_precursor = Dict{UInt32, Float32}()

            @inbounds for i in eachindex(direct_mask)
                if !direct_mask[i]
                    continue
                end

                precursor_idx = UInt32(gdf.precursor_idx[i])
                weight_val = Float32(gdf.weight[i])
                if haskey(best_weight_by_precursor, precursor_idx)
                    if weight_val > best_weight_by_precursor[precursor_idx]
                        best_weight_by_precursor[precursor_idx] = weight_val
                    end
                else
                    best_weight_by_precursor[precursor_idx] = weight_val
                end
            end

            if isempty(best_weight_by_precursor)
                continue
            end

            protein_name_val = gdf.inferred_protein_group[1]
            if ismissing(protein_name_val)
                continue
            end

            protein_name = String(protein_name_val)
            target = Bool(gdf.target[1])
            entrap_id = UInt8(gdf.entrap_id[1])

            run_max_weight = maximum(values(best_weight_by_precursor))
            if !isfinite(run_max_weight) || run_max_weight <= 0.0f0
                continue
            end

            normalized_precursors = Pair{UInt32, Float32}[]
            sizehint!(normalized_precursors, length(best_weight_by_precursor))
            for (precursor_idx, precursor_weight) in best_weight_by_precursor
                push!(
                    normalized_precursors,
                    precursor_idx => Float32(clamp(precursor_weight / run_max_weight, 0.0f0, 1.0f0))
                )
            end
            sort!(normalized_precursors, by = x -> (-x.second, x.first))
            protein_key = (protein_name, target, entrap_id)
            push!(
                get!(protein_run_votes, protein_key, ConsensusRunVote[]),
                (
                    pg_score = pg_score,
                    run_order = Int64(run_order),
                    normalized_precursors = normalized_precursors
                )
            )
        end
    end

    for (protein_key, run_votes) in protein_run_votes
        sort!(run_votes, by = x -> (-x.pg_score, x.run_order))
        n_selected = min(length(run_votes), CONSENSUS_PRECURSOR_MAX_RUNS)

        for selected_rank in 1:n_selected
            run_vote = run_votes[selected_rank]
            run_weight = Float64(run_vote.pg_score) * _consensus_run_decay(selected_rank)
            protein_total_vote[protein_key] = get(protein_total_vote, protein_key, 0.0) + run_weight

            for precursor in run_vote.normalized_precursors
                key = (protein_key[1], protein_key[2], protein_key[3], precursor.first)
                consensus_weight_sums[key] = get(consensus_weight_sums, key, 0.0) + (run_weight * Float64(precursor.second))
            end
        end
    end

    precursor_relative_weight = Dict{Tuple{String, Bool, UInt8, UInt32}, Float32}()
    protein_precursor_values = Dict{Tuple{String, Bool, UInt8}, Vector{Float32}}()
    for ((protein_name, target, entrap_id, precursor_idx), score_sum) in consensus_weight_sums
        protein_key = (protein_name, target, entrap_id)
        total_vote = get(protein_total_vote, protein_key, 0.0)
        total_vote > 0.0 || continue
        relative_weight = Float32(clamp(score_sum / total_vote, 0.0, 1.0))
        precursor_relative_weight[(protein_name, target, entrap_id, precursor_idx)] = relative_weight
        push!(get!(protein_precursor_values, protein_key, Float32[]), relative_weight)
    end

    protein_mean_relative_weight = Dict{Tuple{String, Bool, UInt8}, Float32}()
    protein_profiled_precursor_count = Dict{Tuple{String, Bool, UInt8}, Int32}()
    for (protein_key, relative_weights) in protein_precursor_values
        protein_mean_relative_weight[protein_key] = isempty(relative_weights) ? 0.0f0 : Float32(sum(relative_weights) / length(relative_weights))
        protein_profiled_precursor_count[protein_key] = Int32(length(relative_weights))
    end

    return (
        precursor_relative_weight = precursor_relative_weight,
        protein_mean_relative_weight = protein_mean_relative_weight,
        protein_profiled_precursor_count = protein_profiled_precursor_count
    )
end

function build_precursor_rank_consensus(psm_refs::Vector{PSMFileReference})
    return build_precursor_relative_weight_consensus(psm_refs)
end

"""
    group_psms_by_protein(df::DataFrame)

Transform PSMs into protein groups by aggregating peptides.
Returns a DataFrame with one row per protein group.
"""
function group_psms_by_protein(
    df::DataFrame;
    precursor_consensus::Union{Nothing, NamedTuple} = nothing
)
    if nrow(df) == 0
        # Return empty protein groups DataFrame with expected schema
        if precursor_consensus === nothing
            return DataFrame(
                protein_name = String[],
                target = Bool[],
                entrap_id = UInt8[],
                n_peptides = Int64[],
                peptide_list = String[],
                pg_score = Float32[],
                any_common_peps = Bool[],
                top_pep_weight = Float32[]
            )
        else
            return DataFrame(
                protein_name = String[],
                target = Bool[],
                entrap_id = UInt8[],
                n_peptides = Int64[],
                peptide_list = String[],
                pg_score = Float32[],
                any_common_peps = Bool[],
                top_pep_weight = Float32[],
                consensus_precursor_relative_weight_support = Float32[],
                consensus_precursor_relative_weight_enrichment = Float32[]
            )
        end
    end

    prob_col = _protein_group_probability_column(df)
    consensus_relative_weight = precursor_consensus === nothing ? nothing : precursor_consensus.precursor_relative_weight
    protein_mean_relative_weight = precursor_consensus === nothing ? nothing : precursor_consensus.protein_mean_relative_weight

    # Group by protein
    grouped = groupby(df, [:inferred_protein_group, :target, :entrap_id])
    
    # Aggregate to protein groups
    protein_groups = combine(grouped) do gdf
        quant_mask = gdf.use_for_protein_quant .== true
        quant_peptides, pg_score = _quant_peptides_and_pg_score(gdf, quant_mask, prob_col)
        n_peptides = length(quant_peptides)

        top_pep_weight = 0.0f0
        best_weight_by_peptide = Dict{String, Float32}()
        for i in 1:nrow(gdf)
            if quant_mask[i] != true
                continue
            end
            weight_val = Float32(gdf.weight[i])
            pep = gdf.sequence[i]
            if haskey(best_weight_by_peptide, pep)
                if weight_val > best_weight_by_peptide[pep]
                    best_weight_by_peptide[pep] = weight_val
                end
            else
                best_weight_by_peptide[pep] = weight_val
            end
        end

        if !isempty(best_weight_by_peptide)
            top_pep_weight = maximum(values(best_weight_by_peptide))
        end

        consensus_precursor_relative_weight_support = 0.0f0
        consensus_precursor_relative_weight_enrichment = 0.0f0
        if precursor_consensus !== nothing
            observed_precursors = Set{UInt32}()
            for i in eachindex(quant_mask)
                if quant_mask[i]
                    push!(observed_precursors, UInt32(gdf.precursor_idx[i]))
                end
            end

            if !isempty(observed_precursors)
                protein_key = (
                    String(gdf.inferred_protein_group[1]),
                    Bool(gdf.target[1]),
                    UInt8(gdf.entrap_id[1])
                )
                matched_precursors = 0

                for precursor_idx in observed_precursors
                    key = (protein_key[1], protein_key[2], protein_key[3], precursor_idx)
                    if haskey(consensus_relative_weight, key)
                        consensus_precursor_relative_weight_support += Float32(consensus_relative_weight[key])
                        matched_precursors += 1
                    end
                end

                mean_relative_weight = get(protein_mean_relative_weight, protein_key, 0.0f0)
                if matched_precursors > 0 && mean_relative_weight > 0.0f0
                    expected_random_support =
                        Float32(matched_precursors) * mean_relative_weight
                    consensus_precursor_relative_weight_enrichment =
                        consensus_precursor_relative_weight_support / expected_random_support
                end
            end
        end

        has_common = any(
            quant_mask .&
            (gdf.missed_cleavage .== 0) .&
            (gdf.Mox .== 0)
        )

        if precursor_consensus === nothing
            DataFrame(
                n_peptides = n_peptides,
                peptide_list = join(quant_peptides, ";"),
                pg_score = pg_score,
                any_common_peps = has_common,
                top_pep_weight = top_pep_weight
            )
        else
            DataFrame(
                n_peptides = n_peptides,
                peptide_list = join(quant_peptides, ";"),
                pg_score = pg_score,
                any_common_peps = has_common,
                top_pep_weight = top_pep_weight,
                consensus_precursor_relative_weight_support = consensus_precursor_relative_weight_support,
                consensus_precursor_relative_weight_enrichment = consensus_precursor_relative_weight_enrichment
            )
        end
    end
    
    # Rename the grouping column
    rename!(protein_groups, :inferred_protein_group => :protein_name)
    
    return protein_groups
end

"""
    filter_by_min_peptides(min_peptides::Int)

Filter protein groups that don't meet minimum peptide requirement.
"""
function filter_by_min_peptides(min_peptides::Int)
    desc = "filter_by_min_peptides(min=$min_peptides)"
    
    op = function(df)
        if hasproperty(df, :n_peptides)
            # For protein groups
            filter!(row -> row.n_peptides >= min_peptides, df)
        end
        return df
    end
    
    return desc => op
end

"""
    add_protein_features(protein_catalog::Dict)

Add protein-level features like peptide coverage.
"""
function add_protein_features(protein_catalog::Dict)
    desc = "add_protein_features"
    
    op = function(df)
        if !hasproperty(df, :protein_name)
            return df
        end
        
        n_rows = nrow(df)
        n_possible = Vector{Int64}(undef, n_rows)
        peptide_coverage = Vector{Float32}(undef, n_rows)
        grouped_catalog_cache = Dict{
            @NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},
            Set{String}
        }()
        
        for i in 1:n_rows
            key = (
                protein_name = df.protein_name[i],
                target = df.target[i],
                entrap_id = df.entrap_id[i]
            )

            possible_peptides = if haskey(protein_catalog, key)
                protein_catalog[key]
            elseif haskey(grouped_catalog_cache, key)
                grouped_catalog_cache[key]
            elseif occursin(';', key.protein_name)
                merged_peptides = Set{String}()
                for member in split(key.protein_name, ';')
                    member_key = (
                        protein_name = strip(member),
                        target = key.target,
                        entrap_id = key.entrap_id
                    )
                    if haskey(protein_catalog, member_key)
                        union!(merged_peptides, protein_catalog[member_key])
                    end
                end
                grouped_catalog_cache[key] = merged_peptides
                merged_peptides
            else
                Set{String}()
            end

            n_possible[i] = length(possible_peptides)

            if n_possible[i] > 0
                peptide_coverage[i] = Float32(df.n_peptides[i]) / Float32(n_possible[i])
            else
                peptide_coverage[i] = 0.0f0
            end
        end
        
        df.n_possible_peptides = n_possible
        df.peptide_coverage = peptide_coverage
        
        return df
    end
    
    return desc => op
end

#==========================================================
High-Level Interface
==========================================================#

"""
    perform_protein_inference_pipeline(psm_refs, output_folder, precursors, protein_catalog; kwargs...)

Perform protein inference using a composable pipeline approach.

# Arguments
- `psm_refs`: Vector of PSM file references
- `output_folder`: Directory for protein group output
- `precursors`: Library precursor information
- `protein_catalog`: Pre-computed protein-to-peptide mappings
- `min_peptides`: Minimum peptides per protein group (default: 2)

# Returns
- `pg_refs`: Vector of protein group file references
- `psm_to_pg_mapping`: Dictionary mapping PSM paths to protein group paths
"""
function perform_protein_inference_pipeline(
    psm_refs::Vector{PSMFileReference},
    output_folder::String,
    precursors::LibraryPrecursors,
    protein_catalog::Dict;
    min_peptides::Int = 2
)
    # Ensure output folder exists
    !isdir(output_folder) && mkpath(output_folder)
    
    # Build pre-inference pipeline
    pre_inference_pipeline = TransformPipeline() |>
        add_peptide_metadata(precursors) |>
        validate_peptide_data()
    
    # Process each file
    pg_refs = ProteinGroupFileReference[]
    psm_to_pg_mapping = Dict{String, String}()
    indexed_refs = collect(enumerate(psm_refs))
    
    for (idx, psm_ref) in ProgressBar(indexed_refs)
        if !exists(psm_ref)
            continue
        end
        
        # Step 1: Apply pre-inference pipeline to add necessary columns
        apply_pipeline!(psm_ref, pre_inference_pipeline)
        
        # Step 2: Run inference on the prepared PSMs
        prepared_df = load_dataframe(psm_ref)
        inference_result = apply_inference_to_dataframe(prepared_df, precursors)
        
        # Step 3: Update PSMs with inference results
        psm_update_pipeline = TransformPipeline() |>
            add_inferred_protein_column(inference_result) |> #Protein group assigned to each peptide
            add_quantification_flag(inference_result) #Whether or not the peptide should be used for protein quant and inference (non ambiguous)
        
        # Apply updates to the same PSM file (which now has necessary columns)
        apply_pipeline!(psm_ref, psm_update_pipeline)
        
    end

    precursor_consensus = build_precursor_relative_weight_consensus(psm_refs)

    for (idx, psm_ref) in ProgressBar(indexed_refs)
        if !exists(psm_ref)
            continue
        end

        updated_psms = load_dataframe(psm_ref)
        weight_calibration = estimate_weight_detection_model(updated_psms)
        @user_info "Weight protein coverage calibration file_idx=$(idx) n_unique_peptides=$(weight_calibration.n_unique_peptides) log_threshold=$(weight_calibration.log_threshold) sigma_log=$(weight_calibration.sigma_log) used_fallback=$(weight_calibration.used_fallback) profiled_rank_count=$(weight_calibration.profiled_rank_count) rank_profile_examples=$(weight_calibration.rank_profile_examples) used_rank_profile_fallback=$(weight_calibration.used_rank_profile_fallback)"

        protein_groups_df = group_psms_by_protein(
            updated_psms;
            precursor_consensus = precursor_consensus
        )

        post_inference_pipeline = TransformPipeline() |>
            filter_by_min_peptides(min_peptides) |>
            add_protein_features(protein_catalog) |>
            add_weight_observation_features(weight_calibration) |>
            add_pg_score_interaction_features() |>
            add_consensus_precursor_relative_weight_interaction_feature()

        initial_rows = nrow(protein_groups_df)
        for (desc, op) in post_inference_pipeline.operations
            protein_groups_df = op(protein_groups_df)
            #@debug_l1 "Pipeline operation on protein groups" operation=desc rows_before=initial_rows rows_after=nrow(protein_groups_df)
            initial_rows = nrow(protein_groups_df)
        end

        pg_filename = "protein_groups_$(lpad(idx, 3, '0')).arrow"
        pg_path = joinpath(output_folder, pg_filename)

        if nrow(protein_groups_df) > 0
            protein_groups_df[!, :file_idx] = fill(Int64(idx), nrow(protein_groups_df))
            writeArrow(pg_path, protein_groups_df)
            pg_ref = ProteinGroupFileReference(pg_path)
            push!(pg_refs, pg_ref)
            psm_to_pg_mapping[file_path(psm_ref)] = pg_path
        else
            @user_warn "No protein groups to write after filtering" file_idx=idx pg_path=pg_path
        end
    end
    
    return pg_refs, psm_to_pg_mapping
end
