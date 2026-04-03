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

const PREFIX_SHAPE_CONFIDENCE_PIVOT = 0.5f0
@inline function _shape_confidence(
    shape_strength::Float32,
    shape_confidence_scale::Float32
)::Float32
    shape_strength <= 0.0f0 && return 0.0f0
    tau = max(shape_confidence_scale, eps(Float32))
    return Float32(1.0 - exp(-Float64(shape_strength) / Float64(tau)))
end

@inline function _apply_shape_confidence(
    raw_shape::Float32,
    shape_confidence::Float32
)::Float32
    return Float32(shape_confidence * (raw_shape - PREFIX_SHAPE_CONFIDENCE_PIVOT))
end

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
        all_species = getProteomeIdentifiers(precursors)::AbstractVector{<:AbstractString}
        all_base_pep_ids = getBasePepId(precursors)::AbstractVector{UInt32}
        all_structural_mods = getStructuralMods(precursors)::AbstractVector{Union{Missing, String}}
        all_isotopic_mods = getIsotopicMods(precursors)::AbstractVector{Union{Missing, String}}
        all_pair_ids = getPairIdx(precursors)
        
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

        if !hasproperty(df, :species)
            species = Vector{String}(undef, n_rows)
            for i in 1:n_rows
                species[i] = join(sort(unique(split(coalesce(all_species[precursor_idx[i]], ""), ';'))), ';')
            end
            df.species = species
        end

        if !hasproperty(df, :base_pep_id)
            base_pep_ids = Vector{UInt32}(undef, n_rows)
            for i in 1:n_rows
                base_pep_ids[i] = all_base_pep_ids[precursor_idx[i]]
            end
            df.base_pep_id = base_pep_ids
        end

        if !hasproperty(df, :structural_mods)
            structural_mods = Vector{String}(undef, n_rows)
            for i in 1:n_rows
                structural_mods[i] = coalesce(all_structural_mods[precursor_idx[i]], "")
            end
            df.structural_mods = structural_mods
        end

        if !hasproperty(df, :isotopic_mods)
            isotopic_mods = Vector{String}(undef, n_rows)
            for i in 1:n_rows
                isotopic_mods[i] = coalesce(all_isotopic_mods[precursor_idx[i]], "")
            end
            df.isotopic_mods = isotopic_mods
        end

        if !hasproperty(df, :pair_id)
            pair_ids = Vector{UInt32}(undef, n_rows)
            for i in 1:n_rows
                pair_ids[i] = extract_pair_idx(all_pair_ids, precursor_idx[i])
            end
            df.pair_id = pair_ids
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
        required_cols = [
            :sequence,
            :accession_numbers,
            :is_decoy,
            :entrap_id,
            :precursor_idx,
            :base_pep_id,
            :structural_mods,
            :isotopic_mods
        ]
        
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

"""
    _median_abs_deviation(values)

Return the median absolute deviation of `values`.
"""
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

Estimate per-file abundance calibration for protein coverage surprise features.
Uses unique inferred peptides with valid positive integrated abundance values
when available, falling back to the legacy `weight` column otherwise.
"""
function estimate_weight_detection_model(df::DataFrame)
    max_rank = 12
    min_rank_support = 8
    default_model = (
        log_threshold = 0.0f0,
        rank_drop_profile = Float32[0.0f0],
        rank_scale_profile = Float32[1.0f0],
        profiled_rank_count = 1
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

        weight_val = Float64(_protein_rollup_weight_value(df, i))
        if !isfinite(weight_val) || weight_val <= 0.0
            continue
        end

        target_val = Bool(df.target[i])
        entrap_val = UInt8(df.entrap_id[i])
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
    if !isempty(pooled_residuals)
        mad_val = _median_abs_deviation(pooled_residuals)
        if isfinite(mad_val) && mad_val > 0.0
            sigma_log = clamp(mad_val / 0.6744897501960817, 0.25, 2.5)
        end
    end

    rank_drop_profile = fill(0.0f0, max_rank)
    rank_scale_profile = fill(Float32(sigma_log), max_rank)
    profiled_rank_count = 1
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
    end

    return (
        log_threshold = Float32(log_threshold),
        rank_drop_profile = rank_drop_profile[1:profiled_rank_count],
        rank_scale_profile = rank_scale_profile[1:profiled_rank_count],
        profiled_rank_count = profiled_rank_count
    )
end

"""
    add_weight_observation_features(calibration)

Add weight-based protein coverage features using an empirical top-rank exceedance model.
"""
function add_weight_observation_features(calibration::NamedTuple)
    desc = "add_weight_observation_features"

    op = function(df)
        n_rows = nrow(df)
        expected_excess_from_top = Vector{Float32}(undef, n_rows)
        coverage_log_ratio = Vector{Float32}(undef, n_rows)

        log_threshold = Float64(calibration.log_threshold)
        rank_drop_profile = Float64.(calibration.rank_drop_profile)
        rank_scale_profile = Float64.(calibration.rank_scale_profile)
        profiled_rank_count = max(Int(calibration.profiled_rank_count), 1)

        for i in 1:n_rows
            top_weight = Float64(df.top_pep_weight[i])
            N_total = max(Int(df.n_possible_peptides[i]), 0)
            k_obs = max(Int(df.n_peptides[i]), 0)

            if top_weight <= 0.0
                expected_excess_from_top[i] = 0.0f0
                coverage_log_ratio[i] = 0.0f0
                continue
            end

            log_top_weight = log(top_weight)
            effective_rank_count = min(max(N_total, 1), profiled_rank_count)
            if effective_rank_count <= 1
                expected_excess_from_top[i] = 0.0f0
                coverage_log_ratio[i] = 0.0f0
                continue
            end

            expected_excess = 0.0
            for rank_idx in 2:effective_rank_count
                rank_scale = rank_scale_profile[rank_idx]
                rank_scale = (isfinite(rank_scale) && rank_scale > 0.0) ? rank_scale : 1.0
                expected_log_weight = log_top_weight + rank_drop_profile[rank_idx]
                rank_z = (expected_log_weight - log_threshold) / rank_scale
                rank_excess = if rank_z > 20.0
                    rank_z
                elseif rank_z < -20.0
                    exp(rank_z)
                else
                    log1p(exp(rank_z))
                end
                expected_excess += rank_excess
            end

            observed_additional = min(max(k_obs - 1, 0), effective_rank_count - 1)

            expected_excess_from_top[i] = Float32(expected_excess)
            coverage_log_ratio[i] = Float32(log((observed_additional + 0.5) / (expected_excess + 0.5)))
        end

        df.expected_excess_from_top = expected_excess_from_top
        df.coverage_log_ratio = coverage_log_ratio
        return df
    end

    return desc => op
end

"""
    _protein_group_probability_column(df)

Choose the precursor probability column used for the initial protein roll-up.
Prefer direct run-specific precursor probabilities when present.
"""
function _protein_group_probability_column(df::AbstractDataFrame)
    if hasproperty(df, :prec_prob)
        return :prec_prob
    end
    if hasproperty(df, :MBR_boosted_prec_prob)
        return :MBR_boosted_prec_prob
    end
    error("Missing required precursor probability column for protein roll-up")
end

const PROTEIN_ROLLUP_PROB_EPS = 1.0f-6
const PROTEIN_ROLLUP_PRECURSOR_NONE_PSEUDOCOUNT = 0.01f0
const ProteinRollupPrecursorRow = @NamedTuple{
    precursor_idx::UInt32,
    pair_id::UInt32,
    base_pep_id::UInt32,
    sequence::String,
    charge::UInt8,
    structural_mods::String,
    isotopic_mods::String,
    pep::Float32,
    prob::Float32,
    score::Float32,
    best_weight::Float32
}
const ProteinRollupModifiedPeptideRow = @NamedTuple{
    base_pep_id::UInt32,
    sequence::String,
    structural_mods::String,
    isotopic_mods::String,
    log_none_sum::Float64,
    pep::Float32,
    prob::Float32,
    score::Float32,
    best_weight::Float32,
    precursor_count::Int32
}
const ProteinRollupPeptideRow = @NamedTuple{
    sequence::String,
    log_none_sum::Float64,
    pep::Float32,
    prob::Float32,
    score::Float32,
    best_weight::Float32,
    modified_peptide_count::Int32
}

@inline function _sanitize_rollup_probability(value)::Float32
    prob = Float32(value)
    if !isfinite(prob)
        return PROTEIN_ROLLUP_PROB_EPS
    end
    return clamp(prob, PROTEIN_ROLLUP_PROB_EPS, 1.0f0 - PROTEIN_ROLLUP_PROB_EPS)
end

@inline function _precursor_log_none_for_rollup(prob::Float32)::Float64
    adjusted_none_prob = clamp(
        (1.0f0 - prob) + PROTEIN_ROLLUP_PRECURSOR_NONE_PSEUDOCOUNT,
        PROTEIN_ROLLUP_PROB_EPS,
        1.0f0
    )
    return log(Float64(adjusted_none_prob))
end

"""
    _protein_rollup_quant_mask(df; q_value_threshold=0.01f0)

Return the precursor mask used for protein roll-up. Rows must be quant-eligible
and, when available, still pass both run-specific and global precursor q-value
thresholds.
"""
function _protein_rollup_quant_mask(
    df::AbstractDataFrame;
    q_value_threshold::Float32 = 0.01f0
)
    quant_mask = df.use_for_protein_quant .== true

    if hasproperty(df, :MBR_boosted_qval)
        quant_mask .&= coalesce.(df.MBR_boosted_qval .<= q_value_threshold, false)
    elseif hasproperty(df, :qval)
        quant_mask .&= coalesce.(df.qval .<= q_value_threshold, false)
    end

    if hasproperty(df, :MBR_boosted_global_qval)
        quant_mask .&= coalesce.(df.MBR_boosted_global_qval .<= q_value_threshold, false)
    elseif hasproperty(df, :global_qval)
        quant_mask .&= coalesce.(df.global_qval .<= q_value_threshold, false)
    end

    return quant_mask
end

@inline function _probability_from_log_none_sum(log_none_sum::Float64)::Float32
    if !isfinite(log_none_sum)
        return 1.0f0 - PROTEIN_ROLLUP_PROB_EPS
    end
    return Float32(clamp(-expm1(log_none_sum), Float64(PROTEIN_ROLLUP_PROB_EPS), 1.0 - Float64(PROTEIN_ROLLUP_PROB_EPS)))
end

@inline function _none_probability_from_log_none_sum(log_none_sum::Float64)::Float32
    if !isfinite(log_none_sum)
        return PROTEIN_ROLLUP_PROB_EPS
    end
    return Float32(clamp(exp(log_none_sum), Float64(PROTEIN_ROLLUP_PROB_EPS), 1.0))
end

@inline function _score_from_log_none_sum(log_none_sum::Float64)::Float32
    if isnan(log_none_sum)
        return 0.0f0
    end
    if isinf(log_none_sum)
        return log_none_sum < 0.0 ? Float32(-log(Float64(PROTEIN_ROLLUP_PROB_EPS))) : 0.0f0
    end
    return Float32(max(-log_none_sum, 0.0))
end

@inline function _protein_rollup_weight_value(
    gdf::AbstractDataFrame,
    row_idx::Int
)::Float32
    if hasproperty(gdf, :peak_area) && !ismissing(gdf.peak_area[row_idx])
        peak_area_val = Float32(gdf.peak_area[row_idx])
        if isfinite(peak_area_val) && peak_area_val > 0.0f0
            return peak_area_val
        end
    end
    return Float32(gdf.weight[row_idx])
end

"""
    _build_protein_rollup(gdf, quant_mask, prob_col)

Roll up precursor-level probability evidence to modified peptides, then
peptides, and return the peptide-level protein score plus reusable
intermediate rows.
"""
function _build_protein_rollup(
    gdf::AbstractDataFrame,
    quant_mask::AbstractVector{Bool},
    prob_col::Symbol
)
    if !hasproperty(gdf, prob_col)
        error("Missing required $(prob_col) column for protein roll-up")
    end
    if !hasproperty(gdf, :base_pep_id)
        error("Missing required :base_pep_id column for protein roll-up")
    end

    prob_by_precursor = Dict{UInt32, Float32}()
    best_weight_by_precursor = Dict{UInt32, Float32}()
    base_pep_id_by_precursor = Dict{UInt32, UInt32}()
    pair_id_by_precursor = Dict{UInt32, UInt32}()
    sequence_by_precursor = Dict{UInt32, String}()
    charge_by_precursor = Dict{UInt32, UInt8}()
    structural_mods_by_precursor = Dict{UInt32, String}()
    isotopic_mods_by_precursor = Dict{UInt32, String}()

    @inbounds for i in eachindex(quant_mask)
        quant_mask[i] || continue

        precursor_idx = UInt32(gdf.precursor_idx[i])
        prob_val = _sanitize_rollup_probability(gdf[!, prob_col][i])
        weight_val = _protein_rollup_weight_value(gdf, i)

        if haskey(prob_by_precursor, precursor_idx)
            prob_by_precursor[precursor_idx] = max(prob_by_precursor[precursor_idx], prob_val)
            if weight_val > best_weight_by_precursor[precursor_idx]
                best_weight_by_precursor[precursor_idx] = weight_val
            end
        else
            prob_by_precursor[precursor_idx] = prob_val
            best_weight_by_precursor[precursor_idx] = weight_val
            base_pep_id_by_precursor[precursor_idx] = UInt32(gdf.base_pep_id[i])
            pair_id_by_precursor[precursor_idx] =
                hasproperty(gdf, :pair_id) && !ismissing(gdf.pair_id[i]) ? UInt32(gdf.pair_id[i]) : zero(UInt32)
            sequence_by_precursor[precursor_idx] = String(gdf.sequence[i])
            charge_by_precursor[precursor_idx] =
                hasproperty(gdf, :charge) ? UInt8(gdf.charge[i]) : UInt8(0)
            structural_mods_by_precursor[precursor_idx] =
                hasproperty(gdf, :structural_mods) && !ismissing(gdf.structural_mods[i]) ? String(gdf.structural_mods[i]) : ""
            isotopic_mods_by_precursor[precursor_idx] =
                hasproperty(gdf, :isotopic_mods) && !ismissing(gdf.isotopic_mods[i]) ? String(gdf.isotopic_mods[i]) : ""
        end
    end

    precursor_rows = ProteinRollupPrecursorRow[]
    sizehint!(precursor_rows, length(prob_by_precursor))
    for precursor_idx in keys(prob_by_precursor)
        prob_val = prob_by_precursor[precursor_idx]
        none_prob_val = Float32(1.0f0 - prob_val)
        score_val = Float32(-log(none_prob_val))
        push!(precursor_rows, (
            precursor_idx = precursor_idx,
            pair_id = pair_id_by_precursor[precursor_idx],
            base_pep_id = base_pep_id_by_precursor[precursor_idx],
            sequence = sequence_by_precursor[precursor_idx],
            charge = charge_by_precursor[precursor_idx],
            structural_mods = structural_mods_by_precursor[precursor_idx],
            isotopic_mods = isotopic_mods_by_precursor[precursor_idx],
            pep = none_prob_val,
            prob = prob_val,
            score = score_val,
            best_weight = best_weight_by_precursor[precursor_idx]
        ))
    end

    if isempty(precursor_rows)
        return (
            precursor_rows = precursor_rows,
            modified_peptide_rows = ProteinRollupModifiedPeptideRow[],
            peptide_rows = ProteinRollupPeptideRow[],
            pg_score = 0.0f0,
            n_peptides = 0,
            peptide_list = String[],
            top_pep_weight = 0.0f0
        )
    end

    modified_log_none_sum = Dict{Tuple{UInt32, String, String}, Float64}()
    modified_best_weight = Dict{Tuple{UInt32, String, String}, Float32}()
    modified_sequence = Dict{Tuple{UInt32, String, String}, String}()
    modified_precursor_count = Dict{Tuple{UInt32, String, String}, Int32}()

    for precursor_row in precursor_rows
        mod_key = (
            precursor_row.base_pep_id,
            precursor_row.structural_mods,
            precursor_row.isotopic_mods
        )
        modified_log_none_sum[mod_key] = get(modified_log_none_sum, mod_key, 0.0) + _precursor_log_none_for_rollup(precursor_row.prob)
        modified_best_weight[mod_key] = max(
            get(modified_best_weight, mod_key, 0.0f0),
            precursor_row.best_weight
        )
        modified_sequence[mod_key] = precursor_row.sequence
        modified_precursor_count[mod_key] = get(modified_precursor_count, mod_key, Int32(0)) + Int32(1)
    end

    modified_peptide_rows = ProteinRollupModifiedPeptideRow[]
    sizehint!(modified_peptide_rows, length(modified_log_none_sum))
    for (mod_key, log_none_sum) in modified_log_none_sum
        push!(modified_peptide_rows, (
            base_pep_id = mod_key[1],
            sequence = modified_sequence[mod_key],
            structural_mods = mod_key[2],
            isotopic_mods = mod_key[3],
            log_none_sum = log_none_sum,
            pep = _none_probability_from_log_none_sum(log_none_sum),
            prob = _probability_from_log_none_sum(log_none_sum),
            score = _score_from_log_none_sum(log_none_sum),
            best_weight = modified_best_weight[mod_key],
            precursor_count = modified_precursor_count[mod_key]
        ))
    end
    sort!(modified_peptide_rows; by = row -> row.score, rev = true)

    peptide_log_none_sum = Dict{String, Float64}()
    peptide_best_weight = Dict{String, Float32}()
    peptide_modified_count = Dict{String, Int32}()

    for modified_row in modified_peptide_rows
        peptide_key = modified_row.sequence
        peptide_log_none_sum[peptide_key] = get(peptide_log_none_sum, peptide_key, 0.0) + modified_row.log_none_sum
        peptide_best_weight[peptide_key] = max(
            get(peptide_best_weight, peptide_key, 0.0f0),
            modified_row.best_weight
        )
        peptide_modified_count[peptide_key] = get(peptide_modified_count, peptide_key, Int32(0)) + Int32(1)
    end

    peptide_rows = ProteinRollupPeptideRow[]
    sizehint!(peptide_rows, length(peptide_log_none_sum))
    for (sequence, log_none_sum) in peptide_log_none_sum
        push!(peptide_rows, (
            sequence = sequence,
            log_none_sum = log_none_sum,
            pep = _none_probability_from_log_none_sum(log_none_sum),
            prob = _probability_from_log_none_sum(log_none_sum),
            score = _score_from_log_none_sum(log_none_sum),
            best_weight = peptide_best_weight[sequence],
            modified_peptide_count = peptide_modified_count[sequence]
        ))
    end
    sort!(peptide_rows; by = row -> row.score, rev = true)

    pg_score = isempty(peptide_rows) ? 0.0f0 : Float32(sum(row.score for row in peptide_rows))
    top_pep_weight = isempty(peptide_rows) ? 0.0f0 : maximum(row.best_weight for row in peptide_rows)

    return (
        precursor_rows = precursor_rows,
        modified_peptide_rows = modified_peptide_rows,
        peptide_rows = peptide_rows,
        pg_score = pg_score,
        n_peptides = length(peptide_rows),
        peptide_list = [row.sequence for row in peptide_rows],
        top_pep_weight = top_pep_weight
    )
end

const ConsensusRunVote = @NamedTuple{
    pg_score::Float32,
    run_order::Int64,
    normalized_precursors::Vector{Pair{UInt32, Float32}}
}

@inline function _consensus_runs_to_keep(n_runs::Int)::Int
    return n_runs <= 0 ? 0 : ceil(Int, sqrt(n_runs))
end

"""
    _consensus_vote_precedes(left_pg_score, left_run_order, right_vote)

Return whether a candidate run vote should rank ahead of `right_vote`.
"""
@inline function _consensus_vote_precedes(
    left_pg_score::Float32,
    left_run_order::Int64,
    right_vote::ConsensusRunVote
)::Bool
    return (left_pg_score > right_vote.pg_score) ||
           ((left_pg_score == right_vote.pg_score) && (left_run_order < right_vote.run_order))
end

"""
    _insert_top_consensus_vote!(run_votes, vote)

Insert `vote` into the in-place score-sorted run list for one protein.
"""
function _insert_top_consensus_vote!(
    run_votes::Vector{ConsensusRunVote},
    vote::ConsensusRunVote
)
    insert_idx = length(run_votes) + 1

    @inbounds for i in eachindex(run_votes)
        if _consensus_vote_precedes(vote.pg_score, vote.run_order, run_votes[i])
            insert_idx = i
            break
        end
    end

    insert!(run_votes, insert_idx, vote)

    return run_votes
end

"""
    build_precursor_consensus(psm_refs::Vector{PSMFileReference})

Build a dataset-level precursor relative-weight consensus within each inferred
protein group. Consensus weights are derived from quant precursors after
normalizing each precursor by the max precursor weight in that protein run.
Each run's vote is weighted by the run's current protein `pg_score`. After all
runs are collected for a protein, only the top `ceil(sqrt(n_runs))` runs by
`pg_score` are retained for the consensus.
"""
function build_precursor_consensus(
    psm_refs::Vector{PSMFileReference};
    q_value_threshold::Float32 = 0.01f0
)
    consensus_weight_sums = Dict{Tuple{String, Bool, UInt8, UInt32}, Float64}()
    protein_total_vote = Dict{Tuple{String, Bool, UInt8}, Float64}()
    protein_run_votes = Dict{Tuple{String, Bool, UInt8}, Vector{ConsensusRunVote}}()
    precursor_metadata = Dict{
        Tuple{String, Bool, UInt8, UInt32},
        @NamedTuple{
            pair_id::UInt32,
            sequence::String,
            charge::UInt8,
            structural_mods::String,
            isotopic_mods::String
        }
    }()
    precursor_run_counts = Dict{Tuple{String, Bool, UInt8, UInt32}, Int32}()

    for (run_order, psm_ref) in enumerate(psm_refs)
        if !exists(psm_ref)
            continue
        end

        df = load_dataframe(psm_ref)
        prob_col = _protein_group_probability_column(df)

        for gdf in groupby(df, [:inferred_protein_group, :target, :entrap_id])
            quant_mask = _protein_rollup_quant_mask(gdf; q_value_threshold = q_value_threshold)
            rollup = _build_protein_rollup(gdf, quant_mask, prob_col)

            if isempty(rollup.precursor_rows)
                continue
            end

            protein_name_val = gdf.inferred_protein_group[1]
            if ismissing(protein_name_val)
                continue
            end

            protein_name = String(protein_name_val)
            target = Bool(gdf.target[1])
            entrap_id = UInt8(gdf.entrap_id[1])

            for precursor_row in rollup.precursor_rows
                key = (protein_name, target, entrap_id, precursor_row.precursor_idx)
                precursor_run_counts[key] = get(precursor_run_counts, key, Int32(0)) + Int32(1)
                if !haskey(precursor_metadata, key)
                    precursor_metadata[key] = (
                        pair_id = precursor_row.pair_id,
                        sequence = precursor_row.sequence,
                        charge = precursor_row.charge,
                        structural_mods = precursor_row.structural_mods,
                        isotopic_mods = precursor_row.isotopic_mods
                    )
                end
            end

            run_max_weight = maximum(precursor_row.best_weight for precursor_row in rollup.precursor_rows)
            if !isfinite(run_max_weight) || run_max_weight <= 0.0f0
                continue
            end

            normalized_precursors = Pair{UInt32, Float32}[]
            sizehint!(normalized_precursors, length(rollup.precursor_rows))
            for precursor_row in rollup.precursor_rows
                push!(
                    normalized_precursors,
                    precursor_row.precursor_idx => Float32(clamp(precursor_row.best_weight / run_max_weight, 0.0f0, 1.0f0))
                )
            end
            protein_key = (protein_name, target, entrap_id)
            _insert_top_consensus_vote!(
                get!(
                    () -> sizehint!(ConsensusRunVote[], length(psm_refs)),
                    protein_run_votes,
                    protein_key
                ),
                (
                    pg_score = rollup.pg_score,
                    run_order = Int64(run_order),
                    normalized_precursors = normalized_precursors
                )
            )
        end
    end

    selected_run_votes = Dict{Tuple{String, Bool, UInt8}, Vector{ConsensusRunVote}}()
    for (protein_key, run_votes) in protein_run_votes
        runs_to_keep = _consensus_runs_to_keep(length(run_votes))
        selected_votes = run_votes[1:runs_to_keep]
        selected_run_votes[protein_key] = selected_votes
        for run_vote in selected_votes
            run_weight = Float64(run_vote.pg_score)
            protein_total_vote[protein_key] = get(protein_total_vote, protein_key, 0.0) + run_weight

            for precursor in run_vote.normalized_precursors
                key = (protein_key[1], protein_key[2], protein_key[3], precursor.first)
                consensus_weight_sums[key] = get(consensus_weight_sums, key, 0.0) + (run_weight * Float64(precursor.second))
            end
        end
    end

    relative_weight = Dict{Tuple{String, Bool, UInt8, UInt32}, Float32}()
    protein_precursor_values = Dict{Tuple{String, Bool, UInt8}, Vector{Float32}}()
    for ((protein_name, target, entrap_id, precursor_idx), score_sum) in consensus_weight_sums
        protein_key = (protein_name, target, entrap_id)
        total_vote = get(protein_total_vote, protein_key, 0.0)
        total_vote > 0.0 || continue
        precursor_relative_weight = Float32(clamp(score_sum / total_vote, 0.0, 1.0))
        relative_weight[(protein_name, target, entrap_id, precursor_idx)] = precursor_relative_weight
        push!(get!(protein_precursor_values, protein_key, Float32[]), precursor_relative_weight)
    end

    mean_relative_weight = Dict{Tuple{String, Bool, UInt8}, Float32}()
    profiled_precursor_count = Dict{Tuple{String, Bool, UInt8}, Int32}()
    shape_strength = Dict{Tuple{String, Bool, UInt8}, Float32}()
    precursors_by_protein = Dict{Tuple{String, Bool, UInt8}, Vector{Pair{UInt32, Float32}}}()
    for (protein_key, relative_weights) in protein_precursor_values
        mean_relative_weight[protein_key] = isempty(relative_weights) ? 0.0f0 : Float32(sum(relative_weights) / length(relative_weights))
        profiled_precursor_count[protein_key] = Int32(length(relative_weights))
    end
    for (protein_key, run_votes) in selected_run_votes
        total_pg_score = sum(Float64(run_vote.pg_score) for run_vote in run_votes)
        shape_strength[protein_key] = Float32(log1p(total_pg_score))
    end
    shape_strength_values = collect(values(shape_strength))
    shape_confidence_scale = isempty(shape_strength_values) ? 1.0f0 :
        Float32(max(Statistics.median(shape_strength_values), eps(Float32)))
    for ((protein_name, target, entrap_id, precursor_idx), precursor_relative_weight) in relative_weight
        push!(
            get!(precursors_by_protein, (protein_name, target, entrap_id), Pair{UInt32, Float32}[]),
            precursor_idx => precursor_relative_weight
        )
    end
    for precursor_weights in values(precursors_by_protein)
        sort!(precursor_weights; by = last, rev = true)
    end
    return (
        relative_weight = relative_weight,
        mean_relative_weight = mean_relative_weight,
        profiled_precursor_count = profiled_precursor_count,
        shape_strength = shape_strength,
        shape_confidence_scale = shape_confidence_scale,
        selected_run_votes = selected_run_votes,
        precursors_by_protein = precursors_by_protein,
        precursor_metadata = precursor_metadata,
        precursor_run_counts = precursor_run_counts
    )
end

"""
    _shape_consensus_inputs(precursor_rows, protein_key, precursor_consensus)

Prepare the observed precursor weights and protein key used for shape scoring.
"""
function _shape_consensus_inputs(
    precursor_rows::AbstractVector,
    protein_key::Tuple{String, Bool, UInt8},
    precursor_consensus::NamedTuple
)
    best_weight_by_precursor = Dict{UInt32, Float32}()
    for row in precursor_rows
        best_weight_by_precursor[row.precursor_idx] = row.best_weight
    end

    return (
        protein_key = protein_key,
        best_weight_by_precursor = best_weight_by_precursor
    )
end

"""
    _precursor_consensus_prefix_features(precursor_rows, protein_key, precursor_consensus)

Score how well the run matches the consensus precursor profile up to the weakest
observed consensus precursor.
"""
function _precursor_consensus_prefix_features(
    precursor_rows::AbstractVector,
    protein_key::Tuple{String, Bool, UInt8},
    precursor_consensus::NamedTuple
)
    shape_inputs = _shape_consensus_inputs(precursor_rows, protein_key, precursor_consensus)
    shape_protein_key = shape_inputs.protein_key
    best_weight_by_precursor = shape_inputs.best_weight_by_precursor

    profiled_precursor_count = Int(get(precursor_consensus.profiled_precursor_count, shape_protein_key, Int32(0)))
    shape_strength = hasproperty(precursor_consensus, :shape_strength) ?
        get(precursor_consensus.shape_strength, shape_protein_key, 0.0f0) : 0.0f0
    shape_confidence_scale = hasproperty(precursor_consensus, :shape_confidence_scale) ?
        precursor_consensus.shape_confidence_scale : 1.0f0
    shape_confidence = _shape_confidence(shape_strength, shape_confidence_scale)

    if isempty(best_weight_by_precursor)
        return (
            prefix_shape_raw = 0.0f0,
            prefix_shape = 0.0f0,
            threshold = 0.0f0,
            matched_precursors = 0,
            prefix_consensus_sum = 0.0f0,
            run_prefix_sum = 0.0f0,
            profiled_precursor_count = profiled_precursor_count,
            shape_strength = shape_strength,
            shape_confidence = shape_confidence
        )
    end

    consensus_relative_weight = precursor_consensus.relative_weight
    consensus_precursors = get(
        precursor_consensus.precursors_by_protein,
        shape_protein_key,
        Pair{UInt32, Float32}[]
    )
    isempty(consensus_precursors) && return (
        prefix_shape_raw = 0.0f0,
        prefix_shape = 0.0f0,
        threshold = 0.0f0,
        matched_precursors = 0,
        prefix_consensus_sum = 0.0f0,
        run_prefix_sum = 0.0f0,
        profiled_precursor_count = profiled_precursor_count,
        shape_strength = shape_strength,
        shape_confidence = shape_confidence
    )

    run_max_weight = maximum(values(best_weight_by_precursor))
    (!isfinite(run_max_weight) || run_max_weight <= 0.0f0) && return (
        prefix_shape_raw = 0.0f0,
        prefix_shape = 0.0f0,
        threshold = 0.0f0,
        matched_precursors = 0,
        prefix_consensus_sum = 0.0f0,
        run_prefix_sum = 0.0f0,
        profiled_precursor_count = profiled_precursor_count,
        shape_strength = shape_strength,
        shape_confidence = shape_confidence
    )

    observed_consensus_weights = Float32[]
    matched_precursors = 0
    for precursor_idx in keys(best_weight_by_precursor)
        c = get(consensus_relative_weight, (shape_protein_key[1], shape_protein_key[2], shape_protein_key[3], precursor_idx), 0.0f0)
        if c > 0.0f0
            push!(observed_consensus_weights, c)
            matched_precursors += 1
        end
    end
    isempty(observed_consensus_weights) && return (
        prefix_shape_raw = 0.0f0,
        prefix_shape = 0.0f0,
        threshold = 0.0f0,
        matched_precursors = 0,
        prefix_consensus_sum = 0.0f0,
        run_prefix_sum = 0.0f0,
        profiled_precursor_count = profiled_precursor_count,
        shape_strength = shape_strength,
        shape_confidence = shape_confidence
    )

    threshold = minimum(observed_consensus_weights)
    prefix_consensus_sum = 0.0f0
    run_prefix_sum = 0.0f0

    for precursor in consensus_precursors
        precursor.second < threshold && continue
        run_relative_weight = get(best_weight_by_precursor, precursor.first, 0.0f0) / run_max_weight
        prefix_consensus_sum += precursor.second
        run_prefix_sum += run_relative_weight
    end

    prefix_consensus_sum <= 0.0f0 && return (
        prefix_shape_raw = 0.0f0,
        prefix_shape = 0.0f0,
        threshold = threshold,
        matched_precursors = matched_precursors,
        prefix_consensus_sum = 0.0f0,
        run_prefix_sum = run_prefix_sum,
        profiled_precursor_count = profiled_precursor_count,
        shape_strength = shape_strength,
        shape_confidence = shape_confidence
    )

    if run_prefix_sum <= 0.0f0
        return (
            prefix_shape_raw = 0.0f0,
            prefix_shape = 0.0f0,
            threshold = threshold,
            matched_precursors = matched_precursors,
            prefix_consensus_sum = prefix_consensus_sum,
            run_prefix_sum = run_prefix_sum,
            profiled_precursor_count = profiled_precursor_count,
            shape_strength = shape_strength,
            shape_confidence = shape_confidence
        )
    end

    shape_dot = 0.0f0
    consensus_norm_sq = 0.0f0
    run_norm_sq = 0.0f0
    for precursor in consensus_precursors
        precursor.second < threshold && continue
        run_relative_weight = get(best_weight_by_precursor, precursor.first, 0.0f0) / run_max_weight
        consensus_norm = precursor.second / prefix_consensus_sum
        run_norm = run_relative_weight / run_prefix_sum
        shape_dot += run_norm * consensus_norm
        consensus_norm_sq += consensus_norm * consensus_norm
        run_norm_sq += run_norm * run_norm
    end

    if consensus_norm_sq <= 0.0f0 || run_norm_sq <= 0.0f0
        prefix_shape_raw = 0.0f0
    else
        prefix_shape_raw = Float32(clamp(shape_dot / sqrt(consensus_norm_sq * run_norm_sq), 0.0f0, 1.0f0))
    end
    return (
        prefix_shape_raw = prefix_shape_raw,
        prefix_shape = _apply_shape_confidence(prefix_shape_raw, shape_confidence),
        threshold = threshold,
        matched_precursors = matched_precursors,
        prefix_consensus_sum = prefix_consensus_sum,
        run_prefix_sum = run_prefix_sum,
        profiled_precursor_count = profiled_precursor_count,
        shape_strength = shape_strength,
        shape_confidence = shape_confidence
    )
end

@inline function _protein_group_contains_member(
    protein_name::AbstractString,
    member_name::AbstractString
)::Bool
    protein_name == member_name && return true
    occursin(';', protein_name) || return false
    return any(==(member_name), split(protein_name, ';'))
end

@inline function _precursor_debug_label(
    protein_key::Tuple{String, Bool, UInt8},
    precursor_idx::UInt32,
    precursor_consensus::NamedTuple
)::String
    metadata = get(
        precursor_consensus.precursor_metadata,
        (protein_key[1], protein_key[2], protein_key[3], precursor_idx),
        nothing
    )
    metadata === nothing && return "precursor_idx=$(precursor_idx)"

    mods = String[]
    !isempty(metadata.structural_mods) && push!(mods, metadata.structural_mods)
    !isempty(metadata.isotopic_mods) && push!(mods, metadata.isotopic_mods)
    mod_suffix = isempty(mods) ? "" : " [" * join(mods, " | ") * "]"

    return "$(metadata.sequence)/$(metadata.charge)$(mod_suffix) precursor_idx=$(precursor_idx)"
end

function _collect_precursor_consensus_debug_rows(
    df::AbstractDataFrame,
    file_idx::Int64,
    precursor_consensus::NamedTuple;
    protein_names::AbstractVector{String},
    q_value_threshold::Float32 = 0.01f0,
    file_idx_to_name::Union{Nothing, AbstractDict{Int64, String}} = nothing
)
    isempty(protein_names) && return NamedTuple[]

    prob_col = _protein_group_probability_column(df)
    debug_rows = NamedTuple[]

    for gdf in groupby(df, [:inferred_protein_group, :target, :entrap_id])
        protein_name_val = gdf.inferred_protein_group[1]
        ismissing(protein_name_val) && continue

        protein_name = String(protein_name_val)
        any(name -> _protein_group_contains_member(protein_name, name), protein_names) || continue

        quant_mask = _protein_rollup_quant_mask(gdf; q_value_threshold = q_value_threshold)
        rollup = _build_protein_rollup(gdf, quant_mask, prob_col)
        isempty(rollup.precursor_rows) && continue

        protein_key = (
            protein_name,
            Bool(gdf.target[1]),
            UInt8(gdf.entrap_id[1])
        )
        prefix_features = _precursor_consensus_prefix_features(
            rollup.precursor_rows,
            protein_key,
            precursor_consensus
        )
        shape_inputs = _shape_consensus_inputs(rollup.precursor_rows, protein_key, precursor_consensus)
        best_weight_by_precursor = shape_inputs.best_weight_by_precursor
        shape_protein_key = shape_inputs.protein_key

        run_max_weight = maximum(values(best_weight_by_precursor))
        threshold = prefix_features.threshold
        consensus_precursors = get(
            precursor_consensus.precursors_by_protein,
            shape_protein_key,
            Pair{UInt32, Float32}[]
        )

        observed_precursors = NamedTuple[]
        for precursor_row in sort(collect(rollup.precursor_rows); by = row -> row.best_weight, rev = true)
            consensus_weight = get(
                precursor_consensus.relative_weight,
                (shape_protein_key[1], shape_protein_key[2], shape_protein_key[3], precursor_row.precursor_idx),
                0.0f0
            )
            run_relative_weight = run_max_weight > 0.0f0 ?
                Float32(precursor_row.best_weight / run_max_weight) : 0.0f0
            push!(
                observed_precursors,
                (
                    precursor_idx = precursor_row.precursor_idx,
                    label = _precursor_debug_label(protein_key, precursor_row.precursor_idx, precursor_consensus),
                    best_weight = precursor_row.best_weight,
                    run_relative_weight = run_relative_weight,
                    consensus_weight = consensus_weight,
                    in_prefix = consensus_weight >= threshold
                )
            )
        end

        prefix_precursors = NamedTuple[]
        for precursor in consensus_precursors
            precursor.second < threshold && continue
            push!(
                prefix_precursors,
                (
                    precursor_idx = precursor.first,
                    label = _precursor_debug_label(protein_key, precursor.first, precursor_consensus),
                    consensus_weight = precursor.second,
                    run_relative_weight = run_max_weight > 0.0f0 ?
                        Float32(get(best_weight_by_precursor, precursor.first, 0.0f0) / run_max_weight) : 0.0f0,
                    observed = haskey(best_weight_by_precursor, precursor.first),
                    run_count = Int(get(
                        precursor_consensus.precursor_run_counts,
                        (shape_protein_key[1], shape_protein_key[2], shape_protein_key[3], precursor.first),
                        Int32(0)
                    ))
                )
            )
        end

        push!(
            debug_rows,
            (
                protein_key = protein_key,
                file_idx = file_idx,
                run_name = isnothing(file_idx_to_name) ?
                    "file_$(file_idx)" :
                    get(file_idx_to_name, file_idx, "file_$(file_idx)"),
                n_peptides = rollup.n_peptides,
                peptide_list = join(rollup.peptide_list, ";"),
                pg_score = rollup.pg_score,
                top_pep_weight = rollup.top_pep_weight,
                matched_precursors = prefix_features.matched_precursors,
                threshold = prefix_features.threshold,
                prefix_consensus_sum = prefix_features.prefix_consensus_sum,
                run_prefix_sum = prefix_features.run_prefix_sum,
                prefix_shape_raw = prefix_features.prefix_shape_raw,
                shape_strength = prefix_features.shape_strength,
                shape_confidence = prefix_features.shape_confidence,
                prefix_shape = prefix_features.prefix_shape,
                observed_precursors = observed_precursors,
                prefix_precursors = prefix_precursors
            )
        )
    end

    return debug_rows
end

function _write_precursor_consensus_debug_file(
    qc_folder::String,
    precursor_consensus::NamedTuple,
    debug_rows::AbstractVector{<:NamedTuple};
    protein_name::String,
    file_idx_to_name::Union{Nothing, AbstractDict{Int64, String}} = nothing
)
    isdir(qc_folder) || mkpath(qc_folder)
    debug_path = joinpath(qc_folder, "$(protein_name)_precursor_consensus_debug.txt")

    matching_keys = sort!(
        unique([row.protein_key for row in debug_rows]);
        by = key -> (key[1], key[2], key[3])
    )

    open(debug_path, "w") do io
        println(io, "Protein precursor consensus debug")
        println(io, "requested_protein=$(protein_name)")
        println(io, "shape_confidence_scale=$(precursor_consensus.shape_confidence_scale)")
        println(io, "shape_confidence_pivot=$(PREFIX_SHAPE_CONFIDENCE_PIVOT)")
        println(io)

        if isempty(matching_keys)
            println(io, "No matching protein groups found for $(protein_name).")
            return
        end

        for protein_key in matching_keys
            shape_strength_value = get(precursor_consensus.shape_strength, protein_key, 0.0f0)
            shape_confidence_value = _shape_confidence(
                shape_strength_value,
                precursor_consensus.shape_confidence_scale
            )
            println(io, "================================================================================")
            println(
                io,
                "protein_name=$(repr(protein_key[1])) target=$(protein_key[2]) entrap_id=$(protein_key[3])"
            )
            println(
                io,
                "shape_strength=$(shape_strength_value) shape_confidence=$(shape_confidence_value)"
            )
            println(io)

            println(io, "Selected consensus source runs")
            selected_votes = get(precursor_consensus.selected_run_votes, protein_key, ConsensusRunVote[])
            if isempty(selected_votes)
                println(io, "  (none)")
            else
                for (vote_idx, vote) in enumerate(selected_votes)
                    run_name = isnothing(file_idx_to_name) ?
                        "file_$(vote.run_order)" :
                        get(file_idx_to_name, Int64(vote.run_order), "file_$(vote.run_order)")
                    println(
                        io,
                        "  [$(vote_idx)] run_order=$(vote.run_order) run_name=$(repr(run_name)) pg_score=$(vote.pg_score)"
                    )
                    for precursor in sort(collect(vote.normalized_precursors); by = last, rev = true)
                        println(
                            io,
                            "      normalized_weight=$(precursor.second) " *
                            _precursor_debug_label(protein_key, precursor.first, precursor_consensus)
                        )
                    end
                end
            end
            println(io)

            println(io, "Consensus precursors")
            consensus_precursors = get(
                precursor_consensus.precursors_by_protein,
                protein_key,
                Pair{UInt32, Float32}[]
            )
            if isempty(consensus_precursors)
                println(io, "  (none)")
            else
                for precursor in consensus_precursors
                    run_count = Int(get(
                        precursor_consensus.precursor_run_counts,
                        (protein_key[1], protein_key[2], protein_key[3], precursor.first),
                        Int32(0)
                    ))
                    println(
                        io,
                        "  consensus_weight=$(precursor.second) run_count=$(run_count) " *
                        _precursor_debug_label(protein_key, precursor.first, precursor_consensus)
                    )
                end
            end
            println(io)

            println(io, "Per-run shape diagnostics")
            protein_rows = sort!(
                [row for row in debug_rows if row.protein_key == protein_key];
                by = row -> row.file_idx
            )
            if isempty(protein_rows)
                println(io, "  (none)")
            else
                for row in protein_rows
                    println(io, "--------------------------------------------------------------------------------")
                    println(
                        io,
                        "file_idx=$(row.file_idx) run_name=$(repr(row.run_name)) " *
                        "n_peptides=$(row.n_peptides) peptide_list=$(repr(row.peptide_list)) " *
                        "pg_score=$(row.pg_score) top_pep_weight=$(row.top_pep_weight)"
                    )
                    println(
                        io,
                        "matched_precursors=$(row.matched_precursors) tau=$(row.threshold) " *
                        "prefix_consensus_sum=$(row.prefix_consensus_sum) run_prefix_sum=$(row.run_prefix_sum)"
                    )
                    println(
                        io,
                        "prefix_shape_raw=$(row.prefix_shape_raw) shape_strength=$(row.shape_strength) " *
                        "shape_confidence=$(row.shape_confidence) prefix_shape=$(row.prefix_shape)"
                    )
                    println(io, "  Observed precursors")
                    if isempty(row.observed_precursors)
                        println(io, "    (none)")
                    else
                        for precursor in row.observed_precursors
                            println(
                                io,
                                "    best_weight=$(precursor.best_weight) run_relative_weight=$(precursor.run_relative_weight) " *
                                "consensus_weight=$(precursor.consensus_weight) in_prefix=$(precursor.in_prefix) " *
                                precursor.label
                            )
                        end
                    end
                    println(io, "  Prefix precursors (consensus_weight >= tau)")
                    if isempty(row.prefix_precursors)
                        println(io, "    (none)")
                    else
                        for precursor in row.prefix_precursors
                            println(
                                io,
                                "    consensus_weight=$(precursor.consensus_weight) run_relative_weight=$(precursor.run_relative_weight) " *
                                "observed=$(precursor.observed) run_count=$(precursor.run_count) " *
                                precursor.label
                            )
                        end
                    end
                end
            end
            println(io)
        end
    end

    return debug_path
end

"""
    group_psms_by_protein(df::DataFrame)

Transform PSMs into protein groups by aggregating peptides.
Returns a DataFrame with one row per protein group.
"""
function group_psms_by_protein(
    df::DataFrame;
    precursor_consensus::NamedTuple,
    q_value_threshold::Float32 = 0.01f0
)
    if nrow(df) == 0
        return DataFrame(
            protein_name = String[],
            species = String[],
            target = Bool[],
            entrap_id = UInt8[],
            n_peptides = Int64[],
            peptide_list = String[],
            pg_score = Float32[],
            any_common_peps = Bool[],
            top_pep_weight = Float32[],
            precursor_consensus_prefix_shape = Float32[],
            pg_score_x_precursor_consensus_prefix_shape = Float32[]
        )
    end

    df = df[.!ismissing.(df.inferred_protein_group), :]
    if nrow(df) == 0
        return DataFrame(
            protein_name = String[],
            species = String[],
            target = Bool[],
            entrap_id = UInt8[],
            n_peptides = Int64[],
            peptide_list = String[],
            pg_score = Float32[],
            any_common_peps = Bool[],
            top_pep_weight = Float32[],
            precursor_consensus_prefix_shape = Float32[],
            pg_score_x_precursor_consensus_prefix_shape = Float32[]
        )
    end

    prob_col = _protein_group_probability_column(df)
    # Group by protein
    grouped = groupby(df, [:inferred_protein_group, :target, :entrap_id])
    
    # Aggregate to protein groups
    protein_groups = combine(grouped) do gdf
        quant_mask = _protein_rollup_quant_mask(gdf; q_value_threshold = q_value_threshold)
        rollup = _build_protein_rollup(gdf, quant_mask, prob_col)
        n_peptides = rollup.n_peptides
        pg_score = rollup.pg_score
        top_pep_weight = rollup.top_pep_weight
        species = if hasproperty(gdf, :species)
            join(sort!(unique!(collect(skipmissing(String.(gdf.species))))), ';')
        else
            ""
        end

        precursor_consensus_prefix_shape = 0.0f0
        if !isempty(rollup.precursor_rows)
            protein_key = (
                String(gdf.inferred_protein_group[1]),
                Bool(gdf.target[1]),
                UInt8(gdf.entrap_id[1])
            )
            prefix_features = _precursor_consensus_prefix_features(
                rollup.precursor_rows,
                protein_key,
                precursor_consensus
            )
            precursor_consensus_prefix_shape = prefix_features.prefix_shape
        end

        has_common = any(
            quant_mask .&
            (gdf.missed_cleavage .== 0) .&
            (gdf.Mox .== 0)
        )

        DataFrame(
            species = species,
            n_peptides = n_peptides,
            peptide_list = join(rollup.peptide_list, ";"),
            pg_score = pg_score,
            any_common_peps = has_common,
            top_pep_weight = top_pep_weight,
            precursor_consensus_prefix_shape = precursor_consensus_prefix_shape,
            pg_score_x_precursor_consensus_prefix_shape = pg_score * precursor_consensus_prefix_shape
        )
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
    grouped_catalog_cache = Dict{
        @NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},
        Set{String}
    }()
    
    op = function(df)
        if !hasproperty(df, :protein_name)
            return df
        end
        
        n_rows = nrow(df)
        n_possible = Vector{Int64}(undef, n_rows)
        peptide_coverage = Vector{Float32}(undef, n_rows)
        peptide_coverage_logit = Vector{Float32}(undef, n_rows)
        
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

            observed_n = Int(df.n_peptides[i])
            if observed_n > n_possible[i]
                error(
                    "Protein feature count inconsistency: " *
                    "protein_name=$(repr(key.protein_name)) " *
                    "target=$(key.target) " *
                    "entrap_id=$(key.entrap_id) " *
                    "n_peptides=$(observed_n) " *
                    "n_possible_peptides=$(n_possible[i])"
                )
            end

            if n_possible[i] > 0
                peptide_coverage[i] = Float32(df.n_peptides[i]) / Float32(n_possible[i])
                peptide_coverage_logit[i] = smoothed_coverage_logit(df.n_peptides[i], n_possible[i])
            else
                peptide_coverage[i] = 0.0f0
                peptide_coverage_logit[i] = 0.0f0
            end
        end
        
        df.n_possible_peptides = n_possible
        df.peptide_coverage = peptide_coverage
        df.peptide_coverage_logit = peptide_coverage_logit
        
        return df
    end
    
    return desc => op
end

"""
    _precursor_consensus_report_bucket(n_peptides)

Bucket proteins by support size for the precursor consensus QC report.
"""
@inline function _precursor_consensus_report_bucket(n_peptides::Int)
    if n_peptides <= 1
        return :one_peptide
    elseif n_peptides <= 3
        return :two_to_three_peptides
    end
    return :four_plus_peptides
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
    min_peptides::Int = 2,
    q_value_threshold::Float32 = 0.01f0,
    qc_folder::Union{Nothing, String} = nothing,
    file_idx_to_name::Union{Nothing, AbstractDict{Int64, String}} = nothing,
    consensus_debug_protein_names::Vector{String} = String[]
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

    @user_info "Annotating passing PSM files with inferred protein groups and protein-quant flags"
    
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

    precursor_consensus = build_precursor_consensus(psm_refs; q_value_threshold = q_value_threshold)
    consensus_debug_rows = NamedTuple[]

    @user_info "Building per-run protein group tables and protein scoring features"

    for (idx, psm_ref) in ProgressBar(indexed_refs)
        if !exists(psm_ref)
            continue
        end

        updated_psms = load_dataframe(psm_ref)
        weight_calibration = estimate_weight_detection_model(updated_psms)
        if !isempty(consensus_debug_protein_names) && !isnothing(qc_folder)
            append!(
                consensus_debug_rows,
                _collect_precursor_consensus_debug_rows(
                    updated_psms,
                    Int64(idx),
                    precursor_consensus;
                    protein_names = consensus_debug_protein_names,
                    q_value_threshold = q_value_threshold,
                    file_idx_to_name = file_idx_to_name
                )
            )
        end

        protein_groups_df = group_psms_by_protein(
            updated_psms;
            precursor_consensus = precursor_consensus,
            q_value_threshold = q_value_threshold
        )

        post_inference_pipeline = TransformPipeline() |>
            filter_by_min_peptides(min_peptides) |>
            add_protein_features(protein_catalog) |>
            add_weight_observation_features(weight_calibration)

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

    if !isempty(consensus_debug_protein_names) && !isnothing(qc_folder)
        for protein_name in consensus_debug_protein_names
            _write_precursor_consensus_debug_file(
                qc_folder,
                precursor_consensus,
                filter(row -> _protein_group_contains_member(row.protein_key[1], protein_name), consensus_debug_rows);
                protein_name = protein_name,
                file_idx_to_name = file_idx_to_name
            )
        end
    end

    return pg_refs, psm_to_pg_mapping
end
