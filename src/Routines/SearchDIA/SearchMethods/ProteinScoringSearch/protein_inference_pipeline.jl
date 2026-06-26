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

@inline function _insert_top_log_peak_area!(
    top_log_peak_areas::Vector{Float64},
    log_peak_area::Float64,
    max_rank::Int
)::Vector{Float64}
    max_rank <= 0 && return top_log_peak_areas

    if (length(top_log_peak_areas) == max_rank) && (log_peak_area <= top_log_peak_areas[end])
        return top_log_peak_areas
    end

    insert_idx = length(top_log_peak_areas) + 1
    @inbounds for i in eachindex(top_log_peak_areas)
        if log_peak_area > top_log_peak_areas[i]
            insert_idx = i
            break
        end
    end

    insert!(top_log_peak_areas, insert_idx, log_peak_area)
    if length(top_log_peak_areas) > max_rank
        pop!(top_log_peak_areas)
    end

    return top_log_peak_areas
end

"""
    estimate_peak_area_detection_model(df::DataFrame)

Estimate per-file peak-area calibration for protein coverage surprise features.
Uses unique inferred peptides with valid positive integrated peak areas.
"""
function estimate_peak_area_detection_model(df::DataFrame)
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

    best_peak_area_by_protein_peptide = Dict{Tuple{String, Bool, UInt8, String}, Float64}()

    for i in 1:n_rows
        if df.use_for_protein_quant[i] != true
            continue
        end

        protein_name = df.inferred_protein_group[i]
        if ismissing(protein_name)
            continue
        end

        peak_area_val = Float64(df.peak_area[i])

        target_val = Bool(df.target[i])
        entrap_val = UInt8(df.entrap_id[i])
        key = (String(protein_name), target_val, entrap_val, String(df.sequence[i]))
        if haskey(best_peak_area_by_protein_peptide, key)
            if peak_area_val > best_peak_area_by_protein_peptide[key]
                best_peak_area_by_protein_peptide[key] = peak_area_val
            end
        else
            best_peak_area_by_protein_peptide[key] = peak_area_val
        end
    end

    if isempty(best_peak_area_by_protein_peptide)
        return default_model
    end

    log_peak_areas = Float64[log(peak_area) for peak_area in values(best_peak_area_by_protein_peptide)]
    log_threshold = Float64(Statistics.quantile(log_peak_areas, 0.05))

    protein_to_log_peak_areas = Dict{Tuple{String, Bool, UInt8}, Vector{Float64}}()
    protein_to_top_log_peak_areas = Dict{Tuple{String, Bool, UInt8}, Vector{Float64}}()
    for ((protein_name, target_val, entrap_val, _), peak_area) in best_peak_area_by_protein_peptide
        protein_key = (protein_name, target_val, entrap_val)
        log_peak_area = log(peak_area)
        push!(get!(protein_to_log_peak_areas, protein_key, Float64[]), log_peak_area)
        _insert_top_log_peak_area!(
            get!(() -> sizehint!(Float64[], max_rank), protein_to_top_log_peak_areas, protein_key),
            log_peak_area,
            max_rank
        )
    end

    pooled_residuals = Float64[]
    for log_vals in values(protein_to_log_peak_areas)
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

    for ((_, target_val, _), top_log_peak_areas) in protein_to_top_log_peak_areas
        target_val || continue
        length(top_log_peak_areas) >= 2 || continue

        top_log_peak_area = top_log_peak_areas[1]
        local_max_rank = length(top_log_peak_areas)
        for rank_idx in 2:local_max_rank
            push!(rank_drop_samples[rank_idx], top_log_peak_areas[rank_idx] - top_log_peak_area)
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
    add_peak_area_observation_features(calibration::NamedTuple)

Add peak-area-based coverage surprise features to the protein-group table using
the per-run calibration from `estimate_peak_area_detection_model(...)`.
"""
function add_peak_area_observation_features(calibration::NamedTuple)
    desc = "add_peak_area_observation_features"

    op = function(df)
        n_rows = nrow(df)
        expected_excess_from_top_peak_area = Vector{Float32}(undef, n_rows)
        coverage_log_ratio = Vector{Float32}(undef, n_rows)

        log_threshold = Float64(calibration.log_threshold)
        rank_drop_profile = Float64.(calibration.rank_drop_profile)
        rank_scale_profile = Float64.(calibration.rank_scale_profile)
        profiled_rank_count = max(Int(calibration.profiled_rank_count), 1)

        for i in 1:n_rows
            top_peak_area = Float64(df.top_pep_peak_area[i])
            N_total = max(Int(df.n_possible_peptides[i]), 0)
            k_obs = max(Int(df.n_peptides[i]), 0)

            if top_peak_area <= 0.0
                expected_excess_from_top_peak_area[i] = 0.0f0
                coverage_log_ratio[i] = 0.0f0
                continue
            end

            effective_rank_count = min(max(N_total, 1), profiled_rank_count)
            if effective_rank_count <= 1
                expected_excess_from_top_peak_area[i] = 0.0f0
                coverage_log_ratio[i] = 0.0f0
                continue
            end

            log_top_peak_area = log(top_peak_area)
            expected_excess = 0.0
            for rank_idx in 2:effective_rank_count
                rank_scale = rank_scale_profile[rank_idx]
                rank_scale = (isfinite(rank_scale) && rank_scale > 0.0) ? rank_scale : 1.0
                expected_log_peak_area = log_top_peak_area + rank_drop_profile[rank_idx]
                rank_z = (expected_log_peak_area - log_threshold) / rank_scale
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

            expected_excess_from_top_peak_area[i] = Float32(expected_excess)
            coverage_log_ratio[i] = Float32(log((observed_additional + 0.5) / (expected_excess + 0.5)))
        end

        df.expected_excess_from_top_peak_area = expected_excess_from_top_peak_area
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
const PROTEIN_ROLLUP_PRECURSOR_NONE_PSEUDOCOUNT = 0.001f0

const ProteinRollupPrecursorRow = @NamedTuple{
    precursor_idx::UInt32,
    base_pep_id::UInt32,
    sequence::String,
    charge::UInt8,
    structural_mods::String,
    isotopic_mods::String,
    pep::Float32,
    prob::Float32,
    score::Float32,
    best_peak_area::Float32
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
    best_peak_area::Float32,
    precursor_count::Int32
}
const ProteinRollupPeptideRow = @NamedTuple{
    sequence::String,
    log_none_sum::Float64,
    pep::Float32,
    prob::Float32,
    score::Float32,
    best_peak_area::Float32,
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

function _protein_mbr_summary(gdf::AbstractDataFrame, quant_mask::AbstractVector{Bool})
    if !hasproperty(gdf, :mbr_recovered)
        return (recovered_peptides = Int64(0), only_protein = false)
    end

    recovered_precursors = 0
    retained_precursors = 0
    recovered_peptides = Set{String}()

    @inbounds for i in eachindex(quant_mask)
        quant_mask[i] || continue
        retained_precursors += 1

        if Bool(gdf.mbr_recovered[i])
            recovered_precursors += 1
            push!(recovered_peptides, String(gdf.sequence[i]))
        end
    end

    return (
        recovered_peptides = Int64(length(recovered_peptides)),
        only_protein = retained_precursors > 0 && recovered_precursors == retained_precursors
    )
end

"""
Reusable scratch for `_build_protein_rollup`. One dedup `Dict` per roll-up level
(precursor_idx / mod-key / sequence -> row index) plus parallel arrays for the
fields, all reused across protein groups via `_reset_rollup_scratch!` (so the
14 `Dict`s previously allocated per group become 3 reused dedup `Dict`s + reused
arrays). The dedup `Dict`s are used only for lookup, never iterated, so output
order is the deterministic first-seen (array) order — independent of `Dict`
capacity, hence stable across runs.
"""
mutable struct ProteinRollupScratch
    # precursor level (dedup by precursor_idx)
    prec_pos::Dict{UInt32, Int}
    prec_id::Vector{UInt32}
    prec_prob::Vector{Float32}
    prec_peak::Vector{Float32}
    prec_basepep::Vector{UInt32}
    prec_seq::Vector{String}
    prec_charge::Vector{UInt8}
    prec_smods::Vector{String}
    prec_imods::Vector{String}
    # modified-peptide level (dedup by (base_pep_id, structural_mods, isotopic_mods))
    mod_pos::Dict{Tuple{UInt32, String, String}, Int}
    mod_key::Vector{Tuple{UInt32, String, String}}
    mod_logsum::Vector{Float64}
    mod_peak::Vector{Float32}
    mod_seq::Vector{String}
    mod_count::Vector{Int32}
    # peptide level (dedup by sequence)
    pep_pos::Dict{String, Int}
    pep_seq::Vector{String}
    pep_logsum::Vector{Float64}
    pep_peak::Vector{Float32}
    pep_count::Vector{Int32}
end

ProteinRollupScratch() = ProteinRollupScratch(
    Dict{UInt32, Int}(), UInt32[], Float32[], Float32[], UInt32[], String[], UInt8[], String[], String[],
    Dict{Tuple{UInt32, String, String}, Int}(), Tuple{UInt32, String, String}[], Float64[], Float32[], String[], Int32[],
    Dict{String, Int}(), String[], Float64[], Float32[], Int32[],
)

function _reset_rollup_scratch!(s::ProteinRollupScratch)
    empty!(s.prec_pos); empty!(s.prec_id); empty!(s.prec_prob); empty!(s.prec_peak)
    empty!(s.prec_basepep); empty!(s.prec_seq); empty!(s.prec_charge)
    empty!(s.prec_smods); empty!(s.prec_imods)
    empty!(s.mod_pos); empty!(s.mod_key); empty!(s.mod_logsum)
    empty!(s.mod_peak); empty!(s.mod_seq); empty!(s.mod_count)
    empty!(s.pep_pos); empty!(s.pep_seq); empty!(s.pep_logsum)
    empty!(s.pep_peak); empty!(s.pep_count)
    return s
end

@inline function _rollup_mods_or_empty(gdf, col::Symbol, i::Int)
    (hasproperty(gdf, col) && !ismissing(gdf[!, col][i])) ? String(gdf[!, col][i]) : ""
end

"""
    _build_protein_rollup(gdf, quant_mask, prob_col[, scratch])

Roll up precursor-level probability evidence to modified peptides, then
peptides, and return the peptide-level protein score plus reusable
intermediate rows. Output rows are emitted in first-seen (deterministic) order.
The 3-arg form allocates a fresh `ProteinRollupScratch`; pass a reused one (4-arg)
to avoid per-group allocation.
"""
_build_protein_rollup(gdf::AbstractDataFrame, quant_mask::AbstractVector{Bool}, prob_col::Symbol) =
    _build_protein_rollup(gdf, quant_mask, prob_col, ProteinRollupScratch())

function _build_protein_rollup(
    gdf::AbstractDataFrame,
    quant_mask::AbstractVector{Bool},
    prob_col::Symbol,
    s::ProteinRollupScratch,
)
    if !hasproperty(gdf, prob_col)
        error("Missing required $(prob_col) column for protein roll-up")
    end
    if !hasproperty(gdf, :peak_area)
        error("Missing required :peak_area column for protein roll-up")
    end
    if !hasproperty(gdf, :base_pep_id)
        error("Missing required :base_pep_id column for protein roll-up")
    end

    _reset_rollup_scratch!(s)

    # Precursor level: dedup by precursor_idx into parallel arrays (first-seen order).
    @inbounds for i in eachindex(quant_mask)
        quant_mask[i] || continue

        precursor_idx = UInt32(gdf.precursor_idx[i])
        prob_val = _sanitize_rollup_probability(gdf[!, prob_col][i])
        peak_area_val = Float32(gdf.peak_area[i])

        pos = get(s.prec_pos, precursor_idx, 0)
        if pos == 0
            push!(s.prec_id, precursor_idx)
            push!(s.prec_prob, prob_val)
            push!(s.prec_peak, peak_area_val)
            push!(s.prec_basepep, UInt32(gdf.base_pep_id[i]))
            push!(s.prec_seq, String(gdf.sequence[i]))
            push!(s.prec_charge, hasproperty(gdf, :charge) ? UInt8(gdf.charge[i]) : UInt8(0))
            push!(s.prec_smods, _rollup_mods_or_empty(gdf, :structural_mods, i))
            push!(s.prec_imods, _rollup_mods_or_empty(gdf, :isotopic_mods, i))
            s.prec_pos[precursor_idx] = length(s.prec_id)
        else
            s.prec_prob[pos] = max(s.prec_prob[pos], prob_val)
            if peak_area_val > s.prec_peak[pos]
                s.prec_peak[pos] = peak_area_val
            end
        end
    end

    nprec = length(s.prec_id)
    precursor_rows = ProteinRollupPrecursorRow[]
    sizehint!(precursor_rows, nprec)
    @inbounds for p in 1:nprec
        prob_val = s.prec_prob[p]
        none_prob_val = Float32(1.0f0 - prob_val)
        score_val = Float32(-log(none_prob_val))
        push!(precursor_rows, (
            precursor_idx = s.prec_id[p],
            base_pep_id = s.prec_basepep[p],
            sequence = s.prec_seq[p],
            charge = s.prec_charge[p],
            structural_mods = s.prec_smods[p],
            isotopic_mods = s.prec_imods[p],
            pep = none_prob_val,
            prob = prob_val,
            score = score_val,
            best_peak_area = s.prec_peak[p]
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
            top_pep_peak_area = 0.0f0
        )
    end

    # Modified-peptide level: dedup by (base_pep_id, structural_mods, isotopic_mods).
    for precursor_row in precursor_rows
        mod_key = (
            precursor_row.base_pep_id,
            precursor_row.structural_mods,
            precursor_row.isotopic_mods
        )
        log_none = _precursor_log_none_for_rollup(precursor_row.prob)
        mpos = get(s.mod_pos, mod_key, 0)
        if mpos == 0
            push!(s.mod_key, mod_key)
            push!(s.mod_logsum, log_none)
            push!(s.mod_peak, precursor_row.best_peak_area)
            push!(s.mod_seq, precursor_row.sequence)
            push!(s.mod_count, Int32(1))
            s.mod_pos[mod_key] = length(s.mod_key)
        else
            s.mod_logsum[mpos] += log_none
            if precursor_row.best_peak_area > s.mod_peak[mpos]
                s.mod_peak[mpos] = precursor_row.best_peak_area
            end
            s.mod_seq[mpos] = precursor_row.sequence
            s.mod_count[mpos] += Int32(1)
        end
    end

    nmod = length(s.mod_key)
    modified_peptide_rows = ProteinRollupModifiedPeptideRow[]
    sizehint!(modified_peptide_rows, nmod)
    @inbounds for m in 1:nmod
        log_none_sum = s.mod_logsum[m]
        push!(modified_peptide_rows, (
            base_pep_id = s.mod_key[m][1],
            sequence = s.mod_seq[m],
            structural_mods = s.mod_key[m][2],
            isotopic_mods = s.mod_key[m][3],
            log_none_sum = log_none_sum,
            pep = _none_probability_from_log_none_sum(log_none_sum),
            prob = _probability_from_log_none_sum(log_none_sum),
            score = _score_from_log_none_sum(log_none_sum),
            best_peak_area = s.mod_peak[m],
            precursor_count = s.mod_count[m]
        ))
    end

    # Peptide level: dedup by sequence.
    for modified_row in modified_peptide_rows
        peptide_key = modified_row.sequence
        ppos = get(s.pep_pos, peptide_key, 0)
        if ppos == 0
            push!(s.pep_seq, peptide_key)
            push!(s.pep_logsum, modified_row.log_none_sum)
            push!(s.pep_peak, modified_row.best_peak_area)
            push!(s.pep_count, Int32(1))
            s.pep_pos[peptide_key] = length(s.pep_seq)
        else
            s.pep_logsum[ppos] += modified_row.log_none_sum
            if modified_row.best_peak_area > s.pep_peak[ppos]
                s.pep_peak[ppos] = modified_row.best_peak_area
            end
            s.pep_count[ppos] += Int32(1)
        end
    end

    npep = length(s.pep_seq)
    peptide_rows = ProteinRollupPeptideRow[]
    sizehint!(peptide_rows, npep)
    @inbounds for q in 1:npep
        log_none_sum = s.pep_logsum[q]
        push!(peptide_rows, (
            sequence = s.pep_seq[q],
            log_none_sum = log_none_sum,
            pep = _none_probability_from_log_none_sum(log_none_sum),
            prob = _probability_from_log_none_sum(log_none_sum),
            score = _score_from_log_none_sum(log_none_sum),
            best_peak_area = s.pep_peak[q],
            modified_peptide_count = s.pep_count[q]
        ))
    end
    pg_score = isempty(peptide_rows) ? 0.0f0 : Float32(sum(row.score for row in peptide_rows))
    top_pep_peak_area = isempty(peptide_rows) ? 0.0f0 : maximum(row.best_peak_area for row in peptide_rows)
    peptide_list = sort!([row.sequence for row in peptide_rows])

    return (
        precursor_rows = precursor_rows,
        modified_peptide_rows = modified_peptide_rows,
        peptide_rows = peptide_rows,
        pg_score = pg_score,
        n_peptides = length(peptide_rows),
        peptide_list = peptide_list,
        top_pep_peak_area = top_pep_peak_area
    )
end

const ConsensusRunVote = @NamedTuple{
    pg_score::Float32,
    run_order::Int64,
    normalized_precursors::Vector{Pair{UInt32, Float32}}
}

@inline function _consensus_runs_to_keep(n_runs::Int)::Int
    n_runs <= 1 && return 0
    return min(n_runs - 1, max(5, ceil(Int, log2(n_runs))))
end

@inline function _consensus_candidate_runs_to_keep(n_runs::Int)::Int
    target_runs = _consensus_runs_to_keep(n_runs)
    return min(n_runs, target_runs + 1)
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
    _insert_top_consensus_vote!(run_votes, vote, max_votes)

Insert `vote` into the in-place score-sorted run list for one protein,
keeping only the top `max_votes` entries.
"""
function _insert_top_consensus_vote!(
    run_votes::Vector{ConsensusRunVote},
    vote::ConsensusRunVote,
    max_votes::Int
)
    max_votes <= 0 && return run_votes

    if (length(run_votes) == max_votes) &&
       !_consensus_vote_precedes(vote.pg_score, vote.run_order, run_votes[end])
        return run_votes
    end

    insert_idx = length(run_votes) + 1

    @inbounds for i in eachindex(run_votes)
        if _consensus_vote_precedes(vote.pg_score, vote.run_order, run_votes[i])
            insert_idx = i
            break
        end
    end

    insert!(run_votes, insert_idx, vote)
    if length(run_votes) > max_votes
        pop!(run_votes)
    end

    return run_votes
end

"""
    build_precursor_consensus(psm_refs::Vector{PSMFileReference})

Build a dataset-level precursor relative-weight consensus within each inferred
protein group. Consensus weights are derived from quant precursors after
normalizing each precursor by the max precursor peak area in that protein run.
Each run's vote is weighted by the run's current protein `pg_score`. While
streaming the runs, only the top
`min(n_runs - 1, max(5, ceil(log2(n_runs)))) + 1` runs by `pg_score` are
retained as a cached candidate pool so leave-one-run-out consensus can be
computed by subtraction at scoring time.
"""
function build_precursor_consensus(
    psm_refs::Vector{PSMFileReference};
    q_value_threshold::Float32 = 0.01f0
)
    protein_run_votes = Dict{Tuple{String, Bool, UInt8}, Vector{ConsensusRunVote}}()
    protein_observed_run_count = Dict{Tuple{String, Bool, UInt8}, Int32}()
    max_candidate_runs = _consensus_candidate_runs_to_keep(length(psm_refs))
    rollup_scratch = ProteinRollupScratch()   # reused across all groups (sequential loop)

    for (run_order, psm_ref) in enumerate(psm_refs)
        if !exists(psm_ref)
            continue
        end

        df = load_dataframe(psm_ref)
        prob_col = _protein_group_probability_column(df)

        for gdf in groupby(df, [:inferred_protein_group, :target, :entrap_id])
            quant_mask = _protein_rollup_quant_mask(gdf; q_value_threshold = q_value_threshold)
            rollup = _build_protein_rollup(gdf, quant_mask, prob_col, rollup_scratch)

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

            run_max_peak_area = maximum(precursor_row.best_peak_area for precursor_row in rollup.precursor_rows)
            if !isfinite(run_max_peak_area) || run_max_peak_area <= 0.0f0
                continue
            end

            normalized_precursors = Pair{UInt32, Float32}[]
            sizehint!(normalized_precursors, length(rollup.precursor_rows))
            for precursor_row in rollup.precursor_rows
                push!(
                    normalized_precursors,
                    precursor_row.precursor_idx => Float32(clamp(precursor_row.best_peak_area / run_max_peak_area, 0.0f0, 1.0f0))
                )
            end
            protein_key = (protein_name, target, entrap_id)
            protein_observed_run_count[protein_key] = get(protein_observed_run_count, protein_key, Int32(0)) + Int32(1)
            _insert_top_consensus_vote!(
                get!(
                    () -> sizehint!(ConsensusRunVote[], max_candidate_runs),
                    protein_run_votes,
                    protein_key
                ),
                (
                    pg_score = rollup.pg_score,
                    run_order = Int64(run_order),
                    normalized_precursors = normalized_precursors
                ),
                max_candidate_runs
            )
        end
    end

    consensus_weight_sums = Dict{Tuple{String, Bool, UInt8, UInt32}, Float64}()
    protein_total_vote = Dict{Tuple{String, Bool, UInt8}, Float64}()
    selected_run_votes = Dict{Tuple{String, Bool, UInt8}, Vector{ConsensusRunVote}}()
    consensus_target_run_count = Dict{Tuple{String, Bool, UInt8}, Int32}()
    cached_consensus_weight_sums = Dict{Tuple{String, Bool, UInt8}, Dict{UInt32, Float64}}()
    cached_protein_total_vote = Dict{Tuple{String, Bool, UInt8}, Float64}()
    for (protein_key, run_votes) in protein_run_votes
        observed_runs = Int(get(protein_observed_run_count, protein_key, Int32(length(run_votes))))
        target_runs = _consensus_runs_to_keep(observed_runs)
        candidate_runs = _consensus_candidate_runs_to_keep(observed_runs)
        selected_votes = run_votes[1:min(candidate_runs, length(run_votes))]
        selected_run_votes[protein_key] = selected_votes
        consensus_target_run_count[protein_key] = Int32(target_runs)

        cached_weight_sums = Dict{UInt32, Float64}()
        cached_total_vote = 0.0
        for run_vote in selected_votes
            run_weight = Float64(run_vote.pg_score)
            cached_total_vote += run_weight
            for precursor in run_vote.normalized_precursors
                cached_weight_sums[precursor.first] = get(cached_weight_sums, precursor.first, 0.0) + (run_weight * Float64(precursor.second))
            end
        end
        cached_consensus_weight_sums[protein_key] = cached_weight_sums
        cached_protein_total_vote[protein_key] = cached_total_vote

        default_votes = selected_votes[1:min(target_runs, length(selected_votes))]
        for run_vote in default_votes
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

    profiled_precursor_count = Dict{Tuple{String, Bool, UInt8}, Int32}()
    shape_strength = Dict{Tuple{String, Bool, UInt8}, Float32}()
    for (protein_key, relative_weights) in protein_precursor_values
        profiled_precursor_count[protein_key] = Int32(length(relative_weights))
    end
    for (protein_key, run_votes) in selected_run_votes
        total_pg_score = sum(Float64(run_vote.pg_score) for run_vote in run_votes)
        shape_strength[protein_key] = Float32(log1p(total_pg_score))
    end
    shape_strength_values = collect(values(shape_strength))
    shape_confidence_scale = isempty(shape_strength_values) ? 1.0f0 :
        Float32(max(Statistics.median(shape_strength_values), eps(Float32)))
    return (
        relative_weight = relative_weight,
        profiled_precursor_count = profiled_precursor_count,
        shape_strength = shape_strength,
        shape_confidence_scale = shape_confidence_scale,
        selected_run_votes = selected_run_votes,
        consensus_target_run_count = consensus_target_run_count,
        cached_consensus_weight_sums = cached_consensus_weight_sums,
        cached_protein_total_vote = cached_protein_total_vote
    )
end

"""
    _update_protein_cv_fold_mapping!(protein_to_cv_fold, df, precursors)

Update the protein-to-CV-fold mapping from one annotated passing-PSM table,
using the highest precursor score observed for each inferred protein group.
"""
function _update_protein_cv_fold_mapping!(
    protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
    df::DataFrame,
    precursors::LibraryPrecursors
)
    nrow(df) == 0 && return protein_to_cv_fold

    filtered_df = df[.!ismissing.(df.inferred_protein_group), :]
    nrow(filtered_df) == 0 && return protein_to_cv_fold

    prob_col = _protein_group_probability_column(filtered_df)

    for gdf in groupby(filtered_df, :inferred_protein_group)
        best_score = typemin(Float32)
        best_precursor_idx = zero(UInt32)

        @inbounds for i in 1:nrow(gdf)
            score = Float32(gdf[!, prob_col][i])
            if score > best_score
                best_score = score
                best_precursor_idx = UInt32(gdf.precursor_idx[i])
            end
        end

        best_precursor_idx == zero(UInt32) && continue

        protein_name = String(gdf.inferred_protein_group[1])
        value = (best_score = best_score, cv_fold = UInt8(getCvFold(precursors, best_precursor_idx)))
        if !haskey(protein_to_cv_fold, protein_name)
            insert!(protein_to_cv_fold, protein_name, value)
        elseif best_score > protein_to_cv_fold[protein_name].best_score
            protein_to_cv_fold[protein_name] = value
        end
    end

    return protein_to_cv_fold
end

"""
    _shape_consensus_inputs(precursor_rows, protein_key, precursor_consensus)

Prepare the observed precursor peak areas and protein key used for shape scoring.
"""
function _shape_consensus_inputs(
    precursor_rows::AbstractVector,
    protein_key::Tuple{String, Bool, UInt8},
    precursor_consensus::NamedTuple
)
    best_peak_area_by_precursor = Dict{UInt32, Float32}()
    for row in precursor_rows
        best_peak_area_by_precursor[row.precursor_idx] = row.best_peak_area
    end

    return (
        protein_key = protein_key,
        best_peak_area_by_precursor = best_peak_area_by_precursor
    )
end

"""
    _precursor_consensus_prefix_features(precursor_rows, protein_key, precursor_consensus)

Score how well the run matches the consensus precursor profile up to the weakest
observed consensus precursor using normalized Manhattan similarity. Off-consensus
observed precursors are included as extra dimensions with zero consensus
intensity, but `tau` is still determined only by the matched consensus
precursors.
"""
function _active_consensus_profile(
    protein_key::Tuple{String, Bool, UInt8},
    precursor_consensus::NamedTuple;
    current_run_order::Union{Nothing, Int64} = nothing
)
    candidate_votes = get(precursor_consensus.selected_run_votes, protein_key, ConsensusRunVote[])
    desired_runs = Int(get(precursor_consensus.consensus_target_run_count, protein_key, Int32(length(candidate_votes))))

    if isempty(candidate_votes) || desired_runs <= 0
        return (
            consensus_precursors = Pair{UInt32, Float32}[],
            profiled_precursor_count = 0,
            shape_strength = 0.0f0
        )
    end

    active_votes = ConsensusRunVote[]
    excluded_votes = ConsensusRunVote[]
    sizehint!(active_votes, min(desired_runs, length(candidate_votes)))
    sizehint!(excluded_votes, length(candidate_votes))

    for vote in candidate_votes
        if !isnothing(current_run_order) && (vote.run_order == current_run_order)
            push!(excluded_votes, vote)
        else
            push!(active_votes, vote)
        end
    end

    if length(active_votes) > desired_runs
        append!(excluded_votes, active_votes[(desired_runs + 1):end])
        resize!(active_votes, desired_runs)
    end

    cached_weight_sums = get(precursor_consensus.cached_consensus_weight_sums, protein_key, Dict{UInt32, Float64}())
    cached_total_vote = get(precursor_consensus.cached_protein_total_vote, protein_key, 0.0)

    excluded_total_vote = 0.0
    excluded_weight_sums = Dict{UInt32, Float64}()
    for vote in excluded_votes
        vote_weight = Float64(vote.pg_score)
        excluded_total_vote += vote_weight
        for precursor in vote.normalized_precursors
            excluded_weight_sums[precursor.first] = get(excluded_weight_sums, precursor.first, 0.0) + (vote_weight * Float64(precursor.second))
        end
    end

    active_total_vote = cached_total_vote - excluded_total_vote
    active_total_vote > 0.0 || return (
        consensus_precursors = Pair{UInt32, Float32}[],
        profiled_precursor_count = 0,
        shape_strength = 0.0f0
    )

    consensus_precursors = Pair{UInt32, Float32}[]
    sizehint!(consensus_precursors, length(cached_weight_sums))
    for (precursor_idx, cached_sum) in cached_weight_sums
        adjusted_sum = cached_sum - get(excluded_weight_sums, precursor_idx, 0.0)
        adjusted_sum > 0.0 || continue
        push!(consensus_precursors, precursor_idx => Float32(clamp(adjusted_sum / active_total_vote, 0.0, 1.0)))
    end
    sort!(consensus_precursors; by = last, rev = true)

    return (
        consensus_precursors = consensus_precursors,
        profiled_precursor_count = length(consensus_precursors),
        shape_strength = Float32(log1p(active_total_vote))
    )
end

function _precursor_consensus_prefix_features(
    precursor_rows::AbstractVector,
    protein_key::Tuple{String, Bool, UInt8},
    precursor_consensus::NamedTuple;
    current_run_order::Union{Nothing, Int64} = nothing
)
    shape_inputs = _shape_consensus_inputs(precursor_rows, protein_key, precursor_consensus)
    shape_protein_key = shape_inputs.protein_key
    best_peak_area_by_precursor = shape_inputs.best_peak_area_by_precursor

    active_consensus = _active_consensus_profile(
        shape_protein_key,
        precursor_consensus;
        current_run_order = current_run_order
    )
    profiled_precursor_count = active_consensus.profiled_precursor_count
    shape_strength = active_consensus.shape_strength
    shape_confidence_scale = hasproperty(precursor_consensus, :shape_confidence_scale) ?
        precursor_consensus.shape_confidence_scale : 1.0f0
    shape_confidence = _shape_confidence(shape_strength, shape_confidence_scale)

    if isempty(best_peak_area_by_precursor)
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

    consensus_relative_weight = Dict{Tuple{String, Bool, UInt8, UInt32}, Float32}()
    consensus_precursors = active_consensus.consensus_precursors
    for precursor in consensus_precursors
        consensus_relative_weight[(shape_protein_key[1], shape_protein_key[2], shape_protein_key[3], precursor.first)] = precursor.second
    end
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

    run_max_peak_area = maximum(values(best_peak_area_by_precursor))
    (!isfinite(run_max_peak_area) || run_max_peak_area <= 0.0f0) && return (
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
    off_consensus_precursors = UInt32[]
    matched_precursors = 0
    for precursor_idx in keys(best_peak_area_by_precursor)
        c = get(consensus_relative_weight, (shape_protein_key[1], shape_protein_key[2], shape_protein_key[3], precursor_idx), 0.0f0)
        if c > 0.0f0
            push!(observed_consensus_weights, c)
            matched_precursors += 1
        else
            push!(off_consensus_precursors, precursor_idx)
        end
    end

    threshold = isempty(observed_consensus_weights) ? 0.0f0 : minimum(observed_consensus_weights)
    prefix_consensus_sum = 0.0f0
    run_prefix_sum = 0.0f0
    included_consensus_weights = Float32[]
    included_run_weights = Float32[]

    for precursor in consensus_precursors
        precursor.second < threshold && continue
        run_relative_weight = get(best_peak_area_by_precursor, precursor.first, 0.0f0) / run_max_peak_area
        prefix_consensus_sum += precursor.second
        run_prefix_sum += run_relative_weight
        push!(included_consensus_weights, precursor.second)
        push!(included_run_weights, run_relative_weight)
    end

    for precursor_idx in off_consensus_precursors
        run_relative_weight = get(best_peak_area_by_precursor, precursor_idx, 0.0f0) / run_max_peak_area
        run_prefix_sum += run_relative_weight
        push!(included_consensus_weights, 0.0f0)
        push!(included_run_weights, run_relative_weight)
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

    l1_distance = 0.0f0
    for i in eachindex(included_consensus_weights, included_run_weights)
        consensus_norm = included_consensus_weights[i] / prefix_consensus_sum
        run_norm = included_run_weights[i] / run_prefix_sum
        l1_distance += abs(run_norm - consensus_norm)
    end

    prefix_shape_raw = Float32(clamp(1.0f0 - 0.5f0 * l1_distance, 0.0f0, 1.0f0))
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

"""
    group_psms_by_protein(df::DataFrame)

Transform PSMs into protein groups by aggregating peptides.
Returns a DataFrame with one row per protein group.
"""
function group_psms_by_protein(
    df::DataFrame;
    precursor_consensus::NamedTuple,
    current_run_order::Union{Nothing, Int64} = nothing,
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
            top_pep_peak_area = Float32[],
            precursor_consensus_prefix_shape = Float32[],
            mbr_recovered_peptides = Int64[],
            mbr_only_protein = Bool[]
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
            top_pep_peak_area = Float32[],
            precursor_consensus_prefix_shape = Float32[],
            mbr_recovered_peptides = Int64[],
            mbr_only_protein = Bool[]
        )
    end

    prob_col = _protein_group_probability_column(df)
    # Group by protein
    grouped = groupby(df, [:inferred_protein_group, :target, :entrap_id])

    # Aggregate to protein groups. combine() runs the do-block multi-threaded
    # (threads=true default), so use one scratch per thread. The do-block is pure
    # compute (no yield/I/O), so the task does not migrate threads and threadid()
    # is stable for the call; sized by maxthreadid() (threadid() can exceed
    # nthreads() with an interactive threadpool).
    rollup_scratches = [ProteinRollupScratch() for _ in 1:Threads.maxthreadid()]
    protein_groups = combine(grouped) do gdf
        quant_mask = _protein_rollup_quant_mask(gdf; q_value_threshold = q_value_threshold)
        rollup = _build_protein_rollup(gdf, quant_mask, prob_col, rollup_scratches[Threads.threadid()])
        mbr_summary = _protein_mbr_summary(gdf, quant_mask)
        n_peptides = rollup.n_peptides
        pg_score = rollup.pg_score
        top_pep_peak_area = rollup.top_pep_peak_area
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
                precursor_consensus;
                current_run_order = current_run_order
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
            top_pep_peak_area = top_pep_peak_area,
            precursor_consensus_prefix_shape = precursor_consensus_prefix_shape,
            mbr_recovered_peptides = mbr_summary.recovered_peptides,
            mbr_only_protein = mbr_summary.only_protein
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
        filter!(row -> row.n_peptides >= min_peptides, df)
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
    build_protein_group_tables(psm_refs, output_folder, protein_catalog; precursors, kwargs...)

Build per-run protein-group tables and protein scoring features from PSM tables
that have already been annotated with protein inference results.

# Arguments
- `psm_refs`: Vector of PSM file references
- `output_folder`: Directory for protein group output
- `protein_catalog`: Pre-computed protein-to-peptide mappings
- `precursors`: Library precursors used for protein CV fold assignment
- `min_peptides`: Minimum peptides per protein group (default: 2)

# Returns
- `pg_refs`: Vector of protein group file references
- `psm_to_pg_mapping`: Dictionary mapping PSM paths to protein group paths
- `protein_to_cv_fold`: Protein-to-CV-fold mapping built during protein-group table construction
"""
function build_protein_group_tables(
    psm_refs::Vector{PSMFileReference},
    output_folder::String,
    protein_catalog::Dict;
    precursors::LibraryPrecursors,
    min_peptides::Int = 2,
    q_value_threshold::Float32 = 0.01f0
)
    !isdir(output_folder) && mkpath(output_folder)

    pg_refs = ProteinGroupFileReference[]
    psm_to_pg_mapping = Dict{String, String}()
    protein_to_cv_fold = Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}()
    indexed_refs = collect(enumerate(psm_refs))

    precursor_consensus = build_precursor_consensus(psm_refs; q_value_threshold = q_value_threshold)

    @user_info "Building per-run protein group tables and protein scoring features"

    for (idx, psm_ref) in ProgressBar(indexed_refs)
        if !exists(psm_ref)
            continue
        end

        updated_psms = load_dataframe(psm_ref)
        _update_protein_cv_fold_mapping!(protein_to_cv_fold, updated_psms, precursors)
        peak_area_calibration = estimate_peak_area_detection_model(updated_psms)

        protein_groups_df = group_psms_by_protein(
            updated_psms;
            precursor_consensus = precursor_consensus,
            current_run_order = Int64(idx),
            q_value_threshold = q_value_threshold
        )

        post_inference_pipeline = TransformPipeline() |>
            filter_by_min_peptides(min_peptides) |>
            add_protein_features(protein_catalog) |>
            add_peak_area_observation_features(peak_area_calibration)

        for (desc, op) in post_inference_pipeline.operations
            protein_groups_df = op(protein_groups_df)
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

    return pg_refs, psm_to_pg_mapping, protein_to_cv_fold
end
