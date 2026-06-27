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
PSM score merging and processing
==========================================================#

const MBR_RESCUE_PSMS_DIRNAME = "mbr_rescue_psms"
const MAIN_SEARCH_PSMS_DIRNAME = "main_search_psms"
const MBR_RESCUE_RECOVERED_SUFFIX = ".mbr_recovered.arrow"

function mbr_rescue_fold_path(main_fold_path::AbstractString)
    parts = splitpath(String(main_fold_path))
    idx = findlast(==(MAIN_SEARCH_PSMS_DIRNAME), parts)
    idx === nothing && throw(ArgumentError(
        "Expected a fold path under $MAIN_SEARCH_PSMS_DIRNAME, got $main_fold_path"
    ))
    parts[idx] = MBR_RESCUE_PSMS_DIRNAME
    return joinpath(parts...)
end

function get_existing_mbr_rescue_fold_paths(main_fold_paths::Vector{String})
    rescue_paths = String[]
    for fpath in main_fold_paths
        rescue_path = mbr_rescue_fold_path(fpath)
        isfile(rescue_path) && push!(rescue_paths, rescue_path)
    end
    return rescue_paths
end

function mbr_rescue_recovered_path(rescue_fold_path::AbstractString)
    return String(rescue_fold_path) * MBR_RESCUE_RECOVERED_SUFFIX
end


"""
    get_pep_interpolation(merged_psms_path::String, score_col::Symbol;
                          fdr_scale_factor=1.0f0)

Create an interpolation function mapping scores to posterior error probabilities.

PEP values are estimated across all runs using `get_PEP!` and then fit with a
linear interpolation.  The returned object can be used with
`add_interpolated_column` to populate per-file PEP columns.
"""
function get_pep_interpolation(
    merged_psms_path::String,
    score_col::Symbol;
    fdr_scale_factor::Float32 = 1.0f0,
)
    df = DataFrame(Arrow.Table(merged_psms_path))
    scores  = Vector{Float32}(df[!, score_col])
    targets = Vector{Bool}(df[!, :target])

    pep_vals = Vector{Float32}(undef, length(scores))
    get_PEP!(scores, targets, pep_vals; doSort=true,
             fdr_scale_factor=fdr_scale_factor)

    # Collapse redundant PEP values by taking the minimum score for each PEP
    pep_to_score = Dict{Float32, Float32}()
    order = sortperm(scores)
    for idx in order
        p = pep_vals[idx]
        s = scores[idx]
        if !haskey(pep_to_score, p)
            pep_to_score[p] = s
        end
    end

    # Sort the resulting pairs by score for interpolation
    ordered = sort(collect(pep_to_score), by = last)
    xs_tmp = Float32[last(pair) for pair in ordered]
    ys_tmp = Float32[first(pair) for pair in ordered]

    # Remove duplicate knot values explicitly
    xs = Float32[]
    ys = Float32[]
    prev_x = NaN32
    for (x, y) in zip(xs_tmp, ys_tmp)
        if x != prev_x && !isnan(x)
            push!(xs, x)
            push!(ys, clamp(y, 0.0f0, 1.0f0))
            prev_x = x
        end
    end

    if length(xs) < 2
        @user_warn "PEP interpolation: only $(length(xs)) unique score points — using flat default. Affects only the :pep column; q-values are unchanged."
        xs = Float32[0.0, 1.0]
        ys = Float32[1.0, 0.0]
    end

    Interpolations.deduplicate_knots!(xs)

    return linear_interpolation(xs, ys; extrapolation_bc=Interpolations.Flat())
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
        # select the columns needed to identify globally unique precursors or proteins
        if score_col == :global_prob # precursors
            select!(psms_scores, [:precursor_idx, :target, score_col])
        elseif score_col == :global_pg_score # proteins
            select!(psms_scores, [:protein_name, :target, :entrap_id, score_col])
        end
        psms_scores = unique(psms_scores)
    end

    Q = size(psms_scores, 1)
    M = ceil(Int, (Q - 1) / min_pep_points_per_bin) + 1
    bin_qval, bin_mean_prob = Vector{Float32}(undef, M), Vector{Float32}(undef, M)
    bin_size = 0
    bin_idx = 0
    mean_prob, targets, decoys = 0.0f0, 0, 0
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
        if bin_idx == 0 || bin_size == min_pep_points_per_bin
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
        @user_warn "Q-value interpolation: only $(length(xs)) unique score points — using flat default. Likely indicates a small or low-diversity input."
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
    return linear_interpolation(xs, ys; extrapolation_bc = Line())
end
