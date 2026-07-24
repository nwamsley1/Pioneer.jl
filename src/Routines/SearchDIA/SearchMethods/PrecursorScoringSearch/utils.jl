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
    return _pep_interp_from_vectors(scores, targets; fdr_scale_factor = fdr_scale_factor)
end

# PEP interpolation from pre-extracted columns. Split out (like _qvalue_spline_from_vectors) so the
# merged file is read once for both the q-spline and the PEP interp (T1).
function _pep_interp_from_vectors(
    scores::Vector{Float32},
    targets::Vector{Bool};
    fdr_scale_factor::Float32 = 1.0f0,
)
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

    # Hoist the columns into concrete typed vectors ONCE, then delegate to the vector core.
    # Accessing `psms_scores[!, sym][i]` inside the per-row loops (Q can be millions) re-resolves
    # the column by symbol every iteration and returns an abstract `AbstractVector`, boxing each
    # element read.
    target_col_vec = Vector{Bool}(psms_scores[!, :target])
    score_col_vec = Vector{Float32}(psms_scores[!, score_col])
    return _qvalue_spline_from_vectors(target_col_vec, score_col_vec;
        min_pep_points_per_bin = min_pep_points_per_bin, fdr_scale_factor = fdr_scale_factor)
end

# Q-value spline from pre-extracted columns, in score-descending order (the merged file's sort).
# Split out of get_qvalue_spline so build_qvalue_spline_from_refs can read the merged file ONCE and
# feed both this and the PEP interpolation, instead of each re-reading it (A3 double-read, T1).
function _qvalue_spline_from_vectors(
    target_col_vec::Vector{Bool},
    score_col_vec::Vector{Float32};
    min_pep_points_per_bin::Int = 1000,
    fdr_scale_factor::Float32 = 1.0f0,
)
    Q = length(target_col_vec)
    M = ceil(Int, (Q - 1) / min_pep_points_per_bin) + 1
    bin_qval, bin_mean_prob = Vector{Float32}(undef, M), Vector{Float32}(undef, M)
    bin_size = 0
    bin_idx = 0
    mean_prob, targets, decoys = 0.0f0, 0, 0
    for i in range(1, Q)
        targets += target_col_vec[i]
        decoys += (1 - target_col_vec[i])
    end
    min_q_val = typemax(Float32)
    for i in reverse(range(1, Q))
        bin_size += 1
        targets -= target_col_vec[i]
        decoys -= (1 - target_col_vec[i])
        mean_prob += score_col_vec[i]
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

