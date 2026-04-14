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

# ==========================================================
# QC plotting
# ==========================================================

function plot_target_protection_spectral_angle_qc(
    psm_paths::Vector{String},
    output_path::String;
    max_rows::Int64 = 1_000_000
)
    total_rows = 0
    for path in psm_paths
        isfile(path) || continue
        tbl = Arrow.Table(path)
        if hasproperty(tbl, :target_protection_spectral_angle) && hasproperty(tbl, :target)
            total_rows += length(tbl[:target_protection_spectral_angle])
        end
    end

    if total_rows == 0
        @user_warn "Skipping target protection spectral angle QC plots: no target_protection_spectral_angle rows found"
        return nothing
    end

    stride = max(1, cld(total_rows, max_rows))
    spectral_angles = Float32[]
    target_neighbor_spectral_angles = Float32[]
    targets = Bool[]
    nearest_decoy_is_pair = Bool[]
    sequence_lengths = UInt8[]
    has_sequence_length = true
    has_pair_flag = true
    has_target_neighbor = true
    rows_seen = 0

    for path in psm_paths
        isfile(path) || continue
        tbl = Arrow.Table(path)
        (hasproperty(tbl, :target_protection_spectral_angle) && hasproperty(tbl, :target)) || continue

        has_lengths = hasproperty(tbl, :sequence_length)
        has_sequence_length &= has_lengths
        file_has_pair_flag = hasproperty(tbl, :target_protection_pair_decoy)
        has_pair_flag &= file_has_pair_flag
        file_has_target_neighbor = hasproperty(tbl, :target_neighbor_spectral_angle)
        has_target_neighbor &= file_has_target_neighbor
        n = length(tbl[:target_protection_spectral_angle])
        for i in 1:n
            rows_seen += 1
            ((rows_seen - 1) % stride == 0) || continue

            spectral_angle = Float32(tbl[:target_protection_spectral_angle][i])
            isfinite(spectral_angle) || continue
            target_neighbor_spectral_angle = file_has_target_neighbor ? Float32(tbl[:target_neighbor_spectral_angle][i]) : 0.0f0
            isfinite(target_neighbor_spectral_angle) || (target_neighbor_spectral_angle = 0.0f0)

            push!(spectral_angles, clamp(spectral_angle, 0.0f0, 1.0f0))
            push!(target_neighbor_spectral_angles, clamp(target_neighbor_spectral_angle, 0.0f0, 1.0f0))
            push!(targets, Bool(tbl[:target][i]))
            push!(nearest_decoy_is_pair, file_has_pair_flag ? Bool(tbl[:target_protection_pair_decoy][i]) : false)
            if has_lengths
                push!(sequence_lengths, UInt8(tbl[:sequence_length][i]))
            end
        end
    end

    if isempty(spectral_angles)
        @user_warn "Skipping target protection spectral angle QC plots: sampled rows were empty"
        return nothing
    end

    plots = Plots.Plot[]
    bins = range(0.0, 1.0, length=51)
    if !has_pair_flag
        @user_warn "Target protection paired-decoy QC split unavailable: target_protection_pair_decoy column missing"
    end
    if !has_target_neighbor
        @user_warn "Target-neighbor QC split unavailable: target_neighbor_spectral_angle column missing"
    end
    target_pair_mask = targets .& nearest_decoy_is_pair
    target_other_mask = targets .& .!nearest_decoy_is_pair
    decoy_mask = .!targets
    bin_centers = Float32.((collect(bins[1:end-1]) .+ collect(bins[2:end])) ./ 2)

    function log2_bin_counts(values::AbstractVector{Float32})
        counts = zeros(Int64, length(bin_centers))
        for value in values
            bin_idx = searchsortedlast(bins, value)
            bin_idx = clamp(bin_idx, 1, length(counts))
            counts[bin_idx] += 1
        end
        positive_mask = counts .> 0
        return bin_centers[positive_mask], log2.(Float32.(counts[positive_mask]))
    end

    p_hist = Plots.plot(
        title = "Target protection spectral angle",
        xlabel = "Max decoy spectral angle",
        ylabel = "log2(count)",
        legend = :topright,
        layout = (1, 1)
    )
    if any(target_pair_mask)
        target_centers, target_log_counts = log2_bin_counts(spectral_angles[target_pair_mask])
        Plots.plot!(
            p_hist,
            target_centers,
            target_log_counts,
            marker = :circle,
            linewidth = 2,
            label = "target paired decoy"
        )
    end
    if any(target_other_mask)
        target_centers, target_log_counts = log2_bin_counts(spectral_angles[target_other_mask])
        Plots.plot!(
            p_hist,
            target_centers,
            target_log_counts,
            marker = :square,
            linewidth = 2,
            label = "target other decoy"
        )
    end
    if any(decoy_mask)
        decoy_centers, decoy_log_counts = log2_bin_counts(spectral_angles[decoy_mask])
        Plots.plot!(
            p_hist,
            decoy_centers,
            decoy_log_counts,
            marker = :diamond,
            linewidth = 2,
            label = "decoy"
        )
    end
    if has_target_neighbor && any(targets)
        target_centers, target_log_counts = log2_bin_counts(target_neighbor_spectral_angles[targets])
        Plots.plot!(
            p_hist,
            target_centers,
            target_log_counts,
            marker = :utriangle,
            linewidth = 2,
            label = "target nearest target"
        )
    end
    if has_target_neighbor && any(decoy_mask)
        decoy_centers, decoy_log_counts = log2_bin_counts(target_neighbor_spectral_angles[decoy_mask])
        Plots.plot!(
            p_hist,
            decoy_centers,
            decoy_log_counts,
            marker = :dtriangle,
            linewidth = 2,
            label = "decoy nearest target"
        )
    end
    push!(plots, p_hist)

    if has_sequence_length && length(sequence_lengths) == length(spectral_angles)
        function median_by_length(values::AbstractVector{Float32}, mask::AbstractVector{Bool})
            lengths = sort(unique(sequence_lengths[mask]))
            used_lengths = Int64[]
            medians = Float32[]
            for len in lengths
                len_mask = mask .& (sequence_lengths .== len)
                n = count(len_mask)
                n >= 10 || continue
                push!(used_lengths, Int64(len))
                push!(medians, Float32(median(values[len_mask])))
            end
            return used_lengths, medians
        end

        p_len = Plots.plot(
            title = "Median nearest neighbor spectral angle by length",
            xlabel = "Peptide length",
            ylabel = "Median nearest neighbor spectral angle",
            legend = :topright,
            layout = (1, 1)
        )
        target_pair_lengths, target_pair_medians = median_by_length(spectral_angles, target_pair_mask)
        target_other_lengths, target_other_medians = median_by_length(spectral_angles, target_other_mask)
        decoy_lengths, decoy_medians = median_by_length(spectral_angles, decoy_mask)
        if !isempty(target_pair_lengths)
            Plots.plot!(
                p_len,
                target_pair_lengths,
                target_pair_medians,
                marker = :circle,
                label = "target paired decoy"
            )
        end
        if !isempty(target_other_lengths)
            Plots.plot!(
                p_len,
                target_other_lengths,
                target_other_medians,
                marker = :square,
                label = "target other decoy"
            )
        end
        if !isempty(decoy_lengths)
            Plots.plot!(
                p_len,
                decoy_lengths,
                decoy_medians,
                marker = :diamond,
                label = "decoy"
            )
        end
        if has_target_neighbor
            target_lengths, target_medians = median_by_length(target_neighbor_spectral_angles, targets)
            decoy_target_lengths, decoy_target_medians = median_by_length(target_neighbor_spectral_angles, decoy_mask)
            if !isempty(target_lengths)
                Plots.plot!(
                    p_len,
                    target_lengths,
                    target_medians,
                    marker = :utriangle,
                    label = "target nearest target"
                )
            end
            if !isempty(decoy_target_lengths)
                Plots.plot!(
                    p_len,
                    decoy_target_lengths,
                    decoy_target_medians,
                    marker = :dtriangle,
                    label = "decoy nearest target"
                )
            end
        end
        push!(plots, p_len)
    end

    if isfile(output_path)
        safeRm(output_path, nothing; force=true)
    end
    save_multipage_pdf(plots, output_path)
    @user_info "Wrote target protection spectral angle QC plots to $output_path"
    return output_path
end

# ==========================================================
# Trace Selection
# ==========================================================

"""
    get_best_traces(second_pass_psms_paths::Vector{String}, min_prob::Float32=0.75f0)
    -> Set{@NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}}

Identify best scoring isotope traces for each precursor.

# Process
1. Accumulates scores across files using MBR_boosted_trace_prob if available, otherwise trace_prob
2. Selects highest scoring trace per precursor
3. Returns set of best precursor-isotope combinations
"""
function get_best_traces(
    second_pass_psms_paths::Vector{String},
    min_prob::Float32 = 0.75f0
)

    #The sum of scores for a given precursor trace (precursor_idx and isotopes_captured) accross the 
    #entire experiment (all runs)!
    psms_trace_scores = Dictionary{
            @NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}, Float32}()

    # Track aggregated stats
    total_rows_processed = 0
    files_processed = 0
    
    for (file_idx, file_path) in enumerate(second_pass_psms_paths)

        if splitext(file_path)[end] != ".arrow"
            continue
        end
        
        row_score = zero(Float32)
        psms_table = Arrow.Table(file_path)
        n_rows = length(psms_table[1])
        total_rows_processed += n_rows
        files_processed += 1

        # Use MBR_boosted_trace_prob if available, otherwise trace_prob
        use_mbr_column = hasproperty(psms_table, :MBR_boosted_trace_prob)
        score_column = use_mbr_column ? :MBR_boosted_trace_prob : :trace_prob

        for i in range(1, n_rows)
            psms_key = (precursor_idx = psms_table[:precursor_idx][i],  isotopes_captured = psms_table[:isotopes_captured][i])

            row_score = psms_table[score_column][i]
            if haskey(psms_trace_scores, psms_key)
                psms_trace_scores[psms_key] = psms_trace_scores[psms_key] + row_score
            else
                insert!(
                    psms_trace_scores,
                    psms_key,
                    row_score
                )
            end
        end
    end

    #Convert the dictionary to a DataFrame 
    psms_trace_df = DataFrame(
        (precursor_idx = [key[:precursor_idx] for key in keys(psms_trace_scores)],
        isotopes_captured = [key[:isotopes_captured] for key in keys(psms_trace_scores)],
        score = [val for val in values(psms_trace_scores)])
        );
    #Now retain only the very best trace!
    psms_trace_df[!,:best_trace] .= false;
    gpsms = groupby(psms_trace_df,:precursor_idx)
    for (precursor_idx, psms) in pairs(gpsms)
        psms[argmax(psms[!,:score]),:best_trace] = true
    end
    filter!(x->x.best_trace, psms_trace_df);
    traces_passing = Set([(precursor_idx = x.precursor_idx, isotopes_captured = x.isotopes_captured) for x in eachrow(psms_trace_df)]);
    return traces_passing
end

# ==========================================================
# PSM score merging and processing
# ==========================================================


"""
    get_pep_spline(merged_psms_path::String, score_col::Symbol;
                   min_pep_points_per_bin=5000, n_spline_bins=20) -> UniformSpline

Create posterior error probability spline from merged scores.

Returns spline for PEP calculation based on target/decoy distributions.
"""
function get_pep_spline(
                            merged_psms_path::String,
                            score_col::Symbol;
                            min_pep_points_per_bin = 5000,
                            n_spline_bins = 20
)

    psms_scores = Arrow.Table(merged_psms_path)
    Q = length(psms_scores[1])
    M = ceil(Int, Q / min_pep_points_per_bin)
    bin_target_fraction, bin_mean_prob = Vector{Float32}(undef, M), Vector{Float32}(undef, M)
    bin_size = 0
    bin_idx = 0
    mean_prob, targets = 0.0f0, 0
    for i in range(1, Q)
        bin_size += 1
        targets += psms_scores[:target][i]
        mean_prob += psms_scores[score_col][i]
        if bin_size == min_pep_points_per_bin
            bin_idx += 1
            bin_target_fraction[bin_idx] = targets/bin_size
            bin_mean_prob[bin_idx] = mean_prob/bin_size
            bin_size, targets, mean_prob = zero(Int64), zero(Int64), zero(Float32)
        end
    end
    bin_target_fraction[end] = targets/max(bin_size, 1)
    bin_mean_prob[end] = mean_prob/max(bin_size, 1)
    try 
        if length(bin_target_fraction)<20
            @user_warn "Less than 20 bins to estimate PEP. PEP results suspect..."
        end
        return UniformSpline(bin_target_fraction, bin_mean_prob, 3, 3)
    catch
        @user_warn "Failed to estimate PEP spline"
        return UniformSpline(SVector{4, Float32}([0, 0, 0, 0]), 3, 0.0f0, 1.0f0, 100.0f0)
    end
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
        @user_warn "Insufficient unique points for PEP interpolation, using default"
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
        @user_warn "Insufficient unique points for q-value interpolation, using default"
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
