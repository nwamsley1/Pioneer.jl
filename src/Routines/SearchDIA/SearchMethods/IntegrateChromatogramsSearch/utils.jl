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
    build_chrom_index(prec_col::AbstractVector{UInt32}, n_rows::Int)

Build a lookup from precursor ID → row range in a sorted chromatogram table.

The chromatogram DataFrame is sorted by `(precursor_idx, rt)`, so all rows for
a given precursor are contiguous. This function does a single pass to record
each precursor's `start:stop` row range and tracks the longest chromatogram
(used to size scratch buffers).

Returns `(chrom_index, max_chrom_len)`.
"""
function build_chrom_index(prec_col::AbstractVector{UInt32}, n_rows::Int)
    if n_rows == 0
        return Dict{UInt32, UnitRange{Int}}(), 0
    end
    max_len = 0
    chrom_index = Dict{UInt32, UnitRange{Int}}()
    start_idx = 1
    current_key = prec_col[1]
    for row_idx in 2:n_rows
        next_key = prec_col[row_idx]
        if next_key != current_key
            chrom_index[current_key] = start_idx:(row_idx - 1)
            max_len = max(max_len, row_idx - start_idx)
            start_idx = row_idx
            current_key = next_key
        end
    end
    chrom_index[current_key] = start_idx:n_rows
    max_len = max(max_len, n_rows - start_idx + 1)
    return chrom_index, max_len
end

function build_chrom_index(chromatograms::DataFrame, ::CombineTraces)
    return build_chrom_index(
        chromatograms[!, :precursor_idx]::AbstractVector{UInt32},
        nrow(chromatograms),
    )
end

function build_chrom_index(chromatograms::DataFrame, ::SeperateTraces)
    n_rows = nrow(chromatograms)
    if n_rows == 0
        return Dict{Tuple{UInt32, Tuple{Int8, Int8}}, UnitRange{Int}}(), 0
    end

    prec_col = chromatograms[!, :precursor_idx]::AbstractVector{UInt32}
    iso_col = chromatograms[!, :isotopes_captured]::AbstractVector{Tuple{Int8, Int8}}
    max_len = 0
    chrom_index = Dict{Tuple{UInt32, Tuple{Int8, Int8}}, UnitRange{Int}}()
    start_idx = 1
    current_key = (prec_col[1], iso_col[1])
    for row_idx in 2:n_rows
        next_key = (prec_col[row_idx], iso_col[row_idx])
        if next_key != current_key
            chrom_index[current_key] = start_idx:(row_idx - 1)
            max_len = max(max_len, row_idx - start_idx)
            start_idx = row_idx
            current_key = next_key
        end
    end
    chrom_index[current_key] = start_idx:n_rows
    max_len = max(max_len, n_rows - start_idx + 1)
    return chrom_index, max_len
end

"""
    find_nearest_scan(scan_indices::AbstractVector{UInt32}, target_scan::UInt32)

Return the 1-based position within `scan_indices` whose value is closest to
`target_scan`. Used to map the PSM's apex scan index (global scan number) into
the local chromatogram coordinate for a sorted precursor trace.
"""
function find_nearest_scan(scan_indices::AbstractVector{<:Integer}, target_scan::Integer)
    nearest_idx = 1
    min_diff = typemax(Int64)
    for j in eachindex(scan_indices)
        d = abs(Int64(scan_indices[j]) - Int64(target_scan))
        if d < min_diff
            min_diff = d
            nearest_idx = j
        end
    end
    return nearest_idx
end

function sort_chromatograms_for_integration!(chromatograms::DataFrame, ::CombineTraces)
    nrow(chromatograms) == 0 && return chromatograms
    fast_df_sort!(chromatograms, [:precursor_idx, :rt])
    return chromatograms
end

function sort_chromatograms_for_integration!(chromatograms::DataFrame, ::SeperateTraces)
    nrow(chromatograms) == 0 && return chromatograms
    fast_df_sort!(chromatograms, [:precursor_idx, :isotopes_captured, :rt])
    return chromatograms
end

function sort_chromatograms_for_integration!(chromatograms::DataFrame)
    return sort_chromatograms_for_integration!(chromatograms, CombineTraces(0.25f0))
end

"""
    compute_chromatogram_isotope_sets(trace_type)

Return whether chromatogram rows need a temporary `:isotopes_captured` grouping
column.

All trace modes still need row-level `:precursor_fraction_transmitted`; separate
trace mode additionally needs row-level isotope labels so rows can be grouped by
`(precursor_idx, isotopes_captured)`.
"""
compute_chromatogram_isotope_sets(trace_type::IsotopeTraceType) = seperateTraces(trace_type)

# Explicit debug-only target list for writing a small set of chromatogram plots.
# Plot generation is still gated by `DEBUG_CONSOLE_LEVEL[] >= 1`.
const DEBUG_CHROM_TARGET_PRECURSOR_IDXS = Ref{Set{UInt32}}(Set{UInt32}())
const DEBUG_CHROM_TARGET_PRECURSOR_IDX = DEBUG_CHROM_TARGET_PRECURSOR_IDXS

@inline function debug_should_plot_chromatogram(precursor_idx::Integer)
    return DEBUG_CONSOLE_LEVEL[] >= 1 &&
        UInt32(precursor_idx) in DEBUG_CHROM_TARGET_PRECURSOR_IDXS[]
end

function debug_sanitize_chromatogram_filename(value::AbstractString)
    safe = replace(String(value), r"[^A-Za-z0-9]+" => "_")
    safe = replace(safe, r"^_+|_+$" => "")
    return isempty(safe) ? "run" : safe
end

function debug_chromatogram_plot_path(
    output_root::AbstractString,
    ms_file_idx::Integer,
    file_name::AbstractString,
    precursor_idx::Integer,
)
    safe_name = debug_sanitize_chromatogram_filename(file_name)
    return joinpath(
        String(output_root),
        "qc_plots",
        "chromatogram_integration_debug",
        "precursor_$(UInt32(precursor_idx))_file_$(Int(ms_file_idx))_$(safe_name).png",
    )
end

"""
    debug_write_target_chromatogram_plots(...)

Write targeted chromatogram debug plots only when debug logging is enabled and
the precursor is listed in `DEBUG_CHROM_TARGET_PRECURSOR_IDXS`. Plots rerun the
same integration path used for quantification and show the original expected
PSM apex as the red line.
"""
function debug_write_target_chromatogram_plots(
    chromatograms::DataFrame,
    passing_psms::DataFrame,
    min_fraction_transmitted::Float32,
    λ::Float32,
    output_root::AbstractString,
    ms_file_idx::Integer,
    file_name::AbstractString,
    ;
    debug_plot_data_by_precursor::Union{Nothing,AbstractDict} = nothing,
)
    DEBUG_CONSOLE_LEVEL[] >= 1 || return nothing

    target_precursor_idxs = DEBUG_CHROM_TARGET_PRECURSOR_IDXS[]
    isempty(target_precursor_idxs) && return nothing
    target_rows = findall(
        row -> UInt32(passing_psms[row, :precursor_idx]) in target_precursor_idxs,
        axes(passing_psms, 1),
    )
    isempty(target_rows) && return nothing

    nrow(chromatograms) == 0 && return nothing

    prec_col = chromatograms[!, :precursor_idx]::AbstractVector{UInt32}
    rt_all = chromatograms[!, :rt]::AbstractVector{Float32}
    scan_idx_all = chromatograms[!, :scan_idx]::AbstractVector{UInt32}
    intensity_all = chromatograms[!, :intensity]::AbstractVector{Float32}
    fraction_all = chromatograms[!, :precursor_fraction_transmitted]::AbstractVector{Float32}

    for row_idx in target_rows
        target_precursor_idx = UInt32(passing_psms[row_idx, :precursor_idx])
        chrom_rows = findall(==(target_precursor_idx), prec_col)
        isempty(chrom_rows) && continue

        ordered_rows = chrom_rows[sortperm(@view(rt_all[chrom_rows]))]
        N = length(ordered_rows)
        N <= 0 && continue
        avg_cycle_time = N < 2 ? 1.0f0 : (rt_all[last(ordered_rows)] - rt_all[first(ordered_rows)]) / N
        expected_apex_scan = hasproperty(passing_psms, :scan_idx) ?
            UInt32(passing_psms[row_idx, :scan_idx]) :
            UInt32(0)
        selected_apex_scan = hasproperty(passing_psms, :new_best_scan) ?
            UInt32(passing_psms[row_idx, :new_best_scan]) :
            UInt32(0)
        target_scan = selected_apex_scan != UInt32(0) ?
            selected_apex_scan :
            (
                expected_apex_scan != UInt32(0) ?
                    expected_apex_scan :
                    scan_idx_all[ordered_rows[argmax(@view(intensity_all[ordered_rows]))]]
            )
        apex_scan = find_nearest_scan(scan_idx_all[ordered_rows], target_scan)
        debug_apex_scan = expected_apex_scan != UInt32(0) ?
            find_nearest_scan(scan_idx_all[ordered_rows], expected_apex_scan) :
            apex_scan
        plot_path = debug_chromatogram_plot_path(
            output_root,
            ms_file_idx,
            file_name,
            target_precursor_idx,
        )
        plot_title = "file $(Int(ms_file_idx)) $(file_name) precursor $(target_precursor_idx)"
        plot_data_ref = debug_plot_data_by_precursor === nothing ? nothing : Ref{Any}(nothing)

        ws = WHWorkspace(N)
        state = Chromatogram(zeros(Float32, N), zeros(Float32, N), 0)
        integrate_chrom(
            rt_all[ordered_rows],
            scan_idx_all[ordered_rows],
            intensity_all[ordered_rows],
            fraction_all[ordered_rows],
            apex_scan,
            ws,
            state,
            avg_cycle_time,
            λ,
            min_fraction_transmitted = min_fraction_transmitted,
            debug_plot_path = plot_path,
            debug_plot_title = plot_title,
            debug_plot_data = plot_data_ref,
            debug_apex_scan = debug_apex_scan,
        )
        if debug_plot_data_by_precursor !== nothing
            debug_plot_data_by_precursor[target_precursor_idx] = plot_data_ref[]
        end
    end

    return nothing
end

"""
    select_quant_trace_by_transmission(chromatograms)

For separate-trace quantification, choose one isotope-capture trace per
precursor using the highest precursor transmission fraction. Ties are resolved
by the isotope tuple for deterministic output.
"""
function select_quant_trace_by_transmission(chromatograms::DataFrame)
    prec_col = chromatograms[!, :precursor_idx]::AbstractVector{UInt32}
    iso_col = chromatograms[!, :isotopes_captured]::AbstractVector{Tuple{Int8, Int8}}
    fraction_col = chromatograms[!, :precursor_fraction_transmitted]::AbstractVector{Float32}
    selected = Dict{UInt32, Tuple{Tuple{Int8, Int8}, Float32}}()

    for i in eachindex(prec_col)
        pid = prec_col[i]
        iso = iso_col[i]
        fraction = fraction_col[i]
        if !haskey(selected, pid) || fraction > selected[pid][2] ||
           (fraction == selected[pid][2] && iso < selected[pid][1])
            selected[pid] = (iso, fraction)
        end
    end

    return selected
end

"""
    apply_quant_trace_selection!(psms, selected_quant_trace) -> Vector{Tuple{Int8, Int8}}

Apply the selected separate-trace quantification channel to the PSM table.

Updates `:precursor_fraction_transmitted` in place — that is a real output
column, and in separate-trace mode it must report the fraction of the trace
actually quantified — and RETURNS the per-PSM isotope-capture key rather than
writing it into the frame.

The isotope key is deliberately not a column. `:isotopes_captured` was retired
from the PSM schema, and its only live consumer here is
`chromatogram_index_key(::SeperateTraces, ...)`, which uses it to look rows up
in the chromatogram index. Nothing reads it back out: the output alias in
writeCSVTables maps `:isotopes_captured` to a source column
`:isotopes_captured_traces` that nothing in the pipeline produces, so the rename
guarded on its presence never fires. Materialising a column would therefore add
a per-PSM field to the staged table that only this function writes and only this
function reads — and would make the separate-trace output schema differ from
combined for no gain. A local vector costs two bytes a row and only in separate
mode.

Precursors with no chromatogram row get an unmatchable `(-1, -1)` key. That is
the correct outcome rather than a fallback: there is no extracted chromatogram
to integrate, so the precursor keeps zero area and is not quantifiable.
"""
function apply_quant_trace_selection!(psms::DataFrame, selected_quant_trace::Dict{UInt32, Tuple{Tuple{Int8, Int8}, Float32}})
    prec_col = psms[!, :precursor_idx]::AbstractVector{UInt32}
    fraction_col = psms[!, :precursor_fraction_transmitted]::AbstractVector{Float32}
    isotopes_captured = Vector{Tuple{Int8, Int8}}(undef, length(prec_col))

    @inbounds for i in eachindex(prec_col)
        selected = get(selected_quant_trace, prec_col[i], nothing)
        if selected === nothing
            isotopes_captured[i] = (Int8(-1), Int8(-1))
            continue
        end
        isotopes_captured[i] = selected[1]
        fraction_col[i] = selected[2]
    end

    return isotopes_captured
end

function combined_trace_seed_score_column(psms::DataFrame)
    hasproperty(psms, :lgbm_score) && return :lgbm_score
    throw(ArgumentError(
        "Combined trace seed selection requires :lgbm_score on passing PSMs",
    ))
end

"""
    select_combined_trace_seed_rows_by_score(psms)

For combined-trace chromatogram integration, choose one MainSearch seed row per
precursor. The selected `scan_idx` is the scan whose PSM had the highest
LightGBM score, not the final chromatographic apex. Ties are resolved by lower
`scan_idx`.
"""
function select_combined_trace_seed_rows_by_score(
    precursor_idx::AbstractVector{UInt32},
    scan_idx::AbstractVector{UInt32},
    score::AbstractVector{<:Real},
)
    selected = Dict{UInt32, Tuple{Float32, UInt32, Int}}()

    for row in eachindex(precursor_idx)
        pid = precursor_idx[row]
        scan = scan_idx[row]
        score_val = Float32(score[row])
        current = get(selected, pid, nothing)

        if current === nothing ||
           score_val > current[1] ||
           (score_val == current[1] && scan < current[2])
            selected[pid] = (score_val, scan, Int(row))
        end
    end

    selected_rows = Int[current[3] for current in values(selected)]
    sort!(selected_rows)
    return selected_rows
end

function select_combined_trace_seed_rows_by_score(psms::DataFrame)
    score_col = psms[!, combined_trace_seed_score_column(psms)]
    return select_combined_trace_seed_rows_by_score(
        psms[!, :precursor_idx]::AbstractVector{UInt32},
        psms[!, :scan_idx]::AbstractVector{UInt32},
        score_col,
    )
end

function select_combined_trace_seed_psms_by_score(psms::DataFrame)
    return psms[select_combined_trace_seed_rows_by_score(psms), :]
end

@inline function chromatogram_index_key(
    ::CombineTraces,
    precursor_idx::AbstractVector{UInt32},
    isotopes_captured,
    row_idx::Integer,
)
    return precursor_idx[row_idx]
end

@inline function chromatogram_index_key(
    ::SeperateTraces,
    precursor_idx::AbstractVector{UInt32},
    isotopes_captured::AbstractVector{Tuple{Int8, Int8}},
    row_idx::Integer,
)
    return (precursor_idx[row_idx], isotopes_captured[row_idx])
end

"""
    integrate_precursors(chromatograms, isotope_trace_type, min_fraction_transmitted, precursor_idx,
                         apex_scan_idx, peak_area, new_best_scan, points_integrated,
                         integration_start_scan, integration_stop_scan;
                         isotopes_captured=nothing, λ=1.0f0)

Integrate chromatographic peaks for multiple precursors in parallel.

The input chromatogram table must already be sorted for `isotope_trace_type`
and must contain row-level `:precursor_fraction_transmitted`. Combined mode
uses one trace per precursor; separate mode keys rows by
`(precursor_idx, isotopes_captured)` and receives the PSM-level
`isotopes_captured` vector identifying the chosen quantification trace.

For each precursor, this maps the MainSearch seed scan into the local trace and
calls `integrate_chrom` (WH smoothing -> second-derivative bounds -> baseline
subtraction -> trapezoidal integration). Results are written into
`peak_area`, `new_best_scan`, `points_integrated`, and the integration-boundary
scan indices, along with `peak_area_unsubtracted` -- the area over the same
window before baseline subtraction, which is what the not-quantifiable rule in
`integrate_chrom` tests against.
"""
function integrate_precursors(chromatograms::DataFrame,
                             isotope_trace_type::IsotopeTraceType,
                             min_fraction_transmitted::Float32,
                             precursor_idx::AbstractVector{UInt32},
                             apex_scan_idx::AbstractVector{UInt32},
                             peak_area::AbstractVector{Float32},
                             new_best_scan::AbstractVector{UInt32},
                             points_integrated::AbstractVector{UInt32},
                             integration_start_scan::AbstractVector{UInt32},
                             integration_stop_scan::AbstractVector{UInt32},
                             peak_area_unsubtracted::AbstractVector{Float32};
                             isotopes_captured = nothing,
                             λ::Float32 = 1.0f0,
                             )
    n_pad = Int64(0)
    rt_all = chromatograms[!, :rt]::AbstractVector{Float32}
    scan_idx_all = chromatograms[!, :scan_idx]::AbstractVector{UInt32}
    intensity_all = chromatograms[!, :intensity]::AbstractVector{Float32}
    fraction_all = chromatograms[!, :precursor_fraction_transmitted]::AbstractVector{Float32}

    chrom_index, max_chrom_len = build_chrom_index(chromatograms, isotope_trace_type)
    N = max_chrom_len + (2*n_pad)

    # Partition work into chunks, then distribute chunks to exactly nthreads() tasks
    # (one task per thread — prevents workspace corruption from task switching at yield points)
    _n_threads = Threads.nthreads()
    all_chunks = collect(partitionThreadTasks(length(precursor_idx), 10, _n_threads))
    n_tasks = min(_n_threads, length(all_chunks))

    task_ranges = [Int[] for _ in 1:n_tasks]
    for (chunk_idx, _) in enumerate(all_chunks)
        push!(task_ranges[((chunk_idx - 1) % n_tasks) + 1], chunk_idx)
    end

    n_thread_slots = Threads.maxthreadid()
    ws_by_thread = [WHWorkspace(N) for _ in 1:n_thread_slots]
    state_by_thread = [Chromatogram(zeros(Float32, N), zeros(Float32, N), 0) for _ in 1:n_thread_slots]

    function run_integration_batch!(batch_id::Int, ws::WHWorkspace, state::Chromatogram)
        for chunk_idx in task_ranges[batch_id]
            chunk = all_chunks[chunk_idx]
            for i in chunk
                apex_scan = apex_scan_idx[i]
                chrom_key = chromatogram_index_key(
                    isotope_trace_type,
                    precursor_idx,
                    isotopes_captured,
                    i,
                )

                if !haskey(chrom_index, chrom_key)
                    continue
                end
                chrom_range = chrom_index[chrom_key]

                avg_cycle_time = (rt_all[last(chrom_range)] - rt_all[first(chrom_range)]) / length(chrom_range)

                # Skip entirely zero-valued chromatograms
                if isnothing(findfirst(x -> x > 0.0f0, @view intensity_all[chrom_range]))
                    continue
                end

                # Map global apex scan index to local position in chromatogram
                apex_scan = find_nearest_scan(
                    @view(scan_idx_all[chrom_range]), apex_scan
                )

                peak_area[i], new_best_scan[i], points_integrated[i],
                    integration_start_scan[i], integration_stop_scan[i],
                    _, _, peak_area_unsubtracted[i], _ =
                    integrate_chrom(
                    @view(rt_all[chrom_range]),
                    @view(scan_idx_all[chrom_range]),
                    @view(intensity_all[chrom_range]),
                    @view(fraction_all[chrom_range]),
                    apex_scan,
                    ws,
                    state,
                    avg_cycle_time,
                    λ,
                    min_fraction_transmitted = min_fraction_transmitted,
                    n_pad = n_pad,
                )
                reset!(state)
            end
        end
        return nothing
    end

    Threads.@threads :static for batch_id in 1:n_tasks
        tid = Threads.threadid()
        run_integration_batch!(batch_id, ws_by_thread[tid], state_by_thread[tid])
    end

    # Clamp NaN and negative values to zero for downstream processing
    for i in eachindex(precursor_idx)
        if isnan(peak_area[i]) || peak_area[i] < 0f0
            peak_area[i] = 0f0
        end
    end

    return nothing
end

function integrate_precursors(chromatograms::DataFrame,
                             min_fraction_transmitted::Float32,
                             precursor_idx::AbstractVector{UInt32},
                             apex_scan_idx::AbstractVector{UInt32},
                             peak_area::AbstractVector{Float32},
                             new_best_scan::AbstractVector{UInt32},
                             points_integrated::AbstractVector{UInt32},
                             integration_start_scan::AbstractVector{UInt32},
                             integration_stop_scan::AbstractVector{UInt32},
                             peak_area_unsubtracted::AbstractVector{Float32};
                             λ::Float32 = 1.0f0,
                             )
    return integrate_precursors(
        chromatograms,
        CombineTraces(Float32(min_fraction_transmitted)),
        min_fraction_transmitted,
        precursor_idx,
        apex_scan_idx,
        peak_area,
        new_best_scan,
        points_integrated,
        integration_start_scan,
        integration_stop_scan,
        peak_area_unsubtracted;
        λ = λ,
    )
end

# Rank weights for the 8 smoothed shadow fragments. `_fragment_rank_weights` is deterministic in its
# argument, so the 8-fragment result is a constant; computing it per precursor group allocated one
# 8-element Vector per group for nothing. Read-only -- never mutate.
const MBR_FRAGMENT_RANK_WEIGHTS_8 = _fragment_rank_weights(8)

@inline function _chrom_shadow_intensity(
    shadow_columns,
    row::Integer,
    rank::Integer,
)
    value = Float32(shadow_columns[rank][row])
    return isfinite(value) ? max(value, 0.0f0) : 0.0f0
end

@inline function _chrom_positive_intensity(
    columns,
    row::Integer,
    rank::Integer,
)
    value = Float32(columns[rank][row])
    return isfinite(value) ? max(value, 0.0f0) : 0.0f0
end

@inline function _chrom_spectrum_sqrt_tuple_and_sum(columns, row::Integer)
    values = ntuple(
        rank -> _chrom_positive_intensity(columns, row, rank),
        8,
    )
    total = sum(values)
    total > 0.0f0 ||
        return MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT, 0.0f0
    inverse_total = inv(total)
    return ntuple(
        rank -> sqrt(values[rank] * inverse_total),
        8,
    ), total
end

@inline function _chrom_hellinger_score_from_sqrt(
    observed::NTuple{8, Float32},
    fitted::NTuple{8, Float32},
)
    (_mbr_sqrt_tuple_is_valid(observed) &&
     _mbr_sqrt_tuple_is_valid(fitted)) || return 0.0f0
    bhattacharyya = 0.0f0
    @inbounds for rank in 1:8
        bhattacharyya += observed[rank] * fitted[rank]
    end
    hellinger_squared = max(0.0f0, 1.0f0 - bhattacharyya)
    return Float32(-log2(max(hellinger_squared, 1.0f-10)))
end

@inline function _chrom_manhattan_score(
    observed_columns,
    fitted_columns,
    row::Integer,
)
    observed = ntuple(
        rank -> _chrom_positive_intensity(observed_columns, row, rank),
        8,
    )
    fitted = ntuple(
        rank -> _chrom_positive_intensity(fitted_columns, row, rank),
        8,
    )
    observed_sum = sum(observed)
    observed_sum > 0.0f0 || return 0.0f0
    distance = 0.0f0
    @inbounds for rank in 1:8
        distance += abs(fitted[rank] - observed[rank])
    end
    return Float32(-log2(distance / observed_sum + 1.0f-10))
end

"""
    _ChromCorrWorkspace

Reusable scratch for `_chrom_shadow_fragment_correlation_features`, which runs once per precursor
group. It previously allocated 11 vectors per call (8 fragment traces + their outer `Vector`, the
weight trace, and the 8-element correlation vector); on a 6-file KEAP1 run that was 742,179 groups.

Construct it with the largest group length the caller will ever pass (`_ChromCorrWorkspace(max_pts)`)
so the buffers are allocated exactly once. Each call then `resize!`s to its own `n_points`, which
stays within that capacity: shrinking is a length update and regrowing never reallocates, so the
steady state is genuinely zero-allocation rather than merely amortised.

Nothing here escapes the function — only the `isbits` `_ChromCorrFeatures` result does — so one
workspace shared across every group of a kernel invocation is safe. The kernel's group loop is
sequential; a threaded caller would need one workspace per thread.
"""
struct _ChromCorrWorkspace
    fragment_traces::NTuple{8, Vector{Float32}}
    weight_trace::Vector{Float32}
    correlation_to_weight::Vector{Float32}
end

_ChromCorrWorkspace(max_points::Integer = 0) = _ChromCorrWorkspace(
    ntuple(_ -> Vector{Float32}(undef, max_points), 8),
    Vector{Float32}(undef, max_points),
    Vector{Float32}(undef, 8),
)

function _chrom_shadow_fragment_correlation_features(
    rows::AbstractVector{<:Integer},
    shadow_columns,
    weight_column,
    bitvec_rank_table,
    workspace::_ChromCorrWorkspace,
)
    n_points = length(rows)
    n_points < 2 && return (
        n_correlated = UInt8(0),
        corr_mask = UInt8(0),
        bitvec_rank = UInt16(0),
        strength = 0.0f0,
        effective_n = 0.0f0,
        best_weight = 0.0f0,
    )

    # Resized to EXACTLY n_points rather than viewed: `_frag_pcor` and `maximum` then see a plain
    # `Vector{Float32}` instead of a `SubArray`, which measured ~80% faster in this phase. The
    # workspace is pre-sized to the largest group, so these resizes stay within capacity and do not
    # allocate; every element is overwritten below, so moving data would be harmless anyway.
    fragment_traces = workspace.fragment_traces
    @inbounds for rank in 1:8
        resize!(fragment_traces[rank], n_points)
    end
    weight_trace = workspace.weight_trace
    resize!(weight_trace, n_points)
    @inbounds for point_idx in 1:n_points
        chrom_row = rows[point_idx]
        weight_value = Float32(weight_column[chrom_row])
        weight_trace[point_idx] =
            isfinite(weight_value) ? max(weight_value, 0.0f0) : 0.0f0
        for rank in 1:8
            fragment_traces[rank][point_idx] =
                _chrom_shadow_intensity(shadow_columns, chrom_row, rank)
        end
    end

    has_signal = ntuple(
        rank -> maximum(fragment_traces[rank]) > 0.0f0,
        8,
    )
    correlation_to_weight = workspace.correlation_to_weight
    @inbounds for rank in 1:8
        correlation_to_weight[rank] = has_signal[rank] ?
            _frag_pcor(fragment_traces[rank], weight_trace) :
            0.0f0
    end

    n_correlated = UInt8(0)
    correlation_mask = UInt8(0)
    @inbounds for rank in 1:8
        has_signal[rank] || continue
        if correlation_to_weight[rank] > 0.7f0
            n_correlated += UInt8(1)
            correlation_mask |= UInt8(1) << (rank - 1)
        end
    end
    # `_fragment_rank_weights(8)` is a pure function of the constant 8, but allocated a fresh
    # 8-element Vector on every call -- once per precursor group. Hoisted to a const.
    strength, effective_n =
        _positive_corr_summary(correlation_to_weight, MBR_FRAGMENT_RANK_WEIGHTS_8)

    best_rank = 0
    best_consensus = typemin(Float32)
    @inbounds for rank in 1:8
        has_signal[rank] || continue
        consensus = 0.0f0
        n_pairs = 0
        for other_rank in 1:8
            (other_rank == rank || !has_signal[other_rank]) && continue
            consensus += _frag_pcor(
                fragment_traces[rank],
                fragment_traces[other_rank],
            )
            n_pairs += 1
        end
        average = n_pairs > 0 ?
            consensus / n_pairs :
            typemin(Float32)
        if average > best_consensus
            best_consensus = average
            best_rank = rank
        end
    end
    best_weight = best_rank > 0 ?
        correlation_to_weight[best_rank] :
        0.0f0
    return (
        n_correlated = n_correlated,
        corr_mask = correlation_mask,
        bitvec_rank =
            _bitvec_pattern_rank(bitvec_rank_table, correlation_mask),
        strength = strength,
        effective_n = effective_n,
        best_weight = best_weight,
    )
end

@inline function _chrom_smoothed_shadow_sqrt_tuple_and_sum(
    shadow_columns,
    apex_row::Integer,
    left_row::Integer,
    right_row::Integer,
)
    # Single assignment on purpose. `denominator` is captured by the `ntuple` closure below, and a
    # captured variable that is *reassigned* in the enclosing scope gets boxed as `Core.Box` (type
    # `Any`). That boxing made `signal / denominator` a runtime dispatch, `values` an
    # `NTuple{8, Any}`, and `sum(values)` a dynamic `afoldl` chain -- once per PSM row.
    denominator =
        2.0f0 + (left_row > 0 ? 1.0f0 : 0.0f0) + (right_row > 0 ? 1.0f0 : 0.0f0)
    values = ntuple(rank -> begin
        signal =
            2.0f0 *
            _chrom_shadow_intensity(shadow_columns, apex_row, rank)
        left_row > 0 &&
            (signal += _chrom_shadow_intensity(
                shadow_columns,
                left_row,
                rank,
            ))
        right_row > 0 &&
            (signal += _chrom_shadow_intensity(
                shadow_columns,
                right_row,
                rank,
            ))
        Float32(signal / denominator)
    end, 8)
    total = sum(values)
    total > 0.0f0 ||
        return MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT, 0.0f0
    inverse_total = inv(total)
    return ntuple(
        rank -> sqrt(values[rank] * inverse_total),
        8,
    ), total
end

function _chrom_temporal_fragment_data(
    rows::AbstractVector{<:Integer},
    shadow_columns,
    weight_column,
    scan_column,
    start_scan::UInt32,
    stop_scan::UInt32,
)
    (start_scan == UInt32(0) || stop_scan == UInt32(0)) &&
        return Float32[], MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT
    scan_low, scan_high = minmax(start_scan, stop_scan)

    # Count qualifying scans first, then fill by index: the trace is exactly
    # MBR_TEMPORAL_TRACE_STRIDE (=9) floats per qualifying scan, so `push!` and its geometric
    # regrowth are unnecessary. Predicate must match the fill loop below exactly.
    n_qual = 0
    @inbounds for chrom_row in rows
        scan = UInt32(scan_column[chrom_row])
        scan_low <= scan <= scan_high || continue
        wv = Float32(weight_column[chrom_row])
        (isfinite(wv) ? max(wv, 0.0f0) : 0.0f0) > 0.0f0 || continue
        n_qual += 1
    end
    trace = Vector{Float32}(undef, MBR_TEMPORAL_TRACE_STRIDE * n_qual)
    pos = 0
    weighted_sqrt = zeros(Float32, 8)
    total_weight = 0.0f0

    @inbounds for chrom_row in rows
        scan = UInt32(scan_column[chrom_row])
        scan_low <= scan <= scan_high || continue
        weight_value = Float32(weight_column[chrom_row])
        weight = isfinite(weight_value) ?
            max(weight_value, 0.0f0) :
            0.0f0
        weight > 0.0f0 || continue
        trace[pos + 1] = weight

        spectrum_total = 0.0f0
        for rank in 1:8
            intensity =
                _chrom_shadow_intensity(shadow_columns, chrom_row, rank)
            trace[pos + 1 + rank] = intensity
            spectrum_total += intensity
        end
        pos += MBR_TEMPORAL_TRACE_STRIDE
        spectrum_total > 0.0f0 || continue
        inverse_total = inv(spectrum_total)
        for rank in 1:8
            intensity =
                _chrom_shadow_intensity(shadow_columns, chrom_row, rank)
            weighted_sqrt[rank] +=
                weight * sqrt(intensity * inverse_total)
        end
        total_weight += weight
    end
    # `total_weight` is reassigned in the loop above, so capturing it directly in the `ntuple`
    # closure boxed it as `Core.Box` (type `Any`) -- making the division a runtime dispatch and the
    # result an `NTuple{8, Any}`, once per PSM row. Copying to a single-assignment local removes the
    # capture of the reassigned variable, so nothing is boxed. (`weighted_sqrt` is only mutated
    # through `setindex!`, never reassigned, so capturing it was always fine.)
    #
    # Same bug, same file, as the `denominator` box in
    # `_chrom_smoothed_shadow_sqrt_tuple_and_sum`. Kept as a division rather than multiplying by
    # `inv(total_weight)` so the result stays bit-identical.
    weight_total = total_weight
    temporal_mean = weight_total > 0.0f0 ?
        ntuple(rank -> weighted_sqrt[rank] / weight_total, 8) :
        MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT
    return trace, temporal_mean
end

# Concrete replacement for the values previously stored in a Dict{UInt32, NamedTuple}; a bare
# NamedTuple is an abstract type, so every value was boxed and retrieved by dynamic dispatch.
struct _ChromCorrFeatures
    n_correlated::UInt8
    corr_mask::UInt8
    bitvec_rank::UInt16
    strength::Float32
    effective_n::Float32
    best_weight::Float32
end

# (apex, left, right) chromatogram rows for `selected_scan`, reproducing the old
# `sort!(rows; by = scan_idx)` + adjacent-position logic without the per-row Dict.
#
# `order === nothing` means the stored row order is already scan order (CombineTraces sorts
# [:precursor_idx, :rt], and rt is monotonic in scan_idx), so positions ARE row indices.
# Otherwise `order` is a permutation whose [lo, hi] slice lists the precursor's rows in scan order
# (SeperateTraces sorts [:precursor_idx, :isotopes_captured, :rt], so scan_idx is not monotonic
# across the precursor's range and the rows must be re-ordered exactly as the old code did).
@inline function _chrom_neighbor_rows(
    scan_column, order, lo::Int, hi::Int, selected_scan::UInt32,
)
    lo > hi && return (0, 0, 0)
    @inline row_at(pos::Int) = order === nothing ? pos : Int(order[pos])
    left, right = lo, hi
    @inbounds while left <= right          # binary search over scan order
        mid = (left + right) >>> 1
        if UInt32(scan_column[row_at(mid)]) < selected_scan
            left = mid + 1
        else
            right = mid - 1
        end
    end
    pos = left
    (pos > hi || @inbounds(UInt32(scan_column[row_at(pos)])) != selected_scan) &&
        return (0, 0, 0)
    return (
        row_at(pos),
        pos > lo ? row_at(pos - 1) : 0,
        pos < hi ? row_at(pos + 1) : 0,
    )
end

# Contiguous sub-range of [lo, hi] whose isotopes_captured equals `iso` (middle sort key).
@inline function _chrom_isotope_subrange(iso_column, lo::Int, hi::Int, iso)
    lo > hi && return 1:0
    s = 0; e = -1
    @inbounds for i in lo:hi
        if Tuple{Int8, Int8}(iso_column[i]) == iso
            s == 0 && (s = i)
            e = i
        elseif s != 0
            break                          # contiguous, so stop at the first miss after a hit
        end
    end
    return s == 0 ? (1:0) : (s:e)
end

"""
    add_mbr_integrated_spectra_to_psms!(psms, chromatograms, rt_to_irt)

Attach the receiver chromatogram evidence consumed by post-integration MBR.
The temporal trace is restricted to the exact peak-integration bounds.
"""
function add_mbr_integrated_spectra_to_psms!(
    passing_psms::DataFrame,
    chromatograms::DataFrame,
    rt_to_irt;
    bitvec_rank_table = nothing,
    scan_tic::Union{Nothing, Vector{Float32}} = nothing,
)
    n = nrow(passing_psms)
    for column in MBR_INTEGRATED_FRAGMENT_SQRT_COLUMNS
        passing_psms[!, column] = zeros(Float32, n)
    end
    for column in MBR_INTEGRATED_TEMPORAL_MEAN_SQRT_COLUMNS
        passing_psms[!, column] = zeros(Float32, n)
    end
    passing_psms[!, MBR_INTEGRATED_TEMPORAL_TRACE_COLUMN] =
        [Float32[] for _ in 1:n]
    passing_psms[!, MBR_INTEGRATED_APEX_IRT_COLUMN] = fill(NaN32, n)
    passing_psms[!, MBR_INTEGRATED_WEIGHT_COLUMN] = fill(NaN32, n)
    passing_psms[
        !,
        MBR_INTEGRATED_LOG2_INTENSITY_EXPLAINED_COLUMN,
    ] = fill(NaN32, n)
    passing_psms[
        !,
        MBR_INTEGRATED_FITTED_MANHATTAN_DISTANCE_COLUMN,
    ] = fill(NaN32, n)
    passing_psms[!, MBR_INTEGRATED_FITTED_HELLINGER_COLUMN] =
        fill(NaN32, n)
    passing_psms[
        !,
        MBR_INTEGRATED_SMOOTHED_2D_SHADOW_HELLINGER_COLUMN,
    ] = fill(NaN32, n)
    passing_psms[!, MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_COLUMN] =
        zeros(UInt8, n)
    passing_psms[!, MBR_INTEGRATED_FRAG_CORR_BITVEC_COLUMN] =
        zeros(UInt8, n)
    passing_psms[
        !,
        MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_BITVEC_RANK_COLUMN,
    ] = zeros(UInt16, n)
    passing_psms[!, MBR_INTEGRATED_FRAG_CORR_STRENGTH_COLUMN] =
        zeros(Float32, n)
    passing_psms[!, MBR_INTEGRATED_FRAG_CORR_EFFECTIVE_N_COLUMN] =
        zeros(Float32, n)
    passing_psms[!, MBR_INTEGRATED_FRAG_CORR_BEST_WEIGHT_COLUMN] =
        zeros(Float32, n)
    (n == 0 || nrow(chromatograms) == 0) && return passing_psms

    hasproperty(chromatograms, :shadow_frag1_int) ||
        error("Integrated MBR requires chromatogram shadow fragments")
    scan_tic === nothing &&
        error("Integrated MBR requires the per-scan spectrum intensity vector")
    hasproperty(passing_psms, :integration_start_scan) ||
        error("Integrated MBR requires integration start scans")
    hasproperty(passing_psms, :integration_stop_scan) ||
        error("Integrated MBR requires integration stop scans")

    shadow_columns = ntuple(
        rank -> chromatograms[!, Symbol("shadow_frag$(rank)_int")],
        8,
    )
    fitted_columns = all(
        column -> hasproperty(chromatograms, column),
        FITTED_FRAGMENT_INTENSITY_COLUMNS,
    ) ? ntuple(
        rank -> chromatograms[
            !,
            FITTED_FRAGMENT_INTENSITY_COLUMNS[rank],
        ],
        8,
    ) : nothing
    rt_column = chromatograms.rt
    precursor_column = chromatograms.precursor_idx
    scan_column = chromatograms.scan_idx
    weight_column = chromatograms.intensity

    # Computed here rather than in the kernel: both touch the DataFrames directly.
    use_separate_traces =
        hasproperty(chromatograms, :isotopes_captured) &&
        hasproperty(passing_psms, :isotopes_captured)
    chromatogram_isotopes = use_separate_traces ?
        chromatograms.isotopes_captured :
        nothing

    # FUNCTION BARRIER -- see _add_mbr_integrated_kernel! for why this is not just tidiness.
    _add_mbr_integrated_kernel!(
        passing_psms,
        n,
        nrow(chromatograms),
        rt_column,
        precursor_column,
        scan_column,
        weight_column,
        shadow_columns,
        fitted_columns,
        chromatogram_isotopes,
        use_separate_traces,
        rt_to_irt,
        bitvec_rank_table,
        scan_tic,
        passing_psms.precursor_idx,
        passing_psms.new_best_scan,
        passing_psms.scan_idx,
        passing_psms.integration_start_scan,
        passing_psms.integration_stop_scan,
        passing_psms.peak_area,
        use_separate_traces ? passing_psms.isotopes_captured : nothing,
    )
    return passing_psms
end


# Function barrier for the hot part of add_mbr_integrated_spectra_to_psms!.
#
# WHY THIS EXISTS: `df.col` infers as `AbstractVector` and `ntuple(r -> df[!, c], 8)` as
# `NTuple{8, AbstractVector}`. Used in a loop *inline in the same function*, every element access is
# then a dynamic dispatch that boxes its result -- measured at 673 ms and 240 MB for a 5 M-element
# sum, versus 0.5 ms and 0 B when the same column is reached through a function argument (1300x).
# The previously reported 2.28 GB / 10,146 ms for two linear scans over 8.7 M rows was this.
#
# Passing the columns as ARGUMENTS is what fixes it: Julia specialises on the runtime type, so no
# type assertion is needed and nothing breaks if a column's concrete type ever changes. Tuples are
# covariant, so `shadow_columns` arrives here as a concrete `NTuple{8, Vector{...}}` even though the
# caller could only infer `NTuple{8, AbstractVector}`.
function _add_mbr_integrated_kernel!(
    passing_psms::DataFrame,
    n::Int,
    n_chrom::Int,
    rt_column,
    precursor_column,
    scan_column,
    weight_column,
    shadow_columns,
    fitted_columns,
    chromatogram_isotopes,
    use_separate_traces::Bool,
    rt_to_irt,
    bitvec_rank_table,
    scan_tic,
    psm_pid,
    psm_new_best_scan,
    psm_scan_idx,
    psm_start_scan,
    psm_stop_scan,
    psm_peak_area,
    psm_isotopes,
)
    # DIAGNOSTIC (PIONEER_MBR_PHASE_DIAG=1): per-phase bytes/time inside this function, accumulated
    # across files. Chromatogram Integration allocates 77.7 GB on this branch vs 4.3 GB without MBR;
    # this attributes that to the four index-building phases vs the per-row write loop.
    _diag = get(ENV, "PIONEER_MBR_PHASE_DIAG", "0") == "1"
    _t0 = time(); _a0 = Base.gc_bytes()

    # Chromatograms are sorted [:precursor_idx, :rt] (SeperateTraces: [:precursor_idx,
    # :isotopes_captured, :rt]) by sort_chromatograms_for_integration! at
    # IntegrateChromatogramsSearch.jl:354, BEFORE this runs at :386. So every precursor's rows are
    # already contiguous and time-ordered, which makes the old index structures redundant:
    #   rows_by_pid / rows_by_trace : Dict of one heap Vector{Int} per precursor -> a group table
    #                                 built by count-then-fill (no push!)
    #   per-pid sort!(rows; by=...) : data is already in scan order
    #   neighbors Dict              : held ONE tuple-keyed entry per chromatogram row (measured
    #                                 8,717,899 = exactly the row count) purely to answer
    #                                 "(pid, selected_scan) -> row"; that is a searchsortedfirst
    #                                 inside the precursor's range, and the neighbours are the
    #                                 adjacent indices.
    #   correlation_by_pid          : Dict{UInt32, NamedTuple} -- bare NamedTuple is ABSTRACT, so
    #                                 every value was boxed. Now a concretely-typed Vector indexed
    #                                 by group.

    # pass 1: count precursor groups (no allocation)
    n_groups = 0
    @inbounds let r = 1
        while r <= n_chrom
            pid = UInt32(precursor_column[r])
            n_groups += 1
            while r <= n_chrom && UInt32(precursor_column[r]) == pid
                r += 1
            end
        end
    end
    # pass 2: fill exact-size group table
    group_pid = Vector{UInt32}(undef, n_groups)
    group_lo  = Vector{Int32}(undef, n_groups)
    group_hi  = Vector{Int32}(undef, n_groups)
    @inbounds let r = 1, g = 0
        while r <= n_chrom
            pid = UInt32(precursor_column[r]); s0 = r
            while r <= n_chrom && UInt32(precursor_column[r]) == pid
                r += 1
            end
            g += 1
            group_pid[g] = pid; group_lo[g] = Int32(s0); group_hi[g] = Int32(r - 1)
        end
    end

    if _diag
        MBR_PHASE_DIAG[:rows_by_pid_bytes] += Base.gc_bytes() - _a0
        MBR_PHASE_DIAG[:rows_by_pid_ms] += round(Int, (time() - _t0) * 1000)
        MBR_PHASE_DIAG[:n_chrom_rows] += n_chrom
        MBR_PHASE_DIAG[:n_psm_rows] += n
        MBR_PHASE_DIAG[:n_pids] += n_groups
        _t0 = time(); _a0 = Base.gc_bytes()
    end

    # SeperateTraces stores rows grouped by isotope set, so scan_idx is not monotonic across a
    # precursor. Reproduce the old per-precursor sort with ONE Int32 permutation (4 B/row) instead
    # of a Dict holding a tuple entry per chromatogram row.
    scan_order = if use_separate_traces
        perm = collect(Int32(1):Int32(n_chrom))
        @inbounds for g in 1:n_groups
            lo, hi = Int(group_lo[g]), Int(group_hi[g])
            hi > lo && sort!(view(perm, lo:hi); by = i -> UInt32(scan_column[i]))
        end
        perm
    else
        nothing
    end

    if _diag
        # Read the phase timers FIRST. The monotonicity check below is an O(n_chrom) diagnostic; when
        # it ran before this read it was timed as part of the phase and inflated the reported cost of
        # this block by ~8.4 s -- a diagnostic measuring itself.
        MBR_PHASE_DIAG[:neighbors_bytes] += Base.gc_bytes() - _a0
        MBR_PHASE_DIAG[:neighbors_ms] += round(Int, (time() - _t0) * 1000)
        # Assumption check for the combined-trace fast path: the stored order must be scan order.
        if scan_order === nothing
            bad = 0
            @inbounds for g in 1:n_groups, r in (Int(group_lo[g]) + 1):Int(group_hi[g])
                UInt32(scan_column[r]) <= UInt32(scan_column[r - 1]) && (bad += 1)
            end
            MBR_PHASE_DIAG[:nonmonotonic_rows] += bad
        end
        _t0 = time(); _a0 = Base.gc_bytes()
    end

    # Correlation features are a function of the precursor's whole chromatogram group, so compute
    # one per group into a concretely-typed vector (was Dict{UInt32, NamedTuple}, abstract).
    correlation_by_group = Vector{_ChromCorrFeatures}(undef, n_groups)
    # One set of scratch buffers for every group, pre-sized to the largest group so the per-group
    # `resize!`s never reallocate. Correlation always consumes a group's full row range, so the
    # widest group extent is an exact bound.
    max_group_points = 0
    @inbounds for g in 1:n_groups
        span = Int(group_hi[g]) - Int(group_lo[g]) + 1
        span > max_group_points && (max_group_points = span)
    end
    corr_workspace = _ChromCorrWorkspace(max_group_points)
    @inbounds for g in 1:n_groups
        # The old `rows_by_pid[pid]` vectors were sorted IN PLACE by scan_idx before this ran, so
        # correlation consumed scan-ordered rows. Feed the same order here.
        lo, hi = Int(group_lo[g]), Int(group_hi[g])
        c = _chrom_shadow_fragment_correlation_features(
            scan_order === nothing ? (lo:hi) : view(scan_order, lo:hi),
            shadow_columns,
            weight_column,
            bitvec_rank_table,
            corr_workspace,
        )
        correlation_by_group[g] = _ChromCorrFeatures(
            c.n_correlated, c.corr_mask, c.bitvec_rank,
            c.strength, c.effective_n, c.best_weight,
        )
    end

    if _diag
        MBR_PHASE_DIAG[:correlation_bytes] += Base.gc_bytes() - _a0
        MBR_PHASE_DIAG[:correlation_ms] += round(Int, (time() - _t0) * 1000)
        _t0 = time(); _a0 = Base.gc_bytes()
    end

    # Bind every column touched in the loop ONCE. `df[row, ::Symbol]` and `df.col` each do a
    # Dict{Symbol,Int} lookup and return an abstractly-typed column, so the ~34 accesses per row
    # were dynamic. These bindings are concrete, so the loop body specialises.
    out_frag_sqrt = ntuple(
        r -> getproperty(passing_psms, MBR_INTEGRATED_FRAGMENT_SQRT_COLUMNS[r])::Vector{Float32},
        Val(8),
    )
    out_temporal_mean = ntuple(
        r -> getproperty(passing_psms, MBR_INTEGRATED_TEMPORAL_MEAN_SQRT_COLUMNS[r])::Vector{Float32},
        Val(8),
    )
    out_trace = getproperty(passing_psms, MBR_INTEGRATED_TEMPORAL_TRACE_COLUMN)::Vector{Vector{Float32}}
    out_apex_irt = getproperty(passing_psms, MBR_INTEGRATED_APEX_IRT_COLUMN)::Vector{Float32}
    out_weight = getproperty(passing_psms, MBR_INTEGRATED_WEIGHT_COLUMN)::Vector{Float32}
    out_log2_explained =
        getproperty(passing_psms, MBR_INTEGRATED_LOG2_INTENSITY_EXPLAINED_COLUMN)::Vector{Float32}
    out_fitted_manhattan =
        getproperty(passing_psms, MBR_INTEGRATED_FITTED_MANHATTAN_DISTANCE_COLUMN)::Vector{Float32}
    out_fitted_hellinger =
        getproperty(passing_psms, MBR_INTEGRATED_FITTED_HELLINGER_COLUMN)::Vector{Float32}
    out_smoothed_2d =
        getproperty(passing_psms, MBR_INTEGRATED_SMOOTHED_2D_SHADOW_HELLINGER_COLUMN)::Vector{Float32}
    out_n_corr = getproperty(passing_psms, MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_COLUMN)::Vector{UInt8}
    out_corr_bitvec = getproperty(passing_psms, MBR_INTEGRATED_FRAG_CORR_BITVEC_COLUMN)::Vector{UInt8}
    out_bitvec_rank = getproperty(
        passing_psms, MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_BITVEC_RANK_COLUMN,
    )::Vector{UInt16}
    out_corr_strength =
        getproperty(passing_psms, MBR_INTEGRATED_FRAG_CORR_STRENGTH_COLUMN)::Vector{Float32}
    out_corr_effective_n =
        getproperty(passing_psms, MBR_INTEGRATED_FRAG_CORR_EFFECTIVE_N_COLUMN)::Vector{Float32}
    out_corr_best_weight =
        getproperty(passing_psms, MBR_INTEGRATED_FRAG_CORR_BEST_WEIGHT_COLUMN)::Vector{Float32}


    @inbounds for row in 1:n
        pid = UInt32(psm_pid[row])
        selected_scan = UInt32(psm_new_best_scan[row])
        selected_scan == UInt32(0) &&
            (selected_scan = UInt32(psm_scan_idx[row]))
        gidx = searchsortedfirst(group_pid, pid)
        has_group = gidx <= n_groups && group_pid[gidx] == pid
        neighbor_rows = has_group ?
            _chrom_neighbor_rows(
                scan_column, scan_order, Int(group_lo[gidx]), Int(group_hi[gidx]), selected_scan,
            ) : (0, 0, 0)
        fragment_sqrt, shadow_sum = neighbor_rows[1] == 0 ?
            (MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT, 0.0f0) :
            _chrom_smoothed_shadow_sqrt_tuple_and_sum(
                shadow_columns,
                neighbor_rows[1],
                neighbor_rows[2],
                neighbor_rows[3],
            )
        for rank in 1:8
            out_frag_sqrt[rank][row] = fragment_sqrt[rank]
        end

        temporal_rows = if !has_group
            1:0                                   # empty range, no allocation
        elseif use_separate_traces
            # isotopes_captured is the MIDDLE sort key, so its groups are contiguous inside the
            # precursor's range; narrow with a second search rather than a second Dict.
            isotope_set = Tuple{Int8, Int8}(psm_isotopes[row])
            _chrom_isotope_subrange(
                chromatogram_isotopes, Int(group_lo[gidx]), Int(group_hi[gidx]), isotope_set,
            )
        else
            Int(group_lo[gidx]):Int(group_hi[gidx])
        end
        temporal_trace, temporal_mean =
            _chrom_temporal_fragment_data(
                temporal_rows,
                shadow_columns,
                weight_column,
                scan_column,
                UInt32(psm_start_scan[row]),
                UInt32(psm_stop_scan[row]),
            )
        out_trace[row] = temporal_trace
        for rank in 1:8
            out_temporal_mean[rank][row] = temporal_mean[rank]
        end

        if has_group
            correlation = correlation_by_group[gidx]
            out_n_corr[row] = correlation.n_correlated
            out_corr_bitvec[row] = correlation.corr_mask
            out_bitvec_rank[row] = correlation.bitvec_rank
            out_corr_strength[row] = correlation.strength
            out_corr_effective_n[row] = correlation.effective_n
            out_corr_best_weight[row] = correlation.best_weight
        end

        apex_row = neighbor_rows[1]
        apex_row == 0 && continue
        out_apex_irt[row] = Float32(rt_to_irt(Float32(rt_column[apex_row])))
        out_weight[row] = Float32(psm_peak_area[row])
        # _chrom_neighbor_rows only returns a row whose scan == selected_scan, so the apex row's
        # scan IS selected_scan; index the per-scan vector directly.
        out_log2_explained[row] = log2(
            max(shadow_sum, 1.0f-20) /
            max(@inbounds(scan_tic[selected_scan]), 1.0f-20),
        )
        if fitted_columns !== nothing
            fitted_sqrt, _ =
                _chrom_spectrum_sqrt_tuple_and_sum(
                    fitted_columns,
                    apex_row,
                )
            apex_shadow_sqrt, _ =
                _chrom_spectrum_sqrt_tuple_and_sum(
                    shadow_columns,
                    apex_row,
                )
            out_fitted_manhattan[row] = _chrom_manhattan_score(
                shadow_columns,
                fitted_columns,
                apex_row,
            )
            out_fitted_hellinger[row] =
                _chrom_hellinger_score_from_sqrt(apex_shadow_sqrt, fitted_sqrt)
            out_smoothed_2d[row] =
                _chrom_hellinger_score_from_sqrt(fragment_sqrt, fitted_sqrt)
        end
    end
    if _diag
        MBR_PHASE_DIAG[:perrow_bytes] += Base.gc_bytes() - _a0
        MBR_PHASE_DIAG[:perrow_ms] += round(Int, (time() - _t0) * 1000)
        MBR_PHASE_DIAG[:n_files] += 1
        _mbr_phase_diag_report()
    end
    return passing_psms
end

# Accumulators for the phase diagnostic above. Printed after every file so a crash mid-run still
# leaves usable numbers.
const MBR_PHASE_DIAG = Dict{Symbol, Int}(
    :rows_by_pid_bytes => 0, :rows_by_pid_ms => 0,
    :neighbors_bytes => 0,   :neighbors_ms => 0,
    :correlation_bytes => 0, :correlation_ms => 0,
    :perrow_bytes => 0,      :perrow_ms => 0,
    :n_chrom_rows => 0, :n_psm_rows => 0, :n_pids => 0, :n_neighbors => 0, :n_files => 0,
    :nonmonotonic_rows => 0,
)

function _mbr_phase_diag_report()
    d = MBR_PHASE_DIAG
    gb(x) = round(x / 2^30, digits = 2)
    tot = d[:rows_by_pid_bytes] + d[:neighbors_bytes] + d[:correlation_bytes] + d[:perrow_bytes]
    pct(x) = tot > 0 ? round(100 * x / tot, digits = 1) : 0.0
    @user_info """
    MBR phase diagnostic (cumulative over $(d[:n_files]) file(s)):
      chrom rows=$(d[:n_chrom_rows])  psm rows=$(d[:n_psm_rows])  pids=$(d[:n_pids])  neighbors entries=$(d[:n_neighbors])
      combined-mode non-scan-ordered rows (MUST be 0): $(d[:nonmonotonic_rows])
      rows_by_pid  + rows_by_trace : $(gb(d[:rows_by_pid_bytes])) GB  $(d[:rows_by_pid_ms]) ms  ($(pct(d[:rows_by_pid_bytes]))%)
      neighbors Dict + sort!       : $(gb(d[:neighbors_bytes])) GB  $(d[:neighbors_ms]) ms  ($(pct(d[:neighbors_bytes]))%)
      correlation_by_pid           : $(gb(d[:correlation_bytes])) GB  $(d[:correlation_ms]) ms  ($(pct(d[:correlation_bytes]))%)
      per-row write loop           : $(gb(d[:perrow_bytes])) GB  $(d[:perrow_ms]) ms  ($(pct(d[:perrow_bytes]))%)
      TOTAL in this function       : $(gb(tot)) GB"""
    return nothing
end

#==========================================================
Chromatogram Building Functions
==========================================================#
"""
    extract_chromatograms(spectra::MassSpecData, passing_psms::DataFrame,
                         rt_index::retentionTimeIndex, search_context::SearchContext,
                         params::IntegrateChromatogramSearchParameters,
                         ms_file_idx::Int64) -> DataFrame

Extract chromatograms for all passing PSMs.

Uses parallel processing to build chromatograms across scan ranges.
"""

# Quad-window helpers (originally in selectTransitions/rtIndexTransitionSelection.jl;
# kept here as the sole remaining caller after the fused migration).

function getQuadrupoleBounds(quad_func::QuadTransmissionFunction)
    return getPrecMinBound(quad_func), getPrecMaxBound(quad_func)
end

function getPrecursorWindowRange(
    precs, min_prec_mz::Float32, max_prec_mz::Float32,
    isotope_err_bounds::Tuple{I, I}
) where {I<:Integer}
    start = searchsortedfirst(precs, by = x->last(x),
        min_prec_mz - first(isotope_err_bounds)*C13_C12_MASS_DIFF/2)
    stop = searchsortedlast(precs, by = x->last(x),
        max_prec_mz + last(isotope_err_bounds)*C13_C12_MASS_DIFF/2)
    return start, stop
end

function withinQuadrupoleBounds(
    prec_mz::Float32, prec_charge::UInt8,
    min_prec_mz::Float32, max_prec_mz::Float32,
    isotope_err_bounds::Tuple{I, I}
) where {I<:Integer}
    mz_low = min_prec_mz - first(isotope_err_bounds)*C13_C12_MASS_DIFF/prec_charge
    mz_high = max_prec_mz + last(isotope_err_bounds)*C13_C12_MASS_DIFF/prec_charge
    return mz_low ≤ prec_mz ≤ mz_high
end

"""
    collect_rt_window_precursors!(precs_temp, rt_index, rt_start_idx, rt_stop_idx,
                                   precursors_passing, prec_mzs, prec_charges,
                                   prec_sulfur_counts, iso_splines, quad_func,
                                   precursor_transmission, isotope_err_bounds,
                                   min_fraction_transmitted, precursor_rt_map,
                                   scan_rt, rt_binned_tol, rt_tol_fallback) -> Int

Walk the RT-bin range, applying quad-window, precursors_passing allowlist,
per-precursor RT, isotope_err_bounds, and min_fraction_transmitted filters.
Writes the surviving precursor ids into `precs_temp[1:n]` and returns `n`.

Extracted from the classic `RTIndexedTransitionSelection` path so
`run_fused!(FusedRTIndexed(), ...)` can consume a pre-filtered precursor
list without repeating the filtering (skipped via `check_prec_filters`).
"""
function collect_rt_window_precursors!(
    precs_temp::Vector{UInt32},
    rt_index::retentionTimeIndex,
    rt_start_idx::Int64, rt_stop_idx::Int64,
    precursors_passing::Union{Set{UInt32}, Nothing},
    prec_mzs::AbstractArray{Float32},
    prec_charges::AbstractArray{UInt8},
    prec_sulfur_counts::AbstractArray{UInt8},
    iso_splines::IsotopeSplineModel,
    quad_transmission_func::QuadTransmissionFunction,
    precursor_transmission::Vector{Float32},
    isotope_err_bounds::Tuple{I, I},
    min_fraction_transmitted::Float32,
    precursor_rt_map::Union{Dict{UInt32, Float32}, Nothing},
    scan_rt::Float32,
    rt_binned_tol::Union{RTBinnedTolerance, Nothing},
    rt_tol_fallback::Float32) where {I<:Integer}

    min_prec_mz, max_prec_mz = getQuadrupoleBounds(quad_transmission_func)
    size = 0
    has_rt_filter = precursor_rt_map !== nothing

    for rt_bin_idx in rt_start_idx:rt_stop_idx
        precs = rt_index.rt_bins[rt_bin_idx].prec
        start, stop = getPrecursorWindowRange(precs, min_prec_mz, max_prec_mz, isotope_err_bounds)

        for i in start:stop
            prec_idx = first(precs[i])
            (!isnothing(precursors_passing) && prec_idx ∉ precursors_passing) && continue

            if has_rt_filter
                prec_rt_val = get(precursor_rt_map, prec_idx, NaN32)
                if !isnan(prec_rt_val)
                    prec_tol = rt_binned_tol !== nothing ? get_rt_tol(rt_binned_tol, prec_rt_val) : rt_tol_fallback
                    abs(scan_rt - prec_rt_val) > prec_tol && continue
                end
            end

            prec_charge = prec_charges[prec_idx]
            prec_mz = prec_mzs[prec_idx]
            prec_sulfur_count = prec_sulfur_counts[prec_idx]

            !withinQuadrupoleBounds(prec_mz, prec_charge, min_prec_mz, max_prec_mz, isotope_err_bounds) && continue

            fraction_transmitted = getPrecursorFractionTransmitted!(
                precursor_transmission, iso_splines, (1, 5),
                quad_transmission_func, prec_mz, prec_charge, prec_sulfur_count)
            fraction_transmitted < min_fraction_transmitted && continue

            prec_isotope_set = getPrecursorIsotopeSet(prec_mz, prec_charge, quad_transmission_func)
            last(prec_isotope_set) < 0 && continue

            size += 1
            if size > length(precs_temp)
                append!(precs_temp, Vector{UInt32}(undef, length(precs_temp)))
            end
            precs_temp[size] = prec_idx
        end
    end
    return size
end

function extract_chromatograms(
    spectra::MassSpecData,
    passing_psms::DataFrame,
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64,
    chrom_type::CHROMATOGRAM
)
    if typeof(chrom_type)==typeof(MS2CHROM())
        ms_order_select = 2
    else
        ms_order_select = 1
    end

    thread_tasks = partition_scans(spectra, Threads.nthreads(), ms_order_select = ms_order_select)

    # Pre-create Sets for each thread to avoid concurrent DataFrame access
    # This eliminates race conditions when multiple threads call passing_psms[!, :precursor_idx]
    precursor_set = Set(passing_psms[!, :precursor_idx])  # shared read-only across threads (was N copies)

    # Build precursor RT map for per-precursor symmetric window filtering
    _pids = passing_psms[!, :precursor_idx]::Vector{UInt32}
    _rts = passing_psms[!, :rt]::Vector{Float32}
    precursor_rt_map = Dict{UInt32, Float32}()
    sizehint!(precursor_rt_map, length(_pids))
    for i in eachindex(_pids)
        precursor_rt_map[_pids[i]] = _rts[i]
    end

    # One entry per scan, shared across threads. partition_scans gives each thread a DISJOINT
    # set of scan indices, so the writes never collide.
    scan_tic = zeros(Float32, length(spectra))

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]

            return build_chromatograms(
                spectra,
                last(thread_task),
                precursor_set,
                precursor_rt_map,
                rt_index,
                search_context,
                search_data,
                params,
                ms_file_idx,
                chrom_type;
                scan_tic = scan_tic,
            )
        end
    end
    thread_dfs = fetch.(tasks)
    return vcat(thread_dfs...), scan_tic
end

"""
    build_chromatograms(spectra::MassSpecData, scan_range::Vector{Int64},
                       precursors_passing::Set{UInt32}, rt_index::retentionTimeIndex,
                       search_context::SearchContext, search_data::SearchDataStructures,
                       params::IntegrateChromatogramSearchParameters,
                       ms_file_idx::Int64) -> DataFrame

Build chromatograms for a range of scans with RT bin caching.

# Process
1. Tracks RT bins for efficient transition selection
2. Selects transitions based on RT windows
3. Matches peaks and performs deconvolution
4. Records chromatogram points with weights
"""
function build_chromatograms(
    spectra::MassSpecData,
    scan_range::Vector{Int64},
    precursors_passing::Set{UInt32},
    precursor_rt_map::Dict{UInt32, Float32},
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64,
    ::MS2CHROM;
    scan_tic::Union{Nothing, Vector{Float32}} = nothing,
)
    # Fused-kernel working arrays.
    Hs = getHsFused(search_data)
    weights = getTempWeights(search_data)
    colnorm2 = getColNorm2(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    spectral_scores = getMainSearchSpectralScores(search_data)
    collect_mbr_evidence = params.match_between_runs
    fused_scratch = getFusedScratch(search_data)
    corr_mz = getScanCorrectedMz(search_data)
    obs_low = getScanObsLow(search_data)
    obs_high = getScanObsHigh(search_data)
    isotopes_buf = getIsotopes(search_data)
    prec_trans_buf = getPrecursorTransmission(search_data)
    id_to_col = getIdToCol(search_data)
    unscored_psms = getTuningUnscoredPsms(search_data)  # placeholder — FusedRTIndexed record_match! is no-op

    # Points per scan measured at 23.81 on 6-file Olsen Exploris (366,129 MS2 scans ->
    # 8,717,899 chromatogram rows). The old estimate of 100 over-allocated 4.2x, which cost
    # ~2.1 GB of cumulative allocation per run at 80 B/row.
    #
    # 32 leaves ~34 % headroom so the common case never resizes. Note that a resize is NOT free
    # in cumulative terms: with initial E and k doublings the total allocated is E*(2^(k+1)-1),
    # so E=20 + one doubling costs 60x scans while E=32 with no doubling costs 32x. The doubling
    # below still covers datasets that run richer than this one (library size, RT tolerance and
    # isolation width all move this number, and only one dataset has been measured).
    estimated_points = length(scan_range) * 32
    initial_size = max(estimated_points, 10000)
    chromatograms = collect_mbr_evidence ?
        Vector{MS2MBRChromObject}(undef, initial_size) :
        Vector{MS2ChromObject}(undef, initial_size)

    # RT bin tracking state — RT index is in empirical RT space
    rt_bin_start, rt_bin_stop = 1, 1
    prec_mz = zero(Float32)
    rt_idx = 0
    precs_temp = getPrecIds(search_data)
    prec_temp_size = 0
    irt_tol = getIrtErrors(search_context)[ms_file_idx]
    has_rt_tol = haskey(getRtTolerances(search_context), ms_file_idx)
    rt_binned_tol = has_rt_tol ? getRtTolerance(search_context, ms_file_idx) : nothing
    rt_irt_model = getRtIrtModel(search_context, ms_file_idx)
    nce_model = getNceModel(search_context, ms_file_idx)
    mass_error_model = getMassErrorModel(search_context, ms_file_idx)
    quad_model = getQuadTransmissionModel(search_context, ms_file_idx)
    spec_lib = getSpecLib(search_context)
    precursors = getPrecursors(spec_lib)
    prec_mz_arr = getMz(precursors)
    prec_charge_arr = getCharge(precursors)
    prec_sulfur_arr = getSulfurCount(precursors)
    prec_irt_arr = getIrt(precursors)
    frag_lookup = getFragmentLookupTable(spec_lib)

    kind = FusedRTIndexed(params.prec_estimation, UInt8(params.max_frag_rank))

    for scan_idx in scan_range
        ((scan_idx<1) | (scan_idx > length(spectra))) && continue
        msn = getMsOrder(spectra, scan_idx)
        msn ∉ params.spec_order && continue

        rt = getRetentionTime(spectra, scan_idx)
        # The total ion current is a property of the SCAN, not of any one precursor. Record it
        # once per scan; it used to be copied onto every chromatogram row of the scan.
        #
        # Read the instrument's stored TIC instead of re-summing the peak list. Re-summing was
        # the only full traversal of the intensity column in this loop -- fragment matching
        # touches matched peaks only -- so its cost scales with peaks-per-scan, which is ~10x
        # higher on Astral than on Exploris.
        #
        # The two are not bit-identical, but they agree closely and one-sidedly: over 60,969
        # MS2 scans, sum/TIC has median 0.99814 (mean 0.99719, sd 0.00344), with TIC >= sum in
        # 86% of scans and exact equality in many. The peak list is a near-complete subset of
        # what the TIC counted. Since the value is consumed as log2(shadow_sum / TIC), that
        # shifts the feature by a median 0.0027 log2 units against a range of ~10.
        if collect_mbr_evidence && scan_tic !== nothing
            @inbounds scan_tic[scan_idx] = Float32(getTIC(spectra, scan_idx))
        end

        if rt_binned_tol !== nothing
            rt_tol_local = get_rt_tol(rt_binned_tol, Float32(rt))
        else
            h = 0.1f0
            local_slope = abs((rt_irt_model(rt + h) - rt_irt_model(rt - h)) / (2f0 * h))
            rt_tol_local = irt_tol / max(local_slope, 0.01f0)
        end

        rt_bin_start_new = max(searchsortedfirst(rt_index.rt_bins, rt - rt_tol_local, lt=(r,x)->r.lb<x) - 1, 1)
        rt_bin_stop_new = min(searchsortedlast(rt_index.rt_bins, rt + rt_tol_local, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

        prec_mz_new = getCenterMz(spectra, scan_idx)
        if (rt_bin_start_new != rt_bin_start) || (rt_bin_stop_new != rt_bin_stop) || (prec_mz_new != prec_mz)
            rt_bin_start = rt_bin_start_new
            rt_bin_stop = rt_bin_stop_new
            prec_mz = prec_mz_new
        end

        quad_func = getQuadTransmissionFunction(
            quad_model,
            getCenterMz(spectra, scan_idx),
            getIsolationWidthMz(spectra, scan_idx)
        )

        # 1. Enumerate precursors in the RT window (+ quad + per-prec RT +
        #    min_fraction_transmitted + allowlist filters). Writes to
        #    precs_temp[1:prec_temp_size].
        prec_temp_size = collect_rt_window_precursors!(
            precs_temp, rt_index, rt_bin_start, rt_bin_stop,
            precursors_passing, prec_mz_arr, prec_charge_arr, prec_sulfur_arr,
            getIsoSplines(search_data), quad_func, prec_trans_buf,
            params.isotope_err_bounds, params.min_fraction_transmitted,
            precursor_rt_map, Float32(rt), rt_binned_tol, rt_tol_local)

        if prec_temp_size == 0
            reset!(id_to_col); reset!(Hs)
            continue
        end

        # 2. Pre-correct peak m/z for the scan.
        scan_mz  = getMzArray(spectra, scan_idx)
        scan_int = getIntensityArray(spectra, scan_idx)
        peak_mz_len = prepare_scan_peaks!(corr_mz, obs_low, obs_high,
                                          mass_error_model, scan_mz, scan_int,
                                          Float32(rt))

        # 3. Fused match+build. iRT + iso_err_bounds filters already done
        #    upstream — skipped by FusedRTIndexed's check_prec_filters = false.
        nmatches, nmisses = run_fused!(
            kind,
            Hs, unscored_psms, id_to_col, fused_scratch,
            corr_mz, obs_low, obs_high, peak_mz_len,
            isotopes_buf, prec_trans_buf,
            frag_lookup, nce_model,
            precs_temp, 1:prec_temp_size,
            prec_mz_arr, prec_charge_arr, prec_sulfur_arr, prec_irt_arr,
            getIsoSplines(search_data), quad_func, mass_error_model,
            scan_int, 0f0, Float32(Inf),
            (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx)),
            params.n_frag_isotopes,
            params.isotope_err_bounds)

        if nmatches > 2
            # Resize weight buffers for new columns.
            n_active_cols = n_active(id_to_col)
            if n_active_cols > length(weights)
                new_entries = n_active_cols - length(weights) + 1000
                resize!(weights, length(weights) + new_entries)
                resize!(colnorm2, length(colnorm2) + new_entries)
                collect_mbr_evidence &&
                    resize!(
                        spectral_scores,
                        length(spectral_scores) + new_entries,
                    )
            end

            initialize_weights!(id_to_col, weights, precursor_weights)

            solve_deconvolution!(
                calibrated_chromatogram_deconvolution_solver(search_context, params.deconvolution_solver),
                Hs, residuals, weights, colnorm2,
                getMu(search_data), getObserved(search_data),
                params.max_iter_outer, params.max_diff)
            # Only the 16 rank-1..8 fragment intensities are read back from spectral_scores here
            # (MS2MBRChromObject has no other score field), so the five aggregate metrics
            # getDistanceMetrics computes are not calculated.
            collect_mbr_evidence &&
                getFragmentIntensities!(weights, residuals, Hs, spectral_scores)

            # Record chromatogram points (weighted for matches, zero otherwise).
            for j in 1:prec_temp_size
                rt_idx += 1
                if rt_idx + 1 > length(chromatograms)
                    resize!(chromatograms, length(chromatograms) * 2)
                end
                col = id_to_col[precs_temp[j]]
                score = !collect_mbr_evidence || iszero(col) ?
                    nothing :
                    spectral_scores[col]
                if collect_mbr_evidence
                    chromatograms[rt_idx] = MS2MBRChromObject(
                        Float32(getRetentionTime(spectra, scan_idx)),
                        iszero(col) ? zero(Float32) : weights[col],
                        scan_idx,
                        precs_temp[j],
                        score === nothing ? 0.0f0 : Float32(score.shadow_frag1_int),
                        score === nothing ? 0.0f0 : Float32(score.shadow_frag2_int),
                        score === nothing ? 0.0f0 : Float32(score.shadow_frag3_int),
                        score === nothing ? 0.0f0 : Float32(score.shadow_frag4_int),
                        score === nothing ? 0.0f0 : Float32(score.shadow_frag5_int),
                        score === nothing ? 0.0f0 : Float32(score.shadow_frag6_int),
                        score === nothing ? 0.0f0 : Float32(score.shadow_frag7_int),
                        score === nothing ? 0.0f0 : Float32(score.shadow_frag8_int),
                        score === nothing ? 0.0f0 : Float32(score.fitted_frag1_int),
                        score === nothing ? 0.0f0 : Float32(score.fitted_frag2_int),
                        score === nothing ? 0.0f0 : Float32(score.fitted_frag3_int),
                        score === nothing ? 0.0f0 : Float32(score.fitted_frag4_int),
                        score === nothing ? 0.0f0 : Float32(score.fitted_frag5_int),
                        score === nothing ? 0.0f0 : Float32(score.fitted_frag6_int),
                        score === nothing ? 0.0f0 : Float32(score.fitted_frag7_int),
                        score === nothing ? 0.0f0 : Float32(score.fitted_frag8_int),
                    )
                else
                    chromatograms[rt_idx] = MS2ChromObject(
                        Float32(getRetentionTime(spectra, scan_idx)),
                        iszero(col) ? zero(Float32) : weights[col],
                        scan_idx,
                        precs_temp[j],
                    )
                end
            end

            update_precursor_weights!(id_to_col, weights, precursor_weights)
        else
            # Too few matches → zero-weight points for every enumerated precursor.
            for j in 1:prec_temp_size
                rt_idx += 1
                if rt_idx + 1 > length(chromatograms)
                    resize!(chromatograms, length(chromatograms) * 2)
                end
                if collect_mbr_evidence
                    chromatograms[rt_idx] = MS2MBRChromObject(
                        Float32(getRetentionTime(spectra, scan_idx)),
                        zero(Float32),
                        scan_idx,
                        precs_temp[j],
                        0.0f0, 0.0f0, 0.0f0, 0.0f0,
                        0.0f0, 0.0f0, 0.0f0, 0.0f0,
                        0.0f0, 0.0f0, 0.0f0, 0.0f0,
                        0.0f0, 0.0f0, 0.0f0, 0.0f0,
                    )
                else
                    chromatograms[rt_idx] = MS2ChromObject(
                        Float32(getRetentionTime(spectra, scan_idx)),
                        zero(Float32),
                        scan_idx,
                        precs_temp[j],
                    )
                end
            end
        end

        reset_scan_arrays!(id_to_col, Hs, unscored_psms)
    end

    return DataFrame(@view(chromatograms[1:rt_idx]))
end

# MS1 chromatogram extraction is currently unwired (the ms1_quant knob
# was deleted from the public schema). The body is preserved in a block-
# comment so it can be revived — but it uses classic ion_matches /
# buildDesignMatrix! / PrecursorMatch machinery that the fused cleanup
# PR will delete. If MS1 quant is re-enabled in the future, this needs
# a fused port (a FusedMs1 kind + an MS1-specific collector, or a
# parallel lightweight kernel that doesn't touch SparseArrayFused).
#=
"""
    build_chromatograms(spectra::MassSpecData, scan_range::Vector{Int64},
                       precursors_passing::Set{UInt32}, rt_index::retentionTimeIndex,
                       search_context::SearchContext, search_data::SearchDataStructures,
                       params::IntegrateChromatogramSearchParameters,
                       ms_file_idx::Int64) -> DataFrame

Build chromatograms for a range of scans with RT bin caching.

# Process
1. Tracks RT bins for efficient transition selection
2. Selects transitions based on RT windows
3. Matches peaks and performs deconvolution
4. Records chromatogram points with weights
"""
function build_chromatograms(
    spectra::MassSpecData,
    scan_range::Vector{Int64},
    precursors_passing::Set{UInt32},
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64,
    ::MS1CHROM
)
    # Initialize working arrays
    mem = getMs1MassErrorModel(search_context, ms_file_idx)
    Hs = getHs(search_data)
    weights = getTempWeights(search_data)
    colnorm2 = getColNorm2(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    chromatograms = Vector{MS1ChromObject}(undef, 500000)  # Initial size
    ion_templates = Vector{Isotope{Float32}}(undef, 100000)
    ion_matches = [PrecursorMatch{Float32}() for _ in range(1, 10000)]
    ion_misses = [UnmatchedIon() for _ in range(1, 10000)]
    precursors = getPrecursors(getSpecLib(search_context))
    seqs = [getSequence(precursors)[pid] for pid in precursors_passing]
    pids = [pid for pid in precursors_passing]
    pcharge = [getCharge(precursors)[pid] for pid in precursors_passing]
    pmz = [getMz(precursors)[pid] for pid in precursors_passing]
    isotopes_dict = getIsotopes(seqs, pmz, pids, pcharge, QRoots(5), 5)

    # NEW: Create m/z grouping map for MS1
    mz_grouping = MzGroupingMap(UInt32(100000))  # 5 decimal place precision

    # RT bin tracking state — RT index is in empirical RT space
    rt_bin_start, rt_bin_stop = 1, 1
    ion_idx = 0
    rt_idx = 0
    precs_temp = getPrecIds(search_data)  # Use search_data's prec_ids
    prec_temp_size = 0
    irt_tol = getIrtErrors(search_context)[ms_file_idx]
    # RT-adaptive tolerance for MS1
    has_rt_tol_ms1 = haskey(getRtTolerances(search_context), ms_file_idx)
    rt_binned_tol_ms1 = has_rt_tol_ms1 ? getRtTolerance(search_context, ms_file_idx) : nothing
    rt_irt_model = getRtIrtModel(search_context, ms_file_idx)
    id_to_col = getIdToCol(search_data)
    spectral_scores = getSpectralScores(search_data)
    unscored_psms = getUnscoredPsms(search_data)
    i = 1
    for scan_idx in scan_range

        ((scan_idx<1) | (scan_idx > length(spectra))) && continue
        # Process MS1 scans
        #msn = getMsOrder(spectra, scan_idx)
        #msn ∉ one(UInt8) && continue
        if getMsOrder(spectra, scan_idx) != 1
            continue
        end
        iso_count = Dictionary{UInt32, @NamedTuple{matched_mono::Bool, iso_count::UInt8}}()
        # Calculate RT window — RT index is in empirical RT space
        rt = getRetentionTime(spectra, scan_idx)

        if rt_binned_tol_ms1 !== nothing
            # RT-adaptive: look up RT tolerance directly in RT space
            rt_tol_local = get_rt_tol(rt_binned_tol_ms1, Float32(rt))
        else
            # Fallback: convert global iRT tolerance to RT space via local slope
            h = 0.1f0
            local_slope = abs((rt_irt_model(rt + h) - rt_irt_model(rt - h)) / (2f0 * h))
            rt_tol_local = irt_tol / max(local_slope, 0.01f0)
        end

        rt_bin_start = max(searchsortedfirst(rt_index.rt_bins, rt - rt_tol_local, lt=(r,x)->r.lb<x) - 1, 1)
        rt_bin_stop = min(searchsortedlast(rt_index.rt_bins, rt + rt_tol_local, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

        # Update transitions if window changed
        prec_temp_size = 0
        ion_idx = 0
        for rt_bin_idx in rt_bin_start:rt_bin_stop
            precs = rt_index.rt_bins[rt_bin_idx].prec
            for i in 1:length(precs)
                prec_idx = first(precs[i])
                if prec_idx in precursors_passing
                    prec_temp_size += 1
                    if prec_temp_size > length(precs_temp)
                        append!(precs_temp, Vector{UInt32}(undef, 1000))
                    end
                    precs_temp[prec_temp_size] = prec_idx
                    for iso in isotopes_dict[prec_idx]
                        ion_idx += 1
                        ion_templates[ion_idx] = iso
                    end
                end
            end
        end
        #Probably more efficient way to do this later 
        for i in range(1, ion_idx)
            _ion_ = ion_templates[i]
            pid = getPrecID(_ion_)
            if haskey(iso_count, pid)
                matched_mono = false
                if iso_count[pid].matched_mono
                    matched_mono = true
                elseif getIsoIdx(_ion_)==one(UInt8)
                    matched_mono = true
                end
                iso_count[pid] = (matched_mono = matched_mono, iso_count = iso_count[pid].iso_count + one(UInt8))
            else
                insert!(iso_count, 
                pid,
                (matched_mono = getIsoIdx(_ion_)==one(UInt8), iso_count = one(UInt8))
            )
            end
            iso_count
        end
        sort!(@view(ion_templates[1:ion_idx]), by = x->(getMZ(x)), alg=PartialQuickSort(1:ion_idx))
        # Match peaks
        nmatches, nmisses = matchPeaks!(
            ion_matches,
            ion_misses,
            ion_templates,
            ion_idx,
            getMzArray(spectra, scan_idx),
            getIntensityArray(spectra, scan_idx),
            mem,
            getHighMz(spectra, scan_idx)
        )

        #nmisses -= 1
        sort!(@view(ion_matches[1:nmatches]), alg=QuickSort, lt=ion_match_lt)
        #println("nmatches $nmatches nmisses $nmisses")
        # Process matches
        if nmatches > 2
            i += 1

            # Reset grouping for this scan
            reset!(mz_grouping)

            # Use MS1-specific design matrix construction with m/z grouping
            buildDesignMatrixMS1!(
                Hs,
                ion_matches,
                ion_misses,
                nmatches,
                nmisses,
                mz_grouping,
                precursors
            )

            # Handle array resizing
            if id_to_col.size > length(weights)
                new_entries = id_to_col.size - length(weights) + 1000
                resize!(weights, length(weights) + new_entries)
                resize!(colnorm2, length(colnorm2) + new_entries)
                resize!(spectral_scores, length(spectral_scores) + new_entries)
                # Avoid list comprehension allocation - use direct resize and loop
                old_length = length(unscored_psms)
                resize!(unscored_psms, old_length + new_entries)
                for i in (old_length + 1):length(unscored_psms)
                    unscored_psms[i] = eltype(unscored_psms)()
                end
            end

            initialize_weights!(id_to_col, weights, precursor_weights)

            # Solve deconvolution
            solve_deconvolution!(
                calibrated_chromatogram_deconvolution_solver(search_context, params.deconvolution_solver),
                Hs, residuals, weights, colnorm2,
                getMu(search_data), getObserved(search_data),
                params.max_iter_outer, params.max_diff
            )

            # NEW: Distribute grouped coefficients back to individual precursors
            distribute_ms1_coefficients!(
                precursor_weights,  # Array indexed by precursor ID
                weights,            # Array indexed by column number (group coefficients)
                mz_grouping
            )

            # Record chromatogram points with weights
            for j in 1:prec_temp_size
                rt_idx += 1
                if rt_idx + 1 > length(chromatograms)
                    resize!(chromatograms, length(chromatograms) * 2)  # Exponential growth
                end

                # Get weight from precursor_weights array (now contains distributed group coefficients)
                prec_id = precs_temp[j]
                weight = prec_id <= length(precursor_weights) ? precursor_weights[prec_id] : 0.0f0

                if weight > 0.0f0
                    chromatograms[rt_idx] = MS1ChromObject(
                        Float32(getRetentionTime(spectra, scan_idx)),
                        weight,
                        iso_count[prec_id].matched_mono,
                        iso_count[prec_id].iso_count,
                        scan_idx,
                        prec_id
                    )
                else
                    chromatograms[rt_idx] = MS1ChromObject(
                        Float32(getRetentionTime(spectra, scan_idx)),
                        zero(Float32),
                        false,
                        UInt8(0),
                        scan_idx,
                        precs_temp[j]
                    )
                end
            end

            # Update precursor weights - already handled by distribute_ms1_coefficients!
            # No need to update here since distribution already populated precursor_weights

        else
            for j in 1:prec_temp_size
                rt_idx += 1
                if rt_idx + 1 > length(chromatograms)
                    resize!(chromatograms, length(chromatograms) * 2)  # Exponential growth
                end
                chromatograms[rt_idx] = MS1ChromObject(
                    Float32(getRetentionTime(spectra, scan_idx)),
                    zero(Float32),
                    false,
                    UInt8(0),
                    scan_idx,
                    precs_temp[j]
                )
            end
        end

        reset_scan_arrays!(id_to_col, Hs, unscored_psms)
    end

    return DataFrame(@view(chromatograms[1:rt_idx]))
end
=#

#==========================================================
Chromatogram Building Functions
==========================================================#
"""
    process_final_psms!(psms::DataFrame, search_context::SearchContext,
                       parsed_fname::String, ms_file_idx::Int64)

Process final PSMs after integration.

# Added Columns
- Protein information (accession numbers, indices)
- Peak areas and normalization
- Sequence information
- Modification details
- File metadata
"""
function process_final_psms!(
    psms::DataFrame,
    search_context::SearchContext,
    parsed_fname::String,
    ms_file_idx::Int64
)
    # Drop only genuinely unusable areas. A zero area is an identification that
    # passed FDR but could not be quantified -- either no peak survived baseline
    # subtraction or integration rejected the window (see
    # QUANT_MIN_AREA_SURVIVING_RATIO). Those rows are kept so the identification
    # reaches protein inference and the output tables; downstream treats a zero
    # area as "not quantified in this run" rather than as an abundance of zero.
    n_before = nrow(psms)
    filter!(row -> !isnan(row.peak_area::Float32), psms)
    n_nan = n_before - nrow(psms)
    n_unquantified = count(row -> row.peak_area::Float32 <= 0.0f0, eachrow(psms))
    if n_nan > 0
        @debug_l1 "[$parsed_fname] dropped $n_nan / $n_before PSMs with NaN peak_area"
    end
    if n_unquantified > 0
        @debug_l1 "[$parsed_fname] $n_unquantified / $(nrow(psms)) PSMs identified " *
                  "but not quantified (peak_area == 0)"
    end
    # Add columns
    precursors = getPrecursors(getSpecLib(search_context))
    n = size(psms, 1)
    accession_numbers = Vector{String}(undef, n)
    ms_file_idxs = Vector{UInt16}(undef, n)
    species = Vector{String}(undef, n)
    peak_area = Vector{Union{Missing, Float32}}(undef, n)
    peak_area_normalized = Vector{Union{Missing, Float32}}(undef, n)
    structural_mods = Vector{Union{Missing, String}}(undef, n)
    isotopic_mods = Vector{Union{Missing, String}}(undef, n)
    charge = Vector{UInt8}(undef, n)
    sequence = Vector{String}(undef, n)
    peptide_start_positions = Vector{String}(undef, n)
    file_name = Vector{String}(undef, n)
    psms_precursor_idx = psms[!,:precursor_idx]::Vector{UInt32}

    accession_col = getAccessionNumbers(precursors)
    for i in range(1, n)
        pid = psms_precursor_idx[i]
        accession_numbers[i] = accession_col[pid]
    end
    psms[!, :accession_numbers] = accession_numbers

    # No sort here. MaxLFQSearch sorts the merged PSMs by :inferred_protein_group
    # before its chunked-merge so chunk boundaries align with protein-group
    # boundaries. ProteinInferenceSearch (which runs between this method and
    # MaxLFQ) is what populates :inferred_protein_group, so a sort by that
    # column at this stage isn't possible anyway.

    parsed_fname = getFileIdToName(getMSData(search_context), ms_file_idx)
    proteome_col = getProteomeIdentifiers(precursors)
    structural_mods_col = getStructuralMods(precursors)
    isotopic_mods_col = getIsotopicMods(precursors)
    charge_col = getCharge(precursors)
    sequence_col = getSequence(precursors)
    start_idx_col = getStartIdx(precursors)
    for i in range(1, n)
        pid = psms_precursor_idx[i]
        ms_file_idxs[i] = UInt32(ms_file_idx)
        species[i] = join(sort(unique(split(coalesce(proteome_col[pid], ""),';'))),';')
        peak_area[i] = psms[i,:peak_area]
        peak_area_normalized[i] = zero(Float32)
        structural_mods[i] = structural_mods_col[pid]
        isotopic_mods[i] = isotopic_mods_col[pid]
        charge[i] = charge_col[pid]
        sequence[i] = sequence_col[pid]
        peptide_start_positions[i] = _format_start_idx(start_idx_col[pid])
        file_name[i] = parsed_fname
    end

    psms[!,:ms_file_idx] = ms_file_idxs
    psms[!,:species] = species
    psms[!,:peak_area] = peak_area
    psms[!,:peak_area_normalized] = peak_area_normalized
    psms[!,:structural_mods] = structural_mods
    psms[!,:isotopic_mods] = isotopic_mods
    psms[!,:charge] = charge 
    psms[!,:sequence] = sequence
    psms[!,:peptide_start_positions] = peptide_start_positions
    psms[!,:file_name] = file_name
    
    return nothing
end
