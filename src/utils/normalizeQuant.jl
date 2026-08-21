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
Minimum PSMs per RT bin. A bin median is only meaningful when taken over enough
values; the previous fixed-bin-count scheme let this fall to 1, at which point
`median_rts` was just the sorted per-PSM retention times and the fitted
correction was noise. See dev_docs/quant_spline_guard.
"""
const QUANT_MIN_BIN_OCCUPANCY = 100

"""
Minimum fraction of runs in which a precursor must be quantified to serve as a
normalization anchor. With two or more runs, at least two observations are
always required; a one-run analysis uses its own precursors and therefore
reduces to the identity correction.
"""
const QUANT_MIN_ANCHOR_RUN_FRACTION = 0.5

"""
Minimum matched precursors required to estimate a run-wide loading offset.
Unlike an RT spline, a single robust median does not require multiple bins, but
it still needs enough anchors to avoid applying a noisy constant correction.
"""
const QUANT_MIN_GLOBAL_ANCHORS = 100

"""
    _min_bins_for_spline(spline_n_knots) -> Int

Fewest RT bins that may be handed to `UniformSpline`.

One more than the coefficient count, so the least-squares system is
over-determined rather than exactly determined. The design matrix is full rank as
of the `_n_spline_coeffs` fix, so a square system is no longer a failure -- but
fitting `n` coefficients through exactly `n` points interpolates the bin medians
instead of smoothing them, which is not what this estimator is for.
"""
_min_bins_for_spline(spline_n_knots::Int) = _n_spline_coeffs(spline_n_knots) + 1

"""
    _occupancy_bins(n, min_occupancy, max_bins) -> Vector{UnitRange{Int}}

Contiguous bin ranges over `n` RT-sorted rows, each holding at least
`min_occupancy` rows, and at most `max_bins` bins.

Rather than forming `max_bins` bins and merging the undersized ones pairwise,
this computes how many bins the data can actually fill (`n ÷ min_occupancy`) and
lays them out evenly -- the same outcome with no merge loop to get wrong. Bin
sizes differ by at most one and every row is covered, unlike the previous
`bin_stop = bin_idx * bin_size` scheme, which silently dropped the final
`n mod max_bins` rows.

Returns a single range when `n` supports only one bin, and an empty vector when
`n < min_occupancy`; the caller treats both as "cannot fit a spline here".
"""
function _occupancy_bins(n::Int, min_occupancy::Int, max_bins::Int)
    (n <= 0 || min_occupancy <= 0 || max_bins <= 0) && return UnitRange{Int}[]
    n >= min_occupancy || return UnitRange{Int}[]
    k = min(max_bins, fld(n, min_occupancy))
    k <= 1 && return [1:n]
    base, rem = divrem(n, k)
    bins = Vector{UnitRange{Int}}(undef, k)
    start = 1
    @inbounds for i in 1:k
        len = base + (i <= rem ? 1 : 0)
        bins[i] = start:(start + len - 1)
        start += len
    end
    return bins
end

@inline function _required_anchor_runs(n_runs::Int, min_run_fraction::Real)
    n_runs <= 0 && return 0
    0.0 < min_run_fraction <= 1.0 || throw(ArgumentError(
        "min_run_fraction must be in (0, 1], got $min_run_fraction"
    ))
    n_runs == 1 && return 1
    return min(n_runs, max(2, ceil(Int, n_runs * min_run_fraction)))
end

function _file_precursor_log2_quant(
    psms::DataFrame,
    quant_col_name::Symbol,
    precursor_col_name::Symbol,
)
    hasproperty(psms, precursor_col_name) || throw(ArgumentError(
        "Matched-precursor quant normalization requires column " *
        "$(precursor_col_name)."
    ))
    hasproperty(psms, quant_col_name) || throw(ArgumentError(
        "Quant normalization requires column $(quant_col_name)."
    ))

    precursor_idx = psms[!, precursor_col_name]
    quant = psms[!, quant_col_name]
    values = Dict{UInt32, Float32}()
    sizehint!(values, length(precursor_idx))

    @inbounds for i in eachindex(precursor_idx, quant)
        q_raw = quant[i]
        ismissing(q_raw) && continue
        q = Float64(q_raw)
        (isfinite(q) && q > 0.0) || continue
        pid = UInt32(precursor_idx[i])
        logq = Float32(log2(q))

        # Passing-PSM tables normally contain one row per precursor. If a
        # duplicate is present, keep one run-level value so that run cannot
        # receive extra weight in the cross-run precursor consensus.
        values[pid] = max(get(values, pid, -Inf32), logq)
    end
    return values
end

"""
    getPrecursorQuantConsensus(psms_paths, quant_col_name;
        precursor_col_name=:precursor_idx,
        min_run_fraction=QUANT_MIN_ANCHOR_RUN_FRACTION)

Build a cross-run reference abundance for matched precursors. Each precursor's
reference is the median log2 abundance across runs in which it was quantified.
Only precursors observed in at least `min_run_fraction` of runs (and at least two
runs when multiple runs are present) are retained as normalization anchors.

The run-level de-duplication is deliberate: a precursor contributes at most one
value per run to the consensus.
"""
function getPrecursorQuantConsensus(
    psms_paths::Vector{String},
    quant_col_name::Symbol;
    precursor_col_name::Symbol = :precursor_idx,
    min_run_fraction::Real = QUANT_MIN_ANCHOR_RUN_FRACTION,
)
    required_runs = _required_anchor_runs(length(psms_paths), min_run_fraction)
    required_runs == 0 && return Dict{UInt32, Float32}()

    values_by_precursor = Dict{UInt32, Vector{Float32}}()
    for fpath in psms_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        file_values = _file_precursor_log2_quant(
            psms,
            quant_col_name,
            precursor_col_name,
        )
        for (pid, logq) in file_values
            push!(get!(() -> Float32[], values_by_precursor, pid), logq)
        end
    end

    consensus = Dict{UInt32, Float32}()
    sizehint!(consensus, length(values_by_precursor))
    for (pid, values) in values_by_precursor
        length(values) >= required_runs || continue
        consensus[pid] = Float32(median(values))
    end
    return consensus
end

"""
    getMatchedGlobalOffsets(psms_paths, quant_col_name;
        precursor_col_name=:precursor_idx,
        min_run_fraction=QUANT_MIN_ANCHOR_RUN_FRACTION,
        min_anchors=QUANT_MIN_GLOBAL_ANCHORS)

Estimate one global log2 loading offset per run from matched precursor ratios.
For every run, computes the median of

`log2(file abundance) - precursor consensus`

over eligible anchors, then subtracts the median run offset so that the
experiment-wide correction remains centered at zero. Runs with fewer than
`min_anchors` matched precursors are omitted and therefore left uncorrected.
"""
function getMatchedGlobalOffsets(
    psms_paths::Vector{String},
    quant_col_name::Symbol;
    precursor_col_name::Symbol = :precursor_idx,
    min_run_fraction::Real = QUANT_MIN_ANCHOR_RUN_FRACTION,
    min_anchors::Int = QUANT_MIN_GLOBAL_ANCHORS,
)
    min_anchors > 0 || throw(ArgumentError(
        "min_anchors must be positive, got $min_anchors"
    ))
    precursor_consensus = getPrecursorQuantConsensus(
        psms_paths,
        quant_col_name;
        precursor_col_name = precursor_col_name,
        min_run_fraction = min_run_fraction,
    )

    raw_offsets = Dictionary{String, Float64}()
    raw_offset_values = Float64[]
    for fpath in psms_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        file_values = _file_precursor_log2_quant(
            psms,
            quant_col_name,
            precursor_col_name,
        )
        residuals = Float64[]
        sizehint!(residuals, min(length(file_values), length(precursor_consensus)))
        for (pid, file_logq) in file_values
            reference_logq = get(precursor_consensus, pid, nothing)
            reference_logq === nothing && continue
            push!(residuals, Float64(file_logq - reference_logq))
        end

        if length(residuals) < min_anchors
            @user_warn "Skipping quant normalization for $(basename(fpath)): " *
                       "$(length(residuals)) matched precursor anchors are " *
                       "below the $min_anchors required for a global loading " *
                       "offset. Abundances are left uncorrected."
            continue
        end

        raw_offset = median(residuals)
        insert!(raw_offsets, fpath, raw_offset)
        push!(raw_offset_values, raw_offset)
    end

    isempty(raw_offsets) && return raw_offsets
    experiment_center = median(raw_offset_values)
    offsets = Dictionary{String, Float64}()
    for (fpath, raw_offset) in pairs(raw_offsets)
        insert!(offsets, fpath, raw_offset - experiment_center)
    end
    return offsets
end

"""
    getQuantSplines(psms_paths, quant_col_name; N=100, spline_n_knots=7,
                    min_bin_occupancy=QUANT_MIN_BIN_OCCUPANCY)

Fit marginal RT-dependent log2-intensity splines for each MS file. All valid
observed precursors contribute to the shape estimator. Rows are sorted by
observed iRT, divided into bins of at least `min_bin_occupancy` observations,
and summarized by the median RT and log2 median abundance.

These splines still contain both a run-wide level and an RT shape.
`getRTShapeCorrections` removes the run-wide level before comparing shapes;
the global loading correction comes separately from matched precursors.

Files whose observed precursor count cannot support
`_min_bins_for_spline(spline_n_knots)` such bins get no entry in the returned
dictionary. They can still receive a constant matched global correction.

Returns `(splines_dict, (min_rt, max_rt))` where `splines_dict` maps file path to
fitted spline, and the RT range spans only the files that got one.
"""
function getQuantSplines(
    psms_paths::Vector{String},
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7,
    min_bin_occupancy::Int = QUANT_MIN_BIN_OCCUPANCY,
)
    quant_splines = Dictionary{String, UniformSpline}()
    min_rt, max_rt = typemax(Float32), typemin(Float32)
    min_bins = _min_bins_for_spline(spline_n_knots)
    for fpath in psms_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        hasproperty(psms, :irt_obs) || throw(ArgumentError(
            "RT-dependent quant normalization requires column irt_obs."
        ))
        hasproperty(psms, quant_col_name) || throw(ArgumentError(
            "Quant normalization requires column $(quant_col_name)."
        ))

        observed_rts = Float64[]
        observed_quant = Float64[]
        irt_obs = psms[!, :irt_obs]
        quant = psms[!, quant_col_name]
        sizehint!(observed_rts, length(irt_obs))
        sizehint!(observed_quant, length(quant))
        @inbounds for i in eachindex(irt_obs, quant)
            rt_raw = irt_obs[i]
            quant_raw = quant[i]
            (ismissing(rt_raw) || ismissing(quant_raw)) && continue
            rt = Float64(rt_raw)
            abundance = Float64(quant_raw)
            (isfinite(rt) && isfinite(abundance) && abundance > 0.0) || continue
            push!(observed_rts, rt)
            push!(observed_quant, abundance)
        end

        order = sortperm(observed_rts)
        observed_rts = observed_rts[order]
        observed_quant = observed_quant[order]
        nobservations = length(observed_rts)

        bins = _occupancy_bins(nobservations, min_bin_occupancy, N)
        if length(bins) < min_bins
            @user_warn "Skipping RT-shape estimation for $(basename(fpath)): " *
                       "$nobservations observed precursors support only " *
                       "$(length(bins)) RT bin(s) at >= $min_bin_occupancy " *
                       "precursors each, and a " *
                       "$(_n_spline_coeffs(spline_n_knots))-coefficient spline needs at least " *
                       "$min_bins. Only the matched global correction will be " *
                       "applied when available."
            continue
        end

        if first(observed_rts) < min_rt
            min_rt = first(observed_rts)
        end
        if last(observed_rts) > max_rt
            max_rt = last(observed_rts)
        end

        median_quant = Vector{Float64}(undef, length(bins))
        median_rts = Vector{Float64}(undef, length(bins))
        for (i, b) in enumerate(bins)
            median_rts[i] = median(@view(observed_rts[b]))
            median_quant[i] = log2(median(@view(observed_quant[b])))
        end

        splinefit = UniformSpline(median_quant, median_rts, 3, spline_n_knots)
        insert!(quant_splines, fpath, splinefit)
    end
    return quant_splines, (min_rt, max_rt)
end

"""
    getQuantCorrections(quant_splines, rt_range; N=100)

Compute per-file RT-dependent correction offsets. Evaluates each file's spline
on a grid, computes the cross-file median at each RT point, and returns a
dictionary of interpolated offset functions (file_spline - global_median).

This is the original full marginal correction and is retained as a baseline.
The decomposed `normalizeQuant` pipeline instead calls
`getRTShapeCorrections`, which removes marginal run levels before constructing
the RT component.

Returns an empty dictionary when no file got a spline: there is no cross-file
median to correct toward, and `rt_range` is still the empty-range sentinel
`(Inf, -Inf)`. `applyNormalization!` reads that as "no correction for this file"
and copies the abundances through unchanged (correction factor 1.0).
"""
function getQuantCorrections(
    quant_splines::Dictionary{String, UniformSpline},
    rt_range::Tuple{AbstractFloat, AbstractFloat};
    N = 100)
    isempty(quant_splines) && return Dictionary{String, Any}()
    median_quant = zeros(Float32, (length(quant_splines), N))
    rt_grid = collect(LinRange(first(rt_range), last(rt_range), N))
    for (i, rt) in enumerate(rt_grid)
        j = 1
        for (key, spline) in pairs(quant_splines)
            median_quant[j, i] = spline(rt_grid[i])
            j += 1
        end
    end
    median_quant = reshape(median(median_quant, dims = 1), (N,))
    median_interp = linear_interpolation(rt_grid, median_quant)
    corrections = Dictionary{String, Any}()
    for (key, spline) in pairs(quant_splines)
        offset = zeros(Float32, N)
        for (i, rt) in enumerate(rt_grid)
            offset[i] = spline(rt) - median_interp(rt)
        end
        # Rows outside the spline grid still need a correction. Hold the nearest
        # supported value constant instead of extrapolating an unobserved trend
        # (or throwing on an out-of-range lookup).
        insert!(corrections, key, linear_interpolation(
            rt_grid,
            offset;
            extrapolation_bc = Interpolations.Flat(),
        ))
    end
    return corrections
end

function _median_polish_interactions!(
    values::Matrix{Float64};
    max_iterations::Int = 20,
    tolerance::Float64 = 1.0e-7,
)
    for _ in 1:max_iterations
        row_centers = median(values, dims = 2)
        values .-= row_centers
        column_centers = median(values, dims = 1)
        values .-= column_centers
        max_adjustment = max(
            maximum(abs, row_centers),
            maximum(abs, column_centers),
        )
        max_adjustment <= tolerance && break
    end
    return values
end

"""
    getRTShapeCorrections(quant_splines, rt_range; N=100)

Extract only run-specific RT shape differences from marginal intensity splines.
Spline values on a common RT grid are robustly double-centered with median
polish: per-run levels and the common across-run RT profile are removed. The
remaining interaction has approximately zero median across RT for each run and
zero median across runs at each RT, so it cannot replace the matched global
loading estimate.
"""
function getRTShapeCorrections(
    quant_splines::Dictionary{String, UniformSpline},
    rt_range::Tuple{AbstractFloat, AbstractFloat};
    N::Int = 100,
)
    isempty(quant_splines) && return Dictionary{String, Any}()
    rt_grid = collect(LinRange(first(rt_range), last(rt_range), N))
    file_paths = String[]
    spline_values = Matrix{Float64}(undef, length(quant_splines), N)
    for (row, (fpath, spline)) in enumerate(pairs(quant_splines))
        push!(file_paths, fpath)
        @inbounds for (column, rt) in enumerate(rt_grid)
            spline_values[row, column] = spline(rt)
        end
    end

    _median_polish_interactions!(spline_values)
    corrections = Dictionary{String, Any}()
    for (row, fpath) in enumerate(file_paths)
        insert!(corrections, fpath, linear_interpolation(
            rt_grid,
            copy(@view(spline_values[row, :]));
            extrapolation_bc = Interpolations.Flat(),
        ))
    end
    return corrections
end

"""
    combineQuantCorrections(global_offsets, rt_shape_corrections)

Combine the matched run-wide offset with the centered marginal RT shape.
Only files with a supported matched global offset are normalized. If their RT
shape spline was skipped, they receive the constant global correction alone.
"""
function combineQuantCorrections(
    global_offsets::Dictionary{String, Float64},
    rt_shape_corrections::Dictionary{String, Any},
)
    corrections = Dictionary{String, Any}()
    for (fpath, global_offset) in pairs(global_offsets)
        rt_shape = get(rt_shape_corrections, fpath, nothing)
        if rt_shape === nothing
            let offset = global_offset
                insert!(corrections, fpath, _ -> offset)
            end
        else
            let offset = global_offset, shape = rt_shape
                insert!(corrections, fpath, rt -> offset + shape(rt))
            end
        end
    end
    return corrections
end

"""
    applyNormalization!(psms_paths, quant_col, corrections)

Apply RT-dependent normalization corrections to Arrow files in place.
For each PSM, subtracts the file-specific log2-offset at its observed iRT,
writing a new `{quant_col}_normalized` column.
"""
function applyNormalization!(
    psms_paths::Vector{String},
    quant_col::Symbol,
    corrections::Dictionary{String, Any}
)
    for fpath in psms_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        norm_quant_col = Symbol(string(quant_col) * "_normalized")
        # A file without enough matched anchors receives no correction. Emit the
        # column unchanged rather than omitting it: downstream readers assume
        # every file carries the same schema.
        correction_spline = get(corrections, fpath, nothing)
        if correction_spline === nothing
            psms[!, norm_quant_col] = Float64.(psms[!, quant_col])
            writeArrow(fpath, psms)
            continue
        end
        n = size(psms, 1)
        norm_vals = Vector{Float64}(undef, n)
        for i in range(1, n)
            hc = correction_spline(psms[i, :irt_obs])
            norm_vals[i] = 2^(log2(max(psms[i, quant_col], 0.0)) - hc)
        end
        psms[!, norm_quant_col] = norm_vals
        writeArrow(fpath, psms)
    end
end

"""
    normalizeQuant(second_quant_folder, quant_col_name; N=100, spline_n_knots=7)

End-to-end decomposed quantification normalization. Matched precursor ratios
estimate one global loading offset per run. All observed precursors estimate a
separate marginal RT spline whose run-wide level and common RT profile are
removed by median polish. The constant matched offset and centered RT shape are
then combined and written back to each file.

This uses precursor matching where differential missingness is most damaging
(global loading) while retaining the larger observed population for the local
RT shape.
"""
function normalizeQuant(
    psms_paths::Vector{String},
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7,
    min_anchor_run_fraction::Real = QUANT_MIN_ANCHOR_RUN_FRACTION)

    global_offsets = getMatchedGlobalOffsets(
        psms_paths,
        quant_col_name;
        min_run_fraction = min_anchor_run_fraction,
    )

    shape_paths = String[
        fpath for fpath in psms_paths if haskey(global_offsets, fpath)
    ]
    quant_splines, rt_range = getQuantSplines(
        shape_paths,
        quant_col_name;
        N = N,
        spline_n_knots = spline_n_knots,
    )
    rt_shape_corrections = getRTShapeCorrections(
        quant_splines,
        rt_range;
        N = N,
    )
    quant_corrections = combineQuantCorrections(
        global_offsets,
        rt_shape_corrections,
    )

    applyNormalization!(psms_paths, quant_col_name, quant_corrections)
    return nothing
end

function normalizeQuant(
    second_quant_folder::String,
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7,
    min_anchor_run_fraction::Real = QUANT_MIN_ANCHOR_RUN_FRACTION)

    psms_paths = [fpath for fpath in readdir(second_quant_folder, join=true) if endswith(fpath, ".arrow")]
    normalizeQuant(
        psms_paths,
        quant_col_name;
        N = N,
        spline_n_knots = spline_n_knots,
        min_anchor_run_fraction = min_anchor_run_fraction,
    )
end
