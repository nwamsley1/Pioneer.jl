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
    getQuantSplines(psms_paths, quant_col_name; N=100, spline_n_knots=7,
                    min_bin_occupancy=QUANT_MIN_BIN_OCCUPANCY,
                    min_anchor_run_fraction=QUANT_MIN_ANCHOR_RUN_FRACTION)

Fit matched-precursor RT-dependent quantification splines for each MS file.
First builds a cross-run median log2 abundance for precursors observed in a
configurable fraction of runs. For each file, it then computes

`log2(file abundance) - precursor consensus`

for those matched anchors, bins the residuals by observed iRT into bins of at
least `min_bin_occupancy` anchors (at most `N` bins), and fits a `UniformSpline`
mapping iRT to the median matched-precursor residual.

Files whose matched-anchor count cannot support
`_min_bins_for_spline(spline_n_knots)` such bins get no entry in the returned
dictionary: they are left uncorrected and excluded from the cross-file median.
`applyNormalization!` handles the absent key.

Returns `(splines_dict, (min_rt, max_rt))` where `splines_dict` maps file path to
fitted spline, and the RT range spans only the files that got one.
"""
function getQuantSplines(psms_paths::Vector{String},
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7,
    min_bin_occupancy::Int = QUANT_MIN_BIN_OCCUPANCY,
    min_anchor_run_fraction::Real = QUANT_MIN_ANCHOR_RUN_FRACTION,
    precursor_col_name::Symbol = :precursor_idx)
    precursor_consensus = getPrecursorQuantConsensus(
        psms_paths,
        quant_col_name;
        precursor_col_name = precursor_col_name,
        min_run_fraction = min_anchor_run_fraction,
    )
    quant_splines = Dictionary{String, UniformSpline}()
    min_rt, max_rt = typemax(Float32), typemin(Float32)
    min_bins = _min_bins_for_spline(spline_n_knots)
    for fpath in psms_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        hasproperty(psms, :irt_obs) || throw(ArgumentError(
            "Matched-precursor quant normalization requires column irt_obs."
        ))
        file_values = _file_precursor_log2_quant(
            psms,
            quant_col_name,
            precursor_col_name,
        )

        anchor_rts = Float64[]
        anchor_residuals = Float64[]
        sizehint!(anchor_rts, min(length(file_values), length(precursor_consensus)))
        sizehint!(anchor_residuals, min(length(file_values), length(precursor_consensus)))
        precursor_idx = psms[!, precursor_col_name]
        irt_obs = psms[!, :irt_obs]
        seen = Set{UInt32}()
        sizehint!(seen, length(file_values))
        @inbounds for i in eachindex(precursor_idx, irt_obs)
            pid = UInt32(precursor_idx[i])
            pid in seen && continue
            reference_logq = get(precursor_consensus, pid, nothing)
            reference_logq === nothing && continue
            file_logq = get(file_values, pid, nothing)
            file_logq === nothing && continue
            rt_raw = irt_obs[i]
            ismissing(rt_raw) && continue
            rt = Float64(rt_raw)
            isfinite(rt) || continue
            push!(seen, pid)
            push!(anchor_rts, rt)
            push!(anchor_residuals, Float64(file_logq - reference_logq))
        end

        order = sortperm(anchor_rts)
        anchor_rts = anchor_rts[order]
        anchor_residuals = anchor_residuals[order]
        nanchors = length(anchor_rts)

        bins = _occupancy_bins(nanchors, min_bin_occupancy, N)
        if length(bins) < min_bins
            @user_warn "Skipping quant normalization for $(basename(fpath)): " *
                       "$nanchors matched precursor anchors support only " *
                       "$(length(bins)) RT bin(s) at >= $min_bin_occupancy " *
                       "anchors each, and a " *
                       "$(_n_spline_coeffs(spline_n_knots))-coefficient spline needs at least " *
                       "$min_bins. Abundances are left uncorrected and excluded " *
                       "from the cross-file median."
            continue
        end

        # Only files that get a spline define the correction grid. A skipped file
        # receives no correction, so widening the grid to its RT range would only
        # extrapolate the surviving splines further than their data supports.
        if first(anchor_rts) < min_rt
            min_rt = first(anchor_rts)
        end
        if last(anchor_rts) > max_rt
            max_rt = last(anchor_rts)
        end

        median_quant = Vector{Float64}(undef, length(bins))
        median_rts = Vector{Float64}(undef, length(bins))
        for (i, b) in enumerate(bins)
            median_rts[i] = median(@view(anchor_rts[b]))
            median_quant[i] = median(@view(anchor_residuals[b]))
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
        # Censoring can leave no matched anchors at one or both RT boundaries.
        # Rows outside the anchor-supported range still need a correction; hold
        # the nearest supported value constant instead of extrapolating a trend
        # that was never observed (or throwing on an out-of-range lookup).
        insert!(corrections, key, linear_interpolation(
            rt_grid,
            offset;
            extrapolation_bc = Interpolations.Flat(),
        ))
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
        # getQuantSplines skips files that cannot support a spline, so this key may
        # be absent. Emit the column unchanged rather than omitting it: downstream
        # readers assume every file carries the same schema.
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

End-to-end RT-dependent quantification normalization. Reads all Arrow files,
constructs a cross-run matched-precursor reference, fits per-file residual
splines, computes the cross-file median, and writes corrected abundances back
to each file.

This corrects systematic RT-dependent intensity differences between MS runs
without forcing the marginal observed intensity distributions to match.
"""
function normalizeQuant(
    psms_paths::Vector{String},
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7,
    min_anchor_run_fraction::Real = QUANT_MIN_ANCHOR_RUN_FRACTION)

    quant_splines_dict, rt_range = getQuantSplines(
        psms_paths,
        quant_col_name;
        N = N,
        spline_n_knots = spline_n_knots,
        min_anchor_run_fraction = min_anchor_run_fraction,
    )

    quant_corrections_dict = getQuantCorrections(
        quant_splines_dict, rt_range, N = N)

    applyNormalization!(psms_paths, quant_col_name, quant_corrections_dict)
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
