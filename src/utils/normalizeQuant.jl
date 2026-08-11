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
    _min_bins_for_spline(spline_n_knots) -> Int

Fewest RT bins that may be handed to `UniformSpline`.

Strictly greater than the coefficient count, not equal to it. `UniformSpline`'s
design matrix has `n_knots + 3` columns of which the last can never hold a
non-zero value, so at exactly `n_knots + 3` bins the system is square, `\\`
dispatches to `lu`, and the empty column raises `SingularException`. One more bin
makes it over-determined and `\\` takes the QR path instead.
"""
_min_bins_for_spline(spline_n_knots::Int) = spline_n_knots + 3 + 1

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

"""
    getQuantSplines(psms_paths, quant_col_name; N=100, spline_n_knots=7,
                    min_bin_occupancy=QUANT_MIN_BIN_OCCUPANCY)

Fit RT-dependent quantification splines for each MS file. For each Arrow file,
bins PSMs by observed iRT into bins of at least `min_bin_occupancy` rows (at most
`N` bins), computes median log2-abundance per bin, and fits a `UniformSpline`
mapping iRT -> median log2-abundance.

Files whose PSM count cannot support `_min_bins_for_spline(spline_n_knots)` such
bins get no entry in the returned dictionary: they are left uncorrected and
excluded from the cross-file median. `applyNormalization!` handles the absent key.

Returns `(splines_dict, (min_rt, max_rt))` where `splines_dict` maps file path to
fitted spline, and the RT range spans only the files that got one.
"""
function getQuantSplines(psms_paths::Vector{String},
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7,
    min_bin_occupancy::Int = QUANT_MIN_BIN_OCCUPANCY)
    quant_splines = Dictionary{String, UniformSpline}()
    min_rt, max_rt = typemax(Float32), typemin(Float32)
    min_bins = _min_bins_for_spline(spline_n_knots)
    for fpath in psms_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        fast_df_sort!(psms, (:irt_obs,))
        nprecs = size(psms, 1)

        bins = _occupancy_bins(nprecs, min_bin_occupancy, N)
        if length(bins) < min_bins
            @user_warn "Skipping quant normalization for $(basename(fpath)): " *
                       "$nprecs PSMs support only $(length(bins)) RT bin(s) at " *
                       ">= $min_bin_occupancy PSMs each, and a " *
                       "$(spline_n_knots + 3)-coefficient spline needs at least " *
                       "$min_bins. Abundances are left uncorrected and excluded " *
                       "from the cross-file median."
            continue
        end

        # Only files that get a spline define the correction grid. A skipped file
        # receives no correction, so widening the grid to its RT range would only
        # extrapolate the surviving splines further than their data supports.
        if minimum(psms[!,:irt_obs]) < min_rt
            min_rt = minimum(psms[!,:irt_obs])
        end
        if maximum(psms[!,:irt_obs]) > max_rt
            max_rt = maximum(psms[!,:irt_obs])
        end

        median_quant = Vector{Float64}(undef, length(bins))
        median_rts = Vector{Float64}(undef, length(bins))
        for (i, b) in enumerate(bins)
            median_rts[i] = median(@view(psms[b, :irt_obs]))
            median_quant[i] = log2(median(@view(psms[b, quant_col_name])))
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
"""
function getQuantCorrections(
    quant_splines::Dictionary{String, UniformSpline},
    rt_range::Tuple{AbstractFloat, AbstractFloat};
    N = 100)
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
        insert!(corrections, key, linear_interpolation(rt_grid, offset))
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

End-to-end RT-dependent quantification normalization. Reads all Arrow files
in `second_quant_folder`, fits per-file splines, computes cross-file median,
and writes corrected abundances back to each file.

This corrects for systematic RT-dependent intensity differences between MS runs,
ensuring that the same protein at the same RT has comparable abundance across files.
"""
function normalizeQuant(
    psms_paths::Vector{String},
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7)

    quant_splines_dict, rt_range = getQuantSplines(
        psms_paths, quant_col_name, N = N, spline_n_knots = spline_n_knots)

    quant_corrections_dict = getQuantCorrections(
        quant_splines_dict, rt_range, N = N)

    applyNormalization!(psms_paths, quant_col_name, quant_corrections_dict)
    return nothing
end

function normalizeQuant(
    second_quant_folder::String,
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7)

    psms_paths = [fpath for fpath in readdir(second_quant_folder, join=true) if endswith(fpath, ".arrow")]
    normalizeQuant(psms_paths, quant_col_name; N = N, spline_n_knots = spline_n_knots)
end
