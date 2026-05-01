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
    getQuantSplines(psms_paths, quant_col_name; N=100, spline_n_knots=7)

Fit RT-dependent quantification splines for each MS file. For each Arrow file,
bins PSMs by observed iRT, computes median log2-abundance per bin, and fits a
`UniformSpline` mapping iRT -> median log2-abundance.

Returns `(splines_dict, (min_rt, max_rt))` where `splines_dict` maps file path
to fitted spline, and the RT range spans all files.
"""
function getQuantSplines(psms_paths::Vector{String},
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7)
    quant_splines = Dictionary{String, UniformSpline}()
    min_rt, max_rt = typemax(Float32), typemin(Float32)
    for fpath in psms_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        N_sample = min(N, size(psms, 1))
        fast_df_sort!(psms, (:irt_obs,))
        if minimum(psms[!,:irt_obs]) < min_rt
            min_rt = minimum(psms[!,:irt_obs])
        end
        if maximum(psms[!,:irt_obs]) > max_rt
            max_rt = maximum(psms[!,:irt_obs])
        end
        nprecs = size(psms, 1)
        bin_size = nprecs ÷ N_sample
        median_quant = zeros(Float64, N_sample)
        median_rts = zeros(Float64, N_sample)
        for bin_idx in range(1, N_sample)
            bin_start = (bin_idx - 1) * bin_size + 1
            bin_stop = bin_idx * bin_size
            median_rts[bin_idx] = median(@view(psms[bin_start:bin_stop, :irt_obs]))
            median_quant[bin_idx] = log2(median(@view(psms[bin_start:bin_stop, quant_col_name])))
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
        correction_spline = corrections[fpath]
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
