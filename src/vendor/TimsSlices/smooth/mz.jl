# Copyright (C) 2026 Nathan Wamsley
#
# This file is part of TimsSlices.jl
#
# TimsSlices.jl is free software: you can redistribute it and/or modify
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

# Step 2 of the pipeline (see window.jl): smooth one slice along m/z and hand the result to the peak picker.
#
# Why smooth along m/z at all: Bruker's raw scans are already centroided (one entry per ion per scan, at an
# integer TOF bin), but an ion's bin jitters by +/-1-2 between neighbouring IM scans (E. coli, strongest peaks:
# same bin in the next scan 21%, +/-1 33%, +/-2 16%). After the IM sum, one ion is therefore spread over 2-4
# adjacent bins with a ragged profile, which a peak picker would split into several centroids. A Gaussian of
# `mz_sigma` bins (0.3125 ns of flight time by default: 2.5 bins on the timsTOF Ultra, 1.56 on the Pro; see
# MZ_SIGMA_NS) turns it into one smooth hump with one apex; the weighted-mean position keeps
# the m/z sharp (0.2-0.7 ppm reproducibility across slices, 3-5 ppm without the kernel).
# Measured (2026-09-16/18): no m/z kernel -24% IDs; sigma 3 vs 1 +2.6% (E. coli 50 ng) and +6.7% (human 250 pg)
# at 1% FDR, with ~30% smaller files. Cost: ions closer than ~2 sigma (6 bins, ~25 ppm at m/z 500) merge.

"""
    centroid_slice!(out, sc, lp, kmz, h_mz)

Consume `sc.sp_*` (one slice: TOF bin, IM-summed intensity and raw-entry count, sorted by bin) and append the
slice's centroids to `out`.

The slice is sparse (a few thousand occupied bins out of ~400,000), so it is not convolved as one dense array.
Instead it is cut into **runs**: maximal groups of occupied bins in which consecutive bins are at most
`2h_mz + 2` apart, where `h_mz` is the kernel half-width. Two occupied bins further apart than that have kernel
footprints separated by at least two empty bins, so their smoothed profiles cannot touch and they can be
processed independently. Each run is:

1. laid out as a dense buffer covering its first to last bin plus `h_mz` of padding on each side (so the kernel
   tails of the edge bins fit);
2. convolved by **scattering**: every occupied bin adds `kmz[t] * value` to the `2h_mz + 1` positions around it
   (equivalent to convolution, and cheaper when most positions are empty);
3. given its raw-entry counts at the occupied bins only (counts are not smoothed: `min_scans` asks how many raw
   entries fall inside a centroid's footprint);
4. passed to `centroid_dense!`, which looks for apices only between the first and last occupied bin (the
   padding holds nothing but falling Gaussian tails).
"""
function centroid_slice!(out::FrameSlices, sc::SmoothScratch, lp::LevelParams, kmz::Vector{Float64}, h_mz::Int)
    bins = sc.sp_bin; vals = sc.sp_val; cnts = sc.sp_cnt
    n = length(bins); gap_max = 2h_mz + 2
    ntaps = length(kmz)
    i = 1
    @inbounds while i <= n
        # the run is bins[i:j]: extend while the next occupied bin is within gap_max
        j = i
        while j < n && bins[j + 1] - bins[j] <= gap_max; j += 1; end
        # dense buffer: element 1 is bin t0 = first bin - h_mz; length covers the run plus both paddings
        t0 = Int(bins[i]) - h_mz; L = Int(bins[j]) - Int(bins[i]) + 1 + 2h_mz
        if length(sc.dense) < L
            resize!(sc.dense, L); resize!(sc.dcnt, L)
        end
        dense = sc.dense; dcnt = sc.dcnt
        @simd for a in 1:L
            dense[a] = 0.0; dcnt[a] = Int32(0)
        end
        # scatter each occupied bin p across p-h_mz .. p+h_mz (kernel tap t lands at base + t)
        for q in i:j
            p = Int(bins[q]) - t0 + 1; v = vals[q]
            dcnt[p] = cnts[q]
            base = p - h_mz - 1
            @simd for t in 1:ntaps
                dense[base + t] += kmz[t] * v
            end
        end
        centroid_dense!(out, dense, dcnt, L, t0, lp, h_mz + 1, L - h_mz)
        i = j + 1
    end
    out
end

"Variant A (`centroid = :none`): emit the IM-accumulated per-bin sums as peaks, subject to the `min_scans` cull."
function emit_sparse!(out::FrameSlices, sc::SmoothScratch, min_scans::Int)
    @inbounds for q in eachindex(sc.sp_bin)
        v = sc.sp_val[q]
        sc.sp_cnt[q] >= min_scans || continue
        push_peak!(out, Float64(sc.sp_bin[q]), v)
    end
    out
end
