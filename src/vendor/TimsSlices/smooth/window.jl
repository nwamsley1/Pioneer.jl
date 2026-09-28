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

# The per-window pipeline: slices every `stride` scans of window scans s0:s1-1 of one decoded frame.
#
# A diaPASEF frame is one TIMS ramp: ~930 IM scans, each a TOF spectrum. Its quad windows each cover a range of
# those scans. For every window the converter produces "slices", one every `stride` scans (0.0065 1/K0 by default,
# 8 scans on the timsTOF Ultra ramps; see STRIDE_K0), each a centroided
# spectrum that Pioneer searches as one scan. Per slice:
#
#   1. IM accumulation (im.jl)    Gaussian-weighted sum of the raw scans around the slice centre, per TOF bin
#                                 (sigma 0.004325 1/K0, 5 scans on the Ultra ramps): the ion's signal over its
#                                 mobility peak, in one spectrum.
#   2. m/z smoothing (mz.jl)      Gaussian along the TOF axis (sigma 0.3125 ns = 2.5 bins on the Ultra): merges
#                                 the +/-1-2 bin jitter of an ion across scans into one hump.
#   3. peak picking (centroid.jl) local maxima, footprint walk, weighted-mean position, footprint-sum intensity,
#                                 `min_scans` cull.
#   4. peak cap (`cap_slice!`)    keep the `max_peaks` most intense centroids (1,500 for MS2 by default; MS1
#                                 uncapped).
#
# Steps 1 and 2 together are a separable 2D Gaussian smoothing of the (IM scan x TOF bin) map, evaluated only at
# the slice centres.

"Per-level constants derived from the parameters once per conversion."
struct LevelSetup
    lp::LevelParams
    kim::Vector{Float64}
    kmz::Vector{Float64}
    h_im::Int
    h_mz::Int
end
function LevelSetup(lp::LevelParams)
    kim = im_kernel(lp); kmz = mz_kernel(lp)
    LevelSetup(lp, kim, kmz, length(kim) ÷ 2, length(kmz) ÷ 2)
end

"""
    smooth_window!(out, sc, buf, s0, s1, widx, ls)

Append the slices of window `widx` (scans `s0:s1-1`) of the frame in `buf` to `out`. A slice with no raw
entries in reach, or no surviving centroids, is not recorded.
"""
function smooth_window!(out::FrameSlices, sc::SmoothScratch, buf::FrameBuffer, s0::Integer, s1::Integer, widx::Integer, ls::LevelSetup)
    lp = ls.lp
    for c in s0:lp.stride:s1-1
        im_accumulate!(sc, buf, s0, s1, c, ls.kim, ls.h_im)
        isempty(sc.sp_bin) && continue
        if lp.centroid == :none
            emit_sparse!(out, sc, lp.min_scans)
        else
            centroid_slice!(out, sc, lp, ls.kmz, ls.h_mz)
        end
        lp.max_peaks > 0 && cap_slice!(out, lp.max_peaks, sc.tmp)
        end_slice!(out, c, widx)
    end
    out
end

"""
    cap_slice!(fs, n, tmp)

Keep only the `n` most intense peaks of the slice currently being filled (those after the last `end_slice!`),
preserving their position order. Ties at the threshold are kept (the slice may then exceed `n` slightly).

Finds the n-th largest intensity with a partial sort of a scratch copy (O(m) on average, no full sort), then
compacts the slice in place, keeping every peak at or above it. The peaks stay in m/z order, which the encoder
needs.
"""
function cap_slice!(fs::FrameSlices, n::Int, tmp::Vector{Float64})
    first = Int(fs.ptr[end]); last = length(fs.pos)
    m = last - first + 1
    m <= n && return fs
    length(tmp) < m && resize!(tmp, m)
    copyto!(tmp, 1, fs.val, first, m)
    thr = partialsort!(view(tmp, 1:m), n; rev = true)      # the n-th largest intensity
    w = first
    @inbounds for k in first:last
        if fs.val[k] >= thr
            fs.pos[w] = fs.pos[k]; fs.val[w] = fs.val[k]; w += 1
        end
    end
    resize!(fs.pos, w - 1); resize!(fs.val, w - 1)
    fs
end

"""
    smooth_frame!(out, sc, buf, file, i, setup_ms1, setup_ms2)

All windows of frame row `i` (already decoded into `buf`) -> `out` (reset first).
"""
function smooth_frame!(out::FrameSlices, sc::SmoothScratch, buf::FrameBuffer, f::TdfFile, i::Integer, ls1::LevelSetup, ls2::LevelSetup)
    reset!(out)
    buf.n_peaks == 0 && return out
    ensure_bins!(sc, max(n_bins(f), Int(maximum(view(buf.tof, 1:buf.n_peaks))) + 1))
    ls = is_ms1(f, i) ? ls1 : ls2
    for (widx, w) in enumerate(windows(f, i))
        smooth_window!(out, sc, buf, w.scan_begin, w.scan_end, widx, ls)
    end
    out
end
