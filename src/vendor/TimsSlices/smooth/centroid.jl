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

# Step 3 of the pipeline (see window.jl): peak picking on one smoothed m/z profile. Classic "smooth, find local
# maxima, take the centre of mass" centroiding, applied to the profile the IM and m/z Gaussians produced.
#
# Logic identical to the prototype's centroid_dense!, with two mechanical changes: the left footprint edge is
# tracked during the forward scan (`rise`: where the current strictly increasing positive run began) instead of
# walked back from each apex, and `jlo:jhi` restricts where apices are looked for. The m/z stage passes the
# span of the run's input bins: the smoothed profile is monotone in the padding (a sum of Gaussian tails) and
# cannot hold an apex there, while footprints still walk into the padding.

"""
    centroid_dense!(out, sm, cnt, L, t0, lp, jlo = 2, jhi = L - 1)

Append the centroids of the dense smoothed profile `sm[1:L]` to `out`. Element `k` of `sm` is TOF bin
`t0 + k - 1`; `cnt[k]` is the number of raw (un-smoothed) entries that landed in that bin.

For every bin `j` in `jlo:jhi`:

1. **Apex.** `j` is an apex when `sm[j] > 0`, `sm[j] > sm[j-1]` and `sm[j] >= sm[j+1]`. The mixed `>` / `>=` makes
   a flat top of equal values yield one apex, its leftmost bin, rather than none or several.
2. **Footprint.** The peak's extent is walked outwards from the apex while the profile keeps falling strictly and
   stays positive, at most `max_half` bins each way. It stops at the first bin that rises again (a valley between
   two peaks: each peak keeps its own side) or at zero (the edge of the signal):

   ```
               apex j
                 *
              *     *
           *           *   *        <- the walk stops here: the profile rises again (next peak)
        *                *
     lo ------------------ hi
   ```

   The left edge `lo` is not walked back from each apex: `rise` records where the current strictly increasing
   positive stretch began during the forward scan, which is the same bin (capped at `j - max_half`).
3. **Intensity** = the sum of `sm` over `lo:hi`, i.e. the peak area in the smoothed profile.
4. **Position.** `:wmean` (default): the intensity-weighted mean bin over the footprint, `Σ sm[k]·bin_k / Σ sm[k]`.
   `:gauss`: a parabola through `log sm` at `j-1, j, j+1` (a Gaussian through three points); its vertex offset
   `0.5 (ln_m - ln_p) / (ln_m - 2 ln_0 + ln_p)`, clamped to +/-1 bin, is added to `j`.
5. **Cull.** A centroid whose footprint holds fewer than `min_scans` raw entries (`Σ cnt`) is dropped: it was
   seen too few times across the IM scans and m/z bins that fed it to be more than noise.

Positions are fractional bins; `quantize!` rounds them to the stored fixed-point grid later.
"""
function centroid_dense!(out::FrameSlices, sm::Vector{Float64}, cnt::Vector{Int32}, L::Int, t0::Int, lp::LevelParams,
                         jlo::Int = 2, jhi::Int = L - 1)
    max_half = lp.max_half; gauss = lp.centroid == :gauss; min_scans = lp.min_scans
    jlo = max(2, jlo); jhi = min(L - 1, jhi)
    jlo > jhi && return out
    # `rise` = first bin of the strictly increasing positive stretch that ends at the current bin. Seed it for jlo
    # by walking back at most max_half steps (once per run, not per apex).
    rise = jlo
    @inbounds while rise > 1 && rise > jlo - max_half && sm[rise-1] < sm[rise] && sm[rise-1] > 0; rise -= 1; end
    @inbounds for j in jlo:jhi
        v = sm[j]; vm = sm[j-1]
        # the stretch continues only while the previous bin is positive and lower; otherwise it restarts at j
        j > jlo && !((vm < v) & (vm > 0.0)) && (rise = j)
        # apex test (step 1)
        ((v > 0.0) & (v > vm) & (v >= sm[j+1])) || continue
        # footprint (step 2): left edge from `rise`, right edge walked down
        lo = max(rise, j - max_half)
        hi = j; while hi < L && hi < j + max_half && sm[hi+1] < sm[hi] && sm[hi+1] > 0; hi += 1; end
        # area, weighted bin sum and raw-entry count over the footprint (steps 3-5)
        s = 0.0; ws = 0.0; n_entries = Int32(0)
        for k in lo:hi; s += sm[k]; ws += sm[k] * (t0 + k - 1); n_entries += cnt[k]; end
        n_entries >= min_scans || continue
        if gauss
            # log-parabola vertex from the apex and its two neighbours; den < 0 for any strict local maximum
            l0 = log(v); lm = log(max(vm, 1e-9)); lp_ = log(max(sm[j+1], 1e-9))
            den = lm - 2l0 + lp_
            off = den < 0 ? clamp(0.5 * (lm - lp_) / den, -1.0, 1.0) : 0.0
            tpos = Float64(t0 + j - 1) + off
        else
            tpos = ws / s                                   # intensity-weighted mean bin
        end
        push_peak!(out, tpos, s)
    end
    out
end
