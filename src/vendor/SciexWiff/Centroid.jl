# Centroiding decoded profile bins (notes/format.md §9).
#
# Stored bins are the non-zero TOF bins of a sparse profile on a grid of 8-unit steps; a bin that is
# not stored is zero. A real peak is narrow and about constant in steps (σ ≈ 1 step at every m/z,
# FWHM ≈ 2.3 steps), so peaks are found with a matched filter: the profile is smoothed with a Gaussian
# of `sigma` steps (gaps count as zeros), every local maximum of the smoothed profile is a peak, and the
# peak is reported from the raw bins within `half` steps of that maximum (not past the valley to a
# neighbouring maximum): intensity-weighted mean position, and area.

const BIN_STEP = UInt32(8)

"""
    CentroidParams(; sigma = 2.0, half = 3, min_bins = 2, min_apex = 0.0, min_intensity = 0.0, intensity = :area)

- `sigma`: Gaussian smoothing width in steps (one step = 8 raw TOF units).
- `half`: a peak uses the raw bins within `half` steps of its smoothed maximum.
- `min_bins`: drop peaks with fewer non-zero raw bins in their extent (the vendor drops single-bin peaks: 2).
- `min_apex`: drop peaks whose smoothed maximum is lower (raw intensity units).
- `min_intensity`: drop peaks whose reported intensity is lower.
- `position`: `:wmean` (raw bins in the extent) or `:parabola` (vertex through the smoothed maximum
  and its two neighbours).
- `intensity`: `:area` = 100 · Σ I·Δm/z (the scale the vendor centroider reports) or `:sum` = Σ I.
"""
Base.@kwdef struct CentroidParams
    sigma::Float64 = 2.0
    half::Int = 3
    min_apex::Float64 = 0.0
    min_bins::Int = 2            # drop peaks with fewer non-zero raw bins in their extent
    min_intensity::Float64 = 0.0
    intensity::Symbol = :area
    position::Symbol = :wmean    # :wmean of raw bins, or :parabola through the smoothed maximum
end

"Per-thread scratch and output for `centroid!`."
mutable struct CentroidBuffer
    n::Int
    mz::Vector{Float32}
    intensity::Vector{Float32}
    nbins::Vector{Int32}         # non-zero raw bins behind each centroid
    pos::Vector{Float64}         # centroid position in raw TOF units (m/z = bin_to_mz(pos, a, b))
    dense::Vector{Float32}       # raw profile of the current cluster on the step grid (zeros filled in)
    smooth::Vector{Float32}
    kernel::Vector{Float32}
    kernel_sigma::Float64
end
CentroidBuffer() = CentroidBuffer(0, Float32[], Float32[], Int32[], Float64[], Float32[], Float32[], Float32[], NaN)

function _kernel!(cb::CentroidBuffer, sigma::Float64)
    cb.kernel_sigma == sigma && return cb.kernel
    r = max(1, ceil(Int, 3sigma))
    k = Float32[exp(-0.5 * (j / sigma)^2) for j in -r:r]
    k ./= sum(k)
    cb.kernel = k; cb.kernel_sigma = sigma
    k
end

@inline function _push!(cb::CentroidBuffer, pos::Float64, mz::Float64, val::Float64, nb::Int)
    cb.n += 1
    if cb.n > length(cb.mz)
        resize!(cb.mz, 2cb.n); resize!(cb.intensity, 2cb.n); resize!(cb.nbins, 2cb.n); resize!(cb.pos, 2cb.n)
    end
    @inbounds cb.mz[cb.n] = Float32(mz); cb.intensity[cb.n] = Float32(val); cb.nbins[cb.n] = Int32(nb); cb.pos[cb.n] = pos
    nothing
end

"Centroid one cluster: stored bins lo:hi of `sb`, gaps between them shorter than the kernel reach."
function _centroid_cluster!(cb::CentroidBuffer, sb::ScanBuffer, lo::Int, hi::Int, p::CentroidParams)
    hi - lo + 1 < p.min_bins && return                  # too few stored bins for any peak to qualify
    k = _kernel!(cb, p.sigma); r = length(k) ÷ 2
    B = sb.bin; I = sb.intensity
    s0 = @inbounds Int(B[lo] ÷ BIN_STEP) - r          # grid index of dense[1]
    L = @inbounds Int(B[hi] ÷ BIN_STEP) - s0 + r + 1
    length(cb.dense) < L && (resize!(cb.dense, L); resize!(cb.smooth, L))
    D = cb.dense; S = cb.smooth
    @inbounds fill!(view(D, 1:L), 0f0); @inbounds fill!(view(S, 1:L), 0f0)
    # scatter: each stored bin adds its kernel (the profile is mostly zeros, so this beats a dense convolution)
    @inbounds for j in lo:hi
        q = Int(B[j] ÷ BIN_STEP) - s0 + 1; v = Float32(I[j])
        D[q] = v
        for t in -r:r
            S[q+t] += k[t+r+1] * v
        end
    end
    a, b = sb.cal_a, sb.cal_b
    base = s0 * Float64(BIN_STEP)                       # raw bin of dense[1]
    @inbounds for i in 2:L-1
        (S[i] > S[i-1] && S[i] >= S[i+1] && S[i] >= p.min_apex) || continue
        # extent: up to `half` steps each side, stopping at a valley (the smoothed profile rising again)
        l = i
        while l > 1 && i - l < p.half && S[l-1] <= S[l]
            l -= 1
        end
        u = i
        while u < L && u - i < p.half && S[u+1] <= S[u]
            u += 1
        end
        sw = 0.0; swx = 0.0; area = 0.0; nb = 0
        for q in l:u
            w = Float64(D[q]); w == 0 && continue
            x = base + (q - 1) * Float64(BIN_STEP)
            sw += w; swx += w * x; nb += 1
            p.intensity === :area && (area += w * (bin_to_mz(x + BIN_STEP, a, b) - bin_to_mz(x, a, b)))
        end
        (sw > 0 && nb >= p.min_bins) || continue
        val = p.intensity === :area ? 100area : sw
        val < p.min_intensity && continue
        x = if p.position === :parabola
            den = S[i-1] - 2S[i] + S[i+1]
            base + (i - 1) * Float64(BIN_STEP) + (den < 0 ? 4 * (S[i-1] - S[i+1]) / den : 0.0)
        else
            swx / sw
        end
        _push!(cb, x, bin_to_mz(x, a, b), val, nb)
    end
end

"""
    centroid!(out, sb, params = CentroidParams()) -> out

Centroid the decoded scan in `sb`; results are `out.mz[1:out.n]`, `out.intensity[1:out.n]`, ascending m/z.
"""
function centroid!(cb::CentroidBuffer, sb::ScanBuffer, p::CentroidParams = CentroidParams())
    cb.n = 0
    sb.n == 0 && return cb
    reach = 2 * (length(_kernel!(cb, p.sigma)) ÷ 2) * BIN_STEP   # bins farther apart than this cannot interact
    lo = 1
    @inbounds for j in 2:sb.n
        if sb.bin[j] - sb.bin[j-1] > reach
            _centroid_cluster!(cb, sb, lo, j - 1, p); lo = j
        end
    end
    _centroid_cluster!(cb, sb, lo, sb.n, p)
    cb
end
