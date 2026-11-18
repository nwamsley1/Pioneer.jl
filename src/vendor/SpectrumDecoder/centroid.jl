# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Centroiding of profile points, for scans that store profile data but no peak list (packet type 19, ion-trap MS1).
# Local maximum, then the Gaussian apex from a log-parabola least-squares fit weighted by intensity squared
# (Caruana; Guo, IEEE Signal Processing Magazine 2011).

"""One centroid: apex m/z, fitted height, Gaussian area, FWHM, weighted RMS of log residuals (NaN if the fit fell back
to the apex point), number of profile points used."""
struct Centroid
    mz::Float64
    height::Float32
    area::Float32
    fwhm::Float32
    resid::Float32
    npts::Int16
end

# weighted least squares ln y = c0 + c1 u + c2 u^2 over points a:b (u = (mz - x0)/d), weights y^2 (Guo)
function _fit_peak(mz, y, l::Int, r::Int, i::Int)
    h = 0.5f0 * y[i]
    a = i; b = i
    while a > l && y[a-1] >= h; a -= 1; end
    while b < r && y[b+1] >= h; b += 1; end
    if b - a < 2                                     # need >= 3 points: widen to the apex neighbours
        a = max(l, i - 1); b = min(r, i + 1)
        b - a < 2 && return Centroid(Float64(mz[i]), y[i], NaN32, NaN32, NaN32, Int16(b - a + 1))
    end
    x0 = Float64(mz[i]); d = (Float64(mz[b]) - Float64(mz[a])) / (b - a)
    S0 = S1 = S2 = S3 = S4 = T0 = T1 = T2 = L2 = 0.0
    @inbounds for j in a:b
        yj = Float64(y[j]); yj > 0 || continue
        u = (Float64(mz[j]) - x0) / d; w = yj * yj; ly = Float64(log(Float32(yj))); u2 = u * u
        S0 += w; S1 += w * u; S2 += w * u2; S3 += w * u2 * u; S4 += w * u2 * u2
        T0 += w * ly; T1 += w * u * ly; T2 += w * u2 * ly; L2 += w * ly * ly
    end
    # 3x3 symmetric solve by Cramer's rule: [S0 S1 S2; S1 S2 S3; S2 S3 S4] c = T
    D = S0 * (S2 * S4 - S3 * S3) - S1 * (S1 * S4 - S3 * S2) + S2 * (S1 * S3 - S2 * S2)
    fallback = Centroid(x0, y[i], NaN32, NaN32, NaN32, Int16(b - a + 1))
    (isfinite(D) && abs(D) > 1e-300) || return fallback
    c0 = (T0 * (S2 * S4 - S3 * S3) - S1 * (T1 * S4 - S3 * T2) + S2 * (T1 * S3 - S2 * T2)) / D
    c1 = (S0 * (T1 * S4 - S3 * T2) - T0 * (S1 * S4 - S3 * S2) + S2 * (S1 * T2 - T1 * S2)) / D
    c2 = (S0 * (S2 * T2 - T1 * S3) - S1 * (S1 * T2 - T1 * S2) + T0 * (S1 * S3 - S2 * S2)) / D
    c2 < 0 || return fallback                         # not peak-shaped
    mu = -c1 / (2c2)
    abs(mu) <= (b - a) / 2 + 0.5 || return fallback   # vertex outside the fitted span
    sigma = sqrt(-1 / (2c2)) * d
    height = exp(c0 - c1 * c1 / (4c2))
    # weighted residual sum from the accumulated sums (no second pass): Σw(ly - c·φ)² = L2 - 2c·T + cᵀSc
    rss = L2 - 2 * (c0 * T0 + c1 * T1 + c2 * T2) +
          (c0 * c0 * S0 + 2c0 * c1 * S1 + (c1 * c1 + 2c0 * c2) * S2 + 2c1 * c2 * S3 + c2 * c2 * S4)
    rss = max(rss, 0.0)
    Centroid(x0 + mu * d, Float32(height), Float32(height * sigma * sqrt(2pi)), Float32(2.35482 * sigma),
             Float32(sqrt(rss / S0)), Int16(b - a + 1))
end

"""
    centroid!(out, mz, y; threshold = nothing, snr = 3.0, max_half = 8, gap_factor = 2.5) -> out

Centroids of profile points (`mz` ascending). A point is a candidate apex when it is above the threshold (default
`snr` × median of the positive intensities), not lower than its left neighbour and higher than its right one. The peak
extends at most `max_half` points each way while intensity keeps falling and points stay within `gap_factor` × the
local spacing (so separate runs of a reduced profile are never joined).
"""
function centroid!(out::Vector{Centroid}, mz::AbstractVector, y::AbstractVector;
                   threshold = nothing, snr::Real = 3.0, max_half::Int = 8, gap_factor::Real = 2.5)
    empty!(out); n = length(y); n < 3 && return out
    thr = threshold === nothing ? _median_positive(y) * snr : threshold
    i = 2
    @inbounds while i < n
        yi = y[i]
        if yi > thr && yi >= y[i-1] && yi > y[i+1]
            dref = min(mz[i] - mz[i-1], mz[i+1] - mz[i])
            l = i
            # the apex may share its top with its left neighbour (true apex midway between two points): step past it
            while l > 1 && i - l < max_half && 0 < y[l-1] && (y[l-1] < y[l] || (l == i && y[l-1] == y[l])) &&
                  mz[l] - mz[l-1] < gap_factor * dref
                l -= 1
            end
            r = i
            while r < n && r - i < max_half && 0 < y[r+1] < y[r] && mz[r+1] - mz[r] < gap_factor * dref; r += 1; end
            push!(out, _fit_peak(mz, y, l, r, i))
            i = r + 1
        else
            i += 1
        end
    end
    out
end

_median_positive(y) = _median_positive!(Float32[], y)

function _median_positive!(tmp::Vector{Float32}, y; stride::Int = 1)
    empty!(tmp); for k in 1:stride:length(y); x = y[k]; x > 0 && push!(tmp, x); end
    isempty(tmp) ? 0f0 : Float32(partialsort!(tmp, (length(tmp) + 1) ÷ 2))
end
