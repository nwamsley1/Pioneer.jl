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

abstract type SpectralScores{T<:AbstractFloat} end
struct SpectralScoresMainSearch{T<:AbstractFloat} <: SpectralScores{T}
    gof::T
    max_matched_residual::T
    max_unmatched_residual::T
    fitted_manhattan_distance::T
    fitted_hellinger::T
    spectral_contrast::T
end
function getDistanceMetrics(w::Vector{T},
    r::Vector{T},
    H::AbstractSparseDesignMatrix{Ti,T},
    spectral_scores::Vector{SpectralScoresMainSearch{U}}
   ) where {Ti<:Integer,T,U<:AbstractFloat}

    # Zero residual vector
    @turbo for i in range(1, H.m)
        r[i] = zero(T)
    end

    for n in range(1, H.n_vals)
        if iszero(r[H.rowval[n]])
            r[H.rowval[n]] = -H.x[n]
        end
    end

    for col in range(1, H.n)
        start = H.colptr[col]
        stop = H.colptr[col+1] - 1
        for n in start:stop
            r[H.rowval[n]] += w[col]*H.nzval[n]
        end
    end

    # Single-pass scoring per precursor
    for col in 1:H.n
        # Skip zero-weight columns
        if w[col] <= zero(T)
            spectral_scores[col] = SpectralScoresMainSearch(
                zero(U), zero(U), zero(U), zero(U), zero(U), zero(U)
            )
            continue
        end

        x_sum = zero(T)
        manhattan_distance = zero(T)
        max_matched_residual = zero(T)
        max_unmatched_residual = zero(T)
        sum_of_residuals = zero(T)
        sum_of_fitted_peaks_matched = zero(T)
        sum_of_fitted_peaks_unmatched = zero(T)
        # Hellinger accumulators
        bc_sum = zero(T)         # Bhattacharyya coefficient
        sum_fitted = zero(T)     # sum of fitted_peak (for normalization)
        sum_shadow = zero(T)     # sum of clamped shadow_peak (for normalization)
        dot_product = zero(T)
        pred_norm2 = zero(T)
        shadow_norm2 = zero(T)

        @inbounds @fastmath for i in H.colptr[col]:(H.colptr[col+1]-1)
            x_sum += H.x[i]
            fitted_peak = w[col]*H.nzval[i]
            manhattan_distance += abs(fitted_peak - H.x[i])

            shadow_peak = fitted_peak - r[H.rowval[i]]
            r_abs = abs(r[H.rowval[i]])
            sum_of_residuals += r_abs

            # Hellinger: uses shadow_peak (per-precursor deconvolved observed) vs fitted_peak
            x_i = max(shadow_peak, zero(T))  # clamp negative shadows
            sum_fitted += fitted_peak
            sum_shadow += x_i
            bc_sum += sqrt(fitted_peak * x_i)
            pred_peak = max(fitted_peak, zero(T))
            dot_product += pred_peak * x_i
            pred_norm2 += pred_peak * pred_peak
            shadow_norm2 += x_i * x_i

            if matched_at(H, i)
                sum_of_fitted_peaks_matched += fitted_peak
                if r_abs > max_matched_residual
                    max_matched_residual = r_abs
                end
            else
                sum_of_fitted_peaks_unmatched += fitted_peak
                if r_abs > max_unmatched_residual
                    max_unmatched_residual = r_abs
                end
            end
        end

        sum_of_fitted_peaks = sum_of_fitted_peaks_matched + sum_of_fitted_peaks_unmatched
        # Floor numerator with +1e-10 to keep log2 finite when residuals
        # collapse to zero on perfect/near-perfect single-fragment matches.
        gof = sum_of_fitted_peaks > 0 ? -log2(sum_of_residuals/sum_of_fitted_peaks + T(1e-10)) : zero(T)
        max_matched_residual = sum_of_fitted_peaks_matched > 0 ? -log2(max_matched_residual/sum_of_fitted_peaks_matched + T(1e-10)) : zero(T)
        max_unmatched_residual = sum_of_fitted_peaks > 0 ? -log2(max_unmatched_residual/sum_of_fitted_peaks + 1e-10) : zero(T)
        fitted_manhattan_distance = x_sum > 0 ? -log2(manhattan_distance/x_sum + 1e-10) : zero(T)

        # Hellinger distance: H² = 1 - BC/sqrt(Σμ · Σx)
        hellinger_denom = sqrt(sum_fitted * sum_shadow)
        hellinger_sq = hellinger_denom > 0 ? one(T) - bc_sum / hellinger_denom : one(T)
        fitted_hellinger = -log2(max(hellinger_sq, T(1e-10)))
        spectral_denom = sqrt(pred_norm2 * shadow_norm2)
        spectral_contrast = spectral_denom > zero(T) ? dot_product / spectral_denom : zero(T)

        spectral_scores[col] = SpectralScoresMainSearch(
            U(gof),
            U(max_matched_residual),
            U(max_unmatched_residual),
            U(fitted_manhattan_distance),
            U(fitted_hellinger),
            U(spectral_contrast)
        )
    end
end
