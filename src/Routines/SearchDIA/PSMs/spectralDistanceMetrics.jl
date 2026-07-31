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
struct SpectralScoresMainSearch{T<:AbstractFloat,I<:AbstractFloat} <: SpectralScores{T}
    gof::T
    max_matched_residual::T
    max_unmatched_residual::T
    fitted_manhattan_distance::T
    fitted_hellinger::T
    fitted_frag1_int::I
    fitted_frag2_int::I
    fitted_frag3_int::I
    fitted_frag4_int::I
    fitted_frag5_int::I
    fitted_frag6_int::I
    fitted_frag7_int::I
    fitted_frag8_int::I
    shadow_frag1_int::I
    shadow_frag2_int::I
    shadow_frag3_int::I
    shadow_frag4_int::I
    shadow_frag5_int::I
    shadow_frag6_int::I
    shadow_frag7_int::I
    shadow_frag8_int::I
end
"""
    _accumulate_residuals!(w, r, H)

Fill `r` with the fit residual for every row touched by `H`: `r[row]` becomes the summed fitted
intensity across all columns minus the observed intensity. Both scorers below need this, because
a column's deconvolved share of a peak is `w[col]*H.nzval[i] - r[row]`.
"""
function _accumulate_residuals!(w::Vector{T},
    r::Vector{T},
    H::AbstractSparseDesignMatrix{Ti,T},
   ) where {Ti<:Integer,T}

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
    return nothing
end

function getDistanceMetrics(w::Vector{T},
    r::Vector{T},
    H::AbstractSparseDesignMatrix{Ti,T},
    spectral_scores::Vector{SpectralScoresMainSearch{U,V}}
   ) where {Ti<:Integer,T,U<:AbstractFloat,V<:AbstractFloat}

    _accumulate_residuals!(w, r, H)

    # Single-pass scoring per precursor
    for col in 1:H.n
        # Skip zero-weight columns
        if w[col] <= zero(T)
            spectral_scores[col] = SpectralScoresMainSearch(
                zero(U), zero(U), zero(U), zero(U), zero(U),
                zero(V), zero(V), zero(V), zero(V),
                zero(V), zero(V), zero(V), zero(V),
                zero(V), zero(V), zero(V), zero(V),
                zero(V), zero(V), zero(V), zero(V)
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
        fitted_frag1_int = zero(T)
        fitted_frag2_int = zero(T)
        fitted_frag3_int = zero(T)
        fitted_frag4_int = zero(T)
        fitted_frag5_int = zero(T)
        fitted_frag6_int = zero(T)
        fitted_frag7_int = zero(T)
        fitted_frag8_int = zero(T)
        shadow_frag1_int = zero(T)
        shadow_frag2_int = zero(T)
        shadow_frag3_int = zero(T)
        shadow_frag4_int = zero(T)
        shadow_frag5_int = zero(T)
        shadow_frag6_int = zero(T)
        shadow_frag7_int = zero(T)
        shadow_frag8_int = zero(T)

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
            rank = rank_at(H, i)
            if rank == UInt8(1)
                fitted_frag1_int += fitted_peak
                shadow_frag1_int += x_i
            elseif rank == UInt8(2)
                fitted_frag2_int += fitted_peak
                shadow_frag2_int += x_i
            elseif rank == UInt8(3)
                fitted_frag3_int += fitted_peak
                shadow_frag3_int += x_i
            elseif rank == UInt8(4)
                fitted_frag4_int += fitted_peak
                shadow_frag4_int += x_i
            elseif rank == UInt8(5)
                fitted_frag5_int += fitted_peak
                shadow_frag5_int += x_i
            elseif rank == UInt8(6)
                fitted_frag6_int += fitted_peak
                shadow_frag6_int += x_i
            elseif rank == UInt8(7)
                fitted_frag7_int += fitted_peak
                shadow_frag7_int += x_i
            elseif rank == UInt8(8)
                fitted_frag8_int += fitted_peak
                shadow_frag8_int += x_i
            end

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

        spectral_scores[col] = SpectralScoresMainSearch(
            U(gof),
            U(max_matched_residual),
            U(max_unmatched_residual),
            U(fitted_manhattan_distance),
            U(fitted_hellinger),
            V(fitted_frag1_int),
            V(fitted_frag2_int),
            V(fitted_frag3_int),
            V(fitted_frag4_int),
            V(fitted_frag5_int),
            V(fitted_frag6_int),
            V(fitted_frag7_int),
            V(fitted_frag8_int),
            V(shadow_frag1_int),
            V(shadow_frag2_int),
            V(shadow_frag3_int),
            V(shadow_frag4_int),
            V(shadow_frag5_int),
            V(shadow_frag6_int),
            V(shadow_frag7_int),
            V(shadow_frag8_int)
        )
    end
end

"""
    getFragmentIntensities!(w, r, H, spectral_scores)

Rank-1..8 fitted and deconvolved ("shadow") fragment intensities per column -- the ONLY fields
IntegrateChromatogramsSearch reads back. `MS2MBRChromObject` has no gof, residual, manhattan or
Hellinger field, so `getDistanceMetrics`' five aggregate metrics cannot be consumed there; this
computes the 16 that are and leaves the rest zero.

Dropped relative to getDistanceMetrics, per non-zero: the manhattan and observed-sum accumulators,
the residual magnitude and its matched/unmatched maxima, and the Bhattacharyya `sqrt`. Per column:
five `log2` and one `sqrt`. The residual construction and the fragment accumulation are unchanged,
so the 16 retained values are bit-identical.
"""
function getFragmentIntensities!(w::Vector{T},
    r::Vector{T},
    H::AbstractSparseDesignMatrix{Ti,T},
    spectral_scores::Vector{SpectralScoresMainSearch{U,V}}
   ) where {Ti<:Integer,T,U<:AbstractFloat,V<:AbstractFloat}

    _accumulate_residuals!(w, r, H)

    for col in 1:H.n
        # Skip zero-weight columns
        if w[col] <= zero(T)
            spectral_scores[col] = SpectralScoresMainSearch(
                zero(U), zero(U), zero(U), zero(U), zero(U),
                zero(V), zero(V), zero(V), zero(V),
                zero(V), zero(V), zero(V), zero(V),
                zero(V), zero(V), zero(V), zero(V),
                zero(V), zero(V), zero(V), zero(V)
            )
            continue
        end

        fitted_frag1_int = zero(T)
        fitted_frag2_int = zero(T)
        fitted_frag3_int = zero(T)
        fitted_frag4_int = zero(T)
        fitted_frag5_int = zero(T)
        fitted_frag6_int = zero(T)
        fitted_frag7_int = zero(T)
        fitted_frag8_int = zero(T)
        shadow_frag1_int = zero(T)
        shadow_frag2_int = zero(T)
        shadow_frag3_int = zero(T)
        shadow_frag4_int = zero(T)
        shadow_frag5_int = zero(T)
        shadow_frag6_int = zero(T)
        shadow_frag7_int = zero(T)
        shadow_frag8_int = zero(T)

        @inbounds @fastmath for i in H.colptr[col]:(H.colptr[col+1]-1)
            fitted_peak = w[col]*H.nzval[i]
            shadow_peak = fitted_peak - r[H.rowval[i]]
            x_i = max(shadow_peak, zero(T))  # clamp negative shadows
            rank = rank_at(H, i)
            if rank == UInt8(1)
                fitted_frag1_int += fitted_peak
                shadow_frag1_int += x_i
            elseif rank == UInt8(2)
                fitted_frag2_int += fitted_peak
                shadow_frag2_int += x_i
            elseif rank == UInt8(3)
                fitted_frag3_int += fitted_peak
                shadow_frag3_int += x_i
            elseif rank == UInt8(4)
                fitted_frag4_int += fitted_peak
                shadow_frag4_int += x_i
            elseif rank == UInt8(5)
                fitted_frag5_int += fitted_peak
                shadow_frag5_int += x_i
            elseif rank == UInt8(6)
                fitted_frag6_int += fitted_peak
                shadow_frag6_int += x_i
            elseif rank == UInt8(7)
                fitted_frag7_int += fitted_peak
                shadow_frag7_int += x_i
            elseif rank == UInt8(8)
                fitted_frag8_int += fitted_peak
                shadow_frag8_int += x_i
            end
        end

        spectral_scores[col] = SpectralScoresMainSearch(
            zero(U), zero(U), zero(U), zero(U), zero(U),
            V(fitted_frag1_int),
            V(fitted_frag2_int),
            V(fitted_frag3_int),
            V(fitted_frag4_int),
            V(fitted_frag5_int),
            V(fitted_frag6_int),
            V(fitted_frag7_int),
            V(fitted_frag8_int),
            V(shadow_frag1_int),
            V(shadow_frag2_int),
            V(shadow_frag3_int),
            V(shadow_frag4_int),
            V(shadow_frag5_int),
            V(shadow_frag6_int),
            V(shadow_frag7_int),
            V(shadow_frag8_int)
        )
    end
end
