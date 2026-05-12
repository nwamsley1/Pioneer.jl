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
    scribe::T
    percent_theoretical_ignored::T
    fitted_hellinger::T
    spectral_contrast::T
    fitted_spectral_contrast::T
    matched_ratio::T
end

@inline function _mainsearch_spectral_zero(::Type{U}) where {U<:AbstractFloat}
    return SpectralScoresMainSearch(
        zero(U), zero(U), zero(U), zero(U), zero(U),
        zero(U), zero(U), zero(U), zero(U), zero(U)
    )
end

@inline function _finite_or_zero(x::T) where {T<:AbstractFloat}
    return isfinite(x) ? x : zero(T)
end

function _compute_mainsearch_spectral_metrics(
    w::Vector{T},
    H::AbstractSparseDesignMatrix{Ti,T},
    r::Vector{T},
    col::Integer,
    included_indices::AbstractVector{Int},
) where {Ti<:Integer,T<:AbstractFloat}
    fitted_total = zero(T)
    shadow_total = zero(T)
    @inbounds for i in included_indices
        fitted_peak = max(w[col] * H.nzval[i], zero(T))
        shadow_peak = max(fitted_peak - r[H.rowval[i]], zero(T))
        fitted_total += fitted_peak
        shadow_total += shadow_peak
    end
    fitted_total = max(fitted_total, eps(T))
    shadow_total = max(shadow_total, eps(T))

    x_sum = zero(T)
    manhattan_distance = zero(T)
    max_matched_residual = zero(T)
    max_unmatched_residual = zero(T)
    sum_of_residuals = zero(T)
    sum_of_fitted_peaks_matched = zero(T)
    sum_of_fitted_peaks_unmatched = zero(T)
    sum_of_fitted_peaks_squared = zero(T)
    fitted_dotp = zero(T)
    fitted_dotp_norm1 = zero(T)
    matched_signal = zero(T)
    unmatched_signal = zero(T)
    dot_product = zero(T)
    h2_sum = zero(T)
    x2_sum = zero(T)
    bc_sum = zero(T)
    sum_fitted = zero(T)
    sum_shadow = zero(T)
    scribe_sqrt_fitted_sum = zero(T)
    scribe_sqrt_shadow_sum = zero(T)
    worst_val = typemin(T)
    worst_pos = 0
    worst_pred_signal = zero(T)
    num_matching_peaks = 0
    n = 0

    @inbounds @fastmath for (local_pos, i) in enumerate(included_indices)
        pred_peak = max(H.nzval[i], zero(T))
        observed_peak = max(H.x[i], zero(T))
        fitted_peak = max(w[col] * H.nzval[i], zero(T))
        shadow_peak = max(fitted_peak - r[H.rowval[i]], zero(T))
        r_abs = abs(r[H.rowval[i]])

        x_sum += observed_peak
        manhattan_distance += abs(fitted_peak - observed_peak)
        sum_of_residuals += r_abs
        sum_fitted += fitted_peak
        sum_shadow += shadow_peak
        bc_sum += sqrt(fitted_peak * shadow_peak)
        scribe_sqrt_fitted_sum += sqrt(fitted_peak)
        scribe_sqrt_shadow_sum += sqrt(shadow_peak)

        dot_product += pred_peak * observed_peak
        h2_sum += pred_peak^2
        x2_sum += observed_peak^2
        sum_of_fitted_peaks_squared += fitted_peak^2
        observed_peak > zero(T) && (num_matching_peaks += 1)

        if matched_at(H, i)
            matched_signal += pred_peak
            fitted_dotp += shadow_peak * fitted_peak
            fitted_dotp_norm1 += shadow_peak^2
            sum_of_fitted_peaks_matched += fitted_peak
            max_matched_residual = max(max_matched_residual, r_abs)
        else
            unmatched_signal += pred_peak
            sum_of_fitted_peaks_unmatched += fitted_peak
            max_unmatched_residual = max(max_unmatched_residual, r_abs)
        end

        diff = (shadow_peak / shadow_total) - (fitted_peak / fitted_total)
        if diff > worst_val && shadow_peak > zero(T)
            worst_val = diff
            worst_pos = local_pos
            worst_pred_signal = pred_peak
        end
        n += 1
    end

    sum_of_fitted_peaks = sum_of_fitted_peaks_matched + sum_of_fitted_peaks_unmatched
    gof = sum_of_fitted_peaks > zero(T) ?
        -log2(sum_of_residuals / sum_of_fitted_peaks + T(1e-10)) : zero(T)
    max_matched_residual = sum_of_fitted_peaks_matched > zero(T) ?
        -log2(max_matched_residual / sum_of_fitted_peaks_matched + T(1e-10)) : zero(T)
    max_unmatched_residual = sum_of_fitted_peaks > zero(T) ?
        -log2(max_unmatched_residual / sum_of_fitted_peaks + T(1e-10)) : zero(T)
    fitted_manhattan_distance = x_sum > zero(T) ?
        -log2(manhattan_distance / x_sum + T(1e-10)) : zero(T)

    hellinger_denom = sqrt(sum_fitted * sum_shadow)
    hellinger_sq = hellinger_denom > zero(T) ? one(T) - bc_sum / hellinger_denom : one(T)
    fitted_hellinger = -log2(max(hellinger_sq, T(1e-10)))

    scribe = zero(T)
    if n > 0 && scribe_sqrt_fitted_sum > zero(T) && scribe_sqrt_shadow_sum > zero(T)
        scribe_distance = zero(T)
        @inbounds @fastmath for i in included_indices
            fitted_peak = max(w[col] * H.nzval[i], zero(T))
            shadow_peak = max(fitted_peak - r[H.rowval[i]], zero(T))
            fitted_norm = sqrt(fitted_peak) / scribe_sqrt_fitted_sum
            shadow_norm = sqrt(shadow_peak) / scribe_sqrt_shadow_sum
            scribe_distance += (fitted_norm - shadow_norm)^2
        end
        scribe = -log2(max(scribe_distance / T(n), T(1e-10)))
    end

    spectral_denom = sqrt(h2_sum * x2_sum)
    spectral_contrast = spectral_denom > zero(T) ? dot_product / spectral_denom : zero(T)
    fitted_spectral_denom = sqrt(fitted_dotp_norm1 * sum_of_fitted_peaks_squared)
    fitted_spectral_contrast = fitted_spectral_denom > zero(T) ?
        fitted_dotp / fitted_spectral_denom : zero(T)
    matched_ratio = matched_signal > zero(T) && unmatched_signal > zero(T) ?
        log2(matched_signal / unmatched_signal) : zero(T)

    return (
        gof = _finite_or_zero(gof),
        max_matched_residual = _finite_or_zero(max_matched_residual),
        max_unmatched_residual = _finite_or_zero(max_unmatched_residual),
        fitted_manhattan_distance = _finite_or_zero(fitted_manhattan_distance),
        scribe = _finite_or_zero(scribe),
        fitted_hellinger = _finite_or_zero(fitted_hellinger),
        spectral_contrast = _finite_or_zero(spectral_contrast),
        fitted_spectral_contrast = _finite_or_zero(fitted_spectral_contrast),
        matched_ratio = _finite_or_zero(matched_ratio),
        worst_pos = worst_pos,
        worst_pred_signal = _finite_or_zero(worst_pred_signal),
        num_matching_peaks = num_matching_peaks,
    )
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

    included = Int[]

    # Single-pass scoring per precursor
    for col in 1:H.n
        # Skip zero-weight columns
        if w[col] <= zero(T)
            spectral_scores[col] = _mainsearch_spectral_zero(U)
            continue
        end

        empty!(included)
        @inbounds for i in H.colptr[col]:(H.colptr[col+1]-1)
            push!(included, i)
        end
        isempty(included) && begin
            spectral_scores[col] = _mainsearch_spectral_zero(U)
            continue
        end
        total_pred_signal = zero(T)
        @inbounds for i in included
            total_pred_signal += max(H.nzval[i], zero(T))
        end
        total_pred_signal = max(total_pred_signal, eps(T))

        best = _compute_mainsearch_spectral_metrics(w, H, r, col, included)
        ignored_signal = zero(T)
        while best.num_matching_peaks > 3 && best.worst_pos > 0 && length(included) > 1
            worst_pos = best.worst_pos
            worst_signal = best.worst_pred_signal
            deleteat!(included, worst_pos)
            candidate = _compute_mainsearch_spectral_metrics(w, H, r, col, included)
            if candidate.scribe > best.scribe * T(1.25)
                ignored_signal += worst_signal
                best = candidate
            else
                break
            end
        end
        percent_theoretical_ignored = min(ignored_signal / total_pred_signal, one(T))

        spectral_scores[col] = SpectralScoresMainSearch(
            U(best.gof),
            U(best.max_matched_residual),
            U(best.max_unmatched_residual),
            U(best.fitted_manhattan_distance),
            U(best.scribe),
            U(percent_theoretical_ignored),
            U(best.fitted_hellinger),
            U(best.spectral_contrast),
            U(best.fitted_spectral_contrast),
            U(best.matched_ratio)
        )
    end
end
