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

abstract type ScoredPSM{H,L<:AbstractFloat} <: PSM end

"""
    psm_logfac(N) -> Float64

Stirling-like log-factorial approximation. Returns 0 for N ≤ 0 (the closed
form has a 0·log(0) = NaN at N = 0). Used by `psm_getPoisson`.
"""
@inline function psm_logfac(N::Integer)
    N <= 0 && return 0.0
    Nf = float(N)
    Nf*log(Nf) - Nf + log(Nf*(1 + 4*Nf*(1 + 2*Nf)))/6 + log(π)/2
end

"""
    psm_getPoisson(lam, observed) -> T

Log-likelihood-style Poisson score `log(λ^k · e^-λ) - logfac(k)`.

Floors `λ` at `1e-20` to keep `log(0)` from producing `-Inf` when the
expected match count is degenerate, and short-circuits `logfac(0)` to
avoid `0·log(0) = NaN`. Always finite for any `λ ≥ 0`, `observed ≥ 0`.
"""
@inline function psm_getPoisson(lam::T, observed::Integer) where {T<:AbstractFloat}
    lam_safe = max(lam, T(1e-20))
    log((lam_safe^observed)*exp(-lam_safe)) - psm_logfac(observed)
end


struct MainSearchScoredPSM{H,L<:AbstractFloat} <: ScoredPSM{H,L}
    #Ion Count Statistics
    longest_y::UInt8
    b_count::UInt8
    y_count::UInt8
    total_ions::UInt8       # M0 (monoisotopic) fragment count
    total_ions_iso::UInt8   # Isotope (M1+) fragment count

    #Basic Metrics
    poisson::L
    log2_intensity_explained::L
    error::L

    #Spectral Similarity
    gof::L
    max_matched_residual::L
    max_unmatched_residual::L
    fitted_manhattan_distance::L
    scribe::L
    percent_theoretical_ignored::L
    weight::H

    fitted_hellinger::L

    #Non-scores/Labels
    precursor_idx::UInt32
    ms_file_idx::UInt32
    cycle_idx::UInt32
    scan_idx::UInt32
end
function growScoredPSMs!(scored_psms::Vector{MainSearchScoredPSM{H,L}}, block_size::Int64) where {L,H<:AbstractFloat}
    scored_psms = append!(scored_psms, Vector{MainSearchScoredPSM{H,L}}(undef, block_size))
end
function Score!(scored_psms::Vector{MainSearchScoredPSM{H, L}},
                unscored_PSMs::Vector{ComplexUnscoredPSM{H}},
                spectral_scores::Vector{SpectralScoresMainSearch{L}},
                weight::Vector{H},
                IDtoCOL::AbstractPrecursorMap{UInt16},
                ms_file_idx::Integer,
                cycle_idx::Integer,
                expected_matches::Float64,
                last_val::Int64,
                n_vals::Int64,
                spectrum_intensity::H,
                scan_idx::Int64;
                block_size::Int64 = 10000,
                default_top3_ll::Float32 = Float32(0)
                ) where {L,H<:AbstractFloat}

    getPoisson(lam, observed) = psm_getPoisson(lam, observed)
    start_idx = last_val
    skipped = 0
    skipped_weight = 0
    skipped_frag_count = 0
    for i in range(1, n_vals)

        precursor_idx = UInt32(unscored_PSMs[i].precursor_idx)
        scores_idx = IDtoCOL[precursor_idx]

        if weight[scores_idx] < Float32(1e-6)
            skipped += 1
            skipped_weight += 1
            continue
        end


        if start_idx + i - skipped > length(scored_psms)
            growScoredPSMs!(scored_psms, block_size);
        end

        total_ions = Int64(unscored_PSMs[i].y_count + unscored_PSMs[i].b_count)
        total_ions_iso = Int64(unscored_PSMs[i].isotope_count)

        scored_psms[start_idx + i - skipped] = MainSearchScoredPSM(
            unscored_PSMs[i].longest_y,
            unscored_PSMs[i].b_count,
            unscored_PSMs[i].y_count,
            UInt8(min(total_ions, 255)),
            UInt8(min(total_ions_iso, 255)),
            Float16(getPoisson(expected_matches, total_ions + total_ions_iso)),
            Float16(log2(max(unscored_PSMs[i].b_int + unscored_PSMs[i].y_int, Float32(1e-20))/max(spectrum_intensity, Float32(1e-20)))),
            Float16(log2(max(unscored_PSMs[i].error, Float32(1e-20)))),

            spectral_scores[scores_idx].gof,
            spectral_scores[scores_idx].max_matched_residual,
            spectral_scores[scores_idx].max_unmatched_residual,
            spectral_scores[scores_idx].fitted_manhattan_distance,
            spectral_scores[scores_idx].scribe,
            spectral_scores[scores_idx].percent_theoretical_ignored,
            weight[scores_idx],

            spectral_scores[scores_idx].fitted_hellinger,

            UInt32(unscored_PSMs[i].precursor_idx),
            UInt32(ms_file_idx),
            UInt32(cycle_idx),
            UInt32(scan_idx)
        )
        last_val += 1
    end
    return (last_val=last_val, skipped_weight=skipped_weight, skipped_frag_count=skipped_frag_count, skipped_matched_ratio=0, skipped_topn=0, skipped_spectral_contrast=0)
end
