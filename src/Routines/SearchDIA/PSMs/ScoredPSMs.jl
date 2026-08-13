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

# Counters for diagnostic accounting of skipped PSMs across all Score! calls.
# Reset and read by the SearchDIA driver when PIONEER_LOG_SKIPPED=1.
const SKIPPED_WEIGHT_TOTAL = Threads.Atomic{Int64}(0)
const TOTAL_PSMS_CONSIDERED = Threads.Atomic{Int64}(0)

# Per-thread diagnostic capture of precursors that ever appeared in deconv,
# split by whether they passed the weight filter. Populated only when
# PIONEER_LOG_DROPPED_PRECS=1. The current ms_file_idx is set by the search
# driver before each file's deconv begins.
const DIAG_CURRENT_FILE_IDX = Ref{Int64}(0)
const DIAG_LOG_DROPPED = Ref{Bool}(false)
# Sized at runtime by SearchDIA driver to current Threads.maxthreadid(); we
# default to 1 here so module precompilation works with single-threaded julia.
const DIAG_SEEN_PER_THREAD = Dict{Int, Set{UInt32}}[]
const DIAG_DROPPED_PER_THREAD = Dict{Int, Set{UInt32}}[]
function _diag_ensure_thread_storage()
    n = Threads.maxthreadid()
    while length(DIAG_SEEN_PER_THREAD) < n
        push!(DIAG_SEEN_PER_THREAD, Dict{Int, Set{UInt32}}())
        push!(DIAG_DROPPED_PER_THREAD, Dict{Int, Set{UInt32}}())
    end
end

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
    y_count::UInt8
    total_ions::UInt8       # M0 (monoisotopic) fragment count

    #Basic Metrics
    poisson::L
    log2_intensity_explained::L
    error::L

    #Spectral Similarity
    gof::L
    max_matched_residual::L
    max_unmatched_residual::L
    fitted_manhattan_distance::L
    weight::H

    fitted_hellinger::L

    # Per-rank fitted and shadow intensities from the deconvolved design
    # matrix. These are temporary helper columns for selected-PSM cross-run
    # features and are dropped after the scalar summaries are computed.
    fitted_frag1_int::H
    fitted_frag2_int::H
    fitted_frag3_int::H
    fitted_frag4_int::H
    fitted_frag5_int::H
    fitted_frag6_int::H
    fitted_frag7_int::H
    fitted_frag8_int::H
    shadow_frag1_int::H
    shadow_frag2_int::H
    shadow_frag3_int::H
    shadow_frag4_int::H
    shadow_frag5_int::H
    shadow_frag6_int::H
    shadow_frag7_int::H
    shadow_frag8_int::H

    # Per-rank fragment trace intensities (top 8); sums matched isotope peaks
    # predicted at >=25% of the fragment's most abundant isotope.
    frag1_int::H
    frag2_int::H
    frag3_int::H
    frag4_int::H
    frag5_int::H
    frag6_int::H
    frag7_int::H
    frag8_int::H

    # E7 (Batch E, 2026-05-12): mean |ppm_err| over matched M0 fragments at
    # ranks 1-3. Zero if no top-3 matches in this scan.
    top3_ms2_mass_error_mean::H

    # E6 M0 (Batch E, 2026-05-12): log(b_int + 1) − log(y_int + 1). Real
    # peptides have characteristic b/y intensity ratios; chimeric hits don't.
    log_by_ratio_m0::L

    #Non-scores/Labels
    precursor_idx::UInt32
    ms_file_idx::UInt32
    cycle_idx::UInt32
    scan_idx::UInt32
end
function growScoredPSMs!(scored_psms::Vector{MainSearchScoredPSM{H,L}}, block_size::Int64) where {L,H<:AbstractFloat}
    # Grow in place (geometric buffer growth) without allocating a throwaway
    # `Vector(undef, block_size)` temp each call. New elements are undef, same
    # as the previous append! — they get overwritten by subsequent scores.
    resize!(scored_psms, length(scored_psms) + block_size)
end
function Score!(scored_psms::Vector{MainSearchScoredPSM{H, L}},
                unscored_PSMs::Vector{MainUnscoredPSM{H}},
                spectral_scores::Vector{SpectralScoresMainSearch{S,I}},
                weight::Vector{H},
                IDtoCOL::AbstractPrecursorMap{UInt16},
                ms_file_idx::Int64,
                cycle_idx::Int64,
                expected_matches::Float64,
                last_val::Int64,
                n_vals::Int64,
                spectrum_intensity::H,
                scan_idx::Int64;
                block_size::Int64 = 10000,
                default_top3_ll::Float32 = Float32(0)
                ) where {L,H<:AbstractFloat,S<:AbstractFloat,I<:AbstractFloat}

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
            unscored_PSMs[i].y_count,
            UInt8(min(total_ions, 255)),
            Float16(getPoisson(expected_matches, total_ions + total_ions_iso)),
            Float16(log2(max(unscored_PSMs[i].b_int + unscored_PSMs[i].y_int, Float32(1e-20))/max(spectrum_intensity, Float32(1e-20)))),
            Float16(log2(max(unscored_PSMs[i].error, Float32(1e-20)))),

            L(spectral_scores[scores_idx].gof),
            L(spectral_scores[scores_idx].max_matched_residual),
            L(spectral_scores[scores_idx].max_unmatched_residual),
            L(spectral_scores[scores_idx].fitted_manhattan_distance),
            weight[scores_idx],

            L(spectral_scores[scores_idx].fitted_hellinger),

            H(spectral_scores[scores_idx].fitted_frag1_int),
            H(spectral_scores[scores_idx].fitted_frag2_int),
            H(spectral_scores[scores_idx].fitted_frag3_int),
            H(spectral_scores[scores_idx].fitted_frag4_int),
            H(spectral_scores[scores_idx].fitted_frag5_int),
            H(spectral_scores[scores_idx].fitted_frag6_int),
            H(spectral_scores[scores_idx].fitted_frag7_int),
            H(spectral_scores[scores_idx].fitted_frag8_int),
            H(spectral_scores[scores_idx].shadow_frag1_int),
            H(spectral_scores[scores_idx].shadow_frag2_int),
            H(spectral_scores[scores_idx].shadow_frag3_int),
            H(spectral_scores[scores_idx].shadow_frag4_int),
            H(spectral_scores[scores_idx].shadow_frag5_int),
            H(spectral_scores[scores_idx].shadow_frag6_int),
            H(spectral_scores[scores_idx].shadow_frag7_int),
            H(spectral_scores[scores_idx].shadow_frag8_int),

            unscored_PSMs[i].frag1_int,
            unscored_PSMs[i].frag2_int,
            unscored_PSMs[i].frag3_int,
            unscored_PSMs[i].frag4_int,
            unscored_PSMs[i].frag5_int,
            unscored_PSMs[i].frag6_int,
            unscored_PSMs[i].frag7_int,
            unscored_PSMs[i].frag8_int,

            H(unscored_PSMs[i].top3_ppm_err_count > 0 ?
                unscored_PSMs[i].top3_abs_ppm_err_sum / unscored_PSMs[i].top3_ppm_err_count :
                zero(H)),

            L(log(Float32(unscored_PSMs[i].b_int) + 1f0) -
              log(Float32(unscored_PSMs[i].y_int) + 1f0)),

            UInt32(unscored_PSMs[i].precursor_idx),
            UInt32(ms_file_idx),
            UInt32(cycle_idx),
            UInt32(scan_idx)
        )
        last_val += 1
    end
    Threads.atomic_add!(SKIPPED_WEIGHT_TOTAL, skipped_weight)
    Threads.atomic_add!(TOTAL_PSMS_CONSIDERED, n_vals)
    return (last_val=last_val, skipped_weight=skipped_weight, skipped_frag_count=skipped_frag_count, skipped_matched_ratio=0, skipped_topn=0, skipped_spectral_contrast=0)
end

"""
    TuningScoredPSM{H,L} <: ScoredPSM{H,L}

Slim per-scan row produced by `Score!` for tuning paths
(ParameterTuningSearch, QuadTuningSearch). Same shape as
`MainSearchScoredPSM` minus the MainSearch-only fragment-chromatogram
fields (`rank1_matched`, `top3_matched`, `top5_matched`,
`frag1_int..frag8_int`) since tuning code never consumes them.
"""
struct TuningScoredPSM{H,L<:AbstractFloat} <: ScoredPSM{H,L}
    longest_y::UInt8
    y_count::UInt8
    total_ions::UInt8

    poisson::L
    log2_intensity_explained::L
    error::L

    gof::L
    max_matched_residual::L
    max_unmatched_residual::L
    fitted_manhattan_distance::L
    weight::H

    fitted_hellinger::L

    precursor_idx::UInt32
    ms_file_idx::UInt32
    cycle_idx::UInt32
    scan_idx::UInt32
end

function growScoredPSMs!(scored_psms::Vector{TuningScoredPSM{H,L}}, block_size::Int64) where {L,H<:AbstractFloat}
    resize!(scored_psms, length(scored_psms) + block_size)
end

function Score!(scored_psms::Vector{TuningScoredPSM{H, L}},
                unscored_PSMs::Vector{TuningUnscoredPSM{H}},
                spectral_scores::Vector{SpectralScoresMainSearch{S,I}},
                weight::Vector{H},
                IDtoCOL::AbstractPrecursorMap{UInt16},
                ms_file_idx::Int64,
                cycle_idx::Int64,
                expected_matches::Float64,
                last_val::Int64,
                n_vals::Int64,
                spectrum_intensity::H,
                scan_idx::Int64;
                block_size::Int64 = 10000,
                default_top3_ll::Float32 = Float32(0)
                ) where {L,H<:AbstractFloat,S<:AbstractFloat,I<:AbstractFloat}

    getPoisson(lam, observed) = psm_getPoisson(lam, observed)
    start_idx = last_val
    skipped = 0
    skipped_weight = 0
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

        scored_psms[start_idx + i - skipped] = TuningScoredPSM{H,L}(
            unscored_PSMs[i].longest_y,
            unscored_PSMs[i].y_count,
            UInt8(min(total_ions, 255)),
            L(getPoisson(expected_matches, total_ions + total_ions_iso)),
            L(log2(max(unscored_PSMs[i].b_int + unscored_PSMs[i].y_int, Float32(1e-20))/max(spectrum_intensity, Float32(1e-20)))),
            L(log2(max(unscored_PSMs[i].error, Float32(1e-20)))),

            L(spectral_scores[scores_idx].gof),
            L(spectral_scores[scores_idx].max_matched_residual),
            L(spectral_scores[scores_idx].max_unmatched_residual),
            L(spectral_scores[scores_idx].fitted_manhattan_distance),
            weight[scores_idx],

            L(spectral_scores[scores_idx].fitted_hellinger),

            UInt32(unscored_PSMs[i].precursor_idx),
            UInt32(ms_file_idx),
            UInt32(cycle_idx),
            UInt32(scan_idx)
        )
        last_val += 1
    end
    Threads.atomic_add!(SKIPPED_WEIGHT_TOTAL, skipped_weight)
    Threads.atomic_add!(TOTAL_PSMS_CONSIDERED, n_vals)
    return (last_val=last_val, skipped_weight=skipped_weight, skipped_frag_count=0, skipped_matched_ratio=0, skipped_topn=0, skipped_spectral_contrast=0)
end
