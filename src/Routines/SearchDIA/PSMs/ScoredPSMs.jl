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


struct SimpleScoredPSM{H,L<:AbstractFloat} <: ScoredPSM{H,L}
    #H is "high precision"
    #L is "low precision"

    #Ion Count Statistics
    best_rank::UInt8 #Highest ranking predicted framgent that was observed
    topn::UInt8 #How many of the topN predicted fragments were observed. 

    b_count::UInt8
    y_count::UInt8
    p_count::UInt8
    i_count::UInt8
    #Basic Metrics 
    poisson::L
    error::H

    #Spectral Simmilarity
    scribe::L
    city_block::L
    spectral_contrast::L
    matched_ratio::L
    log2_summed_intensity::L
    entropy_score::L
    percent_theoretical_ignored::L

    #Non-scores/Labels
    precursor_idx::UInt32
    scan_idx::UInt32
end

struct ComplexScoredPSM{H,L<:AbstractFloat} <: ScoredPSM{H,L}
    #H is "high precision"
    #L is "low precision"

    #Ion Count Statistics
    best_rank::UInt8 #Highest ranking predicted framgent that was observed
    best_rank_iso::UInt8
    topn::UInt8 #How many of the topN predicted fragments were observed. 
    topn_iso::UInt8
    longest_y::UInt8
    longest_b::UInt8
    b_count::UInt8
    y_count::UInt8
    p_count::UInt8
    non_cannonical_count::UInt8
    isotope_count::UInt8

    #Basic Metrics 
    poisson::L
    #hyperscore::L
    log2_intensity_explained::L
    error::L

    #Spectral Simmilarity
    spectral_contrast::L
    fitted_spectral_contrast::L
    gof::L
    max_matched_residual::L
    max_unmatched_residual::L 
    fitted_manhattan_distance::L 
    matched_ratio::L 
    percent_theoretical_ignored::L
    scribe::L
    #entropy_score::L
    weight::H

    #Non-scores/Labels
    precursor_idx::UInt32
    ms_file_idx::UInt32
    scan_idx::UInt32
end

struct Ms1ScoredPSM{H,L<:AbstractFloat} <: ScoredPSM{H,L}
    #H is "high precision"
    #L is "low precision"

    m0::Bool #Highest ranking predicted framgent that was observed
    n_iso::UInt8
    big_iso::UInt8 #How many of the topN predicted fragments were observed. 
    m0_error::Union{Missing,L}
    error::L

    #Spectral Simmilarity
    spectral_contrast::L
    fitted_spectral_contrast::L
    gof::L
    max_matched_residual::L
    max_unmatched_residual::L 
    fitted_manhattan_distance::L 
    matched_ratio::L 
    #entropy_score::L
    weight::H

    #Non-scores/Labels
    precursor_idx::UInt32
    ms_file_idx::UInt32
    scan_idx::UInt32
end

function growScoredPSMs!(scored_psms::Vector{SimpleScoredPSM{H,L}}, block_size::Int64) where {L,H<:AbstractFloat}
    scored_psms = append!(scored_psms, Vector{SimpleScoredPSM{H,L}}(undef, block_size))
end

function growScoredPSMs!(scored_psms::Vector{ComplexScoredPSM{H,L}}, block_size::Int64) where {L,H<:AbstractFloat}
    scored_psms = append!(scored_psms, Vector{ComplexScoredPSM{H,L}}(undef, block_size))
end

function growScoredPSMs!(scored_psms::Vector{Ms1ScoredPSM{H,L}}, block_size::Int64) where {L,H<:AbstractFloat}
    scored_psms = append!(scored_psms, Vector{Ms1ScoredPSM{H,L}}(undef, block_size))
end

function Score!(scored_psms::Vector{SimpleScoredPSM{H, L}}, 
                unscored_PSMs::Vector{SimpleUnscoredPSM{H}}, 
                spectral_scores::Vector{SpectralScoresSimple{L}},
                IDtoCOL::ArrayDict{UInt32, UInt16},
                expected_matches::Float64,
                last_val::Int64,
                n_vals::Int64,
                spectrum_intensity::H,
                scan_idx::Int64;
                min_spectral_contrast::H = 0f0,
                min_log2_matched_ratio::H = -Inf32,
                min_frag_count::Int64 = 1,
                max_best_rank::UInt8 = one(UInt8),
                min_topn::Int64 = 1,
                block_size::Int64 = 10000
                ) where {L,H<:AbstractFloat}

    start_idx = last_val
    skipped = 0
    n = 0

    function getPoisson(lam::T, observed::Int64) where {T<:AbstractFloat}
        function logfac(N)
            N*log(N) - N + (log(N*(1 + 4*N*(1 + 2*N))))/6 + log(π)/2
        end
        log((lam^observed)*exp(-lam)) - logfac(observed)
    end

    for i in range(1, n_vals)

        passing_filter = (
            (spectral_scores[i].spectral_contrast) >= min_spectral_contrast
        )&(
            (unscored_PSMs[i].y_count + unscored_PSMs[i].b_count) >= min_frag_count
        )&(
            spectral_scores[i].matched_ratio > min_log2_matched_ratio
        )&(
            (unscored_PSMs[i].topn >= min_topn)
        )&(
            UInt8(unscored_PSMs[i].best_rank) == max_best_rank
        )

        if !passing_filter #Skip this scan
            skipped += 1
            continue
        end
        
        if start_idx + i - skipped > length(scored_psms)
            growScoredPSMs!(scored_psms, block_size);
        end

        precursor_idx = UInt32(unscored_PSMs[i].precursor_idx)
        scores_idx = IDtoCOL[precursor_idx]
        total_ions = Int64(unscored_PSMs[i].y_count + unscored_PSMs[i].b_count + unscored_PSMs[i].p_count + unscored_PSMs[i].i_count)
        scored_psms[start_idx + i - skipped] = SimpleScoredPSM(
            unscored_PSMs[i].best_rank,

            unscored_PSMs[i].topn,

            unscored_PSMs[i].b_count,
            unscored_PSMs[i].y_count,
            unscored_PSMs[i].p_count,
            unscored_PSMs[i].i_count,

            Float16(getPoisson(expected_matches, total_ions)),
            unscored_PSMs[i].error,
            
            spectral_scores[scores_idx].scribe,
            spectral_scores[scores_idx].city_block,
            spectral_scores[scores_idx].spectral_contrast,
            spectral_scores[scores_idx].matched_ratio,
            #Float16(log2((unscored_PSMs[i].intensity)/spectrum_intensity)),
            Float16(log2(unscored_PSMs[i].intensity)),
            spectral_scores[scores_idx].entropy_score,
            spectral_scores[scores_idx].percent_theoretical_ignored,

            UInt32(unscored_PSMs[i].precursor_idx),
            UInt32(scan_idx)
        )
        n += 1
        last_val += 1
    end
    return last_val
end

function Score!(scored_psms::Vector{ComplexScoredPSM{H, L}}, 
                unscored_PSMs::Vector{ComplexUnscoredPSM{H}}, 
                spectral_scores::Vector{SpectralScoresComplex{L}},
                weight::Vector{H}, 
                IDtoCOL::ArrayDict{UInt32, UInt16},
                cycle_idx::Int64,
                expected_matches::Float64,
                last_val::Int64,
                n_vals::Int64,
                spectrum_intensity::H,
                scan_idx::Int64;
                min_spectral_contrast::H = 0f0,
                min_log2_matched_ratio::H = -1f0,
                min_y_count::Int64,
                min_frag_count::Int64 = 4,
                max_best_rank::Int64 = 1,
                min_topn::Int64 = 2,
                block_size::Int64 = 10000
                ) where {L,H<:AbstractFloat}

    #Get Hyperscore. Kong, Leprevost, and Avtonomov https://doi.org/10.1038/nmeth.4256
    #log(Nb!Ny!∑Ib∑Iy)
    function HyperScore(score::ComplexUnscoredPSM{T}) where {T<:AbstractFloat}
        #Ramanujan's approximation for log(!n)
        function logfac(N)
            N*log(N) - N + (log(N*(1 + 4*N*(1 + 2*N))))/6 + log(π)/2
        end
        (abs(logfac(max(1, score.b_count))) + 
        abs(logfac(max(1, score.y_count))) + 
        max(log(score.y_int), 0) + max(log(score.b_int), 0)
        )
    end

    function getPoisson(lam::T, observed::Int64) where {T<:AbstractFloat}
        function logfac(N)
            N*log(N) - N + (log(N*(1 + 4*N*(1 + 2*N))))/6 + log(π)/2
        end
        log((lam^observed)*exp(-lam)) - logfac(observed)
    end
    start_idx = last_val
    skipped = 0
    for i in range(1, n_vals)
        
        passing_filter = (
           (unscored_PSMs[i].y_count) >= min_y_count
        )&(
            (unscored_PSMs[i].y_count + unscored_PSMs[i].b_count + unscored_PSMs[i].isotope_count) >= min_frag_count
        )&(
            (spectral_scores[i].fitted_spectral_contrast) >= min_spectral_contrast
        )&(
            spectral_scores[i].matched_ratio > min_log2_matched_ratio
        )&(
            weight[i] >= 1e-6
        )&(
            (unscored_PSMs[i].topn >= min_topn) |  (unscored_PSMs[i].topn_iso >= min_topn)
        )&(
            (UInt8(unscored_PSMs[i].best_rank) <= max_best_rank) | (UInt8(unscored_PSMs[i].best_rank_iso) <= max_best_rank) 
        )
        #passing_filter = true
        if !passing_filter #Skip this scan
            skipped += 1
            continue
        end

        if start_idx + i - skipped > length(scored_psms)
            growScoredPSMs!(scored_psms, block_size);
        end

        total_ions = Int64(unscored_PSMs[i].y_count + unscored_PSMs[i].b_count)

        precursor_idx = UInt32(unscored_PSMs[i].precursor_idx)
        scores_idx = IDtoCOL[precursor_idx]
        scored_psms[start_idx + i - skipped] = ComplexScoredPSM(
            unscored_PSMs[i].best_rank,
            unscored_PSMs[i].best_rank_iso,
            unscored_PSMs[i].topn,
            unscored_PSMs[i].topn_iso,
            unscored_PSMs[i].longest_y,
            unscored_PSMs[i].longest_b,
            unscored_PSMs[i].b_count,
            unscored_PSMs[i].y_count,
            unscored_PSMs[i].p_count,
            unscored_PSMs[i].non_cannonical_count,
            unscored_PSMs[i].isotope_count,
            Float16(getPoisson(expected_matches, total_ions)),
            #Float16(HyperScore(unscored_PSMs[i])),
            Float16(log2((unscored_PSMs[i].b_int + unscored_PSMs[i].y_int)/spectrum_intensity)),
            Float16(log2(unscored_PSMs[i].error)),
            
            spectral_scores[scores_idx].spectral_contrast,
            spectral_scores[scores_idx].fitted_spectral_contrast,
            spectral_scores[scores_idx].gof,
            spectral_scores[scores_idx].max_matched_residual,
            spectral_scores[scores_idx].max_unmatched_residual,
            spectral_scores[scores_idx].fitted_manhattan_distance,
            spectral_scores[scores_idx].matched_ratio,
            spectral_scores[scores_idx].percent_theoretical_ignored,
            spectral_scores[scores_idx].scribe,
            #spectral_scores[scores_idx].entropy_score,
            weight[scores_idx],

            
            UInt32(unscored_PSMs[i].precursor_idx),
            #UInt32(unscored_PSMs[i].ms_file_idx),
            UInt32(cycle_idx),
            UInt32(scan_idx)
        )
        last_val += 1
    end
    return last_val
end

function Score!(scored_psms::Vector{Ms1ScoredPSM{H, L}}, 
                unscored_PSMs::Vector{Ms1UnscoredPSM{H}}, 
                spectral_scores::Vector{SpectralScoresMs1{L}},
                weight::Vector{H}, 
                IDtoCOL::ArrayDict{UInt32, UInt16},
                cycle_idx::Int64,
                last_val::Int64,
                n_vals::Int64,
                scan_idx::Int64;
                block_size::Int64 = 10000
                ) where {L,H<:AbstractFloat}

    start_idx = last_val
    skipped = 0
    for i in range(1, n_vals)
        
        passing_filter = (
           (unscored_PSMs[i].m0) == true
        )&(
            (unscored_PSMs[i].n_iso) >= 2
        )&(
            weight[i] >= 1e-6
        )
        #passing_filter = true
        if !passing_filter #Skip this scan
            skipped += 1
            continue
        end

        if start_idx + i - skipped > length(scored_psms)
            growScoredPSMs!(scored_psms, block_size);
        end
        precursor_idx = UInt32(unscored_PSMs[i].precursor_idx)
        scores_idx = IDtoCOL[precursor_idx]
        scored_psms[start_idx + i - skipped] = Ms1ScoredPSM(
            unscored_PSMs[i].m0,
            unscored_PSMs[i].n_iso,
            unscored_PSMs[i].big_iso,
            L(coalesce(unscored_PSMs[i].m0_error, zero(H))),
            Float16(log2(unscored_PSMs[i].error + 1e-6)),
            
            spectral_scores[scores_idx].spectral_contrast,
            spectral_scores[scores_idx].fitted_spectral_contrast,
            spectral_scores[scores_idx].gof,
            spectral_scores[scores_idx].max_matched_residual,
            spectral_scores[scores_idx].max_unmatched_residual,
            spectral_scores[scores_idx].fitted_manhattan_distance,
            min(spectral_scores[scores_idx].matched_ratio, Float16(10.0)),
            #spectral_scores[scores_idx].entropy_score,
            weight[scores_idx],

            
            UInt32(unscored_PSMs[i].precursor_idx),
            UInt32(cycle_idx),
            UInt32(scan_idx)
        )
        last_val += 1
    end
    return last_val
end