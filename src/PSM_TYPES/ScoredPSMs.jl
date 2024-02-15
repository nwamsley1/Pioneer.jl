abstract type ScoredPSM{H,L<:AbstractFloat} <: PSM end

using SpecialFunctions

struct SimpleScoredPSM{H,L<:AbstractFloat} <: ScoredPSM{H,L}
    #H is "high precision"
    #L is "low precision"

    #Ion Count Statistics
    best_rank::UInt8 #Highest ranking predicted framgent that was observed
    topn::UInt8 #How many of the topN predicted fragments were observed. 

    b_count::UInt8
    y_count::UInt8

    #Basic Metrics 
    error::H

    #Spectral Simmilarity
    scribe::L
    city_block::L
    spectral_contrast::L
    matched_ratio::L
    log2_intensity_explained::L
    entropy_score::L

    #Non-scores/Labels
    precursor_idx::UInt32
    scan_idx::UInt32
end

struct ComplexScoredPSM{H,L<:AbstractFloat} <: ScoredPSM{H,L}
    #H is "high precision"
    #L is "low precision"

    #Ion Count Statistics
    best_rank::UInt8 #Highest ranking predicted framgent that was observed
    topn::UInt8 #How many of the topN predicted fragments were observed. 
    longest_y::UInt8
    longest_b::UInt8
    b_count::UInt8
    y_count::UInt8
    isotope_count::UInt8
    prec_mz_offset::Float32

    #Basic Metrics 
    poisson::L
    hyperscore::L
    log2_intensity_explained::L
    error::H

    #Spectral Simmilarity
    scribe::L
    scribe_corrected::L
    scribe_fitted::L
    city_block::L
    city_block_fitted::L
    spectral_contrast::L
    spectral_contrast_corrected::L
    matched_ratio::L
    entropy_score::L
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

function Score!(scored_psms::Vector{SimpleScoredPSM{H, L}}, 
                unscored_PSMs::Vector{SimpleUnscoredPSM{H}}, 
                spectral_scores::Vector{SpectralScoresSimple{L}},
                weight::Vector{H}, 
                expected_matches::Float64,
                last_val::Int64,
                n_vals::Int64,
                spectrum_intensity::H,
                scan_idx::Int64;
                min_log2_matched_ratio::H = -Inf32,
                min_frag_count::Int64 = 1,
                min_weight::H = 0f0,
                min_topn::Int64 = 1,
                block_size::Int64 = 10000
                ) where {L,H<:AbstractFloat}

    start_idx = last_val
    skipped = 0
    for i in range(1, n_vals)

        passing_filter = (
            (unscored_PSMs[i].y_count + unscored_PSMs[i].b_count) >= min_frag_count
        )&(
            spectral_scores[i].matched_ratio > min_log2_matched_ratio
        )&(
            weight[i] >= min_weight
        )&(
            UInt8(unscored_PSMs[i].topn) >= min_topn
        )&(
            UInt8(unscored_PSMs[i].best_rank) == 1
        )

        if !passing_filter #Skip this scan
            skipped += 1
            continue
        end

        if start_idx + i - skipped > length(scored_psms)
            growScoredPSMs!(scored_psms, block_size);
        end

        scored_psms[start_idx + i - skipped] = SimpleScoredPSM(
            unscored_PSMs[i].best_rank,
            unscored_PSMs[i].topn,

            unscored_PSMs[i].b_count,
            unscored_PSMs[i].y_count,

            unscored_PSMs[i].error,
            
            spectral_scores[i].scribe,
            spectral_scores[i].city_block,
            spectral_scores[i].spectral_contrast,
            spectral_scores[i].matched_ratio,
            Float16(log2((unscored_PSMs[i].intensity)/spectrum_intensity)),
            spectral_scores[i].entropy_score,
            
            UInt32(unscored_PSMs[i].precursor_idx),
            UInt32(scan_idx)
        )
        last_val += 1
    end
    return last_val
end

function Score!(scored_psms::Vector{ComplexScoredPSM{H, L}}, 
                unscored_PSMs::Vector{ComplexUnscoredPSM{H}}, 
                spectral_scores::Vector{SpectralScoresComplex{L}},
                weight::Vector{H}, 
                expected_matches::Float64,
                last_val::Int64,
                n_vals::Int64,
                spectrum_intensity::H,
                scan_idx::Int64;
                min_log2_matched_ratio::H = -1f0,
                min_frag_count::Int64 = 4,
                min_weight::H = 100f0,
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
            (unscored_PSMs[i].y_count + unscored_PSMs[i].b_count) >= min_frag_count
        )&(
            spectral_scores[i].matched_ratio > min_log2_matched_ratio
        )&(
            weight[i] >= min_weight
        )&(
            UInt8(unscored_PSMs[i].topn) >= min_topn
        )&(
            UInt8(unscored_PSMs[i].best_rank) == 1
        )

        if !passing_filter #Skip this scan
            skipped += 1
            continue
        end

        if start_idx + i - skipped > length(scored_psms)
            growScoredPSMs!(scored_psms, block_size);
        end

        total_ions = Int64(unscored_PSMs[i].y_count + unscored_PSMs[i].b_count)

        scored_psms[start_idx + i - skipped] = ComplexScoredPSM(
            unscored_PSMs[i].best_rank,
            unscored_PSMs[i].topn,
            unscored_PSMs[i].longest_y,
            unscored_PSMs[i].longest_b,
            unscored_PSMs[i].b_count,
            unscored_PSMs[i].y_count,
            unscored_PSMs[i].isotope_count,
            zero(Float32),

            Float16(getPoisson(expected_matches, total_ions)),
            Float16(HyperScore(unscored_PSMs[i])),
            Float16(log2((unscored_PSMs[i].b_int + unscored_PSMs[i].y_int)/spectrum_intensity)),
            unscored_PSMs[i].error,
            
            spectral_scores[i].scribe,
            spectral_scores[i].scribe_corrected,
            spectral_scores[i].scribe_fitted,
            spectral_scores[i].city_block,
            spectral_scores[i].city_block_fitted,
            spectral_scores[i].spectral_contrast,
            spectral_scores[i].spectral_contrast_corrected,
            spectral_scores[i].matched_ratio,
            spectral_scores[i].entropy_score,
            weight[i],

            
            UInt32(unscored_PSMs[i].precursor_idx),
            UInt32(unscored_PSMs[i].ms_file_idx),
            UInt32(scan_idx)
        )
        last_val += 1
    end
    return last_val
end