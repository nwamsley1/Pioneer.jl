include("PSM.jl");

mutable struct LXTandem{T<:AbstractFloat} <: PSM
    best_rank::UInt8 #Highest ranking predicted framgent that was observed
    topn::UInt8 #How many of the topN predicted fragments were observed. 
    longest_y::UInt8
    longest_b::UInt8
    b_count::UInt8
    b_int::T
    y_count::UInt8
    y_int::T
    error::T
    precursor_idx::UInt32
    ms_file_idx::UInt32
end

LXTandem(::Type{Float32}) = LXTandem(UInt8(255), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), Float32(0), zero(UInt8), Float32(0), Float32(0), UInt32(0), UInt32(0))
LXTandem(::Type{Float64}) = LXTandem(UInt8(255), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), Float64(0), zero(UInt8), Float64(0), Float64(0), UInt32(0), UInt32(0))
LXTandem() = LXTandem(Float64)

function ScoreFragmentMatches!(results::Vector{LXTandem{T}}, IDtoCOL::ArrayDict{UInt32, UInt16}, matches::Vector{FragmentMatch{U}}, nmatches::Int64, errdist::Laplace{Float64}) where {T,U<:AbstractFloat}
    for i in range(1, nmatches)
        match = matches[i]
        col = IDtoCOL.vals[getPrecID(match)]
        ModifyFeatures!(results[col], match, errdist)
    end
end

function ModifyFeatures!(score::LXTandem{U}, match::FragmentMatch{T}, errdist::Laplace{Float64}) where {U,T<:Real}
    if getIonType(match) == 'b'
        score.b_count += 1
        score.b_int += getIntensity(match)
        if getFragInd(match) > score.longest_b
            score.longest_b = getFragInd(match)
        end
    elseif getIonType(match) == 'y'
        score.y_count += 1
        score.y_int +=  getIntensity(match)
        if getFragInd(match) > score.longest_y
            score.longest_y = getFragInd(match)
        end
    end

    ppm_err = (getFragMZ(match) - getMatchMZ(match))/(getFragMZ(match)/1e6)
    score.error +=  Distributions.logpdf(errdist, ppm_err*sqrt(getIntensity(match)))

    score.precursor_idx = getPrecID(match)

    if match.predicted_rank < score.best_rank
        score.best_rank = match.predicted_rank
    end

    if match.predicted_rank <= 3
        score.topn += one(UInt8)
    end

end

using SpecialFunctions

struct LibPSM{H,L<:Real} <: PSM
    #H is "high precision"
    #L is "low precision"

    #Ion Count Statistics
    best_rank::UInt8 #Highest ranking predicted framgent that was observed
    topn::UInt8 #How many of the topN predicted fragments were observed. 
    longest_y::UInt8
    longest_b::UInt8
    b_count::UInt8
    y_count::UInt8

    #Basic Metrics 
    poisson::L
    hyperscore::L
    log2_intensity_explained::L
    error::H

    #Spectral Simmilarity
    scribe::L
    scribe_corrected::L
    city_block::L
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

function growScoredPSMs!(scored_psms::Vector{LibPSM{H,L}}, block_size::Int64) where {L,H<:AbstractFloat}
    scored_psms = append!(scored_psms, Vector{LibPSM{H,L}}(undef, block_size))
end

function Score!(scored_psms::Vector{LibPSM{H, L}}, 
                unscored_PSMs::Vector{LXTandem{H}}, 
                spectral_scores::Vector{SpectralScores{L}},
                weight::Vector{H}, 
                expected_matches::Float64,
                n_vals::Int64,
                spectrum_intensity::H; 
                scan_idx::Int64 = 0, 
                min_spectral_contrast::H = 0.6f0, 
                min_matched_ratio::H = -1f0,
                min_frag_count::Int64 = 4,
                min_weight::H = 100f0,
                min_topn::Int64 = 2,
                block_size::Int64 = 10000
                ) where {L,H<:AbstractFloat}

    #Get Hyperscore. Kong, Leprevost, and Avtonomov https://doi.org/10.1038/nmeth.4256
    #log(Nb!Ny!∑Ib∑Iy)
    function HyperScore(score::LXTandem{T}) where {T<:AbstractFloat}
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

    for i in range(1, n_vals)

        passing_filter = ( #Filter Bad PSMs and don't add them to the DataFrame
          spectral_scores[i].spectral_contrast > min_spectral_contrast
        )&(
            (unscored_PSMs[i].y_count + unscored_PSMs[i].b_count) > min_frag_count
        )&(
            spectral_scores[i].matched_ratio > min_matched_ratio
        )&(
            weight[i] > min_weight
        )&(
            UInt8(unscored_PSMs[i].topn) > min_topn
        )&(
            UInt8(unscored_PSMs[i].best_rank) == 1
        )
        passing_filter = true
        if !passing_filter #Skip this scan
            continue
        end
        if i > length(scored_psms)
            growScoredPSMs!(scored_psms, block_size);
        end

        total_ions = Int64(unscored_PSMs[i].y_count + unscored_PSMs[i].b_count)

        scored_psms[i] = LibPSM(
            unscored_PSMs[i].best_rank,
            unscored_PSMs[i].topn,
            unscored_PSMs[i].longest_y,
            unscored_PSMs[i].longest_b,
            unscored_PSMs[i].b_count,
            unscored_PSMs[i].y_count,

            Float16(getPoisson(expected_matches, total_ions)),
            Float16(HyperScore(unscored_PSMs[i])),
            Float16(log2((unscored_PSMs[i].b_int + unscored_PSMs[i].y_int)/spectrum_intensity)),
            unscored_PSMs[i].error,
            
            spectral_scores[i].scribe,
            spectral_scores[i].scribe_corrected,
            spectral_scores[i].city_block,
            spectral_scores[i].spectral_contrast,
            spectral_scores[i].spectral_contrast_corrected,
            spectral_scores[i].matched_ratio,
            spectral_scores[i].entropy_score,
            weight[i],

            UInt32(scan_idx),
            UInt32(unscored_PSMs[i].precursor_idx),
            UInt32(unscored_PSMs[i].ms_file_idx)
        )
    end
end