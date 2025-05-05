abstract type UnscoredPSM{T<:AbstractFloat} <: PSM end

struct SimpleUnscoredPSM{T<:AbstractFloat} <: UnscoredPSM{T}
    best_rank::UInt8 #Highest ranking predicted framgent that was observed
    topn::UInt8 #How many of the topN predicted fragments were observed. 
    b_count::UInt8
    y_count::UInt8
    p_count::UInt8
    i_count::UInt8
    intensity::T
    error::T
    precursor_idx::UInt32
end

SimpleUnscoredPSM{Float32}() = SimpleUnscoredPSM(UInt8(255), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(Float32), zero(Float32), UInt32(0))
#LXTandem(::Type{Float64}) = LXTandem(UInt8(255), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), Float64(0), zero(UInt8), Float64(0), Float64(0), UInt32(0), UInt32(0))
#SimpleUnscoredPSM() = SimpleUnscoredPSM(Float32)

struct ComplexUnscoredPSM{T<:AbstractFloat} <: UnscoredPSM{T}
    best_rank::UInt8 #Highest ranking predicted framgent that was observed
    best_rank_iso::UInt8
    topn::UInt8 #How many of the topN predicted fragments were observed. 
    topn_iso::UInt8
    longest_y::UInt8
    longest_b::UInt8
    isotope_count::UInt8
    b_count::UInt8
    b_int::T
    y_count::UInt8
    y_int::T
    p_count::UInt8
    non_cannonical_count::UInt8
    error::T
    precursor_idx::UInt32
    ms_file_idx::UInt32
end

ComplexUnscoredPSM{Float32}() = ComplexUnscoredPSM(UInt8(255), UInt8(255), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), Float32(0), zero(UInt8), Float32(0), zero(UInt8), zero(UInt8), Float32(0), UInt32(0), UInt32(0))

struct Ms1UnscoredPSM{T<:AbstractFloat} <: UnscoredPSM{T}
    m0::Bool #Highest ranking predicted framgent that was observed
    n_iso::UInt8
    big_iso::UInt8 #How many of the topN predicted fragments were observed. 
    m0_error::Union{Missing,T}
    error::T
    precursor_idx::UInt32
end

Ms1UnscoredPSM{Float32}() = Ms1UnscoredPSM(
    false,
    zero(UInt8),
    zero(UInt8),
    missing,
    zero(Float32),
    zero(UInt32)
)

function ScoreFragmentMatches!(results::Vector{P}, 
                                IDtoCOL::ArrayDict{UInt32, UInt16}, 
                                matches::Vector{M}, 
                                nmatches::Int64, 
                                errdist::MassErrorModel,
                                 m_rank::Int64) where {P<:UnscoredPSM,M<:MatchIon}
    for i in range(1, nmatches)
        match = matches[i]
        prec_id = getPrecID(match)
        col = IDtoCOL[prec_id]
        results[col] = ModifyFeatures!(results[col], prec_id, match, errdist, m_rank)
    end
end

function ModifyFeatures!(score::SimpleUnscoredPSM{T}, prec_id::UInt32, match::FragmentMatch{T}, errdist::MassErrorModel, m_rank::Int64) where {T<:Real}
    
    best_rank = score.best_rank
    topn = score.topn
    b_count = score.b_count
    y_count = score.y_count
    p_count = score.p_count
    i_count = score.i_count
    intensity = score.intensity
    error = score.error
    precursor_idx = prec_id

    if match.is_isotope == false
        if match.predicted_rank < best_rank
            best_rank = match.predicted_rank
        end

        if match.predicted_rank <= m_rank
            topn += one(UInt8)
        end

        if getIonType(match) == one(UInt8)
            b_count += 1
        elseif getIonType(match) == UInt8(2)
            y_count += 1
        elseif getIonType(match) == UInt8(3)
            p_count += 1
        end

    else
        i_count += 1
    end


    ppm_err = (getMZ(match) - getMatchMZ(match))/(getMZ(match)/Float32(1e6))
    error += abs(ppm_err)#Distributions.logpdf(Laplace{Float64}(-1.16443, 4.36538), ppm_err)
    precursor_idx = getPrecID(match)
    intensity += getIntensity(match)

    return SimpleUnscoredPSM{T}(
        best_rank,
        topn,
        b_count,
        y_count, 
        p_count,
        i_count,
        intensity,
        error,
        precursor_idx
    )
end

function ModifyFeatures!(score::ComplexUnscoredPSM{T},  prec_id::UInt32, match::FragmentMatch{T}, errdist::MassErrorModel, m_rank::Int64) where {T<:Real}
    
    best_rank = score.best_rank
    best_rank_iso = score.best_rank_iso
    topn = score.topn
    topn_iso = score.topn_iso
    longest_y = score.longest_y
    longest_b = score.longest_b
    isotope_count = score.isotope_count
    b_count = score.b_count
    b_int = score.b_int
    y_count = score.y_count
    y_int = score.y_int
    p_count = score.p_count
    non_cannonical_count = score.non_cannonical_count
    error = score.error
    precursor_idx = prec_id

    if isIsotope(match)
        isotope_count += 1
        if match.predicted_rank < best_rank_iso
            best_rank_iso = match.predicted_rank
        end
        if match.predicted_rank <= m_rank
            topn_iso += one(UInt8)
        end
    else
        if match.predicted_rank < best_rank
            best_rank = match.predicted_rank
        end
        if match.predicted_rank <= m_rank
            topn += one(UInt8)
        end
        if getIonType(match) == one(UInt8)
            #score.b_count += 1
            #score.b_int += getIntensity(match)
            b_count += 1
            b_int += getIntensity(match)
            if getFragInd(match) > longest_b
                longest_b = getFragInd(match)
            end
        elseif getIonType(match) == UInt8(2)
            y_count += 1
            y_int +=  getIntensity(match)
            if getFragInd(match) > longest_y
                longest_y = getFragInd(match)
            end
        elseif getIonType(match) == UInt8(3)
            p_count += 1
        else
            non_cannonical_count += 1
        end
    end

    ppm_err = (getMZ(match) - getMatchMZ(match))/(getMZ(match)/1e6)
    error += abs(ppm_err)

    precursor_idx = getPrecID(match)

    return ComplexUnscoredPSM{T}(
        best_rank,
        best_rank_iso,
        min(topn, 255),
        min(topn_iso, 255),
        longest_y,
        longest_b,
        min(isotope_count, 255),
        min(b_count,255),
        b_int,
        min(y_count,255),
        y_int,
        p_count,
        non_cannonical_count,
        error,
        precursor_idx,
        score.ms_file_idx)
end

function ModifyFeatures!(score::Ms1UnscoredPSM{T},  prec_id::UInt32, match::PrecursorMatch{T}, errdist::MassErrorModel, m_rank::Int64) where {T<:Real}
    
    m0 = score.m0
    n_iso = score.n_iso
    big_iso = score.big_iso
    m0_error = score.m0_error
    error = score.error
    precursor_idx = prec_id

    n_iso += one(UInt8) 
    if match.iso_idx > big_iso
        big_iso = match.iso_idx
    end
    ppm_err =  (getMZ(match) - match.observed_mz)/(getMZ(match)/1e6)
    error += abs(ppm_err)
    precursor_idx = getPrecID(match)
    if getIsoIdx(match)==one(UInt8)#Is the m0
        m0 = true
        m0_error = ppm_err
    end

    return Ms1UnscoredPSM{T}(
        m0,
        n_iso,
        big_iso,
        m0_error,
        error,
        precursor_idx
)
end
using SpecialFunctions
