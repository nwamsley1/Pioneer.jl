include("PSM.jl");

mutable struct XTandem{T<:Real} <: PSM
    best_rank::UInt8 #Highest ranking predicted framgent that was observed
    topn::UInt8 #How many of the topN predicted fragments were observed. 
    b_count::Int64
    b_int::T
    b_ions::Set{UInt8}
    y_count::Int64
    y_int::T
    y_ions::Set{UInt8}
    intensities::Vector{T}
    error::T
    precursor_idx::UInt32
    ms_file_idx::UInt32
end

#FastXTandem(type::DataType) where {T<:Real} = FastXTandem(0, convert(typeof(T), 0), 0, convert(typeof(T), 0), Set(UInt8[]), convert(typeof(T), 0), UInt32(0), UInt32(0))
#results = UnorderedDictionary{UInt32, FastXTandem}()
#Base.zero(::Type{FragmentMatch}) = FragmentMatch()

XTandem(::Type{Float32}) = XTandem(UInt8(255), zero(UInt8),0, Float32(0), Set(UInt8[]), 0, Float32(0), Set(UInt8[]), Float32[], Float32(0), UInt32(0), UInt32(0))
XTandem(::Type{Float64}) = XTandem(UInt8(255),zero(UInt8),0, Float64(0), Set(UInt8[]), 0, Float64(0), Set(UInt8[]), Float64[], Float64(0), UInt32(0), UInt32(0))
XTandem() = XTandem(Float64)

function ScoreFragmentMatches!(results::UnorderedDictionary{UInt32, XTandem{U}}, matches::Vector{FragmentMatch{T}}, nmatches::Int64, errdist::Laplace{Float64}) where {T,U<:Real}
    #UnorderedDictionary{UInt32, FastXTandem}()
    for i in range(1, nmatches)
        match = matches[i]
        if !isassigned(results, getPrecID(match))
            insert!(results, getPrecID(match), XTandem(U))
            results[getPrecID(match)].ms_file_idx = getMSFileID(match)
        end
        ModifyFeatures!(results[getPrecID(match)], match, match.match_mz, getIntensity(match), errdist)
    end
    #results
end

function ModifyFeatures!(score::XTandem{U}, match::FragmentMatch{T}, mass::Union{Missing, T}, intensity::Union{Missing, T}, errdist::Laplace{Float64}) where {U,T<:Real}
    if getIonType(match) == 'b'
        score.b_count += 1
        score.b_int += intensity
        push!(score.b_ions, getFragInd(match))
    elseif getIonType(match) == 'y'
        score.y_count += 1
        score.y_int += intensity
        push!(score.y_ions, getFragInd(match))
    end
    #ppm_err = ((mass - getFragMZ(match))/(mass/1e6))  - errdist.θ
    #score.error += Distributions.logpdf(errdist, ppm_err) - Distributions.logpdf(Uniform(0, 13.073079713498782), min(abs(ppm_err), 13.073079713498782)) #abs(mass - getFragMZ(match))
    #ppm_err =  (getFragMZ(match) - mass)/(getFragMZ(match)/1e6)
    #score.error += Distributions.logpdf(Laplace((-1.0)*errdist.μ,  errdist.θ), ppm_err) - Distributions.logpdf(Uniform(0, 13.073079713498782), min(abs(ppm_err), 13.073079713498782)) #abs(mass - getFragMZ(match))
    ppm_err = (getFragMZ(match))/(getFragMZ(match)/1e6)
    score.error +=  Distributions.logpdf(errdist, ppm_err*sqrt(match.intensity))
    #Distributions.logpdf(errdist, abs(mass - getFragMZ(match))/(mass/1e6))

    score.precursor_idx = getPrecID(match)
    push!(score.intensities, intensity)
    if match.predicted_rank < score.best_rank
        score.best_rank = match.predicted_rank
    end
    if match.predicted_rank <= 3
        score.topn += one(UInt8)
    end
    #push!(results[prec_id].test, getIonType(Transitions[transition]))
end

using SpecialFunctions

function makePSMsDict(::XTandem{T}) where {T<:Real}
    Dict(

        #Stanard Metrics 
        :hyperscore => Float16[],
        :error => Float16[],
        :y_ladder => UInt8[],
        :b_ladder => UInt8[],
        :poisson => Float16[],
        :total_ions => UInt16[],
        #:entropy => T[],
        :spectrum_peaks => UInt16[],
        :intensity_explained => Float16[],
        :best_rank => UInt8[],
        :topn => UInt8[],

        #Spectral Distrance Metrics
        :spectral_contrast => Float16[],
        :scribe_score => Float16[],
        :city_block => Float16[],
        :matched_ratio => Float16[],
        :entropy_sim => Float16[],

        :weight => Float32[],
        :scan_idx => UInt32[],
        :precursor_idx => UInt32[],
        :ms_file_idx => UInt16[]
    )
end

function Score!(PSMs_dict::Dict, 
                unscored_PSMs::UnorderedDictionary{UInt32, XTandem{T}}, 
                spectrum_peaks::Int, 
                spectrum_intensity::U, 
                expected_matches::Float64,
                #scribe_score::Vector{Float32}, 
                #city_block::Vector{Float32}, 
                #matched_ratio::Vector{Float32}, 
                #spectral_contrast::Vector{Float32}, 
                #entropy_sim::Vector{Float32}, 
                scores::NamedTuple,
                weight::Vector{Float32}, 
                IDtoROW::UnorderedDictionary{UInt32, Tuple{UInt32, UInt8}},
                IDtoROW_weights::UnorderedDictionary{UInt32, UInt32}; 
                scan_idx::Int64 = 0, min_spectral_contrast::Float32 = 0.6, min_frag_count::Int = 4) where {T,U<:Real}

    scribe_score = scores[:scribe]::Vector{Float32}
    city_block = scores[:city_block]::Vector{Float32}
    matched_ratio = scores[:matched_ratio]::Vector{Float32}
    spectral_contrast = scores[:spectral_contrast]::Vector{Float32}
    entropy_sim = scores[:entropy_sim]::Vector{Float32}
    #Get Hyperscore. Kong, Leprevost, and Avtonomov https://doi.org/10.1038/nmeth.4256
    #log(Nb!Ny!∑Ib∑Iy)
    function HyperScore(score::XTandem)
        #Ramanujan's approximation for log(!n)
        function logfac(N)
            N*log(N) - N + (log(N*(1 + 4*N*(1 + 2*N))))/6 + log(π)/2
        end
        (abs(logfac(max(1, score.b_count))) + 
        abs(logfac(max(1, score.y_count))) + 
        max(log(score.y_int), 0) + max(log(score.b_int), 0)
        )
    end

    function getPoisson(lam::T, observed::Int) where {T<:AbstractFloat}
        function logfac(N)
            N*log(N) - N + (log(N*(1 + 4*N*(1 + 2*N))))/6 + log(π)/2
        end
        log((lam^observed)*exp(-lam)) - logfac(observed)
    end

    #Get length of the the longest streak of consecutive y-ions
    function getLongestSet(x_ions_set::Set{UInt8})
        #Must sort y-ions by index
        #Might not be the most efficient. 
        x_ions = sort(collect(x_ions_set));

        #No calculation needed 
        if length(x_ions) <= 1
            return length(x_ions)
        end

        x_streak = 1
        longest_x_streak = 1
        for i in 2:length(x_ions)
            #If this isn't true, the sequence/ladder
            #is broken and needs to start over. 
            if (x_ions[i] - x_ions[i-1]) == 1
                x_streak += 1
            else
                #Need to check if this is the longest
                #streak before resetting the counter
                if x_streak > longest_x_streak
                    longest_x_streak = x_streak
                end
                x_streak = 1
            end
        end
        if x_streak > longest_x_streak
            return x_streak
        end
        return longest_x_streak
    end

    for key in keys(unscored_PSMs)
        if !haskey(IDtoROW, unscored_PSMs[key].precursor_idx)
            continue
        end

        
        index = IDtoROW[unscored_PSMs[key].precursor_idx][1]
        
        if spectral_contrast[index]<min_spectral_contrast
            continue
        end
        
        if (unscored_PSMs[key].y_count + unscored_PSMs[key].b_count) < min_frag_count
            continue
        end
        #if matched_ratio[index]<0.5
        #    continue
        #end
        #index = IDtoROW[unscored_PSMs[key].precursor_idx]
        append!(PSMs_dict[:hyperscore], Float16(HyperScore(unscored_PSMs[key])))
        append!(PSMs_dict[:y_ladder], UInt8(getLongestSet(unscored_PSMs[key].y_ions)))
        append!(PSMs_dict[:b_ladder], UInt8(getLongestSet(unscored_PSMs[key].b_ions)))
        append!(PSMs_dict[:total_ions], UInt16(unscored_PSMs[key].y_count + unscored_PSMs[key].b_count))
        append!(PSMs_dict[:poisson], Float16(getPoisson(expected_matches, unscored_PSMs[key].y_count + unscored_PSMs[key].b_count)))
        
        #append!(PSMs_dict[:error], unscored_PSMs[key].error/(unscored_PSMs[key].y_count + unscored_PSMs[key].b_count))
        append!(PSMs_dict[:error], Float16(log2(abs(unscored_PSMs[key].error))))
        append!(PSMs_dict[:spectrum_peaks], UInt16(spectrum_peaks))
        append!(PSMs_dict[:intensity_explained], Float16((sum(unscored_PSMs[key].b_int) + sum(unscored_PSMs[key].y_int))/spectrum_intensity))
        #append!(PSMs_dict[:entropy], StatsBase.entropy(unscored_PSMs[key].intensities))
        append!(PSMs_dict[:best_rank], UInt8(unscored_PSMs[key].best_rank))
        append!(PSMs_dict[:topn], UInt8(unscored_PSMs[key].topn))

        append!(PSMs_dict[:spectral_contrast], Float16(spectral_contrast[index]))
        append!(PSMs_dict[:scribe_score], Float16(scribe_score[index]))
        append!(PSMs_dict[:city_block], Float16(city_block[index]))
        append!(PSMs_dict[:matched_ratio], Float16(matched_ratio[index]))
        append!(PSMs_dict[:entropy_sim], Float16(entropy_sim[index]))
        append!(PSMs_dict[:weight], Float32(weight[IDtoROW_weights[unscored_PSMs[key].precursor_idx]]))

        append!(PSMs_dict[:scan_idx], UInt32(scan_idx))
        append!(PSMs_dict[:precursor_idx], UInt32(unscored_PSMs[key].precursor_idx))
        append!(PSMs_dict[:ms_file_idx], UInt16(unscored_PSMs[key].ms_file_idx))
    end
end
