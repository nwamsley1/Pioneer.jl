include("PSM.jl");

mutable struct XTandem{T<:Real} <: PSM
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

XTandem(::Type{Float32}) = XTandem(0, Float32(0), Set(UInt8[]), 0, Float32(0), Set(UInt8[]), Float32[], Float32(0), UInt32(0), UInt32(0))
XTandem(::Type{Float64}) = XTandem(0, Float64(0), Set(UInt8[]), 0, Float64(0), Set(UInt8[]), Float64[], Float64(0), UInt32(0), UInt32(0))
XTandem() = XTandem(Float64)

function ScoreFragmentMatches!(results::UnorderedDictionary{UInt32, XTandem{U}}, matches::Vector{FragmentMatch{T}}) where {T,U<:Real}
    #UnorderedDictionary{UInt32, FastXTandem}()
    for match in matches
        if !isassigned(results, getPrecID(match.transition))
            insert!(results, getPrecID(match.transition), XTandem(U))
            results[getPrecID(match.transition)].ms_file_idx = getMSFileID(match)
        end
        ModifyFeatures!(results[getPrecID(match.transition)], match.transition, match.match_mz, match.intensity)
    end
    #results
end

function ModifyFeatures!(score::XTandem{U}, transition::Transition, mass::Union{Missing, T}, intensity::Union{Missing, T}) where {U,T<:Real}
    if getIonType(transition) == 'b'
        score.b_count += 1
        score.b_int += intensity
        push!(score.b_ions, getInd(transition))
    elseif getIonType(transition) == 'y'
        score.y_count += 1
        score.y_int += intensity
        push!(score.y_ions, getInd(transition))
    end
    score.error += abs(mass - getMZ(transition))
    score.precursor_idx = getPrecID(transition)
    push!(score.intensities, intensity)
    #push!(results[prec_id].test, getIonType(Transitions[transition]))
end

using SpecialFunctions

function makePSMsDict(::XTandem{T}) where {T<:Real}
    Dict(
        :hyperscore => T[],
        #:Δhyperscore => Float64[]
        :error => T[],
        :y_ladder => Int8[],
        :b_ladder => Int8[],
        :poisson => T[],
        :total_ions => UInt32[],
        :entropy => T[],
        :spectrum_peaks => UInt32[],
        :intensity_explained => T[],

        :scan_idx => Int64[],
        :precursor_idx => UInt32[],
        :ms_file_idx => UInt32[]
    )
end

function Score!(PSMs_dict::Dict, unscored_PSMs::UnorderedDictionary{UInt32, XTandem{T}}, spectrum_peaks::Int, spectrum_intensity::T, expected_matches::Float64; scan_idx::Int64 = 0) where {T<:Real}
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

    function getPoisson(lam::T, observed::Int)
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
        append!(PSMs_dict[:hyperscore], HyperScore(unscored_PSMs[key]))
        append!(PSMs_dict[:y_ladder], getLongestSet(unscored_PSMs[key].y_ions))
        append!(PSMs_dict[:b_ladder], getLongestSet(unscored_PSMs[key].b_ions))
        append!(PSMs_dict[:total_ions], unscored_PSMs[key].y_count + unscored_PSMs[key].b_count)
        append!(PSMs_dict[:poisson], getPoisson(expected_matches, unscored_PSMs[key].y_count + unscored_PSMs[key].b_count))
        append!(PSMs_dict[:error], unscored_PSMs[key].error/(unscored_PSMs[key].y_count + unscored_PSMs[key].b_count))
        append!(PSMs_dict[:spectrum_peaks], spectrum_peaks)
        append!(PSMs_dict[:intensity_explained], (sum(unscored_PSMs[key].b_int) + sum(unscored_PSMs[key].y_int))/spectrum_intensity)
        append!(PSMs_dict[:entropy], StatsBase.entropy(unscored_PSMs[key].intensities))


        append!(PSMs_dict[:scan_idx], scan_idx)
        append!(PSMs_dict[:precursor_idx], unscored_PSMs[key].precursor_idx)
        append!(PSMs_dict[:ms_file_idx], unscored_PSMs[key].ms_file_idx)
    end
end
