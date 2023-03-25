include("PSM.jl");

mutable struct FastXTandem <: PSM
    b_count::Int64
    b_int::Float32
    y_count::Int64
    y_int::Float32
    last_y::UInt8
    y_start::UInt8
    longest_y::Int8
    error::Float64
end

FastXTandem() = FastXTandem(0, Float32(0), 0, Float32(0), UInt8(0), UInt8(0), UInt8(0), Float64(0))
results = UnorderedDictionary{UInt32, FastXTandem}()
#Base.zero(::Type{FragmentMatch}) = FragmentMatch()



function ScoreFragmentMatches!(results::UnorderedDictionary{UInt32, FastXTandem}, matches::Vector{FragmentMatch})
    #UnorderedDictionary{UInt32, FastXTandem}()
    for match in matches
        if !isassigned(results, getPrecID(match.transition))
            insert!(results, getPrecID(match.transition), FastXTandem())
        end
        ModifyFeatures!(results[getPrecID(match.transition)], match.transition, match.mass, match.intensity)
    end
    #results
end

function ModifyFeatures!(score::FastXTandem, transition::Transition, mass::Union{Missing, Float32}, intensity::Union{Missing, Float32})
    if getIonType(transition) == 'b'
        score.b_count += 1
        score.b_int += intensity
    elseif getIonType(transition) == 'y'
        score.y_count += 1
        score.y_int += intensity
        if score.longest_y == 0
            score.y_start = getInd(transition)
            score.last_y = getInd(transition) - 1
        end
        if (getInd(transition)-score.last_y) == 1
            if (getInd(transition) - score.y_start) > score.longest_y
                score.longest_y = getInd(transition) - score.y_start
            else
                score.y_start = getInd(transition)
            end
        end
    end
    score.error += abs(mass - getMZ(transition))
    #push!(results[prec_id].intensities, (intensities[best_peak]))
    #push!(results[prec_id].test, getIonType(Transitions[transition]))
end

using SpecialFunctions

function makePSMsDict(::FastXTandem)
    Dict(
        :hyperscore => Float64[],
        #:Δhyperscore => Float64[]
        :error => Float64[],
        :y_ladder => Int8[]
    )
end

function Score!(PSMs_dict::Dict, unscored_PSMs::UnorderedDictionary{UInt32, FastXTandem})
    #Get Hyperscore. Kong, Leprevost, and Avtonomov https://doi.org/10.1038/nmeth.4256
    #log(Nb!Ny!∑Ib∑Iy)
    function HyperScore(score::FastXTandem)
        #Ramanujan's approximation for log(!n)
        function logfac(N)
            N*log(N) - N + (log(N*(1 + 4*N*(1 + 2*N))))/6 + log(π)/2
        end
        (abs(logfac(max(1, score.b_count))) + 
        abs(logfac(max(1, score.y_count))) + 
        max(log(score.y_int*score.b_int), 0)
        )
    end
    for key in keys(unscored_PSMs)
        append!(PSMs_dict[:hyperscore], HyperScore(unscored_PSMs[key]))
        append!(PSMs_dict[:error], unscored_PSMs[key].error)
        append!(PSMs_dict[:y_ladder], unscored_PSMs[key].longest_y)
    end
end
