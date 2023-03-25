include("PSM.jl");

mutable struct FastXTandem <: PSM
    b_count::Int64
    b_int::Float32
    y_count::Int64
    y_int::Float32
    y_ions::Set{UInt8}
    error::Float64
    precursor_idx::Int32
end

FastXTandem() = FastXTandem(0, Float32(0), 0, Float32(0), Set(UInt8[]), Float64(0), Int32(0))
results = UnorderedDictionary{UInt32, FastXTandem}()
#Base.zero(::Type{FragmentMatch}) = FragmentMatch()



function ScoreFragmentMatches!(results::UnorderedDictionary{UInt32, FastXTandem}, matches::Vector{FragmentMatch})
    #UnorderedDictionary{UInt32, FastXTandem}()
    for match in matches
        if !isassigned(results, getPrecID(match.transition))
            insert!(results, getPrecID(match.transition), FastXTandem())
        end
        ModifyFeatures!(results[getPrecID(match.transition)], match.transition, match.match_mz, match.intensity)
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
        push!(score.y_ions, getInd(transition))
    end
    score.error += abs(mass - getMZ(transition))
    score.precursor_idx = getPrecID(transition)
    #push!(results[prec_id].intensities, (intensities[best_peak]))
    #push!(results[prec_id].test, getIonType(Transitions[transition]))
end

using SpecialFunctions

function makePSMsDict(::FastXTandem)
    Dict(
        :hyperscore => Float64[],
        #:Δhyperscore => Float64[]
        :error => Float64[],
        :y_ladder => Int8[],
        :scan_idx => Int64[],
        :precursor_idx => Int32[]
    )
end

function Score!(PSMs_dict::Dict, unscored_PSMs::UnorderedDictionary{UInt32, FastXTandem}; scan_idx::Int64 = 0)
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

    #Get length of the the longest streak of consecutive y-ions
    function getLongestY(y_ions_set::Set{UInt8})
        #Must sort y-ions by index
        #Might not be the most efficient. 
        y_ions = sort(collect(y_ions_set));

        #No calculation needed 
        if length(y_ions) <= 1
            return length(y_ions)
        end

        y_streak = 1
        longest_y_streak = 1
        for i in 2:length(y_ions)
            #If this isn't true, the sequence/ladder
            #is broken and needs to start over. 
            if (y_ions[i] - y_ions[i-1]) == 1
                y_streak += 1
            else
                #Need to check if this is the longest
                #streak before resetting the counter
                if y_streak > longest_y_streak
                    longest_y_streak = y_streak
                end
                y_streak = 1
            end
        end
        if y_streak > longest_y_streak
            return y_streak
        end
        return longest_y_streak
    end

    for key in keys(unscored_PSMs)
        append!(PSMs_dict[:hyperscore], HyperScore(unscored_PSMs[key]))
        append!(PSMs_dict[:error], unscored_PSMs[key].error)
        append!(PSMs_dict[:y_ladder], getLongestY(unscored_PSMs[key].y_ions))
        append!(PSMs_dict[:scan_idx], scan_idx)
        append!(PSMs_dict[:precursor_idx], unscored_PSMs[key].precursor_idx)
    end
end
