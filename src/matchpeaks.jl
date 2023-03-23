"""
    FragmentMatch

Type that represents a match between a fragment ion and a mass spectrum peak

### Fields

- transition::Transition -- Represents the fragment ion
- intensity::Float32 -- Intensity of the matching peak
- match_mz::Float32 -- M/Z of the matching empirical peak. 
    NOT the transition mass and may differ from getMZ(transition) by some error
- count::UInt8 -- Number of matches (may want to count the number of matches if the spectrum is researched at 
    different mass offsets, such as in a cross correlation score)
- peak_int::Int64 -- Index of the matching peak in the mass spectrum. 

### Examples

- `FragmentMatch(transition::Transition, intensity::Float32, mass::Float32, count::UInt8, peak_ind::Int64) --
     default internal constructor
- `FragmentMatch()` -- constructor for null/empty precursor

### GetterMethods

- getMZ(f::FragmentMatch) = getMZ(f.transition)
- getLow(f::FragmentMatch) = getLow(f.transition)
- getHigh(f::FragmentMatch) = getHigh(f.transition)
- getPrecID(f::FragmentMatch) = getPrecID(f.transition)
- getCharge(f::FragmentMatch) = getCharge(f.transition)
- getIsotope(f::FragmentMatch) = getIsotope(f.transition)
- getIonType(f::FragmentMatch) = getIonType(f.transition)
- getInd(f::FragmentMatch) = getInd(f.transition)
"""
mutable struct FragmentMatch
    transition::Transition
    intensity::Float32
    match_mz::Float32
    count::UInt8
    peak_ind::Int64
end

FragmentMatch() = FragmentMatch(Transition(), Float32(0), Float32(0), UInt8(0), 0)
getMZ(f::FragmentMatch) = getMZ(f.transition)
getLow(f::FragmentMatch) = getLow(f.transition)
getHigh(f::FragmentMatch) = getHigh(f.transition)
getPrecID(f::FragmentMatch) = getPrecID(f.transition)
getCharge(f::FragmentMatch) = getCharge(f.transition)
getIsotope(f::FragmentMatch) = getIsotope(f.transition)
getIonType(f::FragmentMatch) = getIonType(f.transition)
getInd(f::FragmentMatch) = getInd(f.transition)

"""
    getNearest(transition::Transition, masses::Vector{Union{Missing, Float32}}, peak::Int; δ = 0.01)

Get the mz ratio of an ion

### Input

- `residues::Vector{Residue}`: -- List of amino acid residues in the peptide ion
- `charge::UInt8` -- Charge of the ion
- `modifier::Float32=PROTON*charge + H2O` -- Added to the mass of the ion
- `isotope::UInt8=UInt8(0)` -- Diference in the number of isotopes from the monoisotopic ion. 

### Output

A Float32 representing the mass-to-charge ratio (m/z) of an ion

### Notes

The `modifier` argument ought to depend on the kind of ion. For B ions PROTON*charge is appropriate,
but for 'y' or precursor ions, PROTON*charge + H2O would be appropriate.

### Algorithm 

Sum the amino acid residue masses, add `modifier` + isotope*NEUTRON and then divide the total by the charge. 

### Examples 

#Gets the b6+1 ion MZ
```julia-repl
julia> getIonMZ(getResidues("PEPTIDE")[1:6], UInt8(1), modifier = PROTON)
653.314f0
```
#Gets the y6+1 ion MZ
```julia-repl
julia> getIonMZ(reverse(getResidues("PEPTIDE"))[1:6], UInt8(1))
703.3144f0
```

"""
function getNearest(transition::Transition, masses::Vector{Union{Missing, Float32}}, peak::Int; δ = 0.01)
    smallest_diff = abs(masses[peak]+ δ - getMZ(transition))
    best_peak = peak
    i = 0
    if peak + 1 + i > length(masses)
        return best_peak
    end
    while (masses[peak + 1 + i]+ δ <= getHigh(transition))
        new_diff = abs(masses[peak  + 1 + i] + δ - getMZ(transition))
        if new_diff < smallest_diff
            smallest_diff = new_diff
            best_peak = peak + 1 + i
        end
        i+=1
        if peak + 1 + i > length(masses)
            break
        end
    end
    best_peak
end

function setFragmentMatch!(hits::Vector{FragmentMatch}, match::Int, transition::Transition, mass::Float32, intensity::Float32, peak_ind::Int64)
    function updateFragmentMatch!(match::FragmentMatch, mass::Float32, intensity::Float32, peak_ind::Int64)
        match.intensity += intensity
        match.mass = mass
        match.count += 1
        match.peak_ind = peak_ind
    end
    if isassigned(hits, match)
        updateFragmentMatch!(hits[match], mass, intensity, peak_ind)
    else
        push!(hits, FragmentMatch(transition, intensity, mass, UInt8(1), peak_ind))
    end
end

function matchPeaks!(matches::Vector{FragmentMatch}, Transitions::Vector{Transition}, masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}}, δ::Float64)
    peak, transition, match = 1, 1, 1
    while (peak <= length(masses)) & (transition <= length(Transitions))
        #Is the peak within the tolerance of the transition m/z?
        if (masses[peak]+ δ >= getLow(Transitions[transition]))
            if (masses[peak]+ δ <= getHigh(Transitions[transition]))
                best_peak = getNearest(Transitions[transition], masses, peak, δ=δ)
                setFragmentMatch!(matches, match, Transitions[transition], masses[best_peak], intensities[best_peak], best_peak);
                transition += 1
                match += 1
                continue
            end
            transition += 1
            continue
        end
        peak+=1
    end
end

function matchPeaks(Transitions::Vector{Transition}, masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}}; δs::Vector{Float64} = [0.0])
    hits = Vector{FragmentMatch}()
    for δ in δs
        matchPeaks!(matches, Transitions, masses, intensities, δ)
    end
    hits
end

abstract type Feature end

mutable struct FastXTandem <: Feature
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
function Score(score::FastXTandem)

    function logfac(N)
        N*log(N) - N + (log(N*(1 + 4*N*(1 + 2*N))))/6 + log(π)/2
    end
    
    (abs(logfac(max(1, score.b_count))) + 
     abs(logfac(max(1, score.y_count))) + 
     log(score.y_int*score.b_int)
    )
end

function logfac(N)
    N*log(N) - N + (log(N*(1 + 4*N*(1 + 2*N))))/6 + log(π)/2
end
