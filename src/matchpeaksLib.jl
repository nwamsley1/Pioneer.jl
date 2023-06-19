using Dictionaries
"""
    FragmentMatch

Type that represents a match between a fragment ion and a mass spectrum peak

### Fields

- transition::Transition -- Represents a fragment ion
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
mutable struct FragmentMatch{T<:AbstractFloat}
    transition::Transition
    intensity::T
    match_mz::T
    count::UInt8
    peak_ind::Int64
    scan_idx::UInt32
    ms_file_idx::UInt32
end
export FragmentMatch

FragmentMatch() = FragmentMatch(Transition(), Float64(0), Float64(0), UInt8(0), 0, UInt32(0), UInt32(0))
getMZ(f::FragmentMatch) = getMZ(f.transition)
getLow(f::FragmentMatch) = getLow(f.transition)
getHigh(f::FragmentMatch) = getHigh(f.transition)
getPrecID(f::FragmentMatch) = getPrecID(f.transition)
getCharge(f::FragmentMatch) = getCharge(f.transition)
getIsotope(f::FragmentMatch) = getIsotope(f.transition)
getIonType(f::FragmentMatch) = getIonType(f.transition)
getInd(f::FragmentMatch) = getInd(f.transition)
getPeakInd(f::FragmentMatch) = f.peak_ind
getIntensity(f::FragmentMatch) = f.intensity
getCount(f::FragmentMatch) = f.count
getMSFileID(f::FragmentMatch) = f.ms_file_idx

export getMZ, getLow, getHigh, getPrecID, getCharge, getIsotope
export getIonType, getInd, getPeakInd, getIntensity, getCount, getMSFileID

"""
    getNearest(transition::Transition, masses::Vector{Union{Missing, Float32}}, peak::Int; δ = 0.01)

Finds the `peak` (index of `masses`) nearest in mass to the `transition` but still within the tolerance (getLow(transition)<masses[peak]<getHigh(transition)). 
Starts searching at initially supplied `peak` and increases the index until outside the tolerance. There could be multiple peaks within the tolerance, 
and this function selects the one with the lowest mass error to the fragment ion. 

### Input

- `transition::Transition`: -- Represents a fragment ion
- `masses::Vector{Union{Missing, Float32}}` -- Mass list from a centroided mass spectrum. MUST BE SORTED IN ASCENDING ORDER. 
- `peak::Int` -- An index for a mass in `masses`
- `δ` -- A mass offset that can be applied to each mass in `masses` 

### Output

An Int representing the index of the m/z in `masses` nearest in m/z to that of the fragment m/z. 

### Notes

- It is assumed that masses is sorted in ascending order. 
- In practice, when called from `matchPeaks!`, masses[peak] will already be within the fragment m/z tolerance.
 Usually there will not be another peak in `masses` where this is true, but it is possible for multiple peaks to
 fall within the tolerance. The purpose of this function is to select the best peak (closes in mass) when this happens.  

### Algorithm 

### Examples 

"""
function getNearest(masses::Vector{Union{Missing, T}}, frag_mz::U, mz_high::U, peak::Int; δ::U = 0.0) where {T,U<:AbstractFloat}
    smallest_diff = abs(masses[peak]+ δ - frag_mz)
    best_peak = peak
    i = 0

    #Shorthand to check for BoundsError on `masses`
    boundsCheck() = peak + 1 + i > length(masses)
    if boundsCheck() return best_peak end
    #Iterate through peaks in  `masses` until a peak is encountered that is 
    #greater in m/z than the upper bound of the `transition` tolerance 
    #or until there are no more masses to check. Keep track of the best peak/transition match. 
    while (masses[peak + 1 + i]+ δ <= mz_high)
        new_diff = abs(masses[peak  + 1 + i] + δ - frag_mz)
        if new_diff < smallest_diff
            smallest_diff = new_diff
            best_peak = peak + 1 + i
        end
        i+=1
        if boundsCheck() break end
    end
    best_peak
end
export getNearest

"""
    setFragmentMatch!(matches::Vector{FragmentMatch}, match::Int, transition::Transition, mass::Float32, intensity::Float32, peak_ind::Int64)  

Adds a FragmentMatch to `matches` if `match` is not an index in `matches`, otherwise, updates the match.

### Input

- `hits::Vector{FragmentMatch}`: -- Represents a fragment ion
- `match::Int` -- Index of the match. Must be <=N+1 where N is length(hits) 
- `mass::Float32` -- m/z of the emperical peak matched to the transition
- `intensity::Float32` -- intensity of the emperical peak matched to the transition
- `peak_ind` -- unique index of the emperical peak matched to the transition

### Output

Modifies `matches[match]` if match is <= lenth(matches). Otherwise adds a new FragmentMatch at `matches[match]`

### Notes

- Updating a match in `matches` could be useful if researching the same spectra many times at different mass offsets with `matchPeaks!` 
    This could be done to calculate a cross correlatoin score for example. If a spectrum is only searched once, then matches should
    only be added to `matches` and existing ones never modified. 
- Recording the `peak_ind` could be useful, for example, "chimeric" scoring of a spectrum from a Vector{FragmentMatch} type. The best scoring precursor would
    have fragments matching to known `peak_ind`. The Vector{FragmentMatch} could be rescored but excluding FragmentMatches corresponding
    to those peak_ind's. This would enable a simple chimeric spectra scoring that would not involve completely researching the spectrum (making an additional call to `matchPeaks`). 

### Algorithm 

### Examples 

"""
function setFragmentMatch!(matches::Vector{FragmentMatch{T}}, transition::LibraryFragment{Float64}, mass::T, intensity::T, peak_ind::Int64, scan_idx::UInt32, ms_file_idx::UInt32) where {T<:AbstractFloat}
    #function updateFragmentMatch!(match::FragmentMatch{T}, mass::T, intensity::T, peak_ind::Int64) where {T<:AbstractFloat}
    #    match.intensity += intensity
    #    match.match_mz = mass
    #    match.count += 1
    #    match.peak_ind = peak_ind
    #end
    #if isassigned(matches, match)
    ##    updateFragmentMatch!(matches[match], mass, intensity, peak_ind)
    #else
    push!(matches, FragmentMatch(getIntensity(transition), 
                                 intensity,
                                 Float32(getFragMZ(transition)),
                                 mass,
                                 peak_ind,
                                 getIonIndex(transition),
                                 getFragCharge(transition),
                                 UInt8(0),
                                 true == isyIon(transition) ? 'y' : 'b',
                                 getPrecID(transition),
                                 UInt8(1), 
                                 scan_idx,
                                 ms_file_idx)
        )
    #end
end
export setFragmentMatch!

"""
    function matchPeaks!(matches::Vector{FragmentMatch}, Transitions::Vector{Transition}, masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}}, δ::Float64)  

Finds the best matching peak in a mass spectrum for each transition/fragment ion supplied if the match is within the fragment tolerance. 
    Adds each FragmentMatch to `matches` if not already present. Otherwise, 
    modifies an existing match (see setFragmentMatch!).

### Input

- `matches::Vector{FragmentMatch}`: -- A list representing fragment ions that match peaks in the mass spectrum (`masses`)
- `Transitions::Vector{Transition` -- A list of fragment ions to search for in the spectrum (`masses`). MUST BE SORTED IN ASCENDING ORDER BY `getMZ(transition)`
- `masses::Vector{Union{Missing, Float32}}` -- Mass list from a centroided mass spectrum. MUST BE SORTED IN ASCENDING ORDER.
- `intensities::Vector{Union{Missing, Float32}}` -- The intensity list from a centroided mass spectrum. Must be the same length as `masses`
- `δ::Float64` -- A mass offset that can be applied to each mass in `masses`

### Output

Modifies `matches[match]` if match is <= lenth(matches). Otherwise adds a new FragmentMatch at `matches[match]`

### Notes

- Updating a match in `matches` could be useful if researching the same spectra many times at different mass offsets with `matchPeaks!` 
    This could be done to calculate a cross correlatoin score for example. If a spectrum is only searched once, then matches should
    only be added to `matches` and existing ones never modified. 
- The fragment tolerance is specified by each `Transition`. A `Transition<:Ion` has a field `mz::MzFeature` which specifies the monoisotopic mass
    and also upper and lower bounds (the tolerance). See getLow(ion::Ion) and getHigh(ion::Ion). This is why the user need not supply a fragment tolerance to `matchPeaks!` 
- `masses` and `intensities` contain type unions Union{Missing, Float32}. This method does nothing to check for Missing values, and indeed,
    it is assumed that there are none, and the presence of any Missing values will cause an error. The reason for the type union is an idiosyncracy
    of the Arrow.jl package and how it implements nested data types in Arrow files. 

### Algorithm 

Given a list of fragment ions and a centroided mass spectrum both sorted by m/z, it is efficient to search the spetrum for matches in a "single pass"
through the spectrum. If there are T transitions and P peaks should be O(T+P). If there are multiple peaks within the tolerance for a given 
fragment ion, the peak closest in m/z to the fragment ion is chosen. It is possible to assign the same peak to multiple fragment ions, but 
each fragment ion is only assigned to 0 or 1 peaks. 

### Examples 

"""
function matchPeaks!(matches::Vector{FragmentMatch{T}}, Transitions::Vector{LibraryFragment{Float64}}, masses::Vector{Union{Missing, T}}, intensities::Vector{Union{Missing, T}}, δ::U, scan_idx::UInt32, ms_file_idx::UInt32, min_intensity::T) where {T,U<:AbstractFloat}

    #match is a running count of the number of transitions that have been matched to a peak
    #This is not necessarily the same as `transition` because some (perhaps most)
    #transitions will not match to any peak. 
    peak, transition, match = 1, 1, 1

    function getPPM(transition::LibraryFragment{Float64}, ppm::Float64)
        mz = getFragMZ(transition)
        tol = ppm*mz/1e6
        return mz - tol, mz + tol
    end
    if length(Transitions)<1
        return
    end
    low, high = getPPM(Transitions[transition], 20.0)
    while (peak <= length(masses)) & (transition <= length(Transitions))
        if intensities[peak] <  min_intensity
            peak += 1
            continue
        end
        #Is the peak within the tolerance of the transition m/z?
        if (masses[peak]+ δ >= low)
            if (masses[peak]+ δ <= high)
                #Find the closest matching peak to the transition within the upper and lower bounds (getLow(transition)<=masses[peak]<=getHigh(transition)))
                best_peak = getNearest(masses, getFragMZ(Transitions[transition]), high, peak, δ=δ)
                setFragmentMatch!(matches, Transitions[transition], masses[best_peak], intensities[best_peak], best_peak, scan_idx, ms_file_idx);
                transition += 1
                if transition > length(Transitions)
                    return
                end
                low, high = getPPM(Transitions[transition], 20.0)
                match += 1
                continue
            end
            #Important that this is also within the first if statement. 
            #Need to check the next fragment against the current peak. 
            transition += 1
            if transition > length(Transitions)
                return
            end
            low, high = getPPM(Transitions[transition], 20.0)
            continue
        end
        #No additional matches possible for the current peak. Move on to the next. 
        peak+=1
    end
end
export matchPeaks!

"""
    function matchPeaks!(matches::Vector{FragmentMatch}, Transitions::Vector{Transition}, masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}}, δ::Float64)  

A wrapper for calling `matchPeaks` at different to search spectra at a list of mass offset. 
    Each all to `matchPeaks` Finds the best matching peak in a mass spectrum for each transition/fragment ion supplied if the match is within the fragment tolerance. 
    Adds each FragmentMatch to `matches` if not already present. Otherwise, modifies an existing match (see setFragmentMatch!). (see `matchPeaks` for additional details)

### Input

- `matches::Vector{FragmentMatch}`: -- A list representing fragment ions that match peaks in the mass spectrum (`masses`)
- `Transitions::Vector{Transition` -- A list of fragment ions to search for in the spectrum (`masses`). MUST BE SORTED IN ASCENDING ORDER BY `getMZ(transition)`
- `masses::Vector{Union{Missing, Float32}}` -- Mass list from a centroided mass spectrum. MUST BE SORTED IN ASCENDING ORDER.
- `intensities::Vector{Union{Missing, Float32}}` -- The intensity list from a centroided mass spectrum. Must be the same length as `masses`
- `δ::Float64` -- A mass offset that can be applied to each mass in `masses`

### Output

Modifies `matches[match]` if match is <= lenth(matches). Otherwise adds a new FragmentMatch at `matches[match]`

### Notes

- Searching a mass spectrum many times at different mass offsets could be useful for caculating cross correlation scores. 

### Algorithm 

    See `matchPeaks`

### Examples 

"""
function matchPeaks(Transitions::Vector{LibraryFragment{Float64}}, masses::Vector{Union{Missing, T}}, intensities::Vector{Union{Missing, T}}; δs::Vector{U} = zeros(Float32, (1, )), scan_idx = UInt32(0), ms_file_idx = UInt32(0), min_intensity::Float32 = Float32(0.0)) where {T,U<:AbstractFloat}
    matches = Vector{FragmentMatch{T}}()
    for δ in δs
        matchPeaks!(matches, Transitions, masses, intensities, δ, scan_idx, ms_file_idx, min_intensity)
    end
    matches
end

export FragmentMatch, getNearest, matchPeaks, matchPeaks!