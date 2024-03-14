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
function setMatch!(matches::Vector{M}, 
                    transition::LibraryIon{Float32}, 
                    mass::T, 
                    intensity::T, 
                    peak_ind::Int64, 
                    scan_idx::UInt32, 
                    ms_file_idx::UInt32,
                    i::Int64; 
                    block_size = 10000) where {T<:AbstractFloat,M<:MatchIon{Float32}}
    i += 1
    if i > length(matches)
        append!(matches, [FragmentMatch{Float32}() for _ in range(1, block_size)])
    end

    matches[i] = FragmentMatch(Float32(getIntensity(transition))::Float32, 
                                 intensity::Float32,
                                 Float32(getMZ(transition)),
                                 mass::Float32,
                                 peak_ind::Int64,
                                 getIonPosition(transition)::UInt8,
                                 getFragCharge(transition)::UInt8,
                                 UInt8(0),
                                 true == isyIon(transition) ? 'y' : 'b',
                                 transition.is_isotope,
                                 getPrecID(transition),
                                 UInt8(1), 
                                 scan_idx,
                                 ms_file_idx,
                                 getRank(transition))::FragmentMatch{Float32}
    return i
end

"""
    getNearest(transition::Transition, masses::Vector{Union{Missing, Float32}}, peak::Int; δ = 0.01)

Finds the `peak` (index of `masses`) nearest in mass to the `transition` but still within the tolerance (getLow(transition)<masses[peak]<getHigh(transition)). 
Starts searching at initially supplied `peak` and increases the index until outside the tolerance. There could be multiple peaks within the tolerance, 
and this function selects the one with the lowest mass error to the fragment ion. 

### Input

- `transition::Transition`: -- Represents a fragment ion
- `masses::Vector{Union{Missing, Float32}}` -- Mass list from a centroided mass spectrum. MUST BE SORTED IN ASCENDING ORDER. 
- `peak::Int` -- An index for a mass in `masses`
- `δ` 
-- A mass offset that can be applied to each mass in `masses` 

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
function setNearest!(matches::Vector{M}, 
                     unmatched::Vector{M},
                    masses::AbstractArray{Union{Missing, Float32}}, 
                    intensities::AbstractArray{Union{Missing, Float32}},
                    Ion::DetailedFrag{Float32},
                    max_mz::Float32,
                    mass_err_model::MassErrorModel{Float32},
                    min_max_ppm::Tuple{Float32, Float32},
                    δ::Float32,
                    peak_idx::Int64,
                    scan_idx::UInt32,
                    ms_file_idx::UInt32,
                    matched_idx::Int64,
                    unmatched_idx::Int64
                    )::Tuple{Int64,Int64} where {M<:MatchIon{Float32}}

    #Get maximum and minimum m/z of a theoretical ion that could match a peak given 
    #the intensity. (tolerance is inversely proportional to the sqrt of the intensity)       
    function getMzBounds(mz::Float32, mass_err_model::MassErrorModel{Float32}, intensity::Float32, min_max_ppm::Tuple{Float32, Float32})
        ppm = mass_err_model(intensity)
        ppm = max(
            min(ppm, last(min_max_ppm)), 
            first(min_max_ppm)
            )
        #ppm = Float32(16.1)
        #ppm = max(ppm_tol_param/sqrt(intensity), 5.0f0)
        tol = ppm*mz/Float32(1e6)
        return mz - tol, mz + tol        
    end

    smallest_diff = Float32(Inf)#abs(masses[peak]-δ - frag_mz)
    best_peak, i = 0, 0
    theoretical_mz = getMZ(Ion)
    #Iterate through peaks in  `masses` until a peak is encountered that is 
    #greater in m/z than the upper bound of the `transition` tolerance 
    #or until there are no more masses to check. Keep track of the best peak/transition match. 
    @inbounds while (masses[peak_idx + i]-δ<= max_mz)
        low, high = getMzBounds(theoretical_mz, mass_err_model, intensities[peak_idx + i], min_max_ppm)
        if (masses[peak_idx + i]-δ> low) & (masses[peak_idx + i]-δ< high)
            mz_diff = abs(masses[peak_idx + i]-δ-theoretical_mz)
            if mz_diff < smallest_diff
                smallest_diff = mz_diff
                best_peak = peak_idx + i
            end
        end
        i+=1
        if (peak_idx + 1 + i > length(masses)) break end
    end

    if best_peak > 0
        matched_idx = setMatch!(matches, 
                                Ion, 
                                masses[best_peak]-δ, 
                                intensities[best_peak], 
                                best_peak, 
                                scan_idx, 
                                ms_file_idx,
                                matched_idx);
    else
        unmatched_idx = setMatch!(unmatched,
                    Ion, 
                    zero(Float32), 
                    zero(Float32),
                    unmatched_idx, 
                    scan_idx, 
                    ms_file_idx,
                    unmatched_idx
                    );
    end
    return matched_idx, unmatched_idx
end

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
function matchPeaks!(matches::Vector{M}, #Pre-allocated container for Matched Ions
                    unmatched::Vector{M}, #Pre-allocated container for Unmatched Ions
                    Ions::Vector{I}, #Library Ion Templates to match to the spectrum
                    ion_idx::Int64, #Maximum elemtn in Ions to consider
                    masses::AbstractArray{Union{Missing, Float32}}, intensities::AbstractArray{Union{Missing, Float32}}, 
                    ppm_err::Float32, 
                    mass_err_model::MassErrorModel{Float32},
                    min_max_ppm::Tuple{Float32, Float32},
                    scan_idx::UInt32, 
                    ms_file_idx::UInt32
                    ) where {I<:LibraryIon{Float32},M<:MatchIon{Float32}}

    #matched_idx is a running count of the number of transitions that have been matched to a peak
    #peak is the index of the peak in the spectrum masses[peak]/intensities[peak]
    #ion is the index of the template being compared to the spectrum Ions[ion]
    peak, ion, matched_idx, unmatched_idx = 1, 1, 0, 0

    #Maximum m/z match window for a theoretical fragment ion. Returns a tuple of the minimum and 
    #maximum posible m/z's of a matching peak
    function getMaxErrBounds(theoretical_mz::Float32, max_ppm_tol::Float32)
        mz_tol = Float32(max_ppm_tol*theoretical_mz/1e6)
        return theoretical_mz - mz_tol, theoretical_mz + mz_tol
    end

    if ion_idx<1
        return matched_idx, unmatched_idx
    end

    δ = Float32(ppm_err*(masses[peak]/1e6)) #ppm offset determined in calibration
    corrected_empirical_mz = masses[peak] - δ #corrected empirical mass
    low, high = getMaxErrBounds(getMZ(Ions[ion]), last(min_max_ppm)) #max and min m/z of a matching empirical ion

    while (peak <= length(masses)) & (ion <= ion_idx)
        #Is the peak within the tolerance of the transition m/z?
        if (corrected_empirical_mz>=low) 
            if (corrected_empirical_mz<=high) #Empirical peak is within maximum mass tolerance
                #Set nearest matching peak if there are any.
                #Mass tolerance is intensity dependent,
                #so is it possible for there to be none.
                #In that case, sets an unmatched ion. 
                matched_idx, unmatched_idx = setNearest!(matches,
                                            unmatched,
                                            masses, 
                                            intensities, 
                                            Ions[ion], 
                                            high, 
                                            mass_err_model,
                                            min_max_ppm,
                                            δ,
                                            peak, 
                                            scan_idx,
                                            ms_file_idx,
                                            matched_idx,
                                            unmatched_idx) #Set nearest matching peak if there are any 

                #Progress to the next theoretical ion
                ion += 1
                low, high = getMaxErrBounds(getMZ(Ions[ion]), last(min_max_ppm))
                if ion > ion_idx
                    return matched_idx, unmatched_idx
                end
                continue
            end
            #Progress to next theoretical ion
            unmatched_idx = setMatch!(unmatched,
                                        Ions[ion], 
                                        zero(Float32), 
                                        zero(Float32),
                                        unmatched_idx, 
                                        scan_idx, 
                                        ms_file_idx,
                                        unmatched_idx
                                        );
            ion += 1
            low, high = getMaxErrBounds(getMZ(Ions[ion]), last(min_max_ppm))
            if ion > ion_idx
                return matched_idx, unmatched_idx
            end
            continue
        end
        
        peak += 1 #Progress to next peak
        if peak <= length(masses)
            δ = Float32(ppm_err*(masses[peak]/1e6)) 
            corrected_empirical_mz = masses[peak] - δ
        end
    end

    #Remaining templates with higher m/z than the highest
    #m/z peak in the spectrum are written to 'unmatched'
    while ion <= ion_idx#length(Ions)
        unmatched_idx = setMatch!(unmatched,
                                    Ions[ion], 
                                    zero(Float32), 
                                    zero(Float32),
                                    unmatched_idx, 
                                    scan_idx, 
                                    ms_file_idx,
                                    unmatched_idx
                                    );
        ion += 1
    end

    return matched_idx, unmatched_idx
end
