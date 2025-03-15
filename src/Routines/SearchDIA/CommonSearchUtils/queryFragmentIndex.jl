"""
    exponentialFragmentBinSearch(frag_index_bins::AbstractArray{FragIndexBin},
                               frag_bin_max_idx::UInt32,
                               lower_bound_guess::UInt32,
                               upper_bound_guess::UInt32,
                               frag_mz_min::Float32,
                               frag_mz_max::Float32,
                               step_size::UInt32) -> Tuple{UInt32, UInt32}

Find new lower and upper bounds for fragment bin search using exponential search.

# Arguments
- `frag_index_bins`: Array of fragment index bins
- `frag_bin_max_idx`: Maximum index in fragment bins
- `lower_bound_guess`: Initial lower bound guess
- `upper_bound_guess`: Initial upper bound guess
- `frag_mz_min`: Minimum fragment m/z value
- `frag_mz_max`: Maximum fragment m/z value
- `step_size`: Initial step size for exponential search

# Process
1. Exponentially increases upper bound until fragment m/z range is contained
2. Adjusts lower bound downward if needed using binary reduction
3. Handles edge cases for array bounds

# Returns
Tuple of (new_lower_bound, new_upper_bound) for fragment bin search
"""
function exponentialFragmentBinSearch(frag_index_bins::AbstractArray{FragIndexBin},
                                        frag_bin_max_idx::UInt32,
                                        lower_bound_guess::UInt32,
                                        upper_bound_guess::UInt32,
                                        frag_mz_min::Float32,
                                        frag_mz_max::Float32,
                                        step_size::UInt32)
    #return lower_bound_guess, upper_bound_guess
    initial_lower_bound_guess = lower_bound_guess
    n = zero(UInt8)
    step = one(UInt32)
    #exponential search for new upper and lower bounds 
    while (getHigh(frag_index_bins[upper_bound_guess]) < frag_mz_max) #Need a new upper and lower bound guess 
        lower_bound_guess = upper_bound_guess 
        step = step_size << n 
        upper_bound_guess += step #Exponentially increasing guess 
        if upper_bound_guess > frag_bin_max_idx #If guess exceeds limits
            upper_bound_guess = frag_bin_max_idx #then set to maximum 
            #This was a mistake 
            #if (getHigh(frag_index_bins[lower_bound_guess]) < frag_mz_min)
            #    return upper_bound_guess, upper_bound_guess
            #end
            break
        end
        n += one(UInt8)
    end
    #If this condition is met, it is possible for the next observed 
    #fragment to match to fragment bins below `lower_bound_guess`
    #If this is the case, use the original lower bound. 
    #=
    n = one(UInt8)
    while (getHigh(frag_index_bins[lower_bound_guess]) > frag_mz_min)
        step = step >> one(UInt8)
        #if (step === zero(UInt32))# (step > lower_bound_guess) | (step === zero(UInt32))
        #    if (step > lower_bound_guess)
        #        return initial_lower_bound_guess, upper_bound_guess
        #    else
        #        return lower_bound_guess, upper_bound_guess
        #    end
        #end
        if (step > lower_bound_guess) | (step === zero(UInt32))
            return initial_lower_bound_guess, upper_bound_guess
        end
        lower_bound_guess -= step 
        n += one(UInt8)
    end
    =#
    #println("lower_bound_guess $lower_bound_guess upper_bound_guess $upper_bound_guess")
    return lower_bound_guess, upper_bound_guess
end

"""
    findFirstFragmentBin(frag_index_bins::AbstractArray{FragIndexBin},
                        lower_bound_guess::UInt32,
                        upper_bound_guess::UInt32,
                        frag_min::Float32) -> UInt32

Perform branchless binary search to find first fragment bin that could contain frag_min.

# Arguments
- `frag_index_bins`: Array of fragment index bins
- `lower_bound_guess`: Lower bound for search range
- `upper_bound_guess`: Upper bound for search range
- `frag_min`: Minimum fragment m/z value to find

# Returns
Index of first fragment bin that could contain frag_min

Uses branchless binary search for performance optimization.
"""
function findFirstFragmentBin(frag_index_bins::AbstractArray{FragIndexBin},
                                lower_bound_guess::UInt32,
                                upper_bound_guess::UInt32,
                                frag_min::Float32) #where {T<:AbstractFloat}
    #branchless binary search
    @inbounds @fastmath begin
        len = upper_bound_guess - lower_bound_guess + UInt32(1)
        mid = len>>>0x01
        base = lower_bound_guess
        while len > 1
            base += (getHigh(frag_index_bins[base + mid - UInt32(1)]) < frag_min)*mid
            len -= mid
            mid = len>>>0x01
        end
    end
    return base
end


"""
    searchFragmentBin!(prec_id_to_score::Counter{UInt32, UInt8},
                      fragments::AbstractArray{IndexFragment},
                      frag_id_range::UnitRange{UInt32},
                      window_min::Float32,
                      window_max::Float32)

Search a fragment bin for matches within precursor mass window.

# Arguments
- `prec_id_to_score`: Counter tracking scores per precursor
- `fragments`: Array of indexed fragments
- `frag_id_range`: Range of fragment indices to search
- `window_min`: Minimum precursor m/z window
- `window_max`: Maximum precursor m/z window

# Process
1. Binary search to find first fragment in precursor window
2. Binary search to find last fragment in precursor window
3. Adds match scores for all qualifying fragments

Updates prec_id_to_score in place with new matches.
"""
function searchFragmentBin!(prec_id_to_score::Counter{UInt32, UInt8}, 
                            fragments::AbstractArray{IndexFragment},
                            frag_id_range::UnitRange{UInt32},
                            window_min::Float32, 
                            window_max::Float32)

    #Index of first and last fragments to search 
    lo, hi = first(frag_id_range), last(frag_id_range)
    base = lo
    @inbounds @fastmath begin 

        #Binary search to find the first fragment matching the precursor tolerance
        len = hi - lo + UInt32(1)
        while len > 1
            mid = len >>> 0x01
            base += (getPrecMZ(fragments[base + mid - one(UInt32)]) < window_min)*mid
            len -= mid
        end
        window_start = base

        #Best initial guess is that the next fragment
        #is outside the precursor tolerance. 
        init_guess = min(window_start + one(UInt32),hi) #Get next fragment
        tf = (getPrecMZ(fragments[init_guess])>window_max) #Is the fragment outside the precursor mz tolerance
        hi = tf*init_guess + (one(UInt32)-tf)*hi #Choose the appropriate upper bound

        #Binary search to find the last fragment matcing the precursor tolerance
        len = hi - base  + one(UInt32)
        base = hi
        while len > 1
            mid = len>>>0x01
            base -= (getPrecMZ( fragments[base - mid + one(UInt32)]) > window_max)*mid
            len -= mid
        end
        window_stop = base

        #If these conditions are met, there is no fragment in the query 
        #range that matches the precursor tolerance. 
        if window_start === window_stop
            if (getPrecMZ(fragments[window_start])>window_max)
                return nothing
            end
            if getPrecMZ(fragments[window_stop])<window_min
                return nothing
            end
        end
    end
        
    function addFragmentMatches!(prec_id_to_score::Counter{UInt32, UInt8}, 
                                    fragments::AbstractArray{IndexFragment},
                                    matched_frag_range::UnitRange{UInt32})
        @inline @inbounds for i in matched_frag_range
            frag = fragments[i]
            inc!(prec_id_to_score, getPrecID(frag), getScore(frag))
        end
    end
    
    #For each fragment matching the query, 
    #award its score to its parent ion. 
    @inline addFragmentMatches!(prec_id_to_score, fragments, window_start:window_stop)
    return nothing

end

function queryFragment!(prec_id_to_score::Counter{UInt32, UInt8}, 
                        frag_bin_max_idx::UInt32,
                        lower_bound_guess::UInt32,
                        upper_bound_guess::UInt32,
                        frag_bins::AbstractArray{FragIndexBin},
                        fragments::AbstractArray{IndexFragment},
                        frag_mz_min::Float32, 
                        frag_mz_max::Float32, 
                        prec_mz_min::Float32,
                        prec_mz_max::Float32)# where {T,U<:AbstractFloat}
    #Get new lower and upper bounds for the fragment bin search if necessary 
    lower_bound_guess, upper_bound_guess = exponentialFragmentBinSearch(
        frag_bins,
        frag_bin_max_idx,
        lower_bound_guess,
        upper_bound_guess,
        frag_mz_min,
        frag_mz_max,
        UInt32(2048) #step size
    )
    #First frag_bin matching fragment tolerance
    frag_bin_idx = findFirstFragmentBin(
                                    frag_bins, 
                                    lower_bound_guess,
                                    upper_bound_guess,
                                    frag_mz_min
                                    )
    @inbounds @fastmath begin 
        #No fragment bins contain the fragment m/z
        if iszero(frag_bin_idx)
            return lower_bound_guess, upper_bound_guess
        end

        #Search subsequent frag bins until no more bins or untill a bin is outside the fragment tolerance
        while (frag_bin_idx <= frag_bin_max_idx)
            #Fragment bin is outside the fragment tolerance
            frag_bin = frag_bins[frag_bin_idx]#getFragmentBin(frag_index, frag_bin_idx)
            #This and all subsequent fragment bins cannot match the fragment,
            #so exit the loop 
            if (getLow(frag_bin) > frag_mz_max)
                break
            else
                if frag_bin_max_idx === frag_bin_idx
                    if getHigh(frag_bin) < frag_mz_min
                        break
                    end
                end
                #Range of fragment ions that could match the observed fragment tolerance 
                frag_id_range = getSubBinRange(frag_bin)
                #Search the fragment range for fragments with precursors in the precursor tolerance
                searchFragmentBin!(prec_id_to_score, 
                                    fragments,
                                    frag_id_range, 
                                    prec_mz_min, 
                                    prec_mz_max
                                )
                #Advance to the next fragment bin 
                frag_bin_idx += 1
            end
        end
    end

    #Only reach this point if frag_bin exceeds length(frag_index)
    return lower_bound_guess, upper_bound_guess
end

function searchScan!(prec_id_to_score::Counter{UInt32, UInt8}, 
                    rt_bins::AbstractArray{FragIndexBin},
                    frag_bins::AbstractArray{FragIndexBin},
                    fragments::AbstractArray{IndexFragment},
                    masses::AbstractArray{Union{Missing, U}}, 
                    intensities::AbstractArray{Union{Missing, U}}, 
                    rt_bin_idx::Int64,
                    irt_high::Float32, 
                    mass_err_model::MassErrorModel,
                    quad_transmission_func::QuadTransmissionFunction,
                    isotope_err_bounds::Tuple{UInt8, UInt8}
                    ) where {U<:AbstractFloat}
    prec_min = U(getPrecMinBound(quad_transmission_func) - NEUTRON*first(isotope_err_bounds)/2)
    prec_max = U(getPrecMaxBound(quad_transmission_func) + NEUTRON*last(isotope_err_bounds)/2)

    #@inbounds @fastmath while getLow(rt_bins[rt_bin_idx]) < irt_high
    while getLow(rt_bins[rt_bin_idx]) < irt_high
        #BinRanges
        sub_bin_range = getSubBinRange(rt_bins[rt_bin_idx])
        min_frag_bin, max_frag_bin = first(sub_bin_range), last(sub_bin_range)
        lower_bound_guess, upper_bound_guess = min_frag_bin, min_frag_bin

        for mass in masses
            #Get intensity dependent fragment tolerance.
            corrected_mz = getCorrectedMz(mass_err_model, mass)
            frag_min, frag_max = getMzBoundsReverse(mass_err_model, corrected_mz)
            #For every precursor that could have produced the observed ion
            #award to it the corresponding score
            lower_bound_guess, upper_bound_guess = queryFragment!(prec_id_to_score, 
                                            max_frag_bin,
                                            lower_bound_guess,
                                            upper_bound_guess,
                                            frag_bins,
                                            fragments,
                                            frag_min, 
                                            frag_max, 
                                            prec_min, 
                                            prec_max
                                        )
        end

        rt_bin_idx += 1
        if rt_bin_idx >length(rt_bins)
            rt_bin_idx = length(rt_bins)
            break
        end 
    end

    return nothing#filterPrecursorMatches!(prec_id_to_score, min_score)
end

function filterPrecursorMatches!(prec_id_to_score::Counter{UInt32, UInt8}, min_score::UInt8)
    match_count = countFragMatches(prec_id_to_score, min_score)
    prec_count = getSize(prec_id_to_score) - 1
    sortCounter!(prec_id_to_score);
    return match_count, prec_count
end
