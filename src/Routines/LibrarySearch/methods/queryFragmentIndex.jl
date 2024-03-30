function exponentialFragmentBinSearch(frag_index_bins::AbstractArray{FragIndexBin},
                                        frag_bin_max_idx::UInt32,
                                        lower_bound_guess::UInt32,
                                        upper_bound_guess::UInt32,
                                        frag_mz_max::Float32,
                                        frag_mz_absolute_min::Float32,
                                        step_size::UInt32)
    
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
            if (getHigh(frag_index_bins[lower_bound_guess]) < frag_mz_absolute_min)
                return upper_bound_guess, upper_bound_guess
            end

            break
        end
        n += one(UInt8)
    end
    #If this condition is met, it is possible for the next observed 
    #fragment to match to fragment bins below `lower_bound_guess`
    #If this is the case, use the original lower bound. 
    n = one(UInt8)
    while (getHigh(frag_index_bins[lower_bound_guess]) > frag_mz_absolute_min)
        step = step >> one(UInt8)
        if (step > lower_bound_guess) | (step === zero(UInt32))
            return initial_lower_bound_guess, upper_bound_guess
        end
        lower_bound_guess -= step 
        n += one(UInt8)
    end
    #println("lower_bound_guess $lower_bound_guess upper_bound_guess $upper_bound_guess")
    return lower_bound_guess, upper_bound_guess
end

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
                        frag_mz_absolute_min::Float32,
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
        frag_mz_max,
        frag_mz_absolute_min,
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
                return lower_bound_guess, upper_bound_guess
            else
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
                    ppm_err::Float32, 
                    mass_err_model::MassErrorModel{Float32},
                    min_max_ppm::Tuple{Float32, Float32},
                    prec_mz::Float32, 
                    prec_tol::Float32, 
                    isotope_err_bounds::Tuple{Int64, Int64},
                    min_score::UInt8) where {U<:AbstractFloat}
    
    function getFragTol(mass::Float32, 
                        ppm_err::Float32, 
                        intensity::Float32, 
                        mass_err_model::MassErrorModel{Float32},
                        min_max_ppm::Tuple{Float32, Float32})
        mass -= Float32(ppm_err*mass/1e6)
        ppm = mass_err_model(intensity)
        ppm = max(
                    min(ppm, last(min_max_ppm)), 
                    first(min_max_ppm)
                    )
        tol = ppm*mass/1e6
        return Float32(mass - last(min_max_ppm)*mass/1e6), Float32(mass - tol), Float32(mass + tol)
    end

    function filterPrecursorMatches!(prec_id_to_score::Counter{UInt32, UInt8}, min_score::UInt8)
        match_count = countFragMatches(prec_id_to_score, min_score)
        prec_count = getSize(prec_id_to_score) - 1
        sortCounter!(prec_id_to_score);
        return match_count, prec_count
    end

    prec_min = Float32(prec_mz - prec_tol - NEUTRON*first(isotope_err_bounds)/2)
    prec_max = Float32(prec_mz + prec_tol + NEUTRON*last(isotope_err_bounds)/2)
    @inbounds @fastmath while getLow(rt_bins[rt_bin_idx]) < irt_high

        #BinRanges
        sub_bin_range = getSubBinRange(rt_bins[rt_bin_idx])
        min_frag_bin, max_frag_bin = first(sub_bin_range), last(sub_bin_range)
        lower_bound_guess, upper_bound_guess = min_frag_bin, min_frag_bin

        for (mass, intensity) in zip(masses, intensities)

            #Get intensity dependent fragment tolerance.
            frag_absolute_min, frag_min, frag_max = getFragTol(mass, ppm_err, intensity, mass_err_model, min_max_ppm)

            #For every precursor that could have produced the observed ion
            #award to it the corresponding score
            lower_bound_guess, upper_bound_guess = queryFragment!(prec_id_to_score, 
                                            max_frag_bin,
                                            lower_bound_guess,
                                            upper_bound_guess,
                                            frag_bins,
                                            fragments,
                                            frag_absolute_min,
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

    return filterPrecursorMatches!(prec_id_to_score, min_score)
end

