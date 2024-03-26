function findFirstFragmentBin(frag_index_bins::Arrow.Struct{FragIndexBin, Tuple{Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{UInt32, Vector{UInt32}}}, (:lb, :ub, :first_bin, :last_bin)}, 
                                frag_bin_max_idx::UInt32,
                                lower_bound_guess::UInt32,
                                upper_bound_guess::UInt32,
                                frag_min::Float32, 
                                frag_max::Float32) #where {T<:AbstractFloat}
    #Binary Search
    #Need to add check for last bin. Right now will always include the last bin?
    @inbounds @fastmath begin
        tf =  (getHigh(frag_index_bins[upper_bound_guess]) < frag_min)
        n = zero(UInt8)
        #exponential search for new upper and lower bounds 
        while tf #Need a new upper and lower bound guess 
            lower_bound_guess = upper_bound_guess 
            upper_bound_guess += UInt32(512) << n#UInt32(2048) << n
            if upper_bound_guess > frag_bin_max_idx
                upper_bound_guess = frag_bin_max_idx
                if (getHigh(frag_index_bins[upper_bound_guess]) < frag_min)
                    return lower_bound_guess, lower_bound_guess, upper_bound_guess
                end
                break
            end
            tf = (getHigh(frag_index_bins[upper_bound_guess]) < frag_min)
            n += one(UInt8)
        end
        len = upper_bound_guess - lower_bound_guess + UInt32(1)
        mid = len>>>0x01
        base = lower_bound_guess
        #n = 0
        while len > 1
            #n += 1
            base += (getHigh(frag_index_bins[base + mid - UInt32(1)]) < frag_min)*mid
            len -= mid
            mid = len>>>0x01
        end
        #println("n $n intit_guess $init_guess tf - $tf : diff - ", base - first(frag_bin_range), "diff -2 ", last(frag_bin_range) - first(frag_bin_range))
    end
    #println("lower_bound_guess $lower_bound_guess base $base upper_bound_guess $upper_bound_guess")
    return base, lower_bound_guess, upper_bound_guess#UInt32(max(potential_match, lo + 1) - lo)#UInt32(2*(hi_f - low_f)÷(mid - low_f))
end

function searchPrecursorBin!(prec_id_to_score::Counter{UInt32, UInt8}, 
                            fragments::Arrow.Struct{IndexFragment, Tuple{Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{UInt8, Vector{UInt8}}, Arrow.Primitive{UInt8, Vector{UInt8}}}, (:prec_id, :prec_mz, :score, :charge)},
                            frag_id_range::UnitRange{UInt32},
                            window_min::Float32, 
                            window_max::Float32)

    lo, hi = first(frag_id_range), last(frag_id_range)
    base = lo
    @inbounds @fastmath begin 
        len = hi - lo + UInt32(1)
        gd = one(UInt32)#min(UInt32(2),len)

        while len > 1
            mid = len >>> 0x01
            base += (getPrecMZ(fragments[base + mid - one(UInt32)]) < window_min)*mid
            len -= mid
        end
        window_start = base

        #Next three lines implement an initial guess strategy. 
        #Worth making a branch to only do this if the len is greater than 2? 
        init_guess = min(window_start + gd,hi) #Could exceed bounds. 
        tf = (getPrecMZ(fragments[init_guess])>window_max)
        hi = tf*init_guess + (one(UInt32)-tf)*hi


        len = hi - base  + one(UInt32)
        base = hi
        while len > 1
            mid = len>>>0x01
            base -= (getPrecMZ( fragments[base - mid + one(UInt32)]) > window_max)*mid
            len -= mid
        end
        window_stop = base

        if window_start === window_stop
            if (getPrecMZ(fragments[window_start])>window_max)
                return 
            end
            if getPrecMZ(fragments[window_stop])<window_min
                return 
            end
        end
    end
        
    function addFragmentMatches!(prec_id_to_score::Counter{UInt32, UInt8}, 
                                    fragments::Arrow.Struct{IndexFragment, Tuple{Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{UInt8, Vector{UInt8}}, Arrow.Primitive{UInt8, Vector{UInt8}}}, (:prec_id, :prec_mz, :score, :charge)}, 
                                    matched_frag_range::UnitRange{UInt32})# where {T,U<:AbstractFloat}
        @inline @inbounds for i in matched_frag_range
            frag = fragments[i]
            inc!(prec_id_to_score, getPrecID(frag), getScore(frag))
        end
    end
    
    #Slower, but for clarity, could have used searchsortedfirst and searchsortedlast to get the same result.
    #Could be used for testing purposes. 
    
    @inline addFragmentMatches!(prec_id_to_score, fragments, window_start:window_stop)
    return 

end

function queryFragment!(prec_id_to_score::Counter{UInt32, UInt8}, 
                        frag_bin_max_idx::UInt32,
                        lower_bound_guess::UInt32,
                        upper_bound_guess::UInt32,
                        frag_index::FragmentIndex{Float32}, 
                        frag_min::Float32, frag_max::Float32, 
                        prec_bounds::Tuple{Float32, Float32})# where {T,U<:AbstractFloat}
    
    #First frag_bin matching fragment tolerance
    frag_bin_idx, lower_bound_guess, upper_bound_guess = findFirstFragmentBin(
                                    getFragBins(frag_index), 
                                    frag_bin_max_idx,
                                    lower_bound_guess,
                                    upper_bound_guess,
                                    frag_min, 
                                    frag_max
                                    )

    #No fragment bins contain the fragment m/z
    #println("no frags")
    if iszero(frag_bin_idx)
        return frag_bin_max_idx, lower_bound_guess, upper_bound_guess
    end
    #println("some frags $frag_bin_idx")
    #Search subsequent frag bins until no more bins or untill a bin is outside the fragment tolerance
    @inbounds @fastmath while (frag_bin_idx <= frag_bin_max_idx)
        #Fragment bin is outside the fragment tolerance
        frag_bin = getFragmentBin(frag_index, frag_bin_idx)
        if (getLow(frag_bin) > frag_max)
            #println("last frag bin $frag_bin_idx")
            return UInt32(frag_bin_idx), lower_bound_guess, upper_bound_guess
        else
            frag_id_range = getSubBinRange(frag_bin)
            #Search the fragment range for fragments with precursors in the precursor tolerance
            searchPrecursorBin!(prec_id_to_score, 
                                getFragments(frag_index),
                                frag_id_range, 
                                first(prec_bounds), 
                                last(prec_bounds)
                            )
            #println("searched precursor bin")
            frag_bin_idx += 1
        end
    end
    #Only reach this point if frag_bin exceeds length(frag_index)
    return UInt32(frag_bin_idx - 1), lower_bound_guess, upper_bound_guess
end

function searchScan!(prec_id_to_score::Counter{UInt32, UInt8}, 
                    frag_index::FragmentIndex{Float32}, 
                    masses::AbstractArray{Union{Missing, U}}, intensities::AbstractArray{Union{Missing, U}}, 
                    rt_bin_idx::Int64,
                    irt_high::Float32, 
                    ppm_err::Float32, 
                    mass_err_model::MassErrorModel{Float32},
                    min_max_ppm::Tuple{Float32, Float32},
                    prec_mz::Float32, 
                    prec_tol::Float32, 
                    isotope_err_bounds::Tuple{Int64, Int64}; 
                    min_score::UInt8 = zero(UInt8)) where {U<:AbstractFloat}
    
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
        return Float32(mass - tol), Float32(mass + tol)
    end

    function filterPrecursorMatches!(prec_id_to_score::Counter{UInt32, UInt8}, min_score::UInt8)
        match_count = countFragMatches(prec_id_to_score, min_score)
        prec_count = getSize(prec_id_to_score) - 1
        sortCounter!(prec_id_to_score);
        return match_count, prec_count
    end

    prec_min = Float32(prec_mz - prec_tol - NEUTRON*first(isotope_err_bounds)/2)
    prec_max = Float32(prec_mz + prec_tol + NEUTRON*last(isotope_err_bounds)/2)
    @inbounds @fastmath while getLow(getRTBin(frag_index, rt_bin_idx)) < irt_high

        #BinRanges
        sub_bin_range = getSubBinRange(getRTBin(frag_index, rt_bin_idx))
        min_frag_bin, max_frag_bin = first(sub_bin_range), last(sub_bin_range)
        lower_bound_guess, upper_bound_guess = min_frag_bin, min_frag_bin

        for (mass, intensity) in zip(masses, intensities)
            #mass, intensity = coalesce(mass, zero(U)),  coalesce(intensity, zero(U))
            #Get intensity dependent fragment tolerance
            frag_min, frag_max = getFragTol(mass, ppm_err, intensity, mass_err_model, min_max_ppm)
            #Don't avance the minimum frag bin if it is still possible to encounter a smaller matching fragment 


            min_frag_bin, lower_bound_guess, upper_bound_guess = queryFragment!(prec_id_to_score, 
                                            max_frag_bin,
                                            lower_bound_guess,
                                            upper_bound_guess,
                                            frag_index, 
                                            frag_min, frag_max, 
                                            (prec_min, prec_max)
                                        )
          

            #println("min_frag_bin post $min_frag_bin")
        end
        rt_bin_idx += 1
        #Cannot exceed number of rt bins 
        if rt_bin_idx >length(getRTBins(frag_index))
            rt_bin_idx = length(getRTBins(frag_index))
            break
        end 
    end
    #return 0
    return filterPrecursorMatches!(prec_id_to_score, min_score)
end

