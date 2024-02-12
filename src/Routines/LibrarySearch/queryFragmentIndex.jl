function findFirstFragmentBin(frag_index::Vector{FragBin{Float32}}, frag_min::Float32, frag_max::Float32) #where {T<:AbstractFloat}
    #Binary Search
    lo, hi = 1, length(frag_index)
    potential_match = nothing
    while lo <= hi

        mid = (lo + hi) ÷ 2

        if (frag_min) <= getHighMZ(frag_index[mid])
            if (frag_max) >= getHighMZ(frag_index[mid]) #Frag tolerance overlaps the upper boundary of the frag bin
                potential_match = mid
            end
            hi = mid - 1
        elseif (frag_max) >= getLowMZ(frag_index[mid]) #Frag tolerance overlaps the lower boundary of the frag bin
            if (frag_min) <= getLowMZ(frag_index[mid])
                potential_match = mid
                #return mid
            end
            lo = mid + 1
        end
    end

    return potential_match#, Int64(getPrecBinID(frag_index[potential_match]))
end

function findFirstRTBin(rt_bins::Vector{RTBin{Float32}}, rt_min::Float32, rt_max::Float32)# where {T<:AbstractFloat}
    #Binary Search
    lo, hi = 1, length(rt_bins)
    potential_match = nothing
    while lo <= hi

        mid = (lo + hi) ÷ 2

        if (rt_min) <= getHigh(rt_bins[mid])
            if (rt_max) >= getHigh(rt_bins[mid]) #Frag tolerance overlaps the upper boundary of the frag bin
                potential_match = mid
            end
            hi = mid - 1
        elseif (rt_max) >= getLow(rt_bins[mid]) #Frag tolerance overlaps the lower boundary of the frag bin
            if (rt_min) <= getLow(rt_bins[mid])
                potential_match = mid
                #return mid
            end
            lo = mid + 1
        end
    end

    return potential_match#, Int64(getPrecBinID(frag_index[potential_match]))
end

function searchPrecursorBin!(precs::Counter{UInt32, UInt8}, precursor_bin::PrecursorBin{Float32}, window_min::Float32, window_max::Float32)

    
        N = getLength(precursor_bin)
        lo, hi = 1, N

        @inbounds @fastmath while lo <= hi
            mid = (lo + hi) ÷ 2
            if getPrecMZ(getPrecursor(precursor_bin, mid)) < window_min
                lo = mid + 1
            else
                hi = mid - 1
            end
        end

        window_start = (lo <= N ? lo : return)

        if getPrecMZ(getPrecursor(precursor_bin, window_start)) > window_max
            return 
        end

        lo, hi = window_start, N

        @inbounds @fastmath while lo <= hi
            mid = (lo + hi) ÷ 2
            if getPrecMZ(getPrecursor(precursor_bin, mid)) > window_max
                hi = mid - 1
            else
                lo = mid + 1
            end
        end

        window_stop = hi

    function addFragmentMatches!(precs::Counter{UInt32, UInt8}, precursor_bin::PrecursorBin{Float32}, start::Int, stop::Int)# where {T,U<:AbstractFloat}
        #initialize
        #=precursor_idx = start
        ms1 = MzFeature(MS1[ms1_idx], ppm = Float32(prec_ppm))
        prec = getPrecursor(precursor_bin, precursor_idx)

        
        while true 
            #println(precursor_idx)
            if getPrecMZ(prec) > getLow(ms1)
                if getPrecMZ(prec) < getHigh(ms1) #The precursor matches the MS1 fragment
                    inc!(precs, getPrecID(prec), getIntensity(prec), intensity)
                    #Move on to the next precursor
                    precursor_idx += 1
                    (precursor_idx <= stop) ? prec = getPrecursor(precursor_bin, precursor_idx) : return
                else #The precursor m/z is above the MS1 fragment 
                    #Move on to the next ms1 peak
                    ms1_idx += 1
                    (ms1_idx <= length(MS1)) ? ms1 = MzFeature(MS1[ms1_idx], ppm = Float32(prec_ppm)) : return 
                end
            else #The precursor is below the MS1 fragment
                precursor_idx += 1
                (precursor_idx <= stop) ? prec = getPrecursor(precursor_bin, precursor_idx) : return
            end
        end=#

        @inbounds for precursor_idx in start:stop
            prec = getPrecursor(precursor_bin, precursor_idx)
            inc!(precs, getPrecID(prec), getScore(prec))
        end
    end
    
    #window_start= searchsortedfirst(precursor_bin.precs, window_min, lt=(r,x)->r.prec_mz<x)
    #window_stop =  searchsortedlast(precursor_bin.precs, window_max, lt=(x, r)->r.prec_mz>x)
    addFragmentMatches!(precs, precursor_bin, window_start, window_stop)

    return 

end

function queryFragment!(precs::Counter{UInt32, UInt8}, iRT_low::Float32, iRT_high::Float32, frag_index::FragmentIndex{Float32}, min_frag_bin::Int64, frag_min::Float32, frag_max::Float32, prec_bounds::Tuple{Float32, Float32})# where {T,U<:AbstractFloat}
    
    frag_bin = findFirstFragmentBin(getFragBins(frag_index), frag_min, frag_max)
    #No fragment bins contain the fragment m/z
    if (frag_bin === nothing)
        return min_frag_bin
    #This frag bin has already been searched
    #Should test how often this branch is being hit. 
    #If frequent, may need a different solution. If a fragment mz could match multiple different
    #peaks in the spectra given a fragment tolerarnce, with this arangement it will always be matched to the fragment with the lowest m/Z
    #which could be incorrect. 
    elseif frag_bin <= min_frag_bin
        return min_frag_bin
    end

    i = 1
    #prec_min = prec_mz #prec_tol - Float32(1.51)
    #prec_max = prec_mz + prec_tol + Float32(.51)

    while (frag_bin < length(getFragBins(frag_index))) #getLowMZ(getFragmentBin(frag_index, frag_bin)) <frag_max
    
        #Fragment bin matches the fragment ion
        #println(i)
        i += 1
        if (getLowMZ(getFragmentBin(frag_index, frag_bin)) > frag_max)
            return frag_bin
        else

            #RT Information Goes Here
            rt_bin_idx = frag_index.fragment_bins[frag_bin].sub_bin
            rt_sub_bin = findFirstRTBin(frag_index.rt_bins[rt_bin_idx], iRT_low, iRT_high)

            if rt_sub_bin === nothing
                frag_bin += 1
                continue
            end
            i = 1
            while (getLow(frag_index.rt_bins[rt_bin_idx][rt_sub_bin]) < iRT_high)
                prec_bin = frag_index.rt_bins[rt_bin_idx][rt_sub_bin].sub_bin
                searchPrecursorBin!(precs, getPrecursorBin(frag_index, prec_bin), first(prec_bounds), last(prec_bounds))
                rt_sub_bin += one(UInt32)
                i += 1
                if rt_sub_bin > length(frag_index.rt_bins[rt_bin_idx])
                    break
                end
            end

            frag_bin += 1
        end

    end

    #Only reach this point if frag_bin exceeds length(frag_index)
    return frag_bin - 1
end

function searchScan!(precs::Counter{UInt32, UInt8}, 
                    f_index::FragmentIndex{Float32}, 
                    min_intensity::Float32, 
                    masses::AbstractArray{Union{Missing, U}}, intensities::AbstractArray{Union{Missing, U}}, 
                    iRT_low::Float32, iRT_high::Float32, 
                    ppm_err::Float32, 
                    ppm_tol_param::Float32,
                    min_max_ppm::Tuple{Float32, Float32},
                    prec_mz::Float32, 
                    prec_tol::Float32, 
                    isotope_err_bounds::Tuple{Int64, Int64}; 
                    topN::Int = 20, min_frag_count::Int = 3, min_score::UInt8 = zero(UInt8)) where {U<:AbstractFloat}
    
    function getFragTol(mass::Float32, ppm_err::Float32, 
                        intensity::Float32, ppm_tol_param::Float32,
                        min_max_ppm::Tuple{Float32, Float32})
        mass -= Float32(ppm_err*mass/1e6)
        ppm = max(
                    min(ppm_tol_param/sqrt(intensity), last(min_max_ppm)), 
                    first(min_max_ppm)
                    )
        #ppm = 16.1
        tol = ppm*mass/1e6
        return Float32(mass - tol), Float32(mass + tol)
       #return Float32(mass*(1 - 16.1*mass/1e6)), Float32(mass*(1 + 16.1*mass/1e6))
    end

    function filterPrecursorMatches!(precs::Counter{UInt32, UInt8}, topN::Int, min_score::UInt8)
        #match_count = countFragMatches(precs, min_frag_count, min_ratio)
        match_count = countFragMatches(precs, min_score)
        prec_count = getSize(precs) - 1
        sort!(precs, topN);
        return match_count, prec_count
    end
    #println("TEST")
    min_frag_bin = 0 #Fragment bins with a lower index can be excluded from the search
    prec_min = Float32(prec_mz - prec_tol - NEUTRON*first(isotope_err_bounds)/2)
    prec_max = Float32(prec_mz + prec_tol + NEUTRON*last(isotope_err_bounds)/2)
    for (mass, intensity) in zip(masses, intensities)

        mass, intensity = coalesce(mass, 0.0),  coalesce(intensity, 0.0)

        if intensity<min_intensity
            continue
        end

        FRAGMIN, FRAGMAX = getFragTol(mass, ppm_err, intensity, ppm_tol_param, min_max_ppm)
        #println(typeof(width))
        #println("TEST")
        min_frag_bin = queryFragment!(precs, iRT_low, iRT_high, f_index, min_frag_bin, FRAGMIN, FRAGMAX, (prec_min, prec_max))
    end 

    return filterPrecursorMatches!(precs, topN, min_score)
end

