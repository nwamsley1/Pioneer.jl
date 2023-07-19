function findFirstFragmentBin(frag_index::Vector{FragBin{T}}, frag_min::AbstractFloat, frag_max::AbstractFloat) where {T<:AbstractFloat}
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

function findFirstRTBin(rt_bins::Vector{RTBin{T}}, rt_min::AbstractFloat, rt_max::AbstractFloat) where {T<:AbstractFloat}
    #Binary Search
    lo, hi = 1, length(rt_bins)
    potential_match = nothing
    while lo <= hi

        mid = (lo + hi) ÷ 2

        if (rt_min) <= getLow(rt_bins[mid])
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

function searchPrecursorBin!(precs::Counter{UInt32, UInt8, Float32}, ms1_idx::Int, MS1::Vector{Union{Missing, U}}, prec_ppm::AbstractFloat, intensity::Float32, precursor_bin::PrecursorBin{T}, window_min::Float64, window_max::Float64) where {T,U<:AbstractFloat}
   
    N = getLength(precursor_bin)
    lo, hi = 1, N

    while lo <= hi
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

    while lo <= hi
        mid = (lo + hi) ÷ 2
        if getPrecMZ(getPrecursor(precursor_bin, mid)) > window_max
            hi = mid - 1
        else
            lo = mid + 1
        end
    end

    window_stop = hi

    function addFragmentMatches!(precs::Counter{UInt32, UInt8, Float32}, ms1_idx::Int, MS1::Vector{Union{Missing, U}}, prec_ppm::AbstractFloat, intensity::Float32, precursor_bin::PrecursorBin{T}, start::Int, stop::Int) where {T,U<:AbstractFloat}
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
        for precursor_idx in start:stop
            prec = getPrecursor(precursor_bin, precursor_idx)
            inc!(precs, getPrecID(prec), getIntensity(prec))
        end
    end

    addFragmentMatches!(precs, ms1_idx, MS1, prec_ppm, intensity, precursor_bin, window_start, window_stop)

    return 

end

function queryFragment!(precs::Counter{UInt32, UInt8, Float32}, ms1_idx::Int, MS1::Vector{Union{Missing, U}}, prec_ppm::AbstractFloat, intensity::Float32, frag_index::FragmentIndex{T}, min_frag_bin::Int64, frag_min::Float64, frag_max::Float64, prec_mz::Float32, prec_tol::Float64) where {T,U<:AbstractFloat}
    
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
    prec_min = prec_mz - prec_tol
    prec_max = prec_mz + prec_tol

    while (frag_bin < length(getFragBins(frag_index))) #getLowMZ(getFragmentBin(frag_index, frag_bin)) <frag_max
    
        #Fragment bin matches the fragment ion
        #println(i)
        i += 1
        if (getLowMZ(getFragmentBin(frag_index, frag_bin)) > frag_max)
            return frag_bin
        else
            #RT Information Goes Here
            #rt_bin = findFirstRTBin(frag_index.rt_bins[rt_bin], r)
            searchPrecursorBin!(precs, ms1_idx, MS1, prec_ppm, intensity, getPrecursorBin(frag_index, UInt32(frag_bin)), prec_min, prec_max)
            frag_bin += 1
        end

    end

    #Only reach this point if frag_bin exceeds length(frag_index)
    return frag_bin - 1
end

function searchScan!(precs::Counter{UInt32, UInt8, Float32}, prec_norms::Vector{Float32}, f_index::FragmentIndex{T}, min_intensity::U, masses::Vector{Union{Missing, U}}, intensities::Vector{Union{Missing, U}}, MS1::Vector{Union{Missing, U}}, precursor_window::AbstractFloat, frag_ppm::AbstractFloat, prec_ppm::AbstractFloat, width::AbstractFloat; topN::Int = 20, min_frag_count::Int = 3, min_ratio::Float32 = Float32(0.8)) where {T,U<:AbstractFloat}
    
    getFragTol(mass::U, ppm::AbstractFloat) = mass*(1 - ppm/1e6), mass*(1 + ppm/1e6)

    function filterPrecursorMatches!(precs::Counter{UInt32, UInt8, Float32}, prec_norms::Vector{Float32}, topN::Int, min_frag_count::Int, min_ratio::Float32)
        match_count = countFragMatches(precs, prec_norms, min_frag_count, min_ratio)
        prec_count = getSize(precs) - 1
        sort!(precs, prec_norms, topN);
        return match_count, prec_count
    end
    #println("TEST")
    min_frag_bin = 0 #Fragment bins with a lower index can be excluded from the search
    ms1_idx = 1
    lower_bound = precursor_window - width
    lower_bound += -prec_ppm*lower_bound/(1e6)

    ms1_idx = searchsortedfirst(MS1, lower_bound) #Find the first MS1 peak with m/z greater than the lower_bound of the quadrupole isolation window
    if ms1_idx > length(MS1) #This occurs if every entry in MS1 is below the lower_bound of the quadrupole isolation window
        return 0, 0 #0 fragments mathing to 0 precursors 
    end

    for (mass, intensity) in zip(masses, intensities)

        mass, intensity = coalesce(mass, 0.0),  coalesce(intensity, 0.0)

        if intensity<min_intensity
            continue
        end

        FRAGMIN, FRAGMAX = getFragTol(mass, frag_ppm) 

        min_frag_bin = queryFragment!(precs, ms1_idx, MS1, prec_ppm, intensity, f_index, min_frag_bin, FRAGMIN, FRAGMAX, precursor_window, width)
    end 

    return filterPrecursorMatches!(precs, prec_norms, topN, min_frag_count, min_ratio)
end

