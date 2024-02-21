function findFirstFragmentBin(frag_index::Vector{FragBin{Float32}}, frag_min::Float32, frag_max::Float32) #where {T<:AbstractFloat}
    #Binary Search
    lo, hi = 1, length(frag_index)
    potential_match = nothing
    @inbounds @fastmath while lo <= hi

        mid = (lo + hi) ÷ 2

        if (frag_min) <= getHighMZ(frag_index[mid])
            if (frag_max) >= getHighMZ(frag_index[mid]) #Frag tolerance overlaps the upper boundary of the frag bin
                potential_match = mid
            end
            hi = mid - 1
        elseif (frag_max) >= getLowMZ(frag_index[mid]) #Frag tolerance overlaps the lower boundary of the frag bin
            if (frag_min) <= getLowMZ(frag_index[mid])
                potential_match = mid
            end
            lo = mid + 1
        end
    end

    return potential_match
end

function findFirstRTBin(rt_bins::Vector{RTBin{Float32}}, rt_min::Float32, rt_max::Float32)# where {T<:AbstractFloat}
    #Binary Search
    lo, hi = 1, length(rt_bins)
    potential_match = nothing
    @inbounds @fastmath while lo <= hi

        mid = (lo + hi) ÷ 2

        if (rt_min) <= getHigh(rt_bins[mid])
            if (rt_max) >= getHigh(rt_bins[mid]) #Frag tolerance overlaps the upper boundary of the frag bin
                potential_match = mid
            end
            hi = mid - 1
        elseif (rt_max) >= getLow(rt_bins[mid]) #Frag tolerance overlaps the lower boundary of the frag bin
            if (rt_min) <= getLow(rt_bins[mid])
                potential_match = mid
            end
            lo = mid + 1
        end
    end

    return potential_match
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
        @inline @inbounds for precursor_idx in start:stop
            frag = getPrecursor(precursor_bin, precursor_idx)
            inc!(precs, getPrecID(frag), getScore(frag))
        end
    end
    
    #Slower, but for clarity, could have used searchsortedfirst and searchsortedlast to get the same result.
    #Could be used for testing purposes. 
    #window_start= searchsortedfirst(precursor_bin.precs, window_min, lt=(r,x)->r.prec_mz<x)
    #window_stop =  searchsortedlast(precursor_bin.precs, window_max, lt=(x, r)->r.prec_mz>x)
    addFragmentMatches!(precs, precursor_bin, window_start, window_stop)

    return 

end

function queryFragment!(precs::Counter{UInt32, UInt8}, 
                        min_frag_bin::Int64, 
                        frag_index::FragmentIndex{Float32}, 
                        irt_low::Float32, irt_high::Float32, 
                        frag_min::Float32, frag_max::Float32, 
                        prec_bounds::Tuple{Float32, Float32})# where {T,U<:AbstractFloat}
    
    #First frag_bin matching fragment tolerance
    frag_bin = findFirstFragmentBin(getFragBins(frag_index), frag_min, frag_max) 
    #No fragment bins contain the fragment m/z
    if (frag_bin === nothing)
        return min_frag_bin
    elseif frag_bin <= min_frag_bin
        return min_frag_bin
    end
    #Search subsequent frag bins until no more bins or untill a bin is outside the fragment tolerance
    while (frag_bin < length(getFragBins(frag_index)))
        #println("frag_bin $frag_bin")
        #Fragment bin is outside the fragment tolerance
        if (getLowMZ(getFragmentBin(frag_index, frag_bin)) > frag_max)
            return frag_bin
        else

            rt_bin_idx = frag_index.fragment_bins[frag_bin].sub_bin
            #First rt bin matching the rt tolerance
            rt_sub_bin = findFirstRTBin(frag_index.rt_bins[rt_bin_idx], irt_low, irt_high)

            #No rt bins match the rt tolerance
            if rt_sub_bin === nothing
                frag_bin += 1
                continue
            end
            #Search subsequent rt bins untill no more rt bins or untill an rt bin is outside the retention time tolerance 
            while (getLow(frag_index.rt_bins[rt_bin_idx][rt_sub_bin]) < irt_high)
                prec_bin = frag_index.rt_bins[rt_bin_idx][rt_sub_bin].sub_bin
                #println("prec_bin $prec_bin")
                #Modify precs to include all precursors in the 'prec_bin' within the 'prec_bounds'/precursor tolerance
                #(should be the quadrupole isolation width in most scenarios)
                searchPrecursorBin!(precs, getPrecursorBin(frag_index, prec_bin), first(prec_bounds), last(prec_bounds))
                rt_sub_bin += one(UInt32)
                #No more rt bins
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
                    masses::AbstractArray{Union{Missing, U}}, intensities::AbstractArray{Union{Missing, U}}, 
                    irt_low::Float32, irt_high::Float32, 
                    ppm_err::Float32, 
                    ppm_tol_param::Float32,
                    min_max_ppm::Tuple{Float32, Float32},
                    prec_mz::Float32, 
                    prec_tol::Float32, 
                    isotope_err_bounds::Tuple{Int64, Int64}; 
                    min_score::UInt8 = zero(UInt8)) where {U<:AbstractFloat}
    
    function getFragTol(mass::Float32, ppm_err::Float32, 
                        intensity::Float32, ppm_tol_param::Float32,
                        min_max_ppm::Tuple{Float32, Float32})
        mass -= Float32(ppm_err*mass/1e6)
        ppm = max(
                    min(ppm_tol_param/sqrt(intensity), last(min_max_ppm)), 
                    first(min_max_ppm)
                    )
        tol = ppm*mass/1e6
        return Float32(mass - tol), Float32(mass + tol)
    end

    function filterPrecursorMatches!(precs::Counter{UInt32, UInt8}, min_score::UInt8)
        match_count = countFragMatches(precs, min_score)
        prec_count = getSize(precs) - 1
        sortCounter!(precs);
        return match_count, prec_count
    end
    min_frag_bin = 0
    prec_min = Float32(prec_mz - prec_tol - NEUTRON*first(isotope_err_bounds)/2)
    prec_max = Float32(prec_mz + prec_tol + NEUTRON*last(isotope_err_bounds)/2)
    #println("prec_min $prec_min prec_max $prec_max")
    for (mass, intensity) in zip(masses, intensities)
        mass, intensity = coalesce(mass, 0.0),  coalesce(intensity, 0.0)
        #Get intensity dependent fragment tolerance
        frag_min, frag_max = getFragTol(mass, ppm_err, intensity, ppm_tol_param, min_max_ppm)
        #println("FRAGMIN_MAX $frag_min $frag_max")
        #println("frag_min $frag_min, frag_max $frag_max")
        min_frag_bin = queryFragment!(precs, min_frag_bin, f_index, irt_low, irt_high, frag_min, frag_max, (prec_min, prec_max))
    end 
    #return 0
    return filterPrecursorMatches!(precs, min_score)
end

