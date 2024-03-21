function findFirstFragmentBin(frag_index_bins::Arrow.Struct{FragIndexBin, Tuple{Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{UInt32, Vector{UInt32}}}, (:lb, :ub, :first_bin, :last_bin)}, 
                                frag_bin_range::UnitRange{UInt32},
                                upper_bound_guess::UInt32,
                                frag_min::Float32, 
                                frag_max::Float32) #where {T<:AbstractFloat}
    #Binary Search
    lo, hi = first(frag_bin_range), last(frag_bin_range)
    @fastmath len = hi - lo
    potential_match = zero(UInt32)
    @fastmath mid = (lo + hi)>>>0x01#min(lo + (upper_bound_guess*UInt32(2)), (lo + hi)>>>1)#(hi - lo)÷upper_bound_guess
    @inbounds @fastmath while lo <= hi

        if (frag_min) <= getHigh(frag_index_bins[mid])
            if (frag_max) >= getHigh(frag_index_bins[mid]) #Frag tolerance overlaps the upper boundary of the frag bin
                potential_match = UInt32(mid)
            end
            hi = mid - one(UInt32)
            #mid = (lo + hi)÷2
        elseif (frag_max) >= getLow(frag_index_bins[mid]) #Frag tolerance overlaps the lower boundary of the frag bin
            if (frag_min) <= getLow(frag_index_bins[mid])
                potential_match = UInt32(mid)
            end
            lo = mid + one(UInt32)
            #mid = hi - (hi - lo)÷10
        end
        mid = (lo + hi) >>> 0x01
    end

    return potential_match, one(UInt32)#UInt32(max(potential_match, lo + 1) - lo)#UInt32(2*(hi_f - low_f)÷(mid - low_f))
end
#=
function findFirstFragmentBin(frag_index_bins::Arrow.Struct{FragIndexBin, Tuple{Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{UInt32, Vector{UInt32}}}, (:lb, :ub, :first_bin, :last_bin)}, 
                                frag_bin_range::UnitRange{UInt32},
                                upper_bound_guess::UInt32,
                                frag_min::Float32, 
                                frag_max::Float32) #where {T<:AbstractFloat}
    #Binary Search
    lo, hi = first(frag_bin_range), last(frag_bin_range)
    @fastmath len = hi - lo
    potential_match = zero(UInt32)
    @fastmath mid = (lo + hi)>>>0x01#min(lo + (upper_bound_guess*UInt32(2)), (lo + hi)>>>1)#(hi - lo)÷upper_bound_guess
    @inbounds @fastmath while lo <= hi

        if (frag_min) <= getHigh(frag_index_bins[mid])
            if (frag_max) >= getHigh(frag_index_bins[mid]) #Frag tolerance overlaps the upper boundary of the frag bin
                potential_match = UInt32(mid)
            end
            hi = mid - one(UInt32)
            #mid = (lo + hi)÷2
        elseif (frag_max) >= getLow(frag_index_bins[mid]) #Frag tolerance overlaps the lower boundary of the frag bin
            if (frag_min) <= getLow(frag_index_bins[mid])
                potential_match = UInt32(mid)
            end
            lo = mid + one(UInt32)
            #mid = hi - (hi - lo)÷10
        end
        mid = (lo + hi) >>> 0x01
    end

    return potential_match, one(UInt32)#UInt32(max(potential_match, lo + 1) - lo)#UInt32(2*(hi_f - low_f)÷(mid - low_f))
end

function findFirstFragmentBin(frag_index_bins::Arrow.Struct{FragIndexBin, Tuple{Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{UInt32, Vector{UInt32}}}, (:lb, :ub, :first_bin, :last_bin)}, 
                                frag_bin_range::UnitRange{UInt32},
                                upper_bound_guess::UInt32,
                                frag_min::Float32, 
                                frag_max::Float32) #where {T<:AbstractFloat}
    #Binary Search
    lo, hi = first(frag_bin_range), last(frag_bin_range)
    @fastmath len = hi - lo
    #potential_match = zero(UInt32)
    #@fastmath mid = (lo + hi)>>>0x01#min(lo + (upper_bound_guess*UInt32(2)), (lo + hi)>>>1)#(hi - lo)÷upper_bound_guess
    @fastmath mid = len>>>0x01
    @inbounds @fastmath while (len > 0x01)
        #lo += (getHigh(frag_index_bins[mid-one(UInt32)]) >= frag_min)*mid 
        lo += (getHigh(frag_index_bins[mid + lo - one(UInt32)]) >= frag_min)*(getLow(frag_index_bins[mid + lo - one(UInt32)]) < frag_max)*mid 
        len -= mid 
        mid = len>>>0x01
    end

    @inbounds @fastmath begin
        if (mid == one(UInt32)) & (lo == first(frag_bin_range))
            return zero(UInt32), one(UInt32)
        else
            return lo, one(UInt32)#UInt32(max(potential_match, lo + 1) - lo)#UInt32(2*(hi_f - low_f)÷(mid - low_f))
        end
    end

end
=#
function searchPrecursorBin!(prec_id_to_score::Counter{UInt32, UInt8}, 
                            fragments::Arrow.Struct{IndexFragment, Tuple{Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{UInt8, Vector{UInt8}}, Arrow.Primitive{UInt8, Vector{UInt8}}}, (:prec_id, :prec_mz, :score, :charge)},
                            frag_id_range::UnitRange{UInt32},
                            window_min::Float32, 
                            window_max::Float32)

        lo, hi = first(frag_id_range), last(frag_id_range)
        base = lo
        @inbounds @fastmath begin 
            len = hi - lo + UInt32(1)
            if len > 4 #If query range is sufficiently large, do a branchless binary search
                while len > 1
                    mid = len >>> 0x01
                    base += (getPrecMZ(fragments[base + mid - one(UInt32)]) < window_min)*mid
                    len -= mid
                end
                window_start = base
                len = hi - base + UInt32(1)
                base = hi
                while len > 1
                    mid = len >>> 0x01
                    base -= (getPrecMZ(fragments[base - mid + one(UInt32)]) > window_max)*mid
                    len -= mid
                end
                window_stop = base

                if (getPrecMZ(fragments[window_start])<window_min) | (getPrecMZ(fragments[window_stop])>window_max)
                    return 
                end
            
            else #Query range is small, do linear search
                window_start, window_stop = lo, hi
                #println("a window_start $window_start window_stop $window_stop")
                matched = (getPrecMZ(fragments[window_start])<window_min)
                while (matched) & (window_start < hi)
                    window_start += matched
                    matched = (getPrecMZ(fragments[window_start])<window_min)
                end
                matched = (getPrecMZ(fragments[window_stop])>window_max)
                while matched & (window_stop > window_start)
                    window_stop -= matched
                    matched = (getPrecMZ(fragments[window_stop])>window_max)
                end
                #println("b window_start $window_start window_stop $window_stop")
                if window_start > window_stop
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

#=
function searchPrecursorBin!(prec_id_to_score::Counter{UInt32, UInt8}, 
                            fragments::Arrow.Struct{IndexFragment, Tuple{Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{UInt8, Vector{UInt8}}, Arrow.Primitive{UInt8, Vector{UInt8}}}, (:prec_id, :prec_mz, :score, :charge)},
                            frag_id_range::UnitRange{UInt32},
                            window_min::Float32, 
                            window_max::Float32)

    
        #N = getLength(precursor_bin)
        N =  last(frag_id_range)
        #lo, hi = 1, N
        lo, hi = first(frag_id_range), N
        @inbounds @fastmath begin 
            while lo <= hi
                mid = (lo + hi) >>> 0x01
                if getPrecMZ(fragments[mid]) < window_min
                    lo = mid + one(UInt32)
                else
                    hi = mid - one(UInt32)
                end
            end

            window_start = (lo <= N ? lo : return)


            if getPrecMZ(fragments[window_start]) > window_max
                return 
            end
            
            lo, hi = window_start, N

            while lo <= hi
                mid = (lo + hi) >>> 0x01
                if getPrecMZ(fragments[mid]) > window_max
                    hi = mid - UInt32(1)
                else
                    lo = mid + UInt32(1)
                end
            end

            window_stop = hi
        end
        #=
        window_stop = window_start
        @inbounds @fastmath while getPrecMZ(fragments[window_stop]) < window_max
            window_stop += 1
        end
        window_stop = max(window_stop - 1, window_start)
        =#
        
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
    
    addFragmentMatches!(prec_id_to_score, fragments, window_start:window_stop)
    return 

end
=#
function queryFragment!(prec_id_to_score::Counter{UInt32, UInt8}, 
                        frag_bin_range::UnitRange{UInt32},
                        upper_bound_guess::UInt32,
                        frag_index::FragmentIndex{Float32}, 
                        frag_min::Float32, frag_max::Float32, 
                        prec_bounds::Tuple{Float32, Float32})# where {T,U<:AbstractFloat}
    
    #First frag_bin matching fragment tolerance
    frag_bin_idx, upper_bound_guess = findFirstFragmentBin(getFragBins(frag_index), 
                                    frag_bin_range, #Range of Fragment Bins to Search
                                    upper_bound_guess,
                                    frag_min, frag_max)

    #No fragment bins contain the fragment m/z
    #println("no frags")
    if iszero(frag_bin_idx)
        return first(frag_bin_range), upper_bound_guess
    end
    #println("some frags $frag_bin_idx")
    #Search subsequent frag bins until no more bins or untill a bin is outside the fragment tolerance
    @inbounds while (frag_bin_idx <= last(frag_bin_range))
        #Fragment bin is outside the fragment tolerance
        frag_bin = getFragmentBin(frag_index, frag_bin_idx)
        if (getLow(frag_bin) > frag_max)
            #println("last frag bin $frag_bin_idx")
            return UInt32(frag_bin_idx), upper_bound_guess
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
    return UInt32(frag_bin_idx - 1), upper_bound_guess
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
        return Float32(mass - tol), Float32(mass + tol), Float32(mass - last(min_max_ppm)*mass/1e6)
    end

    function filterPrecursorMatches!(prec_id_to_score::Counter{UInt32, UInt8}, min_score::UInt8)
        match_count = countFragMatches(prec_id_to_score, min_score)
        prec_count = getSize(prec_id_to_score) - 1
        sortCounter!(prec_id_to_score);
        return match_count, prec_count
    end

    prec_min = Float32(prec_mz - prec_tol - NEUTRON*first(isotope_err_bounds)/2)
    prec_max = Float32(prec_mz + prec_tol + NEUTRON*last(isotope_err_bounds)/2)
    while getLow(getRTBin(frag_index, rt_bin_idx)) < irt_high
        #println("rt_bin_idx $rt_bin_idx irt_high $irt_high getRTBin(frag_index, rt_bin_idx) ", getRTBin(frag_index, rt_bin_idx))
        sub_bin_range = getSubBinRange(getRTBin(frag_index, rt_bin_idx))
        min_frag_bin, max_frag_bin = first(sub_bin_range), last(sub_bin_range)
        old_min_frag_bin = min_frag_bin
        upper_bound_guess = UInt32((min_frag_bin + max_frag_bin)÷2)#UInt32(min_frag_bin + 1)
        @inbounds for (mass, intensity) in zip(masses, intensities)
            mass, intensity = coalesce(mass, zero(U)),  coalesce(intensity, zero(U))
            #Get intensity dependent fragment tolerance
            frag_min, frag_max, MIN_ALL_FRAGS = getFragTol(mass, ppm_err, intensity, mass_err_model, min_max_ppm)
            #Don't avance the minimum frag bin if it is still possible to encounter a smaller matching fragment 
            if getHigh(getFragmentBin(frag_index, (min_frag_bin))) >= MIN_ALL_FRAGS
                min_frag_bin = old_min_frag_bin
            end
            #println("frag_min $frag_min, frag_max $frag_max, MIN_ALL_FRAGS $MIN_ALL_FRAGS")
            old_min_frag_bin = min_frag_bin

            min_frag_bin, upper_bound_guess = queryFragment!(prec_id_to_score, 
                                            min_frag_bin:max_frag_bin, 
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

