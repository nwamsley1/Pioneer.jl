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

function searchPrecursorBin!(precs::Dictionary{UInt32, UInt8}, precursor_bin::PrecursorBin{T}, window_min::AbstractFloat, window_max::AbstractFloat) where {T<:AbstractFloat}
   
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

    window_start = (lo <= N ? lo : return nothing, nothing)

    if getPrecMZ(getPrecursor(precursor_bin, window_start)) > window_max
        return nothing, nothing
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

    function addFragmentMatches!(precs::Dictionary{UInt32, UInt8}, precursor_bin::PrecursorBin{T}, window_min::AbstractFloat, window_max::AbstractFloat, start::Int, stop::Int) where {T<:AbstractFloat}
        for precursor_idx in start:stop
            
            precursor = getPrecursor(precursor_bin, precursor_idx)
            prec_mz = getPrecMZ(precursor)
            prec_id = getPrecID(precursor)
            charge = getPrecCharge(precursor)

            prec_mz_min = window_min - 1.0*NEUTRON/charge#upper_tol[charge]
            prec_mz_max = window_max + 3.0*NEUTRON/charge#lower_tol[charge]
            if (prec_mz_min <= prec_mz) & (prec_mz_max >= prec_mz)
                if haskey(precs, prec_id)
                    precs[prec_id] += UInt8(1)
                else
                    insert!(precs, prec_id, UInt8(1))
                end
            end
        end

    end

    addFragmentMatches!(precs, precursor_bin, window_min, window_max, window_start, window_stop)

    return window_start, window_stop

end

#const upper_tol = [(3.0*NEUTRON), (3.0*NEUTRON)/2, (3.0*NEUTRON)/3, (3.0*NEUTRON)/4]
#const lower_tol = [(1*NEUTRON), (1*NEUTRON)/2, (1*NEUTRON)/3, (1*NEUTRON)/4]

function queryFragment!(precs::Dictionary{UInt32, UInt8}, frag_index::FragmentIndex{T}, min_frag_bin::Int64, frag_min::AbstractFloat, frag_max::AbstractFloat, prec_mz::AbstractFloat, prec_tol::AbstractFloat) where {T<:AbstractFloat}
    
    frag_bin = findFirstFragmentBin(getFragBins(frag_index), frag_min, frag_max)
    #No fragment bins contain the fragment m/z
    if (frag_bin === nothing)
        return min_frag_bin
    #This frag bin has already been searched
    elseif frag_bin <= min_frag_bin
        return min_frag_bin
    end

    i = 1
    while (frag_bin < length(getFragBins(frag_index))) #getLowMZ(getFragmentBin(frag_index, frag_bin)) <frag_max
        #Fragment bin matches the fragment ion
        #println(i)
        i += 1
        if (getLowMZ(getFragmentBin(frag_index, frag_bin)) > frag_max)
            return frag_bin
        else
            prec_min = prec_mz - prec_tol - 3.0*NEUTRON
            prec_max = prec_mz + prec_tol + 1.0*NEUTRON
            searchPrecursorBin!(precs, getPrecursorBin(frag_index, UInt32(frag_bin)), prec_min, prec_max)
            frag_bin += 1
        end

    end

    #Only reach this point if frag_bin exceeds length(frag_index)
    return frag_bin - 1
end

function searchScan!(precs::Dictionary{UInt32, UInt8}, f_index::FragmentIndex{T}, massess::Vector{Union{Missing, U}}, intensities::Vector{Union{Missing, U}}, precursor_window::AbstractFloat, ppm::AbstractFloat, width::AbstractFloat; topN::Int = 20, min_frag_count::Int = 3) where {T,U<:AbstractFloat}
    
    getFragTol(mass::U, ppm::AbstractFloat) = mass*(1 - ppm/1e6), mass*(1 + ppm/1e6)

    function filterPrecursorMatches!(precs::Dictionary{UInt32, UInt8}, topN::Int, min_frag_count::Int) where {T<:AbstractFloat}
        #Do not consider peptides wither fewer than 
        match_count = sum(precs)
        prec_count = length(precs)
        filter!(count->(count>=min_frag_count), precs)
        sort!(precs, rev = true)
        #println(precs)
        #Iterator of Peptide ID's for the `topN` scoring peptides
        return Iterators.take(keys(precs), min(topN, length(keys(precs)))), prec_count, match_count
    end

    min_frag_bin = 0

    for (mass, intensity) in zip(massess, intensities)

        mass, intensity = coalesce(mass, 0.0),  coalesce(intensity, 0.0)

        FRAGMIN, FRAGMAX = getFragTol(mass, ppm) 

        min_frag_bin = queryFragment!(precs, f_index, min_frag_bin, FRAGMIN, FRAGMAX, precursor_window, width)
    end 
    #println("PRECS $precs")
    return filterPrecursorMatches!(precs, topN, min_frag_count)
end
