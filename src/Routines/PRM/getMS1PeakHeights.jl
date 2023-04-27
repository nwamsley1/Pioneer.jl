#This function name needs to be more generic?
function rangeQuerySorted(sorted_array::Vector{T}, l_bnd::T, u_bnd::T) where T <: Number
    start = searchsortedfirst(sorted_array, l_bnd ,lt=(t,x)->t<x)
    stop = searchsortedlast(sorted_array, u_bnd,lt=(x, t)->t>x)
    return start, stop
end

function getMS1Peaks!(MS1_MAX_HEIGHTS::UnorderedDictionary{UInt32, T}, 
                        precursors::Dictionary{UInt32, Precursor}, 
                        MS1::Vector{Union{Missing, T}}, 
                        INTENSITIES::Vector{Union{Missing, T}}, 
                        precursor_rts::Vector{R}, 
                        precursor_idxs::Vector{UInt32}, 
                        precursor_ms_file_idxs::Vector{UInt32}, 
                        rt::R, rt_tol::R, left_mz_tol::T, right_mz_tol::T, ms_file_idx::UInt32) where {T <: Number, R <: Number}
    
    #Get precursors for which the best scan RT is within `rt_tol` of the current scan `rt`
    start::Int, stop::Int = rangeQuerySorted(precursor_rts,rt - rt_tol, rt + rt_tol)

    if (stop-start) >= 0
    #Check to see if the MS1 height for each precursor is greater than the maximum previously observed. 
        for (i, psm_idx) in enumerate(start:stop)

            if precursor_ms_file_idxs[psm_idx]!=ms_file_idx
                continue
            end

            precursor_idx = precursor_idxs[psm_idx]
            
            if !isassigned(MS1_MAX_HEIGHTS, precursor_idx) #If this precursor has not been encountered before. 
                insert!(MS1_MAX_HEIGHTS, precursor_idx, Float32(0))
            end

            mz = getMZ(precursors[precursor_idx]) #Get the precursor MZ

            idx = binaryGetNearest(MS1, mz, mz-left_mz_tol, mz+right_mz_tol) #Get the peak index of the peak nearest in mz to the precursor. 

            if idx == 0
                continue
            end

            if INTENSITIES[idx]>=MS1_MAX_HEIGHTS[precursor_idx] #Replace maximum observed MS1 height for the precursor if appropriate. 
                MS1_MAX_HEIGHTS[precursor_idx] = INTENSITIES[idx]
            end
        end
    end
end

function getMS1PeakHeights!(ptable::PrecursorDatabase, retentionTimes::AbstractArray, 
                            masses::AbstractArray,
                            intensities::AbstractArray, 
                            msOrders::AbstractArray,
                            ms1_max_heights::UnorderedDictionary{UInt32, Float32}, 
                            precursor_rts::Vector{Float32}, precursor_idxs::Vector{UInt32}, precursor_ms_file_idxs::Vector{UInt32},
                            rt_tol::Float32, left_mz_tol::Float32, right_mz_tol::Float32, ms_file_idx::UInt32)
    #println("tunction")
    #i = 1

    for scan_idx in eachindex(retentionTimes)
        if msOrders[scan_idx]!=Int32(1) #Skip non MS1 scans. 
            continue
        end
        #i += 1
        #println("test?")
        @inline getMS1Peaks!(getPrecIDToPrecursor(ptable), 
                    masses[scan_idx], 
                    intensities[scan_idx], 
                    ms1_max_heights, 
                    precursor_rts, 
                    precursor_idxs, 
                    precursor_ms_file_idxs,                    #Precursor retention times and id's (sorted by rt)                                      #PrecursorTable
                    retentionTimes[scan_idx],                               #RT of current scan
                    rt_tol, left_mz_tol, right_mz_tol, ms_file_idx #Only consider precursors where the best scan is within the `rt_col` of the current scan
                    );
    end
    #println("i ", i)
end