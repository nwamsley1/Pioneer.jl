function getPrecursorsInRTWindow(retentionTimes::Vector{Float32}, l_bnd::Float32, u_bnd::Float32)
    #print("l_bnd ", l_bnd)
    #print("u_bnd ", u_bnd)
    start = searchsortedfirst(retentionTimes, l_bnd ,lt=(t,x)->t<x)
    stop = searchsortedlast(retentionTimes, u_bnd,lt=(x, t)->t>x)
    #println("test inner")
    #println("start1 ", start)
    #println("stop1 ", stop)
    return start, stop
end

function getMS1Peaks!(precursors::Dictionary{UInt32, Precursor}, 
                        MS1::Vector{Union{Missing, Float32}}, 
                        INTENSITIES::Vector{Union{Missing, Float32}}, 
                        MS1_MAX_HEIGHTS::UnorderedDictionary{UInt32, Float32}, 
                        precursor_rts::Vector{Float32}, 
                        precursor_idxs::Vector{UInt32}, 
                        precursor_ms_file_idxs::Vector{UInt32}, 
                        rt::Float32, rt_tol::Float32, left_mz_tol::Float32, right_mz_tol::Float32, ms_file_idx::UInt32)
    
    #Get precursors for which the best scan RT is within `rt_tol` of the current scan `rt`
    #precursor_idxs =  @view(bestPSMs[!,:precursor_idx][getPrecursorsInRTWindow(bestPSMs[!,:retentionTime],rt - rt_tol, rt + rt_tol)])
    #println("test")
    start, stop = getPrecursorsInRTWindow(precursor_rts,rt - rt_tol, rt + rt_tol)
    #println("start ", start)
    #println("stop ", stop)
    #stop = 1
    #start = 1
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

            if coalesce(INTENSITIES[idx], Float32(0))>=MS1_MAX_HEIGHTS[precursor_idx] #Replace maximum observed MS1 height for the precursor if appropriate. 
                MS1_MAX_HEIGHTS[precursor_idx] = coalesce(INTENSITIES[idx], Float32(0))
            end
        end
    end
end

function getMS1PeakHeights!(ptable::PrecursorTable, retentionTimes::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}, 
                            masses::Arrow.List{Union{Missing, Vector{Union{Missing, Float32}}}, Int32, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}},
                            intensities::Arrow.List{Union{Missing, Vector{Union{Missing, Float32}}}, Int32, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}},
                            msOrders::Arrow.Primitive{Union{Missing, Int32}, Vector{Int32}},
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