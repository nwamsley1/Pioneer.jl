function binaryGetNearest(arr::Vector{Union{Missing, Float32}}, query::Float32, low_tol::Float32, high_tol::Float32)

    #Check special cases (is the answer on the boundary or is the array empty?)
    n = length(arr)
    if n == 0 return 0 end
    if query < arr[1] - low_tol return 0 end
    if query > arr[n] + high_tol return 0 end

    function getNearest(arr::Vector{Union{Missing, Float32}}, lo::Int, hi::Int, query::Float32)
        if hi - lo>1
            smallest_distance = abs(query - arr[lo])
            best_idx = 1
            for (i, mass) in enumerate(@view(arr[lo:hi]))
                if abs(query - mass) < smallest_distance
                    smallest_distance = abs(query - mass)
                    best_idx = i
                end
            end
            return best_idx
        else
            return 1
        end
    end
    lo, hi = 1, n
    while lo <= hi
        mid = (lo + hi) ÷ 2
        if arr[mid] < (query - low_tol)
             lo = mid + 1
        elseif arr[mid] > (query + high_tol)
            hi = mid - 1
        else
            return getNearest(arr,lo, hi, query) + lo - 1
        end
    end

    return 0

end


function getPrecursors(window_center::Float32, precursorList::Vector{Precursor}, params)
    l_bnd, u_bnd = window_center - params[:lower_tol], window_center + params[:upper_tol]
    start, stop = searchsortedfirst(precursorList, l_bnd,lt=(t,x)->getMZ(t)<x), searchsortedlast(precursorList, u_bnd,lt=(x,t)->getMZ(t)>x)
    return @view(precursorList[start:stop])
end

#=
#searchsortedlast(bestPSMs[!,:retentionTime],  x, lt=(t, x)->t<x)
searchsortedfirst(bestPSMs[!,:retentionTime], 50 ,lt=(t,x)->t<x)
searchsortedlast(bestPSMs[!,:retentionTime], 50 ,lt=(x, t)->t>x)
searchsortedfirst(bestPSMs[!,:retentionTime],  x, ltx->50>x)

bestPSMs[170:180,:]
=#

#=
function getMS1PeakHeights(retentionTimes::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}, 
                            masses::Arrow.List{Union{Missing, Vector{Union{Missing, Float32}}}, Int32, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}},
                            intensities::Arrow.List{Union{Missing, Vector{Union{Missing, Float32}}}, Int32, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}},
                            msOrders::Arrow.Primitive{Union{Missing, Int32}, Vector{Int32}},
                            ms1_max_heights::UnorderedDictionary{UInt32, Float32}, 
                            precursor_rts::Vector{Float32}, precursor_idxs::Vector{UInt32}, 
                            precursorsDict::UnorderedDictionary{UInt32, SimplePrecursor}, 
                            rt_tol::Float32, left_mz_tol::Float32, right_mz_tol::Float32)
    #println("tunction")
    #i = 1
    for scan_idx in eachindex(retentionTimes)
        if msOrders[scan_idx]!=Int32(1) #Skip non MS1 scans. 
            continue
        end
        #i += 1
        #println("test?")
    #Get precursors for which the best scan RT is within `rt_tol` of the current scan `rt`
    #precursor_idxs =  @view(bestPSMs[!,:precursor_idx][getPrecursorsInRTWindow(bestPSMs[!,:retentionTime],rt - rt_tol, rt + rt_tol)])
    #println("test")

        start, stop = getPrecursorsInRTWindow(precursor_rts,retentionTimes[scan_idx] - rt_tol, retentionTimes[scan_idx] + rt_tol)

        if (stop-start) > 0
        #Check to see if the MS1 height for each precursor is greater than the maximum previously observed. 
            for (i, psm_idx) in enumerate(start:stop)
                
                precursor_idx = precursor_idxs[psm_idx]

                if !isassigned(ms1_max_heights, precursor_idx) #If this precursor has not been encountered before. 
                    insert!(ms1_max_heights, precursor_idx, Float32(0))
                end

                mz = getMZ(precursorsDict[precursor_idx]) #Get the precursor MZ

                idx = binaryGetNearest(masses[scan_idx], mz, mz-left_mz_tol, mz+right_mz_tol) #Get the peak index of the peak nearest in mz to the precursor. 

                if idx == 0
                    continue
                end
                
                if coalesce(intensities[scan_idx][idx], Float32(0))>=ms1_max_heights[precursor_idx] #Replace maximum observed MS1 height for the precursor if appropriate. 
                    ms1_max_heights[precursor_idx] = coalesce(intensities[scan_idx][idx], Float32(0))
                end
                
            end
        end
    end
    #println("i ", i)
end
=#