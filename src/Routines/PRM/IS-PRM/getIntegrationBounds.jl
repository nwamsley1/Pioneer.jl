function getTargetedScans(array::AbstractArray, lower_bound::Float32, upper_bound::Float32)

    findall(precursor_mz->coalesce(
                        (precursor_mz<(upper_bound)) & 
                        (precursor_mz>(lower_bound)), false), array)

end


getTargetedScans(prec::Ion, array::AbstractArray) = getTargetedScans(array, getLow(prec), getHigh(prec))

function getTargetedScans(precursor_pair::LightHeavyPrecursorPair, 
                            precursors::Dictionary{UInt32, Precursor}, 
                            precursorMZs::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}})

    return (
            light_scan_ids = getTargetedScans(precursors[getLightPrecID(precursor_pair)], precursorMZs),
            heavy_scan_ids = getTargetedScans(precursors[getHeavyPrecID(precursor_pair)], precursorMZs)
    )
end

function setIntegrationBounds!(ptable::ISPRMPrecursorTable, 
                               run_idx::Int64,
                               scan_adresses::Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}},
                               precursorMZs::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}})
    
    for (key, lh_pair) in pairs(getIDToLightHeavyPair(ptable))
        targeted_scans = getTargetedScans(lh_pair, getPrecursors(ptable), precursorMZs)
        #setTargetedScans!(value, getPrecursors(ptable), precursorMZs)
        integration_bounds = getIntegrationBounds(
                                                    getScanCycleUnion(
                                                        scan_adresses[targeted_scans[:light_scan_ids]],
                                                        scan_adresses[targeted_scans[:heavy_scan_ids]]
                                                    )
                                                )
        if (integration_bounds[:upper_bound] - integration_bounds[:lower_bound]) > 3
            setIntegrationBounds!(lh_pair, run_idx, integration_bounds)
        end
    end
end

function getScanAdresses(scan_order::AbstractArray)
    scan_adresses = Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}}(undef,length(scan_order)) 
    ms1 = 0
    msn = 0
    for scan in enumerate(scan_order)
        if scan[2] == 1
            ms1 += 1
            msn = 0
            scan_adresses[scan[1]] = (scan_index = scan[1], ms1 = ms1, msn = msn)
        else
            msn+=1
            scan_adresses[scan[1]] = (scan_index = scan[1], ms1 = ms1, msn = msn)
        end
    end
    scan_adresses
end

function getScanCycleUnion(scan_adresses_1::Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}}, 
                           scan_adresses_2::Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}}
                          )
    sort(collect(Set(map(x->x.ms1, scan_adresses_1))
                ∩Set(map(x->x.ms1, scan_adresses_2))
                )
        )
end

function getIntegrationBounds(scan_indices::Vector{Int64}; max_gap_size::Int = 10)
    if length(scan_indices)<=1
        return (lower_bound = 1, upper_bound = 1)
    end
    start = 1
    stop = 1
    last = scan_indices[stop]
    best_start = 1
    best_stop = 1
    gap_size = 0
    for index in enumerate(diff(scan_indices))
        if gap_size >= max_gap_size
            gap_size = 0
            start = index[1]
            stop = index[1]
            last = scan_indices[index[1]]
            continue
        end

        if index[2] == 1
            last = scan_indices[index[2]]
            stop = index[1]
            if (stop-start)>(best_stop-best_start)
                best_start = start
                best_stop = stop
            end
        else
            gap_size = scan_indices[index[1]] - last
        end
    end
    scan_indices[best_start:(best_stop +1)]
    #(lower_bound = best_start, upper_bound = best_stop+1)
end



