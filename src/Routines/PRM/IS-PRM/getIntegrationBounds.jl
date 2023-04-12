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
scan_adresses = getScanAdresses(GAPDH_VGVNGFGR[:msOrder])

#test.id_to_pepGroup[0x0000005c]
#test.id_to_pep[0x00000027]
#test.lh_pair_id_to_light_heavy_pair[0x0000005e]
#scan_adresses = getScanAdresses(GAPDH_VGVNGFGR[:msOrder])
#ms1_indices = Set(map(x->x.ms1, heavy_adresses))∩Set(map(x->x.ms1, light_adresses))
#test.lh_pair_id_to_light_heavy_pair[0x0000005c]
#function getScanAdresses(scan_order::Arrow.Primitive{Union{Missing, Int32}, Vector{Int32}})
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
#adresses = getScanAdresses(GAPDH_VGVNGFGR[:msOrder]) #P#test.lh_pair_id_to_light_heavy_pair[0x0000005c].light_scan_idxs)
#h = adresses[test.lh_pair_id_to_light_heavy_pair[0x0000005c].heavy_scan_idxs]
#l = adresses[test.lh_pair_id_to_light_heavy_pair[0x0000005c].light_scan_idxs]
#h = getScanAdresses(#test.lh_pair_id_to_light_heavy_pair[0x0000005c].heavy_scan_idxs)
export getScanAdresses
function getScanCycleUnion(scan_adresses_1::Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}}, 
                           scan_adresses_2::Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}}
                          )
    sort(collect(Set(map(x->x.ms1, scan_adresses_1))
                ∩Set(map(x->x.ms1, scan_adresses_2))
                )
        )
end
export getScanCycleUnion
sunion = getScanCycleUnion(l, h)
#light_adresses = getScanAdresses(table.msOrder)[getSub(lightMZ, Float32(10.0), table.precursorMZ)]
#heavy_adresses = getScanAdresses(table.msOrder)[getSub(heavyMZ, Float32(10.0), table.precursorMZ)]
#ms1_indices = getScanCycleUnion(light_adresses, heavy_adresses)
test_indices = Vector{Int64}([1, 2, 10, 11, 12, 13, 60, 61, 62, 63, 64, 65, 66, 67, 68, 100])
test_indices = Vector{Int64}([1, 2, 10, 11, 12, 13, 60, 61, 62, 63, 64, 65, 66, 67, 68, 100])
##getIntegrationBounds(ms1_indices)
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
    (lower_bound = best_start, upper_bound = best_stop+1)
end
export getIntegrationBounds
getIntegrationBounds(sunion)
#getIntegrationBounds(ms1_indices)

#limit this by integration boundaires
#light_scans = [x.scan_index for x in light_adresses if x.ms1∈ms1_indices]
#heavy_scans = [x.scan_index for x in light_adresses if x.ms1∈ms1_indices]
#Function that given a sorted list of transitions gets all the hits
#of_eltype(Float32, table.intensities[1]))
#function getHits(mass_list::Vector{Float32}, ppm::Float64, masses::MappedArray{Float32, 1, Vector{Union{Missing, Float32}}, MappedArrays.var"#7#9"{Float32}, MappedArrays.var"#8#10"{Union{Missing, Float32}}}, 
#    intensities::MappedArray{Float32, 1, Vector{Union{Missing, Float32}}, MappedArrays.var"#7#9"{Float32}, MappedArrays.var"#8#10"{Union{Missing, Float32}}})
#function getHits(mass_list::Vector{Float32}, ppm::Float32, masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}})



