using Tables, Arrow
include("precursor.jl")
#dict = Dict("customer age" => [15, 20, 25],
#                   "first name" => ["Rohit", "Rahul", "Akshat"])
#DataFrame(dict)
#table.precursorMZ
#def getSub(x, low)
#    Set(findall(x->x<478, skipmissing(table.precursorMZ)));
#b = Set(findall(x->x>477, skipmissing(table.precursorMZ)));
lightMZ = getMZ(Precursor(getResidues("VGVNGFGR"), UInt8(2)))
heavyMZ = getMZ(Precursor(getResidues("VGVNGFGR[+10.008269]"), UInt8(2)))


function getSub(mean::Float32, ppm::Float32, array::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}})
    findall(x->coalesce(abs(mean-x)<((mean/1000000.0)*ppm), false), array)
end
export getSub
#ms1_indices = Set(map(x->x.ms1, heavy_adresses))∩Set(map(x->x.ms1, light_adresses))

function getScanAdresses(scan_order::Arrow.Primitive{Union{Missing, Int32}, Vector{Int32}})
    scan_adresses = Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}}(undef,length(scan_order)) 
    ms1 = 0
    msn = 0
    for scan in enumerate(scan_order)
        if scan[2] == 1
            ms1 += 1
            scan_adresses[scan[1]] = (scan_index = scan[1], ms1 = ms1, msn = msn)
            msn=0
        else
            scan_adresses[scan[1]] = (scan_index = scan[1], ms1 = ms1, msn = msn)
            msn+=1
        end
    end
    scan_adresses
end
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
#light_adresses = getScanAdresses(table.msOrder)[getSub(lightMZ, Float32(10.0), table.precursorMZ)]
#heavy_adresses = getScanAdresses(table.msOrder)[getSub(heavyMZ, Float32(10.0), table.precursorMZ)]
#ms1_indices = getScanCycleUnion(light_adresses, heavy_adresses)
test_indices = Vector{Int64}([1, 2, 10, 11, 12, 13, 60, 61, 62, 63, 64, 65, 66, 67, 68, 100])
test_indices = Vector{Int64}([1, 2, 10, 11, 12, 13, 60, 61, 62, 63, 64, 65, 66, 67, 68, 100])
##getIntegrationBounds(ms1_indices)
function getIntegrationBounds(scan_indices::Vector{Int64}; max_gap_size::Int = 10)
    if length(scan_indices)==1
        return (1, 1)
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
    (best_start, best_stop+1)
end
export getIntegrationBounds
#getIntegrationBounds(ms1_indices)

#limit this by integration boundaires
#light_scans = [x.scan_index for x in light_adresses if x.ms1∈ms1_indices]
#heavy_scans = [x.scan_index for x in light_adresses if x.ms1∈ms1_indices]
#Function that given a sorted list of transitions gets all the hits
#of_eltype(Float32, table.intensities[1]))
#function getHits(mass_list::Vector{Float32}, ppm::Float64, masses::MappedArray{Float32, 1, Vector{Union{Missing, Float32}}, MappedArrays.var"#7#9"{Float32}, MappedArrays.var"#8#10"{Union{Missing, Float32}}}, 
#    intensities::MappedArray{Float32, 1, Vector{Union{Missing, Float32}}, MappedArrays.var"#7#9"{Float32}, MappedArrays.var"#8#10"{Union{Missing, Float32}}})
#function getHits(mass_list::Vector{Float32}, ppm::Float32, masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}})
function getHits(MzFeatures::Vector{NamedTuple{(:low, :mass, :high), Tuple{Float32, Float32, Float32}}}, 
    ppm::Float32, masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}})
    #MzFeatures = getMzFeatures(mass_list, ppm)
    feature = 1
    peak = 1
    #println("test");
    while (peak <= length(masses)) & (feature <= length(MzFeatures))
        if coalesce(masses[peak], MzFeatures[feature].low) > MzFeatures[feature].low
            if coalesce(masses[peak], MzFeatures[feature].high) < MzFeatures[feature].high
                #"masses[peak] is in the range(low, high), "
                #There could be multiple peaks in the tolerance and we want to 
                #choose the one that is closest in mass
                smallest_diff = abs(coalesce(masses[peak],MzFeatures[feature].mass) - MzFeatures[feature].mass)
                i = 0
                @inbounds while coalesce(masses[peak+1+i], MzFeatures[feature].high) < MzFeatures[feature].high
                    new_diff = abs(coalesce(masses[peak+1+i],-MzFeatures[feature].mass)  - MzFeatures[feature].mass)
                    if new_diff < smallest_diff
                        smallest_diff = new_diff
                        peak = peak+1+i
                        i = 0
                    end
                    i+=1
                end
                #println(MzFeatures[feature].mass, " has mass ", masses[peak], " and intensity ", intensities[peak])
                feature += 1
                continue
            end
            feature += 1
            continue
        end
        peak+=1
        #println(peak)
        #println(feature)
    end
end
export getHits
function getMzFeatures(mass_list::Vector{Float32}, ppm::Float32)
     (map(x->(low=x-(x/Float32(1000000.0))*ppm, mass = x, high=x+(x/Float32(1000000.0))*ppm), mass_list))
end
export getMzFeatures
test_masses =  Vector{Union{Missing, Float32}}(
    [151.67221f0, missing, 894.0937f0, 894.0938f0, 894.0939f0])


test_intensities =  Vector{Union{Missing, Float32}}([missing for x in test_masses])
test_ppm = Float32(20.0)
test_mz = getMzFeatures(Vector{Float32}([151.67221f0, 700.0, 894.0938f0]), Float32(20.0))
test_mz = getMzFeatures(Vector{Float32}(sort(rand(1000)*1000)), Float32(20.0))
getHits(test_mz, test_ppm, test_masses, test_intensities)