using Tables
#include("precursor.jl")
#dict = Dict("customer age" => [15, 20, 25],
#                   "first name" => ["Rohit", "Rahul", "Akshat"])
#DataFrame(dict)
#table.precursorMZ
#def getSub(x, low)
#    Set(findall(x->x<478, skipmissing(table.precursorMZ)));
#b = Set(findall(x->x>477, skipmissing(table.precursorMZ)));
lightMZ = getMZ(Precursor(getResidues("VGVNGFGR"), UInt8(2)))
heavyMZ = getMZ(Precursor(getResidues("VGVNGFGR[+10.008269]"), UInt8(2)))

#getTransitions(getMZ(Precursor("VGVNGFGR", charge = UInt8(2))))
#function getSub(mean::Float32, array::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}; ppm::Float32)
#    findall(x->coalesce(abs(mean-x)<((mean/1000000.0)*ppm), false), array)
#end

#For testing purposes easier to accept AbstractRray rather 
#than concrete type even if less performant
function getSub(mean::Float32, array::AbstractArray; ppm::Float32)
    findall(x->coalesce(abs(mean-x)<((mean/1000000.0)*ppm), false), array)
end

export getSub
#ms1_indices = Set(map(x->x.ms1, heavy_adresses))∩Set(map(x->x.ms1, light_adresses))

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


function getHits(Transitions::Vector{Transition}, 
    masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}})
    results = Dict{Int64, NamedTuple{(:mass, :error, :index), Tuple{Vector{Float32}, Vector{Float32}, Vector{Int}}}}()

    peak, transition = 1, 1
    while (peak <= length(masses)) & (transition <= length(Transitions))
        #Is the peak within the tolerance of the transition m/z?
        if (masses[peak] >= getLow(Transitions[transition])) & (masses[peak] <= getHigh(Transitions[transition]))

            #There could be multiple peaks in the tolerance and we want to 
            #choose the one that is closest in m/z to the the current transition m/z. 
            smallest_diff = abs(masses[peak] - getMZ(Transitions[transition]))
            best_peak = peak
            i = 0
            @inbounds while masses[peak + 1 + i] <= getHigh(Transitions[transition])
                new_diff = abs(masses[peak  + 1 + i] - getMZ(Transitions[transition]))
                if new_diff < smallest_diff
                    smallest_diff = new_diff
                    best_peak = peak + 1 + i
                end
                i+=1
            end
            prec_id = getPrecID(Transitions[transition])
            if !haskey(results, prec_id)
                results[prec_id] = (mass = Float32[], error = Float32[], index = Int[])
            end
            push!(results[prec_id].mass, getMZ(Transitions[transition]))
            push!(results[prec_id].error, 1e6*(masses[best_peak] - getMZ(Transitions[transition])))/getMZ(Transitions[transition])
            if best_peak == 70
                println(transition)
                println(Transitions[transition])
                println(masses[best_peak])
            end
            #push!(results[prec_id].type, getIonType(Transitions[transition]))
            push!(results[prec_id].index, best_peak)
            transition += 1
            continue
        end
        peak+=1
    end
    results
end
function getNearest(Transitions::Vector{Transition}, masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}}, peak::Int, transition::Int; δ = 0.01)
    smallest_diff = abs(masses[peak]+ δ - getMZ(Transitions[transition]))
    best_peak = peak
    i = 0
    @inbounds while masses[peak + 1 + i]+ δ <= getHigh(Transitions[transition])
        new_diff = abs(masses[peak  + 1 + i] + δ - getMZ(Transitions[transition]))
        if new_diff < smallest_diff
            smallest_diff = new_diff
            best_peak = peak + 1 + i
        end
        i+=1
    end
    best_peak
end

test = Dict{Int64, NamedTuple{(:b_count, :b_int, :y_count, :y_int, :error, :intensities, :test), Tuple{Vector{Int64},Vector{Float32}, Vector{Int64}, Vector{Float32}, Vector{Float64}, Vector{Float32}, Vector{Char}}}}()
for i in 1:10
    test[i] = (b_count = [0], b_int = [Float32(0)], y_count = [0], y_int = [Float32(0)], error = [Float64(0)], intensities = Float32[], test = Char[])
end
δs = map(x->x*0.0001, -75:75)
δs = map(x->x*0.0001, 0:0)
results = UnorderedDictionary{UInt32, FastXTandem}()
abstract type Feature end

mutable struct FastXTandem <: Feature
    b_count::Int64
    b_int::Float32
    y_count::Int64
    y_int::Float32
    last_y::UInt8
    y_start::UInt8
    longest_y::Int8
    error::Float64
end

FastXTandem() = FastXTandem(0, Float32(0), 0, Float32(0), UInt8(0), UInt8(0), UInt8(0), Float64(0))

function ModifyFeatures!(score::FastXTandem, transition::Transition, mass::Union{Missing, Float32}, intensity::Union{Missing, Float32})
    if getIonType(transition) == 'b'
        score.b_count += 1
        score.b_int += intensity
    elseif getIonType(transition) == 'y'
        score.y_count += 1
        score.y_int += intensity
        if score.longest_y == 0
            score.y_start = getInd(transition)
            score.last_y = getInd(transition) - 1
        end
        if (getInd(transition)-score.last_y) == 1
            if (getInd(transition) - score.y_start) > score.longest_y
                score.longest_y = getInd(transition) - score.y_start
            else
                score.y_start = getInd(transition)
            end
        end
    end
    score.error += abs(mass - getMZ(transition))
    #push!(results[prec_id].intensities, (intensities[best_peak]))
    #push!(results[prec_id].test, getIonType(Transitions[transition]))
end

using SpecialFunctions
function Score(score::FastXTandem)
    (lgamma(score.b_count[1]) + 
     lgamma(score.y_count[1])+ 
     log(score.y_int[1]) + 
     log(score.b_int[1])
    )
end
function getHits!(results, δs, Transitions::Vector{Transition}, 
    masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}})
    for δ in δs
        peak, transition = 1, 1
        while (peak <= length(masses)) & (transition <= length(Transitions))
            #Is the peak within the tolerance of the transition m/z?
            if (masses[peak]+ δ >= getLow(Transitions[transition])) & (masses[peak]+ δ <= getHigh(Transitions[transition]))
                best_peak = getNearest(Transitions, masses, intensities, peak, transition, δ=δ)
                prec_id = getPrecID(Transitions[transition])
                if !isassigned(results, prec_id)
                    insert!(results, prec_id, FastXTandem())
                end
                ModifyFeatures!(results[prec_id], Transitions[transition], masses[best_peak], intensities[best_peak]);
                transition += 1
                continue
            end
            peak+=1
        end
    end
end

Base.zero(::Type{FragmentMatch}) = FragmentMatch()
function getHits(δs, Transitions::Vector{Transition}, 
    masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}})
    #hits = zeros(FragmentMatch, length(Transitions))
    hits =  [FragmentMatch() for x in 1:length(Transitions)]
    for δ in δs
        peak, transition = 1, 1
        while (peak <= length(masses)) & (transition <= length(Transitions))
            #Is the peak within the tolerance of the transition m/z?
            if (masses[peak]+ δ >= getLow(Transitions[transition])) & (masses[peak]+ δ <= getHigh(Transitions[transition]))
                best_peak = getNearest(Transitions, masses, intensities, peak, transition, δ=δ)
                prec_id = getPrecID(Transitions[transition])
                hits[transition].transition = Transitions[transition]
                hits[transition].intensity += intensities[best_peak]
                hits[transition].count += 1
                hits[transition].scan_id = best_peak
                transition += 1
                continue
            end
            peak+=1
        end
    end
    hits
end



test_peps = ["MGKVKVGVNG", "FGRIGRLVTR", "AAFNSGKVDI", "VAINDPFIDL", "NYMVYMFQYD", "STHGKFHGTV", "KAENGKLVIN"]
using Base.Iterators
test_features = [x for x in flatten([ getTransitions(Precursor(pep[2], prec_id = UInt32(pep[1])), charge = UInt8(2), y_start = 2, b_start = 2) for pep in enumerate(test_peps)])]
#test_features = getTransitions(Precursor("PEPTIDE", prec_id = UInt32(1)), charge = UInt8(2), y_start = 2, b_start = 2)
#append!(test_features, getTransitions(Precursor("PEPTIDE", prec_id = UInt32(2)), charge = UInt8(1), y_start = 2, b_start = 2))
sort!(test_features, by = transition -> getMZ(transition))
lower_bound = 400.0
upper_bound = 1000.0
n = 1000
test_masses = sort(Vector{Union{Missing, Float32}}(lower_bound .+ rand(Float32, n) .* (upper_bound - lower_bound)))
for t in test_features
    push!(test_masses, getMZ(t))
end
test_masses = sort(test_masses)
export getHits
function getMzFeatures(mass_list::Vector{Float32}, ppm::Float32)
     (map(x->(low=x-(x/Float32(1000000.0))*ppm, mass = x, high=x+(x/Float32(1000000.0))*ppm), mass_list))
end
export getMzFeatures
test_masses =  Vector{Union{Missing, Float32}}(
    [151.67221f0, missing, 894.0937f0, 894.0938f0, 894.0939f0])

lower_bound = 400.0
upper_bound = 1000.0
n = 1000
test_masses = sort(Vector{Union{Missing, Float32}}(lower_bound .+ rand(Float32, n) .* (upper_bound - lower_bound)))
test_features = getMzFeatures(sort(Vector{Float32}(lower_bound .+ rand(Float32, 100) .* (upper_bound - lower_bound))),
Float32(20.0))
#@btime getHitsOld(test_features, Float32(20.0), test_masses, test_masses)


test_intensities =  Vector{Union{Missing, Float32}}([missing for x in test_masses])
test_ppm = Float32(20.0)
test_mz = getMzFeatures(Vector{Float32}([151.67221f0, 700.0, 894.0938f0]), Float32(20.0))
test_mz = getMzFeatures(Vector{Float32}(sort(rand(1000)*1000)), Float32(20.0))
#getHits(test_mz, test_ppm, test_masses, test_intensities)

#=
function integrate2(MzFeatures, masses, intensities)
    n = length(masses)
    m = length(MzFeatures)
    integrated_intensities = zeros(Float32, m)
    peak_idx = 1
    for i in 1:n
        mz = masses[i]
        while peak_idx <= m && MzFeatures[peak_idx].mass < mz
            peak_idx += 1
        end
        if peak_idx > m
            break
        end
        if MzFeatures[peak_idx].mass == mz
            intensity = intensities[i]
            while peak_idx <= m && MzFeatures[peak_idx].mass == mz
                integrated_intensities[peak_idx] += intensity
                peak_idx += 1
            end
        end
    end
    return integrated_intensities
end

function integrate2_turbo(MzFeatures, masses, intensities)
    n = length(masses)
    m = length(MzFeatures)
    integrated_intensities = zeros(Float32, m)
    peak_idx = 1
    @turbo for i=1:n_peaks
        mz = peak_list[i].mz
        intensity = intensities[i]
        peak_idx = findfirst(x -> x.mass >= mz, MzFeatures)
        if peak_idx === nothing
            continue
        end
        while peak_idx <= m && (MzFeatures[peak_idx]).mass == mz
            integrated_intensities[peak_idx] = add_fast(intensity, integrated_intensities[peak_idx])
            peak_idx += 1
        end
    end
    return integrated_intensities
end
=#