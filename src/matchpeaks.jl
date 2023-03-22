GAPDH_VGVNGFGR = Arrow.Table("./data/parquet/GAPDH_VGVNGFGR.arrow")
NRF2_Survey = Arrow.Table("./data/parquet/Nrf2_SQ_052122_survey.arrow")
vg_mass = GAPDH_VGVNGFGR.masses[36]
vg_intensity = GAPDH_VGVNGFGR.intensities[36]
test_peps = ["VGVNGFGR", "VGVNGFGR[+10.008269]"]
test_features = [x for x in flatten([ getTransitions(Precursor(pep[2], prec_id = UInt32(pep[1])), charge = UInt8(2), y_start = 2, b_start = 2, ppm = Float32(40.0)) for pep in enumerate(test_peps)])]
#getTransitions(Precursor("VGVNGFGR", charge = UInt8(2)))
#getTransitions(Precursor("VGVNGFGR[+10.008269]", charge = UInt8(2)))
test_peps = ["MGKVKVGVNG", "FGRIGRLVTR", "AAFNSGKVDI", "VAINDPFIDL", "NYMVYMFQYD", "STHGKFHGTV", "KAENGKLVIN"]
using Base.Iterators

#31
heavy = getTransitions(Precursor("VGVNGFGR[+10.008269]", prec_id = UInt32(2)), charge = UInt8(1), y_start = 2, b_start = 2, ppm = Float32(10.0))
sort!(heavy, by = transition -> getMZ(transition))
results = UnorderedDictionary{UInt32, FastXTandem}()
getHits(δs, heavy, GAPDH_VGVNGFGR.masses[31], GAPDH_VGVNGFGR.intensities[31])
results = UnorderedDictionary{UInt32, FastXTandem}()
map(Score, ScoreFragmentMatches(results, getHits(δs, heavy, GAPDH_VGVNGFGR.masses[120], GAPDH_VGVNGFGR.intensities[120])))
results = UnorderedDictionary{UInt32, FastXTandem}()
map(Score, ScoreFragmentMatches(results, getHits(δs, heavy, GAPDH_VGVNGFGR.masses[31], GAPDH_VGVNGFGR.intensities[31])))
map(Score, score)
#36
light = getTransitions(Precursor("VGVNGFGR", prec_id = UInt32(1)), charge = UInt8(1), y_start = 2, b_start = 2, ppm = Float32(20.0))
sort!(light, by = transition -> getMZ(transition))
results = UnorderedDictionary{UInt32, FastXTandem}()
getHits(δs, light, GAPDH_VGVNGFGR.masses[36], GAPDH_VGVNGFGR.intensities[36])
test_features = [x for x in flatten([ getTransitions(Precursor(pep[2], prec_id = UInt32(pep[1])), charge = UInt8(2), y_start = 2, b_start = 2) for pep in enumerate(test_peps)])]
#test_features = getTransitions(Precursor("PEPTIDE", prec_id = UInt32(1)), charge = UInt8(2), y_start = 2, b_start = 2)
#append!(test_features, getTransitions(Precursor("PEPTIDE", prec_id = UInt32(2)), charge = UInt8(1), y_start = 2, b_start = 2))
sort!(test_features, by = transition -> getMZ(transition))
lower_bound = 400.0
upper_bound = 1000.0
n = 300
test_masses = sort(Vector{Union{Missing, Float32}}(lower_bound .+ rand(Float32, n) .* (upper_bound - lower_bound)))
for t in test_features[1:20]
    push!(test_masses, getMZ(t))
end
test_masses = sort(test_masses)
δs = map(x->x*0.0001, 0:0)
#δs = map(x->x*0.0001, -75:75)
light_heavy = [x for x in flatten([light, heavy])]
sort!(light_heavy, by = transition -> getMZ(transition))
scores = []
for i in 1:length(NRF2_Survey.scanType)
    println("i ", i)
    if coalesce(NRF2_Survey.msOrder[i] == 2, false)
        results = UnorderedDictionary{UInt32, FastXTandem}()
        score = map(Score, ScoreFragmentMatches(results, getHits(δs, heavy, NRF2_Survey.masses[i], NRF2_Survey.intensities[i])))
        if length(score)>0
            println(score)
            push!(scores, score)
        end
    end
end

for i in 1:length(NRF2_Survey.scanType)
    if GAPDH_VGVNGFGR.msOrder[i] == 2
        results = UnorderedDictionary{UInt32, FastXTandem}()
        println(GAPDH_VGVNGFGR.precursorMZ[i])
        println(map(Score, ScoreFragmentMatches(results, getHits(δs, light_heavy, GAPDH_VGVNGFGR.masses[i], GAPDH_VGVNGFGR.intensities[i]))))
    end
end
results = UnorderedDictionary{UInt32, FastXTandem}()
map(Score, ScoreFragmentMatches(results, getHits(δs, transitions, vg_mass, vg_intensity)))

mutable struct FragmentMatch
    transition::Transition
    intensity::Float32
    mass::Float32
    count::UInt8
    peak_ind::Int64
end

FragmentMatch() = FragmentMatch(Transition(), Float32(0), Float32(0), UInt8(0), 0)

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
results = UnorderedDictionary{UInt32, FastXTandem}()
#Base.zero(::Type{FragmentMatch}) = FragmentMatch()

function getNearest(transition::Transition, masses::Vector{Union{Missing, Float32}}, peak::Int; δ = 0.01)
    smallest_diff = abs(masses[peak]+ δ - getMZ(transition))
    best_peak = peak
    i = 0
    @inbounds while masses[peak + 1 + i]+ δ <= getHigh(transition)
        new_diff = abs(masses[peak  + 1 + i] + δ - getMZ(transition))
        if new_diff < smallest_diff
            smallest_diff = new_diff
            best_peak = peak + 1 + i
        end
        i+=1
    end
    best_peak
end

function setFragmentMatch!(match::FragmentMatch, transition::Transition, mass::Float32, intensity::Float32, peak_ind::Int64)
    match.transition = transition
    match.intensity += intensity
    match.mass = mass
    match.count += 1
    match.peak_ind = peak_ind
end


function getHits(δs, Transitions::Vector{Transition}, 
    masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}})
    #hits = zeros(FragmentMatch, length(Transitions))
    #hits =  [FragmentMatch() for transition=1:length(Transitions)]
    hits = Vector{FragmentMatch}()
    for δ in δs
        peak, transition, match = 1, 1, 1
        while (peak <= length(masses)) & (transition <= length(Transitions))
            #Is the peak within the tolerance of the transition m/z?
            if (masses[peak]+ δ >= getLow(Transitions[transition]))
                if (masses[peak]+ δ <= getHigh(Transitions[transition]))
                    best_peak = getNearest(Transitions[transition], masses, peak, δ=δ)
                    if !isassigned(hits, match)
                         push!(hits, FragmentMatch())
                    end
                    setFragmentMatch!(hits[match], Transitions[transition], masses[best_peak], intensities[best_peak], best_peak)
                    transition += 1
                    match += 1
                    continue
                end
                transition += 1
                continue
            end
            peak+=1
        end
    end
    hits
end

function ScoreFragmentMatches(results::UnorderedDictionary{UInt32, FastXTandem}, matches::Vector{FragmentMatch})
    #UnorderedDictionary{UInt32, FastXTandem}()
    for match in matches
        if !isassigned(results, getPrecID(match.transition))
            insert!(results, getPrecID(match.transition), FastXTandem())
        end
        ModifyFeatures!(results[getPrecID(match.transition)], match.transition, match.mass, match.intensity)
    end
    results
end

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
    (lgamma(score.b_count[1]*score.y_count[1]) + 
     log(score.y_int[1]*score.b_int[1])
    )
end

