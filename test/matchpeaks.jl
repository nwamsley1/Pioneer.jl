#using Titus
using Test
using Tables, Arrow, Base.Iterators
#=println(pwd())
include("../src/integrate.jl")
GAPDH_VGVNGFGR = Arrow.Table("../data/parquet/GAPDH_VGVNGFGR.arrow")
GSTM1_RPWFAGNK = Arrow.Table("../data/parquet/GSTM1_RPWFAGNK.arrow")
GSTM4_VAVWGNK = Arrow.Table("../data/parquet/GSTM4_VAVWGNK.arrow")
NRF2_Survey = Arrow.Table("../data/parquet/Nrf2_SQ_052122_survey.arrow")

heavy_VGVNGFGR = getTransitions(Precursor("VGVNGFGR[+10.008269]", prec_id = UInt32(2)), charge = UInt8(1), y_start = 2, b_start = 2, ppm = Float32(10.0))
sort!(heavy_VGVNGFGR, by = transition -> getMZ(transition))
hits = getHits(δs, heavy_VGVNGFGR, GAPDH_VGVNGFGR.masses[31], GAPDH_VGVNGFGR.intensities[31])
results = UnorderedDictionary{UInt32, FastXTandem}()
map(Score, ScoreFragmentMatches(results, hits))

light_VGVNGFGR = getTransitions(Precursor("VGVNGFGR", prec_id = UInt32(1)), charge = UInt8(1), y_start = 2, b_start = 2, ppm = Float32(20.0))
sort!(light_VGVNGFGR, by = transition -> getMZ(transition))
hits = getHits(δs, light_VGVNGFGR, GAPDH_VGVNGFGR.masses[36], GAPDH_VGVNGFGR.intensities[36])
results = UnorderedDictionary{UInt32, FastXTandem}()
map(Score, ScoreFragmentMatches(results, hits))

light_heavy = [x for x in flatten([light_VGVNGFGR, heavy_VGVNGFGR])]
sort!(light_heavy, by = transition -> getMZ(transition))

scores = []
rts = []
is = []
for i in 1:length(NRF2_Survey.scanType)
    if coalesce(NRF2_Survey.msOrder[i] == 2, false)
        results = UnorderedDictionary{UInt32, FastXTandem}();
        ScoreFragmentMatches!(results, getHits(δs, light, NRF2_Survey.masses[i], NRF2_Survey.intensities[i]));
        score = map(Score, results);
        if length(score)>0
            if !isnan(score[1])
                push!(scores, score)
                push!(is, i)
                push!(rts, NRF2_Survey.retentionTime[i])
            end
        end
    end
end
rts[sortperm(scores)][end - 10:end]
scores[sortperm(scores)][end - 10:end]
is[sortperm(scores)][end - 10:end]

for i in 1:length(NRF2_Survey.scanType)
    if GAPDH_VGVNGFGR.msOrder[i] == 2
        results = UnorderedDictionary{UInt32, FastXTandem}()
        println(GAPDH_VGVNGFGR.precursorMZ[i])
        println(map(Score, ScoreFragmentMatches(results, getHits(δs, light_heavy, GAPDH_VGVNGFGR.masses[i], GAPDH_VGVNGFGR.intensities[i]))))
    end
end

test_peps = ["MGKVKVGVNG", "FGRIGRLVTR", "AAFNSGKVDI", "VAINDPFIDL", "NYMVYMFQYD", "STHGKFHGTV", "KAENGKLVIN"]
test_features = [x for x in flatten([ getTransitions(Precursor(pep[2], prec_id = UInt32(pep[1])), charge = UInt8(2), y_start = 2, b_start = 2) for pep in enumerate(test_peps)])]
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
=#
function Tol(a, b, ppm = 2)
    abs(a-b)<=(ppm*minimum((a, b))/1000000)
end

@testset "matchpeaks.jl" begin
    ##########
    #FragmentMatch
    ##########

    test_fragment_match = FragmentMatch()

    @test Tol(getMZ(test_fragment_match), 0.0)
    ##########
    #getNearest
    ##########

    #Transition with a mass of 324.15536f0. 
    #At 20 ppm the lower and upper bounds are 324.1489f0 and 324.16183f0)
    test_t = getTransitions(Precursor("PEPTIDE"))[1]
    #The first three elements in `test_masses` are within the mass tolerance.
    #Element 2 is the closest in m/z to `test_t`. Element 4 is outside the tolerance
    test_masses = Vector{Union{Missing, Float64}}([324.1490f0, 324.15636f0, 324.15736f0, 324.16185f0])
    @test getNearest(test_t, test_masses, 1) == 2

    #This time we will supply an argument peak = 3. This means the search will
    #start at the third element. In this case Element 3 is the closest in m/z to `test_t`
    #and element 4 is outside the m/z tolerance. This is the more typical case. 
    @test getNearest(test_t, test_masses, 3) == 3

    ##########
    #setFragmentMatch
    ##########    

    #Use the empty test_fragment_match from earlier
    @test Tol(test_fragment_match.count, 0)
    @test Tol(test_fragment_match.match_mz, 0)
    @test Tol(test_fragment_match.intensity, 0)
    @test Tol(test_fragment_match.peak_ind, 0)

    setFragmentMatch!([test_fragment_match], 1, test_t, test_masses[2], Float64(1), 2, UInt32(0), UInt32(0))
    @test Tol(test_fragment_match.count, 1)
    @test Tol(test_fragment_match.match_mz, test_masses[2])
    @test Tol(test_fragment_match.count, 1)
    @test Tol(test_fragment_match.peak_ind, 2)

    for i in 1:10
        setFragmentMatch!([test_fragment_match], 1, test_t, test_masses[2], Float64(1), 2, UInt32(0), UInt32(0))
    end

    @test Tol(test_fragment_match.count, 11)

    ##########
    #matchPeaks!
    ##########
    test_t = sort!(getTransitions(Precursor("PEPTIDE"))[[1, 1, 2, 3]], by = transition->getMZ(transition))
    test_masses = Vector{Union{Missing, Float64}}([324.1490f0, 324.15636f0, 324.15736f0, 324.16185f0, 420, 421, 500, 538, 538.27636f0, 538.27739f0, 600])
    test_intensities = Vector{Union{Missing, Float64}}(map(x->1000, test_masses))
    #=map(x->x.mz, test_t)
    4-element Vector{MzFeature}:
    MzFeature(324.15536f0, 324.1489f0, 324.16183f0) #Maps to second mass
    MzFeature(324.15536f0, 324.1489f0, 324.16183f0) #Maps to second mass
    MzFeature(425.20306f0, 425.19455f0, 425.21158f0) #Doesn't map to a mass
    MzFeature(538.2871f0, 538.27637f0, 538.29785f0)=# #Maps to second to last mass (10) .

    #test_t has a duplicate transition which should match to the same peaks. 
    #Multiple transitions can map to the same peak but multiple peaks cannot map to the same transition
    #Transitions 1 and 2 should match to peak 2, and transition 4 should map to peak 10. Transition 3 should not map to any peak. 
    test_matches = Vector{FragmentMatch{Float64}}()
    matchPeaks!(test_matches, test_t, test_masses, test_intensities, 0.0, UInt32(0), UInt32(0))
    @test map(x->getPeakInd(x), test_matches) == [2, 2, 10]

    #Try again but make sure we can match to the first and last peak in the spectrum. 
    test_matches = Vector{FragmentMatch{Float64}}()
    matchPeaks!(test_matches, test_t, test_masses[2:10], test_intensities, 0.0, UInt32(0), UInt32(0))
    @test map(x->getPeakInd(x), test_matches) == [1, 1, length(test_masses[2:10])]

    ##########
    #matchPeaks
    ##########
    #Make sure we get the same result as matchPeaks! when there is only one offset = 0. 
    test_matches = matchPeaks(test_t, test_masses, test_intensities)
    @test map(x->getPeakInd(x), test_matches) == [2, 2, 10]

    #Since δs is of length 3 and the offset is very small, we should match each transition three times. 
    test_matches = matchPeaks(test_t, test_masses, test_intensities, δs = map(x->x*0.0001, -1:1))
    @test map(x->getPeakInd(x), test_matches) == [2, 2, 10]
    @test map(x->getIntensity(x), test_matches) == Float64[3000, 3000, 3000]
    @test map(x->getCount(x), test_matches) == UInt8[3, 3, 3]
    
    #Since δs is of length 2, we should only match each transition twice.  
    test_matches = matchPeaks(test_t, test_masses, test_intensities, δs = map(x->x*0.0001, [-1, 1]))
    @test map(x->getPeakInd(x), test_matches) == [2, 2, 10]
    @test map(x->getIntensity(x), test_matches) == Float64[2000, 2000, 2000]
    @test map(x->getCount(x), test_matches) == UInt8[2, 2, 2]

    #Now since the offset is large, we shoudl not match any peaks. 
    test_matches = matchPeaks(test_t, test_masses, test_intensities, δs = map(x->x*1.0, [-1, 1]))
    @test length(test_matches) == 0

end