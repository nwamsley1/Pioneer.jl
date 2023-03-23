using Titus
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
@testset "Titus.jl" begin
    
end