using Test
using Tables, Arrow, Base.Iterators

function Tol(a, b, ppm = 2)
    abs(a-b)<=(ppm*minimum((a, b))/1000000)
end

@testset "FastXTandem.jl" begin
    #Same first example as in test/matchpeaks.jl
    test_t = sort!(getTransitions(Precursor("PEPTIDE"))[[1, 1, 2, 3]], by = transition->getMZ(transition))
    test_masses = Vector{Union{Missing, Float32}}([324.1490f0, 324.15636f0, 324.15736f0, 324.16185f0, 420, 421, 500, 538, 538.27636f0, 538.27739f0, 600])
    test_intensities = Vector{Union{Missing, Float32}}(map(x->1000, test_masses))
    #=map(x->x.mz, test_t)
    4-element Vector{MzFeature}:
    MzFeature(324.15536f0, 324.1489f0, 324.16183f0) #Maps to second mass
    MzFeature(324.15536f0, 324.1489f0, 324.16183f0) #Maps to second mass
    MzFeature(425.20306f0, 425.19455f0, 425.21158f0) #Doesn't map to a mass
    MzFeature(538.2871f0, 538.27637f0, 538.29785f0)=# #Maps to second to last mass (10) .

    #test_t has a duplicate transition which should match to the same peaks. 
    #Multiple transitions can map to the same peak but multiple peaks cannot map to the same transition
    #Transitions 1 and 2 should match to peak 2, and transition 4 should map to peak 10. Transition 3 should not map to any peak. 
    test_matches = Vector{FragmentMatch}()
    matchPeaks!(test_matches, test_t, test_masses, test_intensities, 0.0)
    unscored_PSMs = UnorderedDictionary{UInt32, FastXTandem}()
    ScoreFragmentMatches!(unscored_PSMs, test_matches)
    scored_PSMs = makePSMsDict(FastXTandem())
    Score!(scored_PSMs, unscored_PSMs)
    scored_PSMs
end