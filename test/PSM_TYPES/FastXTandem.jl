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
    matchPeaks!(test_matches, test_t, test_masses, test_intensities, 0.0, UInt32(0), UInt32(0))
    unscored_PSMs = UnorderedDictionary{UInt32, FastXTandem}()
    ScoreFragmentMatches!(unscored_PSMs, test_matches)
    scored_PSMs = makePSMsDict(FastXTandem())
    Score!(scored_PSMs, unscored_PSMs)

    #No matched y ions so should be zero
    @test scored_PSMs[:y_ladder][1] == 0
    #log(3!) + log(0!) + log(3000*0) = log(6) + 1 = 1.791... 
    @test Tol(scored_PSMs[:hyperscore][1], 1.79189)
    #abs(324.15536f0 -  324.15636f0)*2 + abs(538.2871f0 - 538.27739f0) = 0.01171875f0
    @test Tol(scored_PSMs[:error][1], 0.01171875)

    test_t = sort!(getTransitions(Precursor("PEPTIDE")), by = transition->getMZ(transition))
    #=8-element Vector{Transition}:
 Transition(MzFeature(324.15536f0, 324.1489f0, 324.16183f0), 0x00000000, 'b', 0x03, 0x01, 0x00)
 Transition(MzFeature(376.17142f0, 376.16388f0, 376.17896f0), 0x00000000, 'y', 0x03, 0x01, 0x00)
 Transition(MzFeature(425.20306f0, 425.19455f0, 425.21158f0), 0x00000000, 'b', 0x04, 0x01, 0x00)
 Transition(MzFeature(477.2191f0, 477.20953f0, 477.22864f0), 0x00000000, 'y', 0x04, 0x01, 0x00)
 Transition(MzFeature(538.2871f0, 538.27637f0, 538.29785f0), 0x00000000, 'b', 0x05, 0x01, 0x00)
 Transition(MzFeature(574.2718f0, 574.2603f0, 574.28326f0), 0x00000000, 'y', 0x05, 0x01, 0x00)
 Transition(MzFeature(653.314f0, 653.30096f0, 653.3271f0), 0x00000000, 'b', 0x06, 0x01, 0x00)
 Transition(MzFeature(703.3144f0, 703.30035f0, 703.3284f0), 0x00000000, 'y', 0x06, 0x01, 0x00)
    =#
    #Should have one peak to match every transition 
    test_masses = Vector{Union{Missing, Float32}}([324.15536f0, 376.17142f0, 425.20306f0, 477.2191f0, 
    538.2871f0, 574.2718f0, 653.314f0, 703.3144f0])
    test_intensities = Vector{Union{Missing, Float32}}(map(x->1000, test_masses))
    test_matches = Vector{FragmentMatch}()
    matchPeaks!(test_matches, test_t, test_masses, test_intensities, 0.0, UInt32(0), UInt32(0))
    unscored_PSMs = UnorderedDictionary{UInt32, FastXTandem}()
    ScoreFragmentMatches!(unscored_PSMs, test_matches)
    scored_PSMs = makePSMsDict(FastXTandem())
    Score!(scored_PSMs, unscored_PSMs)

    #No matched y ions so should be zero
    @test scored_PSMs[:y_ladder][1] == 4
    #log(4!) + log(4!) + log(4000*4000) = 22.944206940899946
    @test Tol(scored_PSMs[:hyperscore][1], 22.944206940899946)
    #masses were exact so should be about zero. 
    @test Tol(scored_PSMs[:error][1], 0)

    #Make sure y-ladder counting still works when there is a gap
    test_t = sort!(getTransitions(Precursor("PEPTIDE"))[[1, 2, 3, 4, 5, 7, 8]], by = transition->getMZ(transition))
    #=7-element Vector{Transition}:
7-element Vector{Transition}:
 Transition(MzFeature(324.15536f0, 324.1489f0, 324.16183f0), 0x00000000, 'b', 0x03, 0x01, 0x00)
 Transition(MzFeature(376.17142f0, 376.16388f0, 376.17896f0), 0x00000000, 'y', 0x03, 0x01, 0x00)
 Transition(MzFeature(425.20306f0, 425.19455f0, 425.21158f0), 0x00000000, 'b', 0x04, 0x01, 0x00)
 Transition(MzFeature(538.2871f0, 538.27637f0, 538.29785f0), 0x00000000, 'b', 0x05, 0x01, 0x00)
 Transition(MzFeature(574.2718f0, 574.2603f0, 574.28326f0), 0x00000000, 'y', 0x05, 0x01, 0x00)
 Transition(MzFeature(653.314f0, 653.30096f0, 653.3271f0), 0x00000000, 'b', 0x06, 0x01, 0x00)
 Transition(MzFeature(703.3144f0, 703.30035f0, 703.3284f0), 0x00000000, 'y', 0x06, 0x01, 0x00)
    =#
    #Should have one peak to match every transition 
    test_masses = Vector{Union{Missing, Float32}}([324.15536f0, 376.17142f0, 425.20306f0, 
    538.2871f0, 574.2718f0, 653.314f0, 703.3144f0])
    test_intensities = Vector{Union{Missing, Float32}}(map(x->1000, test_masses))
    test_matches = Vector{FragmentMatch}()
    matchPeaks!(test_matches, test_t, test_masses, test_intensities, 0.0, UInt32(0), UInt32(0))
    unscored_PSMs = UnorderedDictionary{UInt32, FastXTandem}()
    ScoreFragmentMatches!(unscored_PSMs, test_matches)
    scored_PSMs = makePSMsDict(FastXTandem())
    Score!(scored_PSMs, unscored_PSMs)

    @test scored_PSMs[:y_ladder][1] == 2
    #log(4!) + log(3!) + log(4000*3000) = 21.270230507328275
    @test Tol(scored_PSMs[:hyperscore][1], 21.270230507328275)
    #masses were exact so should be about zero. 
    @test Tol(scored_PSMs[:error][1], 0)


    #When charge state is 2, need to make sure we are counting the Y-ladder correctly
    #and not double counting +1 and +2 charge states or isotopes.  
    #So if we see the y4+1 and the y4+2 ion, we should only add one to the count of y-ions,
    #but we should add the intensities of both the y4+1 and y4+2 ions to the y intensity sum. 
    #Likewise, if we see both +0 and +1 isotopes for y4+1, we should sum both intensities but only
    #add one to the count of y-ions. 
    test_t = sort!(getTransitions(Precursor("PEPTIDE"), UInt8[1,2], UInt8[0,1]), by = transition->getMZ(transition))
    test_masses = Vector{Union{Missing, Float32}}(map(x->getMZ(x), test_t))
    test_intensities = Vector{Union{Missing, Float32}}(map(x->1000, test_masses))
    test_matches = Vector{FragmentMatch}()
    matchPeaks!(test_matches, test_t, test_masses, test_intensities, 0.0, UInt32(0), UInt32(0))
    unscored_PSMs = UnorderedDictionary{UInt32, FastXTandem}()
    ScoreFragmentMatches!(unscored_PSMs, test_matches)
    scored_PSMs = makePSMsDict(FastXTandem())
    Score!(scored_PSMs, unscored_PSMs)

    #should hit all 4 y ions
    @test scored_PSMs[:y_ladder][1] == 4
    #log(16!) + log(16!) + log(16000*16000) = 80.70440821460518
    @test Tol(scored_PSMs[:hyperscore][1], 80.70440821460518)
    #masses were exact so should be about zero. 
    @test Tol(scored_PSMs[:error][1], 0)
end