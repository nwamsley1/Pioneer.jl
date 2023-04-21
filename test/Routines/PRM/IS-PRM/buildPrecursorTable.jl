@testset "buildPrecursorTable.jl" begin
    """
    data/peptide_lists/TRANSITION_LIST_TEST1.txt

    protein_name,sequence,precursor_charge,precursor_isotope,transition_names
    PROT_A,PEPTIDER[Harg],2,0,y3+1;y7+1;y6+1;b3+1;b5+2 
    PROT_A,PEPTIDER,2,0,y3+1;y7+1;y6+1;b3+1;b5+2        #Same as above but no heavy modification
    PROT_B,PEPTIDER[Harg],2,0,y5+1;y4+1;y8+1;b4+1;b5+1  #Same as above but different transitions?
    PROT_C,PEPTIDER[Harg],2,0,y3+1;y7+1;y6+1;b3+1;b5+2 
    PROT_A,AMINEACIDK[Hlys],2,0,y5+1;y7+1;y9+1;b3+1;b5+1
    PROT_B,PEPTIC[Carb]KDEK[Hlys],2,0,y4+1;y6+1;y8+1;b3+1;b6+2;b7+1
    PROT_B,PEPTICKDEK[Hlys],2,0,y4+1;y6+1;y8+1;b3+1;b6+2;b7+1 #Same as previous but no [Carb] modification
    PROT_C,D[Hglu]RAGRACE,2,0,y3+1;y7+1;y6+1;b3+1;b5+2          
    PROT_C,D[Hglu]RAGRACE,3,0,y3+1;y7+1;y6+1;b3+1;b5+2  #Same as above but the charge state is different. 


    julia> testPtable.id_to_prot
    3-element UnorderedDictionary{UInt32, Protein}
     0x00000002 │ Protein("PROT_B", Set(UInt32[0x00000003, 0x00000001]))
     0x00000003 │ Protein("PROT_C", Set(UInt32[0x00000004, 0x00000001]))
     0x00000001 │ Protein("PROT_A", Set(UInt32[0x00000002, 0x00000001]))

    julia> testPtable.prot_to_id
    3-element UnorderedDictionary{String, UInt32}
    "PROT_C" │ 0x00000003
    "PROT_B" │ 0x00000002
    "PROT_A" │ 0x00000001

    julia> testPtable.id_to_pepGroup
    4-element UnorderedDictionary{UInt32, PeptideGroup}
    0x00000004 │ PeptideGroup(Set(UInt32[0x00000003]), Set(UInt32[0x0000000a, 0x00000009]), "DRAGRACE")
    0x00000002 │ PeptideGroup(Set(UInt32[0x00000001]), Set(UInt32[0x00000004, 0x00000003]), "AMINEACIDK")
    0x00000003 │ PeptideGroup(Set(UInt32[0x00000002]), Set(UInt32[0x00000005, 0x00000006, 0x00000007, 0x00000008]), "PEPTICKDEK")
    0x00000001 │ PeptideGroup(Set(UInt32[0x00000002, 0x00000003, 0x00000001]), Set(UInt32[0x00000002, 0x00000001]), "PEPTIDER")
    
    julia> testPtable.pepGroup_to_id
    4-element UnorderedDictionary{String, UInt32}
       "PEPTIDER" │ 0x00000001
     "AMINEACIDK" │ 0x00000002
       "DRAGRACE" │ 0x00000004
     "PEPTICKDEK" │ 0x00000003

    julia> testPtable.id_to_pep
    10-element UnorderedDictionary{UInt32, Peptide}
    0x00000005 │ Peptide("PEPTIC[Carb]KDEK", 0x00000003, Set(UInt32[0x00000005]))
    0x00000004 │ Peptide("AMINEACIDK[Hlys]", 0x00000002, Set(UInt32[0x00000004]))
    0x00000006 │ Peptide("PEPTIC[Carb]KDEK[Hlys]", 0x00000003, Set(UInt32[0x00000006]))
    0x00000007 │ Peptide("PEPTICKDEK", 0x00000003, Set(UInt32[0x00000007]))
    0x00000002 │ Peptide("PEPTIDER[Harg]", 0x00000001, Set(UInt32[0x00000002]))
    0x0000000a │ Peptide("D[Hglu]RAGRACE", 0x00000004, Set(UInt32[0x0000000a, 0x0000000c]))
    0x00000009 │ Peptide("DRAGRACE", 0x00000004, Set(UInt32[0x0000000b, 0x00000009]))
    0x00000008 │ Peptide("PEPTICKDEK[Hlys]", 0x00000003, Set(UInt32[0x00000008]))
    0x00000003 │ Peptide("AMINEACIDK", 0x00000002, Set(UInt32[0x00000003]))
    0x00000001 │ Peptide("PEPTIDER", 0x00000001, Set(UInt32[0x00000001]))
    
    julia> testPtable.id_to_prec
    12-element Dictionary{UInt32, Precursor}
    0x0000000b │ Precursor(Residue[Residue(115.02694f0), Residue(156.1011f0), Residue(71.03711f0), Residue(57.02146f0), Residue(156.1011f0), Residue(71.03711f0), Residue(103.00919f0), Residue(129.04259f0)], MzFeature(293.13635f0, 293.1305f0, 293.1422f0), 0x03, 0x00, 0x00000009, 0x0000000b)
    0x0000000c │ Precursor(Residue[Residue(121.02694f0), Residue(156.1011f0), Residue(71.03711f0), Residue(57.02146f0), Residue(156.1011f0), Residue(71.03711f0), Residue(103.00919f0), Residue(129.04259f0)], MzFeature(295.13635f0, 295.13046f0, 295.14224f0), 0x03, 0x00, 0x0000000a, 0x0000000c)
    0x00000009 │ Precursor(Residue[Residue(115.02694f0), Residue(156.1011f0), Residue(71.03711f0), Residue(57.02146f0), Residue(156.1011f0), Residue(71.03711f0), Residue(103.00919f0), Residue(129.04259f0)], MzFeature(439.2009f0, 439.1921f0, 439.2097f0), 0x02, 0x00, 0x00000009, 0x00000009)
    0x0000000a │ Precursor(Residue[Residue(121.02694f0), Residue(156.1011f0), Residue(71.03711f0), Residue(57.02146f0), Residue(156.1011f0), Residue(71.03711f0), Residue(103.00919f0), Residue(129.04259f0)], MzFeature(442.2009f0, 442.19205f0, 442.20975f0), 0x02, 0x00, 0x0000000a, 0x0000000a)
    0x00000001 │ Precursor(Residue[Residue(97.05276f0), Residue(129.04259f0), Residue(97.05276f0), Residue(101.04768f0), Residue(113.08406f0), Residue(115.02694f0), Residue(129.04259f0), Residue(156.1011f0)], MzFeature(478.7378f0, 478.7282f0, 478.74738f0), 0x02, 0x00, 0x00000001, 0x00000001)
    0x00000002 │ Precursor(Residue[Residue(97.05276f0), Residue(129.04259f0), Residue(97.05276f0), Residue(101.04768f0), Residue(113.08406f0), Residue(115.02694f0), Residue(129.04259f0), Residue(166.10938f0)], MzFeature(483.74194f0, 483.73227f0, 483.75162f0), 0x02, 0x00, 0x00000002, 0x00000002)
    0x00000003 │ Precursor(Residue[Residue(71.03711f0), Residue(131.0405f0), Residue(113.08406f0), Residue(114.04293f0), Residue(129.04259f0), Residue(71.03711f0), Residue(103.00919f0), Residue(113.08406f0), Residue(115.02694f0), Residue(128.09496f0)], MzFeature(554.26227f0, 554.25116f0, 554.2734f0), 0x0…
    0x00000004 │ Precursor(Residue[Residue(71.03711f0), Residue(131.0405f0), Residue(113.08406f0), Residue(114.04293f0), Residue(129.04259f0), Residue(71.03711f0), Residue(103.00919f0), Residue(113.08406f0), Residue(115.02694f0), Residue(136.10916f0)], MzFeature(558.2694f0, 558.25824f0, 558.2806f0), 0x02…
    0x00000007 │ Precursor(Residue[Residue(97.05276f0), Residue(129.04259f0), Residue(97.05276f0), Residue(101.04768f0), Residue(113.08406f0), Residue(103.00919f0), Residue(128.09496f0), Residue(115.02694f0), Residue(129.04259f0), Residue(128.09496f0)], MzFeature(580.2868f0, 580.2752f0, 580.2984f0), 0x02…
    0x00000008 │ Precursor(Residue[Residue(97.05276f0), Residue(129.04259f0), Residue(97.05276f0), Residue(101.04768f0), Residue(113.08406f0), Residue(103.00919f0), Residue(128.09496f0), Residue(115.02694f0), Residue(129.04259f0), Residue(136.10916f0)], MzFeature(584.29395f0, 584.2823f0, 584.3056f0), 0x0…
    0x00000005 │ Precursor(Residue[Residue(97.05276f0), Residue(129.04259f0), Residue(97.05276f0), Residue(101.04768f0), Residue(113.08406f0), Residue(160.03065f0), Residue(128.09496f0), Residue(115.02694f0), Residue(129.04259f0), Residue(128.09496f0)], MzFeature(608.79755f0, 608.7854f0, 608.8097f0), 0x0…
    0x00000006 │ Precursor(Residue[Residue(97.05276f0), Residue(129.04259f0), Residue(97.05276f0), Residue(101.04768f0), Residue(113.08406f0), Residue(160.03065f0), Residue(128.09496f0), Residue(115.02694f0), Residue(129.04259f0), Residue(136.10916f0)], MzFeature(612.8046f0, 612.79236f0, 612.8169f0), 0x0…

    julia> testPtable.sorted_prec_ids
    12-element Vector{UInt32}:
    0x0000000b
    0x0000000c
    0x00000009
    0x0000000a
    0x00000001
    0x00000002
    0x00000003
    0x00000004
    0x00000007
    0x00000008
    0x00000005
    0x00000006

    julia> testPtable.prec_id_to_transitions
    12-element Dictionary{UInt32, Vector{Transition}}
    0x00000001 │ Transition[Transition(MzFeature(419.18848f0, 419.18008f0, 419.19687f0), 0x00000001, 'y', 0x03, 0x01, 0x00), Transition(MzFeature(859.4155f0, 859.3983f0, 859.43274f0), 0x00000001, 'y', 0x07, 0x01, 0x00), Transition(MzFeature(730.3729f0, 730.35834f0, 730.3875f0), 0x00000001, 'y', 0x06, 0x0…
    0x00000002 │ Transition[Transition(MzFeature(429.19678f0, 429.1882f0, 429.20535f0), 0x00000002, 'y', 0x03, 0x01, 0x00), Transition(MzFeature(869.4238f0, 869.40643f0, 869.4412f0), 0x00000002, 'y', 0x07, 0x01, 0x00), Transition(MzFeature(740.3812f0, 740.3664f0, 740.39606f0), 0x00000002, 'y', 0x06, 0x01…
    0x00000003 │ Transition[Transition(MzFeature(549.27f0, 549.25903f0, 549.281f0), 0x00000003, 'y', 0x05, 0x01, 0x00), Transition(MzFeature(792.3555f0, 792.33966f0, 792.3714f0), 0x00000003, 'y', 0x07, 0x01, 0x00), Transition(MzFeature(1036.4801f0, 1036.4594f0, 1036.5009f0), 0x00000003, 'y', 0x09, 0x01, …
    0x00000004 │ Transition[Transition(MzFeature(557.2843f0, 557.27313f0, 557.2955f0), 0x00000004, 'y', 0x05, 0x01, 0x00), Transition(MzFeature(800.3698f0, 800.3538f0, 800.3858f0), 0x00000004, 'y', 0x07, 0x01, 0x00), Transition(MzFeature(1044.4944f0, 1044.4735f0, 1044.5153f0), 0x00000004, 'y', 0x09, 0x01…
    0x00000005 │ Transition[Transition(MzFeature(519.2773f0, 519.2669f0, 519.28766f0), 0x00000005, 'y', 0x04, 0x01, 0x00), Transition(MzFeature(792.39197f0, 792.3761f0, 792.40784f0), 0x00000005, 'y', 0x06, 0x01, 0x00), Transition(MzFeature(990.4924f0, 990.47253f0, 990.5122f0), 0x00000005, 'y', 0x08, 0x01…
    0x00000006 │ Transition[Transition(MzFeature(527.29144f0, 527.2809f0, 527.302f0), 0x00000006, 'y', 0x04, 0x01, 0x00), Transition(MzFeature(800.4061f0, 800.39014f0, 800.4221f0), 0x00000006, 'y', 0x06, 0x01, 0x00), Transition(MzFeature(998.50653f0, 998.4866f0, 998.5265f0), 0x00000006, 'y', 0x08, 0x01, …
    0x00000007 │ Transition[Transition(MzFeature(519.2773f0, 519.2669f0, 519.28766f0), 0x00000007, 'y', 0x04, 0x01, 0x00), Transition(MzFeature(735.37054f0, 735.35583f0, 735.38525f0), 0x00000007, 'y', 0x06, 0x01, 0x00), Transition(MzFeature(933.47095f0, 933.4523f0, 933.4896f0), 0x00000007, 'y', 0x08, 0x0…
    0x00000008 │ Transition[Transition(MzFeature(527.29144f0, 527.2809f0, 527.302f0), 0x00000008, 'y', 0x04, 0x01, 0x00), Transition(MzFeature(743.3847f0, 743.3698f0, 743.3996f0), 0x00000008, 'y', 0x06, 0x01, 0x00), Transition(MzFeature(941.4851f0, 941.46625f0, 941.50397f0), 0x00000008, 'y', 0x08, 0x01, …
    0x00000009 │ Transition[Transition(MzFeature(322.10675f0, 322.1003f0, 322.1132f0), 0x00000009, 'y', 0x03, 0x01, 0x00), Transition(MzFeature(762.36755f0, 762.3523f0, 762.3828f0), 0x00000009, 'y', 0x07, 0x01, 0x00), Transition(MzFeature(606.2664f0, 606.2543f0, 606.27856f0), 0x00000009, 'y', 0x06, 0x01,…
    0x0000000a │ Transition[Transition(MzFeature(322.10675f0, 322.1003f0, 322.1132f0), 0x0000000a, 'y', 0x03, 0x01, 0x00), Transition(MzFeature(762.36755f0, 762.3523f0, 762.3828f0), 0x0000000a, 'y', 0x07, 0x01, 0x00), Transition(MzFeature(606.2664f0, 606.2543f0, 606.27856f0), 0x0000000a, 'y', 0x06, 0x01,…
    0x0000000b │ Transition[Transition(MzFeature(322.10675f0, 322.1003f0, 322.1132f0), 0x0000000b, 'y', 0x03, 0x01, 0x00), Transition(MzFeature(762.36755f0, 762.3523f0, 762.3828f0), 0x0000000b, 'y', 0x07, 0x01, 0x00), Transition(MzFeature(606.2664f0, 606.2543f0, 606.27856f0), 0x0000000b, 'y', 0x06, 0x01,…
    0x0000000c │ Transition[Transition(MzFeature(322.10675f0, 322.1003f0, 322.1132f0), 0x0000000c, 'y', 0x03, 0x01, 0x00), Transition(MzFeature(762.36755f0, 762.3523f0, 762.3828f0), 0x0000000c, 'y', 0x07, 0x01, 0x00), Transition(MzFeature(606.2664f0, 606.2543f0, 606.27856f0), 0x0000000c, 'y', 0x06, 0x01,…


    julia> testPtable.lh_pair_id_to_light_heavy_pair
    6-element UnorderedDictionary{UInt32, LightHeavyPrecursorPair}
     0x00000005 │ LightHeavyPrecursorPair(0x00000009, 0x0000000a, "DRAGRACE", "D[Hglu]RAGRACE", 0x00000009, 0x0000000a, {})
     0x00000004 │ LightHeavyPrecursorPair(0x00000007, 0x00000008, "PEPTICKDEK", "PEPTICKDEK[Hlys]", 0x00000007, 0x00000008, {})
     0x00000006 │ LightHeavyPrecursorPair(0x00000009, 0x0000000a, "DRAGRACE", "D[Hglu]RAGRACE", 0x0000000b, 0x0000000c, {})
     0x00000002 │ LightHeavyPrecursorPair(0x00000003, 0x00000004, "AMINEACIDK", "AMINEACIDK[Hlys]", 0x00000003, 0x00000004, {})
     0x00000003 │ LightHeavyPrecursorPair(0x00000005, 0x00000006, "PEPTIC[Carb]KDEK", "PEPTIC[Carb]KDEK[Hlys]", 0x00000005, 0x00000006, {})
     0x00000001 │ LightHeavyPrecursorPair(0x00000001, 0x00000002, "PEPTIDER", "PEPTIDER[Harg]", 0x00000001, 0x00000002, {})

    julia> testPtable.pep_sequence_to_pep_id
    10-element UnorderedDictionary{String, UInt32}
        "AMINEACIDK[Hlys]" │ 0x00000004
            "D[Hglu]RAGRACE" │ 0x0000000a
                "PEPTIDER" │ 0x00000001
        "PEPTICKDEK[Hlys]" │ 0x00000008
                "AMINEACIDK" │ 0x00000003
                "DRAGRACE" │ 0x00000009
    "PEPTIC[Carb]KDEK[Hlys]" │ 0x00000006
        "PEPTIC[Carb]KDEK" │ 0x00000005
            "PEPTIDER[Harg]" │ 0x00000002
                "PEPTICKDEK" │ 0x00000007

    julia> testPtable.simple_precursor_set
    Set{SimplePrecursor} with 6 elements:
        SimplePrecursor("DRAGRACE", 0x03, 0x00, 0x00000009)
        SimplePrecursor("PEPTICKDEK", 0x02, 0x00, 0x00000007)
        SimplePrecursor("AMINEACIDK", 0x02, 0x00, 0x00000003)
        SimplePrecursor("PEPTIDER", 0x02, 0x00, 0x00000001)
        SimplePrecursor("PEPTIC[Carb]KDEK", 0x02, 0x00, 0x00000005)
        SimplePrecursor("DRAGRACE", 0x02, 0x00, 0x00000009)

        julia> testPtable.prec_id_to_lh_pair_id
        12-element UnorderedDictionary{UInt32, UInt32}
         0x00000005 │ 0x00000003
         0x00000004 │ 0x00000002
         0x00000006 │ 0x00000003
         0x00000007 │ 0x00000004
         0x00000002 │ 0x00000001
         0x0000000a │ 0x00000005
         0x0000000b │ 0x00000006
         0x00000009 │ 0x00000005
         0x0000000c │ 0x00000006
         0x00000008 │ 0x00000004
         0x00000003 │ 0x00000002
         0x00000001 │ 0x00000001
        
    mods_dict = Dict("Carb" => Float32(57.021464),
    "Harg" => Float32(10.008269),
    "Hlys" => Float32(8.014199),
    "Hglu" => Float32(6))
    testPtable = ISPRMPrecursorTable()
    buildPrecursorTable!(testPtable, mods_dict, "./data/peptide_lists/TRANSITION_LIST_TEST1.txt")

    """

    mods_dict = Dict("Carb" => Float32(57.021464),
                    "Harg" => Float32(10.008269),
                    "Hlys" => Float32(8.014199),
                    "Hglu" => Float32(6)
                    )

    testPtable = ISPRMPrecursorTable()
    buildPrecursorTable!(testPtable, mods_dict, "./data/peptide_lists/TRANSITION_LIST_TEST1.txt")
                
    @test length(getIDToProt(testPtable)) == 3
    @test length(getProtToID(testPtable)) == 3
    @test length(getIDToPepGroup(testPtable)) == 4
    @test length(getPepGroupToID(testPtable)) == 4
    @test length(getIDToPep(testPtable)) == 10
    @test length(getIDToPrec(testPtable)) == 12
    @test length(getPrecursorIDs(testPtable)) == 12
    @test length(getPrecIDToTransitions(testPtable)) == 12
    @test length(getIDToLightHeavyPair(testPtable)) == 6
    @test length(getPepSeqToPepID(testPtable)) == 10
    @test length(getSimplePrecursors(testPtable)) == 6
    @test length(getPrecIDToLHPairID(testPtable)) == 12

    @test getPepGroupIDs(getProt(testPtable, UInt32(1))) == Set(UInt32[1, 2])
    @test getPepGroupIDs(getProt(testPtable, UInt32(2))) == Set(UInt32[1, 3])
    @test getPepGroupIDs(getProt(testPtable, UInt32(3))) == Set(UInt32[1, 4])

    @test getProtToID(testPtable)["PROT_A"] == UInt32(1)
    @test getProtToID(testPtable)["PROT_B"] == UInt32(2)
    @test getProtToID(testPtable)["PROT_C"] == UInt32(3)

    @test getProtIDs(getPepGroup(testPtable, UInt32(1))) == Set(UInt32[1, 2, 3])
    @test getProtIDs(getPepGroup(testPtable, UInt32(2))) == Set(UInt32[1])
    @test getProtIDs(getPepGroup(testPtable, UInt32(3))) == Set(UInt32[2])
    @test getProtIDs(getPepGroup(testPtable, UInt32(4))) == Set(UInt32[3])

    @test getPepIDs(getPepGroup(testPtable, UInt32(1))) == Set(UInt32[1, 2])
    @test getPepIDs(getPepGroup(testPtable, UInt32(2))) == Set(UInt32[3, 4])
    @test getPepIDs(getPepGroup(testPtable, UInt32(3))) == Set(UInt32[5, 6, 7, 8])
    @test getPepIDs(getPepGroup(testPtable, UInt32(4))) == Set(UInt32[9, 10])

    @test getSeq(getPepGroup(testPtable, UInt32(1))) == "PEPTIDER"
    @test getSeq(getPepGroup(testPtable, UInt32(2))) == "AMINEACIDK"
    @test getSeq(getPepGroup(testPtable, UInt32(3))) == "PEPTICKDEK"
    @test getSeq(getPepGroup(testPtable, UInt32(4))) == "DRAGRACE"

    @test getPrecIDs(getPep(testPtable, UInt32(1))) == Set(UInt32[1])
    @test getPrecIDs(getPep(testPtable, UInt32(2))) == Set(UInt32[2])
    @test getPrecIDs(getPep(testPtable, UInt32(3))) == Set(UInt32[3])
    @test getPrecIDs(getPep(testPtable, UInt32(4))) == Set(UInt32[4])
    @test getPrecIDs(getPep(testPtable, UInt32(5))) == Set(UInt32[5])
    @test getPrecIDs(getPep(testPtable, UInt32(6))) == Set(UInt32[6])
    @test getPrecIDs(getPep(testPtable, UInt32(7))) == Set(UInt32[7])
    @test getPrecIDs(getPep(testPtable, UInt32(8))) == Set(UInt32[8])
    @test getPrecIDs(getPep(testPtable, UInt32(9))) == Set(UInt32[9, 11])
    @test getPrecIDs(getPep(testPtable, UInt32(10))) == Set(UInt32[10, 12])

    @test getPepGroupID(getPep(testPtable, UInt32(1))) == UInt32(1)
    @test getPepGroupID(getPep(testPtable, UInt32(2))) == UInt32(1)
    @test getPepGroupID(getPep(testPtable, UInt32(3))) == UInt32(2)
    @test getPepGroupID(getPep(testPtable, UInt32(4))) == UInt32(2)
    @test getPepGroupID(getPep(testPtable, UInt32(5))) == UInt32(3)
    @test getPepGroupID(getPep(testPtable, UInt32(6))) == UInt32(3)
    @test getPepGroupID(getPep(testPtable, UInt32(7))) == UInt32(3)
    @test getPepGroupID(getPep(testPtable, UInt32(8))) == UInt32(3)
    @test getPepGroupID(getPep(testPtable, UInt32(9))) == UInt32(4)
    @test getPepGroupID(getPep(testPtable, UInt32(10))) == UInt32(4)

    @test getSeq(getPep(testPtable, UInt32(1))) == "PEPTIDER"
    @test getSeq(getPep(testPtable, UInt32(2))) == "PEPTIDER[Harg]"
    @test getSeq(getPep(testPtable, UInt32(3))) == "AMINEACIDK"
    @test getSeq(getPep(testPtable, UInt32(4))) == "AMINEACIDK[Hlys]"
    @test getSeq(getPep(testPtable, UInt32(5))) == "PEPTIC[Carb]KDEK"
    @test getSeq(getPep(testPtable, UInt32(6))) == "PEPTIC[Carb]KDEK[Hlys]"
    @test getSeq(getPep(testPtable, UInt32(7))) == "PEPTICKDEK"
    @test getSeq(getPep(testPtable, UInt32(8))) == "PEPTICKDEK[Hlys]"
    @test getSeq(getPep(testPtable, UInt32(9))) == "DRAGRACE"
    @test getSeq(getPep(testPtable, UInt32(10))) == "D[Hglu]RAGRACE"











end