
@testset "buildPrecursorTable.jl" begin
    """
    data/peptide_lists/PROT_PEPTIDE_TEST1.txt file looks like this. 

    PROT_A	PEPTIDE
    PROT_B	PEPTIDE	
    PROT_A	PEPTIDE
    PROT_C	PEPTIDE
    PROT_A	AMINEACID
    PROT_B	PEPTICKDEK
    PROT_C	DRAGRACER

    There are 3 proteins and 4 peptide groups. With a variable mod on 
    C-terminal K and R, there are 6 distinct peptides. Given three charge
    states for each peptide, there are 3*6 = 18 precursors. A precursor
    table `ptable` to represent this list should look like:

    ptable.id_to_prot
    3-element UnorderedDictionary{UInt32, Protein}
    0x00000002 │ Protein("PROT_B", Set(UInt32[0x00000003, 0x00000001]))
    0x00000003 │ Protein("PROT_C", Set(UInt32[0x00000004, 0x00000001]))
    0x00000001 │ Protein("PROT_A", Set(UInt32[0x00000002, 0x00000001]))

    `ptable.id_to_prot` is saying that "PROT_A" corresponds to peptide groups 1 and 2, 
    "PROT_B" corresponds to peptide groups 1 and 3, and "PROT_C" corresponds
    to peptide groups 1, and 4. This agrees with the order of the peptide list. 

    ptable.prot_to_id
    3-element UnorderedDictionary{String, UInt32}
    "PROT_C" │ 0x00000003
    "PROT_B" │ 0x00000002
    "PROT_A" │ 0x00000001
   
    ptable.id_to_pepGroup
    4-element UnorderedDictionary{UInt32, PeptideGroup}
    0x00000004 │ PeptideGroup(Set(UInt32[0x00000003]), Set(UInt32[0x00000005, 0x00000006]), "DRAGRACER")
    0x00000002 │ PeptideGroup(Set(UInt32[0x00000001]), Set(UInt32[0x00000002]), "AMINEACID")
    0x00000003 │ PeptideGroup(Set(UInt32[0x00000002]), Set(UInt32[0x00000004, 0x00000003]), "PEPTICKDEK")
    0x00000001 │ PeptideGroup(Set(UInt32[0x00000002, 0x00000003, 0x00000001]), Set(UInt32[0x00000001]), "PEPTIDE")
    
    `ptable.id_to_pepGroup` is saying that peptide group #1 is represented 
    by proteins 1, 2, and 3. That is because each of those proteins
    contains the peptide. Peptide groups 2, 3, and 4 correspond to protein
    groups 1, 2, and 3 respectively. Those peptide groups are unique to single
    proteins

    x.pepGroup_to_id
    4-element UnorderedDictionary{String, UInt32}
    "AMINEACID" │ 0x00000002
      "PEPTIDE" │ 0x00000001
   "PEPTICKDEK" │ 0x00000003
    "DRAGRACER" │ 0x00000004

    x.id_to_pep
    6-element UnorderedDictionary{UInt32, Peptide}
    0x00000005 │ Peptide("DRAGRAC[Carb]ER", 0x00000004, Set{UInt32}())
    0x00000004 │ Peptide("PEPTIC[Carb]KDEK[Hlys]", 0x00000003, Set{UInt32}())
    0x00000006 │ Peptide("DRAGRAC[Carb]ER[Harg]", 0x00000004, Set{UInt32}())
    0x00000002 │ Peptide("AMINEAC[Carb]ID", 0x00000002, Set{UInt32}())
    0x00000003 │ Peptide("PEPTIC[Carb]KDEK", 0x00000003, Set{UInt32}())
    0x00000001 │ Peptide("PEPTIDE", 0x00000001, Set{UInt32}())
    """

    test_mods::Dict{String, Float32} = 
    Dict{String, Float32}(
        "Carb" => Float32(57.021464),
        "Harg" => Float32(10),
        "Hlys" => Float32(8),
    )
    fixed_mods = [(p=r"C", r="C[Carb]")]
    var_mods = [(p=r"(K$)", r="[Hlys]"), (p=r"(R$)", r="[Harg]")]

    testPtable = PrecursorTable()
    buildPrecursorTable!(testPtable, fixed_mods, var_mods, 2, "../data/peptide_lists/PROT_PEPTIDE_TEST1.txt")
    #=
    testPtable = PrecursorTable()
    buildPrecursorTable!(testPtable, fixed_mods, var_mods, 2, "./data/peptide_lists/PROT_PEPTIDE_TEST1.txt")
    testPtable.id_to_pepGroup
    =#
    addPrecursors!(testPtable, UInt8[2, 3, 4], UInt8[0], test_mods)
    #Test sizes
    @test length(getIDToProt(testPtable)) == 3
    @test length(getProtToID(testPtable)) == 3
    @test length(getIDToPepGroup(testPtable)) == 4
    @test length(getPepGroupToID(testPtable)) == 4
    @test length(getIDToPep(testPtable)) == 6
    @test length(getPrecursorIDs(testPtable)) == 18

    @test getPepGroupIDs(getIDToProt(testPtable)[UInt32(1)]) == Set(UInt32[0x00000002, 0x00000001])
    @test getPepGroupIDs(getIDToProt(testPtable)[UInt32(2)]) == Set(UInt32[0x00000003, 0x00000001])
    @test getPepGroupIDs(getIDToProt(testPtable)[UInt32(3)]) == Set(UInt32[0x00000004, 0x00000001])

    @test getProtToID(testPtable)["PROT_A"] == UInt32(1)
    @test getProtToID(testPtable)["PROT_B"] == UInt32(2)
    @test getProtToID(testPtable)["PROT_C"] == UInt32(3)

    @test getProtIDs(getIDToPepGroup(testPtable)[UInt32(1)]) == Set(UInt32[0x00000002, 0x00000003, 0x00000001])
    @test getProtIDs(getIDToPepGroup(testPtable)[UInt32(2)]) == Set(UInt32[0x00000001])
    @test getProtIDs(getIDToPepGroup(testPtable)[UInt32(3)]) == Set(UInt32[0x00000002])
    @test getProtIDs(getIDToPepGroup(testPtable)[UInt32(4)]) == Set(UInt32[0x00000003])

    @test getPepIDs(getIDToPepGroup(testPtable)[UInt32(1)]) == Set(UInt32[0x00000001])
    @test getPepIDs(getIDToPepGroup(testPtable)[UInt32(2)]) == Set(UInt32[0x00000002])
    @test getPepIDs(getIDToPepGroup(testPtable)[UInt32(3)]) == Set(UInt32[0x00000004, 0x00000003])
    @test getPepIDs(getIDToPepGroup(testPtable)[UInt32(4)]) == Set(UInt32[0x00000005, 0x00000006])

    @test getPepGroupToID(testPtable)["PEPTIDE"] == UInt32(1)
    @test getPepGroupToID(testPtable)["AMINEACID"] == UInt32(2)
    @test getPepGroupToID(testPtable)["PEPTICKDEK"] == UInt32(3)
    @test getPepGroupToID(testPtable)["DRAGRACER"] == UInt32(4)

    @test getSeq(getIDToPep(testPtable)[UInt32(1)]) == "PEPTIDE"
    @test getSeq(getIDToPep(testPtable)[UInt32(2)]) == "AMINEAC[Carb]ID"
    @test getSeq(getIDToPep(testPtable)[UInt32(3)]) == "PEPTIC[Carb]KDEK" 
    @test getSeq(getIDToPep(testPtable)[UInt32(4)]) == "PEPTIC[Carb]KDEK[Hlys]" 
    @test getSeq(getIDToPep(testPtable)[UInt32(5)]) == "DRAGRAC[Carb]ER"
    @test getSeq(getIDToPep(testPtable)[UInt32(6)]) == "DRAGRAC[Carb]ER[Harg]" 
    @test getPepGroupID(getIDToPep(testPtable)[UInt32(1)]) == 0x00000001
    @test getPepGroupID(getIDToPep(testPtable)[UInt32(2)]) == 0x00000002
    @test getPepGroupID(getIDToPep(testPtable)[UInt32(3)]) == 0x00000003
    @test getPepGroupID(getIDToPep(testPtable)[UInt32(4)]) == 0x00000003
    @test getPepGroupID(getIDToPep(testPtable)[UInt32(5)]) == 0x00000004 
    @test getPepGroupID(getIDToPep(testPtable)[UInt32(6)]) == 0x00000004 
    @test issorted([getPrecursor(testPtable, prec_id) for prec_id in getPrecursorIDs(testPtable)], by=x->getMZ(x))

    """
        Apply to real world data. We have a list of 260 peptides from 90 proteins 
    from the NRF2 SIL peptide array. 
    """
    testPtable = PrecursorTable()
    fixed_mods = [(p=r"C", r="C[Carb]")]
    var_mods = [(p=r"(K$)", r="[Hlys]"), (p=r"(R$)", r="[Harg]")]
    buildPrecursorTable!(testPtable, fixed_mods, var_mods, 2, "../data/NRF2_SIL.txt")
    addPrecursors!(testPtable, UInt8[2, 3, 4], UInt8[0], test_mods)
    @test length(getIDToPepGroup(testPtable)) == 260
    @test length(getPrecursorIDs(testPtable)) == 260*2*3 #2 because each protein ends in a variably modifiable K or R. 3 because 3 charge states. 

    @test getProtNamesFromPepSeq(testPtable, "LAIEAGFR") == Set(["AKR1C1","AKR1C2","AKR1C3"]) 
    @test Set(getPepSeqsFromProt(testPtable, "CD8A")) == Set(["AAEGLDTQR", "TWNLGETVELK", "LGDTFVLTLSDFR"])
    @test Set(getPepSeqsFromProt(testPtable, "ZNF746")) == Set(["LLSLEGR", "GQPLPTPPAPPDPFK"])

end
