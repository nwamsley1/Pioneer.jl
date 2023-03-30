
using Test, Combinatorics, Dictionaries 
function Tol(a, b, ppm = 2)
    abs(a-b)<=(ppm*minimum((a, b))/1000000)
end

@testset "getPrecursors.jl" begin

    ##########
    #Protein
    ##########
    tprot = Protein("TEST_PROTEIN")
    for i in 1:10
        addPepGroup!(tprot, UInt32(i))
    end
    @test tprot.pep_group_ids == Set([UInt32(i) for i in 1:10])

    ##########
    #fixedMods
    ##########
    fixed_mods = [(p=r"C", r="C[Carb]")]
    @test fixedMods("CARACHTER", fixed_mods) == "C[Carb]ARAC[Carb]HTER"
    fixed_mods = [(p=r"C", r="C[Carb]"), (p=r"A", r="A[mod A]")]
    @test fixedMods("RATRACE", fixed_mods) == "RA[mod A]TRA[mod A]C[Carb]E"

    ##########
    #matchVarMods/appyMods!
    ##########
    test_peps_dict = UnorderedDictionary{UInt32, Peptide}()
    group_id = UInt32(1)
    pep_id = UInt32(1)
    var_mods = [(p=r"(E)", r="[modE]"), (p=r"(C)", r="[modC]")]
    test_pep = "PEPEPCPEPEP"
    @test Set(matchVarMods(var_mods, test_pep)) == Set([(2:2, "[modE]"),
                                                    (4:4, "[modE]"),
                                                    (6:6, "[modC]"),
                                                    (8:8, "[modE]"),
                                                    (10:10, "[modE]")])
    #Two or fewer variable modifications
    applyMods!(test_peps_dict, var_mods, test_pep, group_id, pep_id, n=2);
    @test Set(map(pep->getSeq(pep), test_peps_dict)) == Set([
                                                            "PEPEPCPEPEP",
                                                            "PEPEPC[modC]PEPEP",
                                                            "PE[modE]PEPCPEPEP",
                                                            "PEPE[modE]PCPEPEP",
                                                            "PEPEPCPE[modE]PEP",
                                                            "PEPEPCPEPE[modE]P",
                                                            "PE[modE]PE[modE]PCPEPEP",
                                                            "PE[modE]PEPCPE[modE]PEP",
                                                            "PE[modE]PEPCPEPE[modE]P",
                                                            "PEPE[modE]PCPE[modE]PEP",
                                                            "PEPE[modE]PCPEPE[modE]P",
                                                            "PEPEPCPE[modE]PE[modE]P",
                                                            "PE[modE]PEPC[modC]PEPEP",
                                                            "PEPE[modE]PC[modC]PEPEP",
                                                            "PEPEPC[modC]PE[modE]PEP",
                                                            "PEPEPC[modC]PEPE[modE]P"
                                                            ])
    #One or fewer variable modifications
    test_peps_dict = UnorderedDictionary{UInt32, Peptide}()
    group_id = UInt32(1)
    pep_id = UInt32(1)
    applyMods!(test_peps_dict, var_mods, test_pep, group_id, pep_id, n=1);
    @test Set(map(pep->getSeq(pep), test_peps_dict)) == Set([
                                                            "PEPEPCPEPEP",
                                                            "PEPEPC[modC]PEPEP",
                                                            "PE[modE]PEPCPEPEP",
                                                            "PEPE[modE]PCPEPEP",
                                                            "PEPEPCPE[modE]PEP",
                                                            "PEPEPCPEPE[modE]P"
                                                            ])
    #Try some more complicated regex patterns 
    var_mods = [
                (p=r"(K$)", r="[C-term K]"), #C-terminal K
                (p=r"(?<=[ED])P(?=[RK])", r="[special P]"), #P preceeded by and E or a D and followed by an R or K. Use positive-lookahead and lookbehind
               ]
    test_peps_dict = UnorderedDictionary{UInt32, Peptide}()
    group_id = UInt32(1)
    pep_id = UInt32(1)
    test_pep = "PEPRKPK"
    #Matches the proline at position 3 and not the one at position 1
    #Matches the lysine at position 7 and not the one at position 5
    @test Set(matchVarMods(var_mods, test_pep)) == Set([(3:3, "[special P]"),
                                                    (7:7, "[C-term K]")])

    applyMods!(test_peps_dict, var_mods, test_pep, group_id, pep_id, n=2);

    @test Set(map(pep->getSeq(pep), test_peps_dict)) == Set([
                                                            "PEPRKPK",
                                                            "PEP[special P]RKPK",
                                                            "PEPRKPK[C-term K]",
                                                            "PEP[special P]RKPK[C-term K]",
                                                            ])

    #What happens when we try and apply multiple variable modifications to a single residue?
    #Right now modifications are appended in sequence like so K[mod 1][mod 2]. 
    #This behaviour might need to change in the future. 
    var_mods = [
                (p=r"(K$)", r="[C-term K]"), #C-terminal K
                (p=r"K", r="[methyl K]"), #P preceeded by and E or a D and followed by an R or K. Use positive-lookahead and lookbehind
               ]
    test_peps_dict = UnorderedDictionary{UInt32, Peptide}()
    group_id = UInt32(1)
    pep_id = UInt32(1)
    test_pep = "KARENK"
    #Matches the proline at position 3 and not the one at position 1
    #Matches the lysine at position 7 and not the one at position 5
    @test Set(matchVarMods(var_mods, test_pep)) == Set([(1:1, "[methyl K]"),
                                                        (6:6, "[methyl K]"),
                                                        (6:6, "[C-term K]")])

    applyMods!(test_peps_dict, var_mods, test_pep, group_id, pep_id, n=2);
    
    @test Set(map(pep->getSeq(pep), test_peps_dict)) == Set([
                                                            "KARENK",
                                                            "K[methyl K]ARENK",
                                                            "KARENK[methyl K]",
                                                            "KARENK[C-term K]",
                                                            "K[methyl K]ARENK[methyl K]",
                                                            "K[methyl K]ARENK[C-term K]",
                                                            "KARENK[C-term K][methyl K]",
                                                            ])
    ##########
    #buildPrecursorTable!
    ##########
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
    0x00000004 │ PeptideGroup(Set(UInt32[0x00000003]))
    0x00000002 │ PeptideGroup(Set(UInt32[0x00000001]))
    0x00000003 │ PeptideGroup(Set(UInt32[0x00000002]))
    0x00000001 │ PeptideGroup(Set(UInt32[0x00000002, 0x00000003, 0x00000001]))
    
    `ptable.id_to_pepGroup` is saying that peptide group #1 is represented 
    by proteins 1, 2, and 3. That is because each of those proteins
    contains the peptide. Peptide groups 2, 3, and 4 correspond to protein
    groups 1, 2, and 3 respectively. Those peptide groups are unique to single
    proteins

    x.pepGroup_to_id
    4-element UnorderedDictionary{String, UInt32}
     "DRAGRAC[Carb]ER" │ 0x00000004
     "AMINEAC[Carb]ID" │ 0x00000002
             "PEPTIDE" │ 0x00000001
    "PEPTIC[Carb]KDEK" │ 0x00000003

    x.id_to_pep
    6-element UnorderedDictionary{UInt32, Peptide}
    0x00000005 │ Peptide("DRAGRAC[Carb]ER", 0x00000004)
    0x00000004 │ Peptide("PEPTIC[Carb]KDEK[Hlys]", 0x00000003)
    0x00000006 │ Peptide("DRAGRAC[Carb]ER[Harg]", 0x00000004)
    0x00000002 │ Peptide("AMINEAC[Carb]ID", 0x00000002)
    0x00000003 │ Peptide("PEPTIC[Carb]KDEK", 0x00000003)
    0x00000001 │ Peptide("PEPTIDE", 0x00000001)
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
    getPrecursors!(testPtable, UInt8[2, 3, 4], UInt8[0], test_mods)

    #Test sizes
    @test length(getIDToProt(testPtable)) == 3
    @test length(getProtToID(testPtable)) == 3
    @test length(getIDToPepGroup(testPtable)) == 4
    @test length(getPepGroupToID(testPtable)) == 4
    @test length(getIDToPep(testPtable)) == 6
    @test length(getPrecursors(testPtable)) == 18

    @test getPepGroupIDs(getIDToProt(testPtable)[UInt32(1)]) == Set(UInt32[0x00000002, 0x00000001])
    @test getPepGroupIDs(getIDToProt(testPtable)[UInt32(2)]) == Set(UInt32[0x00000003, 0x00000001])
    @test getPepGroupIDs(getIDToProt(testPtable)[UInt32(3)]) == Set(UInt32[0x00000004, 0x00000001])

    @test getProtToID(testPtable)["PROT_A"] == UInt32(1)
    @test getProtToID(testPtable)["PROT_B"] == UInt32(2)
    @test getProtToID(testPtable)["PROT_C"] == UInt32(3)

    @test getProtIDs(getIDToPepGroup(testPtable)[UInt32(1)]) == Set(UInt32[0x00000003, 0x00000002, 0x00000001])
    @test getProtIDs(getIDToPepGroup(testPtable)[UInt32(1)]) == Set(UInt32[0x00000003, 0x00000002, 0x00000001])
    @test getProtIDs(getIDToPepGroup(testPtable)[UInt32(1)]) == Set(UInt32[0x00000003, 0x00000002, 0x00000001])
    @test getProtIDs(getIDToPepGroup(testPtable)[UInt32(1)]) == Set(UInt32[0x00000003, 0x00000002, 0x00000001])

    @test getPepGroupToID(testPtable)["PEPTIDE"] == UInt32(1)
    @test getPepGroupToID(testPtable)["AMINEAC[Carb]ID"] == UInt32(2)
    @test getPepGroupToID(testPtable)["PEPTIC[Carb]KDEK"] == UInt32(3)
    @test getPepGroupToID(testPtable)["DRAGRAC[Carb]ER"] == UInt32(4)

    @test getIDToPep(testPtable)[UInt32(1)] == Peptide("PEPTIDE", 0x00000001)
    @test getIDToPep(testPtable)[UInt32(2)] == Peptide("AMINEAC[Carb]ID", 0x00000002)
    @test getIDToPep(testPtable)[UInt32(3)] == Peptide("PEPTIC[Carb]KDEK", 0x00000003)
    @test getIDToPep(testPtable)[UInt32(4)] == Peptide("PEPTIC[Carb]KDEK[Hlys]", 0x00000003)
    @test getIDToPep(testPtable)[UInt32(5)] == Peptide("DRAGRAC[Carb]ER", 0x00000004)
    @test getIDToPep(testPtable)[UInt32(6)] == Peptide("DRAGRAC[Carb]ER[Harg]", 0x00000004)

    @test issorted(getPrecursors(testPtable), by=x->getMZ(x))

    """
        Apply to real world data. We have a list of 260 peptides from 90 proteins 
    from the NRF2 SIL peptide array. 
    """
    testPtable = PrecursorTable()
    fixed_mods = [(p=r"C", r="C[Carb]")]
    var_mods = [(p=r"(K$)", r="[Hlys]"), (p=r"(R$)", r="[Harg]")]
    buildPrecursorTable!(testPtable, fixed_mods, var_mods, 2, "../data/NRF2_SIL.txt")
    getPrecursors!(testPtable, UInt8[2, 3, 4], UInt8[0], test_mods)
    @test length(getIDToPepGroup(testPtable)) == 260
    @test length(getPrecursors(testPtable)) == 260*2*3 #2 because each protein ends in a variably modifiable K or R. 3 because 3 charge states. 

    @test getProtNamesFromPepSeq(testPtable, "LAIEAGFR") == Set(["AKR1C1","AKR1C2","AKR1C3"]) 
    @test Set(getPepSeqsFromProt(testPtable, "CD8A")) == Set(["AAEGLDTQR", "TWNLGETVELK", "LGDTFVLTLSDFR"])
    @test Set(getPepSeqsFromProt(testPtable, "ZNF746")) == Set(["LLSLEGR", "GQPLPTPPAPPDPFK"])

end
