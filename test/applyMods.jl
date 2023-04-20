
using Test, Combinatorics, Dictionaries 
function Tol(a, b, ppm = 2)
    abs(a-b)<=(ppm*minimum((a, b))/1000000)
end

@testset "applyMods.jl" begin

    ##########
    #Protein
    ##########
    tprot = Protein("TEST_PROTEIN")
    for i in 1:10
        addPepGroupID!(tprot, UInt32(i))
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
end