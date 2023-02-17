using Titus
using Test

@testset "Titus.jl" begin

    
    @test 1==1# Write your tests here.

    #########
    #Tests for 'AA' Struct 
    #########
    @test getAA(AA('A')) == 'A'
    @test getMass(AA('A')) == 71.03711
    "Selenocysteine"
    @test getMass(AA('U')) == 150.95363

    @test_throws ErrorException("The character Z cannot be interpreted as an amino acid!") getAA('Z')
    @test_throws ErrorException("The character % cannot be interpreted as an amino acid!") getAA('%')

    #########
    #Tests for 'Precursor' Struct 
    #########
    @test getAA(Residue("K")) == AA('K')
    @test getAA(Residue('K')) == AA('K')
    @test getMod(Residue("K")) == ("K", nothing)
    @test getMod(Residue('K')) == ("K", nothing)


    @test getAA(Residue("K[+8.014199]")) == AA('K')
    @test getMass(Residue("K[+8.014199]")) == 8.014199 + getMass(AA('K'))
    @test getMod(Residue("K[+8.014199]")) == ("K", 8.014199)

    @test getAA(Residue("K", 8.014199)) == AA('K')
    @test getMass(Residue("K", 8.014199)) == 8.014199 + getMass(AA('K'))
    @test getMod(Residue("K", 8.014199)) == ("K", 8.014199)

    #"C[Carb]" should be a build in modification, which users could add to 
    @test getAA(Residue("C[Carb]")) == AA('C')
    @test getMass(Residue("C[Carb]")) == 57.021464 + getMass(AA('C'))
    @test getMod(Residue("C[Carb]")) == ("C[Carb]", 57.021464)

    @test_throws ErrorException("C[asdf] is not recognized as a valid modification!") getMod(Residue("C[asdf]")) 
    
    #########
    #Tests for 'Frag' Struct 
    #########
    @test_throws ErrorException("Did not specify a fragment charge!", Frag("PEP"))
    @test length(Frag("PEP", 'b', 2)) == 3
    @test getCharge(Frag("PEP", 'b', 2)) == 2
    @test getType(Frag("PEP", 'b', 2)) == 'b'
    @test getMZ(Frag("PEP", 'b', 2)) == 162.581336

    @test_throws ErrorException("Did not specify a fragment charge!", Frag("PEP"))
    @test length(Frag("PEP", 'b', 1)) == 3
    @test getCharge(Frag("PEP", 'b', 1)) == 1
    @test getType(Frag("PEP", 'b', 1)) == 'b'
    @test getMZ(Frag("PEP", 'b', 1)) == 324.155397

    @test length(Frag([Residue('P') , Residue('E'), Residue('P')], 'b', 2)) == length(Frag("PEP", 'b', 2))
    @test getCharge(Frag([Residue('P') , Residue('E'), Residue('P')], 'b', 2)) == getMass(Charge("PEP", 'b', 2))
    @test getType(Frag([Residue('P') , Residue('E'), Residue('P')], 'b', 2)) == getType(Frag("PEP", 'b', 2))
    @test getMZ(Frag([Residue('P') , Residue('E'), Residue('P')], 'b', 2)) == getMass(Frag("PEP", 'b', 2))

    @test length(Frag("TIDE",'y', 1)) == 4
    @test getMZ(Frag("TIDE",'y', 1)) == 477.219119
    @test getMZ(Frag("TIDE", 'y', 2)) == 239.113198
    @test getMZ(Frag([Residue('T') , Residue('I'), Residue('D'), Residue('E')], 'y', 1)) == getMZ(Frag("TIDE",'y', 1))

    @test getMZ(Frag("TIDEK[+8.014199]",'y', 1)) == 613.328281
    @test getMZ(Frag("TIDEK[+8.014199]", 'y', 2)) == 307.167779
    @test getMZ(Frag([Residue('T') , Residue('I'), Residue('D'), Residue('E'), Residue("K[+8.014199]")], 'y', 1)) == getMZ(Frag("TIDEK[+8.014199]", 'y', 2))

    #Should give [M+1] by default
    @test getMZ(Frag("PEPTIDE", 'p', 2)) ==  401.188771

    #Can Specify [M+0], [M+1], or [M+2]
    @test getMZ(Frag("PEPTIDE", 'p', 2, isotope = 2)) == 401.690036
    @test getMZ(Frag("PEPTIDE", 'p', 2, isotope = 1)) == 401.188771
    @test getMZ(Frag("PEPTIDE", 'p', 2, isotope = 0)) == 400.687258

    @test getMZ(Frag("PEPTIDE", 'p', 3)) == 267.79494
    @test getMZ(Frag("PEPTIDE", 'p', 3, isotope = 2)) == 268.129116
    @test getMZ(Frag("PEPTIDE", 'p', 3, isotope = 1)) == 267.79494
    @test getMZ(Frag("PEPTIDE", 'p', 3, isotope = 0)) == 267.460597

    @test getMZ(Frag("PEPTIDEK[+8.014199]", 'p', 2)) == 469.243364
    @test getMZ(Frag("PEPTIDEK[+8.014199]", 'p', 2, isotope = 2)) == 469.74462
    @test getMZ(Frag("PEPTIDEK[+8.014199]", 'p', 2, isotope = 1)) == 469.243364
    @test getMZ(Frag("PEPTIDEK[+8.014199]", 'p', 2, isotope = 0)) == 468.741839

    @test getMZ(Frag([Residue('P'), Residue('E'), Residue('P'), Residue('T') , Residue('I'), Residue('D'), Residue('E'), Residue("K[+8.014199]")], 'p', 2)) == getMZ(Frag("PEPTIDEK[+8.014199]", 'p', 2, isotope = 0))

    #########
    #Tests for 'FragFilter' Struct 
    #########

    #########
    #Tests for 'Precursor' Struct 
    #########
    @test getMZ(Precursor("PEPTIDE", 2))
    @test getMZ(Precursor("PEPTIDE", 3))
    @test getMZ(Precursor("PEPTIDE", 4))
    @test length(Precursor("PEPTIDE", 4))


    @test getSequence(Precursor("PEPTIDE", 4))
    @test getResidues(map(v -> getAA(Residue(v)), "PEPTIDE")) = [AA('P'), AA('E'), AA('P'), AA('T'), AA('I'), AA('D'), AA('E')]

    test_prec = Precursor("PEPTIDE", 3)

    frag!(test_prec)
    @test length(getFrags(test_prec))=21

    #Get the length b3 ion
    @test length(getFrags(test_prec)) = 3
    test_frags = getFrags(test_prec)
    @test length(test_frags[findall(frag -> getType(frag) == 'b' && getCharge(frag) == 1 && length(frag) == 3, test_frags)]) = 3
    @test getMZ(test_frags[findall(frag -> getType(frag) == 'b' && getCharge(frag) == 1 && length(frag) == 3, test_frags)]) = 324.155397

    #Sort the fragments from high mass to low mass
    #Use constructor to get this effect
    #sort!(test_frags, by = frag -> getMZ(frag), rev = true)
    sort!(test_frags, rev = true)
    #y6+1
    @test getMZ(first(getFrags(test_prec))) == 703.314476
    @test getType(first(getFrags(test_prec))) == 'y'
    @test length(first(getFrags(test_prec))) == 6
    @test getCharge(first(getFrags(test_prec))) == 1
    #b3+2
    @test getMZ(last(getFrags(test_prec))) == 162.581336
    @test getMZ(last(getFrags(test_prec))) == 'b'
    @test length(last(getFrags(test_prec))) == 3
    @test getCharge(last(getFrags(test_prec))) == 2
    #Sort low to high
    sort!(test_frags, rev = false)
    #b2+2
    @test getMZ(first(getFrags(test_prec))) == 162.581336
    #y6+1
    @test getMZ(last(getFrags(test_prec))) == 703.314476

    #Get b3-b4 ions and all y ions 2 or greater
    #FragFilter({'b'=> ("low" => 3, "")}
    #    'b' => ("low" => 3, "high" => 4, "max_charge" => 2), 
    #    'y' => ("low" => 3, "high" => 4, "max_charge" => 2),
    #    )
    test_frag_filter = FragFilter(type = 'b', low = 3, high = 4, max_charge = 1) + FragFilter(type = 'y', low = 3, high = ∞, max_charge = 2) 
    frag!(test_prec, test_frag_filter)
    sort!(test_frags)
    @test length(getFrags(test_prec)) = 2 + 9
    @test map(frag -> getMZ(frag), getFrags(test_prec)) == [
                                                            188.589358, #y3+2
                                                            239.113198, #y4+2
                                                            287.639580, #y5+2
                                                            324.155397, #b3+1
                                                            352.160876, #y6+2
                                                            376.171441, #y3+1
                                                            425.203075, #b4+1
                                                            477.219119, #y4+1
                                                            574.271883, #y5+1
                                                            703.314476  #y6+1
                                                            ]

    
    #Compare to Skyline

end
