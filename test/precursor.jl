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
    @test length(test_frags[findall(frag -> getType(frag) == 'b' && getCharge(frag) == 1, test_frags)]) = 3
    @test getMass(test_frags[findall(frag -> getType(frag) == 'b' && getCharge(frag) == 1, test_frags)]) = 3


    frag!(test_prec)
    @test length(getFrags(test_prec))=21
    
    #Compare to Skyline

end
