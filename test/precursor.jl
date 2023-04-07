#using Titus
using Test
const default_mods::Dict{String, Float32} = 
Dict{String, Float32}(
    "Carb" => Float32(57.021464)
)
function Tol(a, b, ppm = 2)
    abs(a-b)<=(ppm*minimum((a, b))/1000000)
end
@testset "precursor.jl" begin

    
    @test 1==1# Write your tests here.

    #########
    #Tests for 'AA' Struct 
    #########
    @test getAA(AA('A')) == 'A'
    @test Tol(getMass(AA('A')), 71.03711)
    #"Selenocysteine"
    @test Tol(getMass(AA('U')), 150.95363)

    @test_throws ErrorException("The character Z cannot be interpreted as an amino acid!") AA('Z')
    @test_throws ErrorException("The character % cannot be interpreted as an amino acid!") AA('%')

    #########
    #Tests for 'Mod' Struct
    #########    
    @test Tol(getMass(Mod("K[+8.014199]")), 8.014199)
    @test Tol(getMass(Mod("K")), 0.0)
    @test Tol(getMass(Mod('K')), 0.0)
    @test Tol(getMass(Mod("C[Carb]", default_mods)), 57.021464)
    #########
    #Tests for 'Residue' Struct
    #########
    #Should throw an error

    @test Tol(getMass(Residue("K[+8.014199]")), 8.014199 + getMass(AA('K')))
    #@test getMod(Residue("K[+8.014199]")) == Mod("K[+8.014199]", 8.014199)

    @test Tol(getMass(Residue("K", Float32(8.014199))), 8.014199 + getMass(AA('K')))
    #@test getMod(Residue("K", Float32(8.014199))) == Mod("K", 8.014199)

    #"C[Carb]" should be a build in modification, which users could add to 
    @test Tol(getMass(Residue("C[Carb]", default_mods)), 57.021464 + getMass(AA('C')))

    @test_throws ErrorException("C[asdf] could not be parsed as given") getMod(Residue("C[asdf]")) 
    
    # #########
    # #Tests for getIonMZ
    # #########
    TIDE = Array{Residue, 1}([Residue('T') , Residue('I'), Residue('D'), Residue('E')])
    @test Tol(getIonMZ(TIDE, 'y', UInt8(1)), 477.219119)
    @test Tol(getIonMZ(TIDE, 'y', UInt8(2)), 239.113198)

    PEP = Array{Residue, 1}([Residue('P') , Residue('E'), Residue('P')])
    @test Tol(getIonMZ(PEP, 'b', UInt8(1)), 324.155397)
    @test Tol(getIonMZ(PEP, 'b', UInt8(2)), 162.581336)

    TIDEK_mod = Array{Residue, 1}([Residue('T') , Residue('I'), Residue('D'), Residue('E'), Residue("K[+8.014199]")])
    @test Tol(getIonMZ(TIDEK_mod, 'y', UInt8(1)), 613.328281)
    @test Tol(getIonMZ(TIDEK_mod, 'y', UInt8(2)), 307.167779)

    PEPTIDE = Array{Residue, 1}([Residue('P') , Residue('E'), Residue('P'),
                                 Residue('T') , Residue('I'), Residue('D'), 
                                 Residue('E')])
    
    # #Can Specify [M+0], [M+1], or [M+2]
    @test Tol(getIonMZ(PEPTIDE,'p',UInt8(2),isotope = UInt8(0)), 400.687258)
    @test Tol(getIonMZ(PEPTIDE,'p',UInt8(2),isotope = UInt8(1)), 401.188771)
    @test Tol(getIonMZ(PEPTIDE,'p',UInt8(2),isotope = UInt8(2)), 401.690036)

    @test Tol(getIonMZ(PEPTIDE,'p',UInt8(3),isotope = UInt8(0)), 267.460597)
    @test Tol(getIonMZ(PEPTIDE,'p',UInt8(3),isotope = UInt8(1)), 267.79494)
    @test Tol(getIonMZ(PEPTIDE,'p',UInt8(3),isotope = UInt8(2)), 268.129116)

    PEPTIDEK_mod = Array{Residue, 1}([Residue('P') , Residue('E'), Residue('P'),
                                      Residue('T') , Residue('I'), Residue('D'), 
                                      Residue('E'), Residue("K[+8.014199]")])

    @test Tol(getIonMZ(PEPTIDEK_mod,'p',UInt8(2),isotope = UInt8(0)), 468.741839)
    @test Tol(getIonMZ(PEPTIDEK_mod,'p',UInt8(2),isotope = UInt8(1)), 469.243364)
    @test Tol(getIonMZ(PEPTIDEK_mod,'p',UInt8(2),isotope = UInt8(2)), 469.74462)                                      
    #########
    #Tests for 'MzFeature' Struct 
    #########

    #########
    #Tests for 'Precursor' Struct 
    #########
    PEPTIDE_1 = Precursor("PEPTIDE", charge = UInt8(2))
    PEPTIDE_2 = Precursor("PEPTIDE", charge = UInt8(3))
    PEPTIDE_3 = Precursor("C[+57.021464]PEC[+57.021464]PTIDE", charge = UInt8(3))
    @test Tol(getMZ(PEPTIDE_1), 400.687258)
    @test Tol(getMZ(PEPTIDE_2), 267.460597)
    @test Tol(getMZ(PEPTIDE_3), 374.147696)
    @test Tol(getHigh(PEPTIDE_1), 400.687258*(1 + 20/1000000))
    @test Tol(getHigh(PEPTIDE_2), 267.460597*(1 + 20/1000000))
    @test Tol(getHigh(PEPTIDE_3), 374.147696*(1 + 20/1000000))
    @test Tol(getLow(PEPTIDE_1), 400.687258*(1 - 20/1000000))
    @test Tol(getLow(PEPTIDE_2), 267.460597*(1 - 20/1000000))
    @test Tol(getLow(PEPTIDE_3), 374.147696*(1 - 20/1000000))
    @test getResidues(PEPTIDE_2) == getResidues("PEPTIDE")
    @test getResidues(PEPTIDE_3) == getResidues("C[+57.021464]PEC[+57.021464]PTIDE")
    @test getCharge(PEPTIDE_1) == UInt8(2)
    @test getCharge(PEPTIDE_3) == UInt8(3)
    #########
    #Tests for 'Transition' Struct 
    #########
    #PEPTIDE_1_res = getResidues(PEPTIDE_1)
    PEPTIDE_3_res = getResidues(PEPTIDE_3)
    #b2+1
    @test Tol(getMZ(Transition(PEPTIDE_3_res, ion_type = 'b', ind = UInt8(2))), 258.090688)
    #b4+2
    @test Tol(getMZ(Transition(PEPTIDE_3_res, ion_type = 'b', ind = UInt8(4), charge = UInt8(2))), 274.085603)
    #y2+1
    @test Tol(getMZ(Transition(PEPTIDE_3_res, ion_type = 'y', ind = UInt8(2))), 263.087377)
    #y5+2
    @test Tol(getMZ(Transition(PEPTIDE_3_res, ion_type = 'y', ind = UInt8(5), charge = UInt8(2))), 287.63958)
    #########
    #Tests for 'getFragIons' function and methods
    #########
    PEPTIDE_frags_charge1 = sort(Float32[
    227.102633, # b2+1
    324.155397, # b3+1
    425.203075, # b4+1
    538.287139, # b5+1
    653.314082, # b6+1]
    263.087377, # y2+1
    376.171441, # y3+1
    477.219119, # y4+1
    574.271883, # y5+1
    703.314476, # y6+1
    ])
    compare_frags = zip(sort(getFragIons(PEPTIDE_1, charge = UInt8(1), 
                                        b_start = 2, y_start = 2
                                        )), 
                        PEPTIDE_frags_charge1
                        )
    @test all([Tol(pair[1], pair[2]) for pair in compare_frags])
    PEPTIDE_frags_charge2 = sort(Float32[
        162.581336, # b3+1
        213.105176, # b4+1
        269.647208, # b5+1
        327.160679, # b6+1
        188.589358, # y3+1
        239.113198, # y4+1
        287.63958, # y5+1
        352.160876, # y6+1
        ])
        compare_frags = zip(sort(getFragIons(PEPTIDE_1, charge = UInt8(2), 
                            b_start = 2, y_start = 2
                            )), 
                            PEPTIDE_frags_charge2
                            )
      @test all([Tol(pair[1], pair[2]) for pair in compare_frags])

      
      compare_frags = zip(sort([getMZ(x) for x in getTransitions(PEPTIDE_1, charge = UInt8(2), y_start = 2, b_start = 2)]), 
      PEPTIDE_frags_charge2
      )

      @test all([Tol(pair[1], pair[2]) for pair in compare_frags])

end
