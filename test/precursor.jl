using Titus
using Test

function Tol(a, b, ppm = 1)
    abs(a-b)<(ppm*minimum((a, b))/10000)
end
@testset "Titus.jl" begin

    
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
    #Tests for 'Residue' Struct 
    #########
    #Should throw an error
    @test getAA(Residue("K")) == AA('K')
    @test getAA(Residue('K')) == AA('K')
    @test getMod(Residue("K")) == Mod()
    @test getMod(Residue('K')) == Mod()


    @test getAA(Residue("K[+8.014199]")) == AA('K')
    @test Tol(getMass(Residue("K[+8.014199]")), 8.014199 + getMass(AA('K')))
    #@test getMod(Residue("K[+8.014199]")) == Mod("K[+8.014199]", 8.014199)

    @test getAA(Residue("K", Float32(8.014199))) == AA('K')
    @test Tol(getMass(Residue("K", Float32(8.014199))), 8.014199 + getMass(AA('K')))
    #@test getMod(Residue("K", Float32(8.014199))) == Mod("K", 8.014199)

    #"C[Carb]" should be a build in modification, which users could add to 
    @test getAA(Residue("C[Carb]", default_mods)) == AA('C')
    @test Tol(getMass(Residue("C[Carb]", default_mods)), 57.021464 + getMass(AA('C')))
    @test getMod(Residue("C[Carb]", default_mods)) == Mod("C[Carb]", default_mods)

    @test_throws ErrorException("C[asdf] could not be parsed as given") getMod(Residue("C[asdf]")) 
    
    # #########
    # #Tests for 'Frag' Struct 
    # #########
    TIDE = Array{Residue, 1}([Residue('T') , Residue('I'), Residue('D'), Residue('E')])
    @test Tol(getMZ(Frag(TIDE, 'y', Int32(1))), 477.219119)
    @test Tol(getMZ(Frag(TIDE, 'y', Int32(2))), 239.113198)

    PEP = Array{Residue, 1}([Residue('P') , Residue('E'), Residue('P')])
    @test Tol(getMZ(Frag(PEP, 'b', Int32(1))), 324.155397)
    @test Tol(getMZ(Frag(PEP, 'b', Int32(2))), 162.581336)

    TIDEK_mod = Array{Residue, 1}([Residue('T') , Residue('I'), Residue('D'), Residue('E'), Residue("K[+8.014199]")])
    @test Tol(getMZ(Frag(TIDEK_mod, 'y', Int32(1))), 613.328281)
    @test Tol(getMZ(Frag(TIDEK_mod, 'y', Int32(2))), 307.167779)

    PEPTIDE = Array{Residue, 1}([Residue('P') , Residue('E'), Residue('P'),
                                 Residue('T') , Residue('I'), Residue('D'), 
                                 Residue('E')])
    
    # #Can Specify [M+0], [M+1], or [M+2]
    @test Tol(getMZ(Frag(PEPTIDE,'p',Int32(2),Int32(0))), 400.687258)
    @test Tol(getMZ(Frag(PEPTIDE,'p',Int32(2),Int32(1))), 401.188771)
    @test Tol(getMZ(Frag(PEPTIDE,'p',Int32(2),Int32(2))), 401.690036)

    @test Tol(getMZ(Frag(PEPTIDE,'p',Int32(3),Int32(0))), 267.460597)
    @test Tol(getMZ(Frag(PEPTIDE,'p',Int32(3),Int32(1))), 267.79494)
    @test Tol(getMZ(Frag(PEPTIDE,'p',Int32(3),Int32(2))), 268.129116)

    PEPTIDEK_mod = Array{Residue, 1}([Residue('P') , Residue('E'), Residue('P'),
                                      Residue('T') , Residue('I'), Residue('D'), 
                                      Residue('E'), Residue("K[+8.014199]")])

    @test Tol(getMZ(Frag(PEPTIDEK_mod,'p',Int32(2),Int32(0))), 468.741839)
    @test Tol(getMZ(Frag(PEPTIDEK_mod,'p',Int32(2),Int32(1))), 469.243364)
    @test Tol(getMZ(Frag(PEPTIDEK_mod,'p',Int32(2),Int32(2))), 469.74462)                                      

    #@time for i in 1:100000000
    #    Frag(TIDEK_mod, 'y', Int32(2))
    #end 

    # #########
    # #Tests for 'FragFilter' Struct 
    # #########

    # #########
    # #Tests for 'Precursor' Struct 
    # #########
    # @test getMZ(Precursor("PEPTIDE", 2))
    # @test getMZ(Precursor("PEPTIDE", 3))
    # @test getMZ(Precursor("PEPTIDE", 4))
    # @test length(Precursor("PEPTIDE", 4))


#     @test getSequence(Precursor("PEPTIDE", 4))
#     @test getResidues(map(v -> getAA(Residue(v)), "PEPTIDE")) = [AA('P'), AA('E'), AA('P'), AA('T'), AA('I'), AA('D'), AA('E')]

#     test_prec = Precursor("PEPTIDE", 3)

#     frag!(test_prec)
#     @test length(getFrags(test_prec))=21

#     #Get the length b3 ion
#     @test length(getFrags(test_prec)) = 3
#     test_frags = getFrags(test_prec)
#     @test length(test_frags[findall(frag -> getType(frag) == 'b' && getCharge(frag) == 1 && length(frag) == 3, test_frags)]) = 3
#     @test getMZ(test_frags[findall(frag -> getType(frag) == 'b' && getCharge(frag) == 1 && length(frag) == 3, test_frags)]) = 324.155397

#     #Sort the fragments from high mass to low mass
#     #Use constructor to get this effect
#     #sort!(test_frags, by = frag -> getMZ(frag), rev = true)
#     sort!(test_frags, rev = true)
#     #y6+1
#     @test getMZ(first(getFrags(test_prec))) == 703.314476
#     @test getType(first(getFrags(test_prec))) == 'y'
#     @test length(first(getFrags(test_prec))) == 6
#     @test getCharge(first(getFrags(test_prec))) == 1
#     #b3+2
#     @test getMZ(last(getFrags(test_prec))) == 162.581336
#     @test getMZ(last(getFrags(test_prec))) == 'b'
#     @test length(last(getFrags(test_prec))) == 3
#     @test getCharge(last(getFrags(test_prec))) == 2
#     #Sort low to high
#     sort!(test_frags, rev = false)
#     #b2+2
#     @test getMZ(first(getFrags(test_prec))) == 162.581336
#     #y6+1
#     @test getMZ(last(getFrags(test_prec))) == 703.314476

#     #Get b3-b4 ions and all y ions 2 or greater
#     #FragFilter({'b'=> ("low" => 3, "")}
#     #    'b' => ("low" => 3, "high" => 4, "max_charge" => 2), 
#     #    'y' => ("low" => 3, "high" => 4, "max_charge" => 2),
#     #    )
#     test_frag_filter = FragFilter(type = 'b', low = 3, high = 4, max_charge = 1) + FragFilter(type = 'y', low = 3, high = ∞, max_charge = 2) 
#     frag!(test_prec, test_frag_filter)
#     sort!(test_frags)
#     @test length(getFrags(test_prec)) = 2 + 9
    
frags = Array{NamedTuple, 1}([
            (ion_type = 'y', ind = Int32(3), charge = Int32(2)),
            (ion_type = 'y', ind = Int32(4), charge = Int32(2)),
            (ion_type = 'y', ind = Int32(5), charge = Int32(2)),
            (ion_type = 'b', ind = Int32(3), charge = Int32(1)),
            (ion_type = 'y', ind = Int32(6), charge = Int32(2)),
            (ion_type = 'y', ind = Int32(3), charge = Int32(1)),
            (ion_type = 'b', ind = Int32(4), charge = Int32(1)),
            (ion_type = 'y', ind = Int32(4), charge = Int32(1)),
            (ion_type = 'y', ind = Int32(5), charge = Int32(1)),
            (ion_type = 'y', ind = Int32(6), charge = Int32(1))])
    
    
test_frags_a = frag!(Peptide("PEPTIDE", Int32(3), default_mods), frags, default_mods)
known_frags_a = [
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
pep = Peptide("PEPTIDE", Int32(3), default_mods)
residues_ = getResidues(pep, default_mods)
#N = 100.000
#test = @timed for i in 1:N
    #frag!(pep, frags, default_mods)
#    map(frag -> getFrag(residues_, frag), frags)
#end 
#println!(test, "total")
#println!(test/N, "per frag")
#N = 249966409
#N = 100000
#test = @timed Threads.@threads for i in 1:N
#    #Frag(TIDEK_mod, 'y', Int32(2))
#    map(frag -> getFrag(residues_, frag), frags)
#end 
#println(test.time, " total")
#println(test.time/(10*N), "per frag")
#println("for ", N, " frags")
@time test = map(pep-> getResidues(pep, default_mods), Vector{String}(["C[Carb]TIDEK[+8.014199]" for x in 1:1000000]))
#function test(pep, frags, default_mods::Dict{String, Float32})
#    #map(frag -> getFrag(residues_, frag), frags)
#    Threads.@threads for i in 1:N
#        frag!(pep, frags, default_mods)
#    end
#end
@time for charge in [Int32(1), Int32(2)]
        for ion_type in ['b','y']
            for pep in test
                getFragIons(pep, 
                            Float32(100.0), 
                            Int32(1), 
                            PROTON*Int32(1),
                            ion_type, 
                            charge,
                            Int32(0))

                getFragIons(reverse(pep), 
                            Float32(100.0), 
                            Int32(1), 
                            PROTON*Int32(1) + H2O,
                            ion_type, 
                            charge,
                            Int32(0))
            end
        end
    end
#@time Threads.@threads testfun(test, default_mods) #end

@test all(map(f -> Tol(f...), zip(known_frags_a, map(frag -> getMZ(frag), test_frags_a))))
@test all(map(f -> Tol(f...), zip(reverse(known_frags_a), map(frag -> getMZ(frag), test_frags_a)))) == false
    
#     #Compare to Skyline

end
