

@testset "getBestPSMs.jl" begin


using DataFrames


"""
    data/peptide_lists/PROT_PEPTIDE_TEST1.txt file looks like this.

   PROT_A	PEPTIDE
    PROT_B	PEPTIDE	
    PROT_A	PEPTIDE
    PROT_C	PEPTIDE
    PROT_A	AMINEACID
    PROT_B	PEPTICKDEK
    PROT_C	DRAGRACER

    See test/Routines/PRM/IS-PRM-Survey/buildPrecursorTable.jl

    This test set makes a fake table of psms, runs "getBestPSMs", and ensures the result is correct. 
"""

#Make the PrecursorDatabase
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
addPrecursors!(testPtable, UInt8[2, 3, 4], UInt8[0], test_mods)

#Fake retention times for MS files
MS_RT = Dict{UInt32, Vector{Float32}}(
    UInt32(1) => Float32[1, 2, 3, 4, 5, 6, 7],
    UInt32(2) => Float32[1, 2, 3]
)

#Fake list of PSMs
PSMs = Dict{Symbol, Vector}(
:hyperscore    => Float64[1.0, 5, 10.0, 10.0, 12.0, 5.0, 1.0, 10, 100, 1000, 100, 100, 200],
:precursor_idx => UInt32[1, 2, 1, 1, 4, 5, 6, 7, 8, 9, 16, 17, 18],
:total_ions    => UInt32[1, 8, 8, 8, 3, 5, 5, 10, 10, 10, 3, 5, 5],
:ms_file_idx   => UInt32[1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1],
:error         => [.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1, 1, 1, 1, 1, 1, 1],
:scan_idx      => [1, 2, 3, 4, 1, 2, 3, 1, 2, 3, 1, 2, 3],
:y_ladder      => Int8[3, 7, 7, 7, 3, 3, 3, 8, 8, 8, 8, 8, 8]
)

#Correct "best_psms" when the minimum_fragment_count is 5
best_psms = DataFrame(Dict{Symbol, Vector}(
    :ms_file_idx => UInt32[1, 2, 1, 1, 2],
    :pep_idx => UInt32[4, 1, 1, 5, 6],
    :error => Float64[0.01, 1.0, 1.0, 0.01, 1.0],
    :hyperscore => Float64[5, 100, 200, 10, 1000],
    :precursor_idx => UInt32[5, 17, 18, 1, 9],
    :scan_idx =>  Int64[2, 2, 3, 3, 3],
    :total_ions => UInt32[5, 5, 5, 8, 10],
    :y_ladder => Int8[3, 8, 8, 7, 8],
    :precursor_charge => UInt8[3, 3, 4, 2, 4],
    :precursor_isotope => UInt8[0, 0, 0, 0, 0],
    :retention_time => Float32[2.0, 2.0, 3.0, 3.0, 3.0],
    :sequence => String["PEPTIC[Carb]KDEK[Hlys]", "PEPTIDE","PEPTIDE","DRAGRAC[Carb]ER","DRAGRAC[Carb]ER[Harg]"],
    :precursor_mz => Float32[408.867, 267.461, 200.847, 545.762,  275.885],
    :protein_names => String["PROT_B", "PROT_A|PROT_B|PROT_C", "PROT_A|PROT_B|PROT_C", "PROT_C","PROT_C"],
))
out = getBestPSMs(PSMs, testPtable, MS_RT, UInt8(5)) # Get best psms with a min_fragment_count of 5
for col in names(out) #Test 
    if any(x -> isa(x, AbstractFloat), out[!, col])
        @test abs(sum(out[!, col] .- best_psms[!, col])) < 0.001
    else
        @test out[!, col] == best_psms[!, col]
    end
end

#Repeat above but with a minimum fragment count of 2. 
best_psms = DataFrame(Dict{Symbol, Vector}(
    :ms_file_idx => UInt32[1, 2, 1, 1, 2],
    :pep_idx => UInt32[4, 1, 1, 5, 6],
    :error => Float64[0.01, 1.0, 1.0, 0.01, 1.0],
    :hyperscore => Float64[12.0, 100.0, 200.0, 10.0, 1000.0],
    :precursor_idx => UInt32[4, 16, 18, 1, 9],
    :scan_idx =>  Int64[1, 1, 3, 3, 3],
    :total_ions => UInt32[3, 3, 5, 8, 10],
    :y_ladder => Int8[3, 8, 8, 7, 8],
    :precursor_charge => UInt8[2, 2, 4, 2, 4],
    :precursor_isotope => UInt8[0, 0, 0, 0, 0],
    :retention_time => Float32[1.0, 1.0, 3.0, 3.0, 3.0],
    :sequence => String["PEPTIC[Carb]KDEK[Hlys]", "PEPTIDE","PEPTIDE","DRAGRAC[Carb]ER","DRAGRAC[Carb]ER[Harg]"],
    :precursor_mz => Float32[612.798 , 400.687, 200.847, 545.762, 275.885],
    :protein_names => String["PROT_B", "PROT_A|PROT_B|PROT_C", "PROT_A|PROT_B|PROT_C", "PROT_C","PROT_C"],
))

out = getBestPSMs(PSMs, testPtable, MS_RT, UInt8(3))
for col in names(out)
    if any(x -> isa(x, AbstractFloat), out[!, col])
        @test abs(sum(out[!, col] .- best_psms[!, col])) < 0.001
    else
        @test out[!, col] == best_psms[!, col]
    end
end

end
