function Tol(a, b, ppm = 2)
    abs(a-b)<=(ppm*minimum((a, b))/1000000)
end
@testset "getMS1PeakHeights.jl" begin

    ##########
    #rangeQuerySorted(sorted_array::Vector{T}, l_bnd::T, u_bnd::T) where T <: Number
    ##########

    rts = Float64[10, 20, 20, 20.1]
    @test rangeQuerySorted(rts, 20.001, 20.1) == (4, 4)
    @test rangeQuerySorted(rts, 19.9, 20.1) == (2, 4)
    @test rangeQuerySorted(rts, 5.0, 20.1) == (1, 4)

    #rangeQuerySorted(sorted_array::Vector{T}, l_bnd::T, u_bnd::T) where T <: Number
    #One of the arguments is an Int but that is not compatible with Vector{Float32} and Float32
    @test_throws MethodError rangeQuerySorted(rts, 5, 20.0)

    
    rts = [10, 11, 12, 13, 14]
    @test rangeQuerySorted(rts, 11, 13) == (2, 4)
    @test rangeQuerySorted(rts, 10, 13) == (1, 4)

    #In this case nothing in the array is in bounds, so gives the last and next to last element
    @test rangeQuerySorted(rts, 15, 15) == (length(rts) + 1, length(rts))
    rts = [10, 10, 10, 11, 12, 13, 14]

    #Nothing in the array is in bounds but the query is less than the minimum value
    @test rangeQuerySorted(rts, 9, 9) == (1, 0)

    ##########
    #getMS1Peaks!
    ##########

#=
```
getMS1Peaks!(precursors::Dictionary{UInt32, Precursor}, 
                        MS1::Vector{Union{Missing, Float32}}, 
                        INTENSITIES::Vector{Union{Missing, Float32}}, 
                        MS1_MAX_HEIGHTS::UnorderedDictionary{UInt32, Float32}, 
                        precursor_rts::Vector{Float32}, 
                        precursor_idxs::Vector{UInt32}, 
                        precursor_ms_file_idxs::Vector{UInt32}, 
                        rt::Float32, rt_tol::Float32, left_mz_tol::Float32, right_mz_tol::Float32, ms_file_idx::UInt32)
```
=#
using Dictionaries
precursors = Dictionary(
    UInt32[1, 2, 3, 4],
    [Precursor(seq) for seq in ["PEPTIDE","DICER","DRAGRACER","AMINEACID"]]
)

#=
precursor mz's
400.68726
400.68726
318.14453
318.14453
517.25146
490.2148
=#

precursor_rts = Float32[10.0, 10.0, 11.0, 11.0, 11.0, 25.0]
precursor_idxs = UInt32[1, 1, 2, 2, 3, 4]
precursor_ms_file_idxs = UInt32[1, 1, 1, 1, 1, 2]

MS1 = allowmissing(Float32[400.68726,500.0,
                            318.14453,318.14453,517.25146,490.2148])
INTENSITIES = allowmissing(Float32[1000,2000,
                                    3000, 4000, 1500, 500])

MS1_MAX_HEIGHTS = UnorderedDictionary{UInt32, Float32}()

#Match the peak for "PEPTIDE"
getMS1Peaks!(MS1_MAX_HEIGHTS,
            precursors,
            MS1,INTENSITIES,
            precursor_rts,precursor_idxs,precursor_ms_file_idxs,
            Float32(10.0), Float32(0.1), Float32(0.1), Float32(0.1), UInt32(1))

@test length(MS1_MAX_HEIGHTS) == 1
@test Tol(MS1_MAX_HEIGHTS[UInt32(1)], 1000.0f0)

#Increase the intensitiy of the matching peak
getMS1Peaks!(MS1_MAX_HEIGHTS,
            precursors,
            MS1,allowmissing(Float32[2000, 2000, 3000, 4000, 1500, 500]),
            precursor_rts,precursor_idxs,precursor_ms_file_idxs,
            Float32(10.0), Float32(0.1), Float32(0.1), Float32(0.1), UInt32(1))

@test length(MS1_MAX_HEIGHTS) == 1
@test Tol(MS1_MAX_HEIGHTS[UInt32(1)], 2000.0f0)

#Increase again but use the wrong ms_file_idx
getMS1Peaks!(MS1_MAX_HEIGHTS,
            precursors,
            MS1,allowmissing(Float32[3000, 2000, 3000, 4000, 1500, 500]),
            precursor_rts,precursor_idxs,precursor_ms_file_idxs,
            Float32(10.0), Float32(0.1), Float32(0.1), Float32(0.1), UInt32(2))

@test length(MS1_MAX_HEIGHTS) == 1
@test Tol(MS1_MAX_HEIGHTS[UInt32(1)], 2000.0f0)

#Match the peak for "DICER" AND "DRAGRACER" but with the wrong ms_file_idx
getMS1Peaks!(MS1_MAX_HEIGHTS,
            precursors,
            MS1,INTENSITIES,
            precursor_rts,precursor_idxs,precursor_ms_file_idxs,
            Float32(11.0), Float32(0.1), Float32(0.1), Float32(0.1), UInt32(2))

@test !isassigned(MS1_MAX_HEIGHTS, UInt32(2))

#Now with the correct ms_file_idx
getMS1Peaks!(MS1_MAX_HEIGHTS,
            precursors,
            MS1,INTENSITIES,
            precursor_rts,precursor_idxs,precursor_ms_file_idxs,
            Float32(11.0), Float32(0.1), Float32(0.1), Float32(0.1), UInt32(1))

@test length(MS1_MAX_HEIGHTS) == 3
@test Tol(MS1_MAX_HEIGHTS[UInt32(2)], 3000.0f0)
@test Tol(MS1_MAX_HEIGHTS[UInt32(3)], 1500.0f0)

#Rt that does not correspond to any precursor
getMS1Peaks!(MS1_MAX_HEIGHTS,
            precursors,
            MS1,INTENSITIES,
            precursor_rts,precursor_idxs,precursor_ms_file_idxs,
            Float32(25.2), Float32(0.1), Float32(0.1), Float32(0.1), UInt32(2))
            
@test length(MS1_MAX_HEIGHTS) == 3
@test Tol(MS1_MAX_HEIGHTS[UInt32(1)], 2000.0f0)
@test Tol(MS1_MAX_HEIGHTS[UInt32(2)], 3000.0f0)
@test Tol(MS1_MAX_HEIGHTS[UInt32(3)], 1500.0f0)

#RT for 
getMS1Peaks!(MS1_MAX_HEIGHTS,
            precursors,
            MS1,INTENSITIES,
            precursor_rts,precursor_idxs,precursor_ms_file_idxs,
            Float32(25.2), Float32(0.3), Float32(0.1), Float32(0.1), UInt32(2))
            
@test length(MS1_MAX_HEIGHTS) == 4
@test Tol(MS1_MAX_HEIGHTS[UInt32(4)], 500.0f0)

getMS1Peaks!(MS1_MAX_HEIGHTS,
            precursors,
            allowmissing(Float32[10000]),allowmissing(Float32[10000]),
            precursor_rts,precursor_idxs,precursor_ms_file_idxs,
            Float32(0), Float32(0.3), Float32(0.1), Float32(0.1), UInt32(2))
            
@test length(MS1_MAX_HEIGHTS) == 4
@test Tol(MS1_MAX_HEIGHTS[UInt32(4)], 500.0f0)
 
##########
#getMS1PeakHeights!
##########
#=
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

#Make a few "fake" scans that have MS1 peaks for some of the peptides in testPtable. 
#=for (key, value) in pairs(testPtable.id_to_prec)
    println("key ", key, " mz ", getMZ(value))
end
key 18 mz 200.84726
key 12 mz 259.8664
key 17 mz 267.4606
key 3 mz 273.3847
key 9 mz 275.8847
key 15 mz 304.9024
key 6 mz 306.9024
key 11 mz 346.15274
key 2 mz 364.17715
key 8 mz 367.5105
key 16 mz 400.68726
key 14 mz 406.20078
key 5 mz 408.86743
key 10 mz 518.7255
key 1 mz 545.76215
key 7 mz 550.76215
key 13 mz 608.79755
key 4 mz 612.79755
=#

#Fake retention times

precursor_rts = Float32[1, 1, 1, 2, 2.5, 2.5, 3, 4, 5, 6, 7, 8]
precursor_idxs = keys(testPtable.id_to_prec)
precursor_ms_file_idxs = UInt32[1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
MS1_MAX_HEIGHTS = UnorderedDictionary{UInt32, Float32}()
masses = [[],
        ]
retention_times = [1, 2, 2.5, 3, 4, 5, 6, 7, 8]
getMS1PeakHeights!(MS1_MAX_HEIGHTS,
                    testPtable,
                    ,
                    ,
                    ,
                    ,
                    precursor_rts,precursor_idxs,precursor_ms_file_idxs,
                    0.1, 0.001, 0.001, UInt32(1)
)

=#
end