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

"""
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
"""
using Dicsionaries

precursors = Dictionary(
    UInt32[1, 2, 3, 4],
    [Precursor(seq) for seq in ["PEPTIDE","DICER","DRAGRACER","AMINEACID"]]
)
"""
precursor mz's
400.68726
400.68726
318.14453
318.14453
517.25146
490.2148
"""
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
 
##########
#getMS1PeakHeights!
##########


end