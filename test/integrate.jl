using Titus
using Test
using Tables, Arrow
println(pwd())
include("../src/integrate.jl")
GAPDH_VGVNGFGR = Arrow.Table("../data/parquet/GAPDH_VGVNGFGR.arrow")
GSTM1_RPWFAGNK = Arrow.Table("../data/parquet/GSTM1_RPWFAGNK.arrow")
GSTM4_VAVWGNK = Arrow.Table("../data/parquet/GSTM4_VAVWGNK.arrow")
NRF2_Survey = Arrow.Table("../data/parquet/Nrf2_SQ_052122_survey.arrow")
#Avoid compiliation time from first call
#@time table = Arrow.Table("../data/parquet/GAPDH_VGVNGFGR.arrow")


@testset "integrate.jl" begin

    ########
    #getSub
    ########

    #Get scans indices where the precursor mass is within a given
    #tolerance of the expected mass for a peptide.
    @test getSub(Float32(300.0), [200.0, 300.0, 300.1, 300.00001, 10000.0], ppm = Float32(20.0)) == [2, 4]
    @test getSub(Float32(300.0), [200.0, 300.0, 300.1, 300.00001, 10000.0], ppm = Float32(0)) == []
    @test getSub(Float32(300.0), [300.0], ppm = Float32(20.0)) == [1]


    ########
    #getScanAdresses
    ########

    #msOrder is a list of scan orders (MS1, MS2, MS3, etc.)
    #We want a list of "scan adresses" from a list of scan orders
    #A scan adress records the each cycle (defined by an MS1 scan)
    #and each MSN scan is numbered starting from 1 within a cycle
    #See the example below 

    test_msOrder = [1, 2, 2, 2, 1, 2, 1, 1, 2, 2];
    test_adresses = Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}}([
     (scan_index = 1, ms1 = 1, msn = 0),
     (scan_index = 2, ms1 = 1, msn = 1),
     (scan_index = 3, ms1 = 1, msn = 2),
     (scan_index = 4, ms1 = 1, msn = 3),
     (scan_index = 5, ms1 = 2, msn = 0),
     (scan_index = 6, ms1 = 2, msn = 1),
     (scan_index = 7, ms1 = 3, msn = 0),
     (scan_index = 8, ms1 = 4, msn = 0),
     (scan_index = 9, ms1 = 4, msn = 1),
     (scan_index = 10, ms1 = 4, msn = 2)]);

    @test getScanAdresses(test_msOrder) == test_adresses


    lightVGVNGFGR = getMZ(Precursor(getResidues("VGVNGFGR"), UInt8(2)))
    heavyVGVNGFGR = getMZ(Precursor(getResidues("VGVNGFGR[+10.008269]"), UInt8(2)))

    ########
    #getScanCycleUnion
    ########

    test_adresses_a = Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}}([
        (scan_index = 3, ms1 = 1, msn = 2), # 1
        (scan_index = 4, ms1 = 1, msn = 3), # 1
        (scan_index = 6, ms1 = 2, msn = 1), #Also in list 2
        (scan_index = 7, ms1 = 3, msn = 1), #Cycle 3 in other list
        (scan_index = 9, ms1 = 5, msn = 1)]);

    test_adresses_b = Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}}([
        (scan_index = 2, ms1 = 1, msn = 1),
        (scan_index = 6, ms1 = 2, msn = 1), #2
        (scan_index = 8, ms1 = 4, msn = 1), #Cycle 4 not in other list
        (scan_index = 9, ms1 = 4, msn = 10),
        (scan_index = 10, ms1 = 5, msn = 20)]);

    @test getScanCycleUnion(test_adresses_a, test_adresses_b) == [1, 2, 5]
    #Given two lists of can adresses, we want to know which MS cycles
    #have at least one scan from both lists. Then we want to keep those
    #scans and discard the rest.

    ########
    #getIntegrationBounds
    ########

    #Given a list of integers the maximum length subarray
    #such that each element is one greater than the one before it 
    #starting with the second element in the subarray. 
    #A maximum gap size should be allowed. If an element in the 
    #array is exactly one greater than the previous element, call
    #it an anchor. The maximum difference between two anchors must
    #not exceed the maximum gap size for the subarray to be valid
    #See the examples below

    test_indices_a = Vector{Int64}([1, 2, 3, 12, 13])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (1, 5)
    test_indices_a = Vector{Int64}([1, 2, 3, 12, 100])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (1, 3)
    test_indices_a = Vector{Int64}([1, 2, 3, 12])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (1, 3)
    test_indices_a = Vector{Int64}([1, 2, 3, 12, 100, 101, 102, 103, 104])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (5, 9)
    test_indices_a = Vector{Int64}([1])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (1, 1)
    test_indices_a = Vector{Int64}([1,2])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (1, 2)
    test_indices_a = Vector{Int64}([1, 2, 10, 11, 12, 13, 60, 61, 62, 63, 64, 65, 66, 67, 68, 100])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (7, 15)

    ########
    #getIntegrationBounds
    ########

    test_masses = Vector{Float32}([151.67221f0, 700.0, 894.0938f0])
    #test_mz = Vector{Float32}([])

end