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




    test_indices_a = Vector{Int64}([1, 2, 10, 11, 12, 13, 60, 61, 62, 63, 64, 65, 66, 67, 68, 100])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (7, 15)
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

end