# ZT mode selection: the file's acquisition metadata alone decides (the user declares ZT when converting).

using Test, Arrow
using Pioneer
using Pioneer: zt_mode, zt_search_k, getAcquisitionMetadata, BasicMassSpecData

function tiny_arrow(path; metadata = nothing)
    tbl = (mz_array = [Union{Missing,Float32}[100f0]], intensity_array = [Union{Missing,Float32}[1f0]],
           scanHeader = [""], scanNumber = Int32[1], packetType = Int32[0], retentionTime = Float32[0.1],
           lowMz = Float32[100], highMz = Float32[1500], TIC = Float32[1], centerMz = Union{Missing,Float32}[missing],
           isolationWidthMz = Union{Missing,Float32}[missing], collisionEnergyField = Union{Missing,Float32}[missing],
           collisionEnergyEvField = Float32[0], msOrder = UInt8[1])
    metadata === nothing ? Arrow.write(path, tbl) : Arrow.write(path, tbl; metadata = metadata)
    BasicMassSpecData(path)
end

@testset "zt_mode" begin
    d = mktempdir()
    zt = tiny_arrow(joinpath(d, "zt.arrow"); metadata = Dict("acquisition_type" => "zt_scan_dia", "q1_bin_step_mz" => "1.0221"))
    sw = tiny_arrow(joinpath(d, "swath.arrow"); metadata = Dict("acquisition_type" => "swath"))
    plain = tiny_arrow(joinpath(d, "plain.arrow"))
    @test getAcquisitionMetadata(zt)["q1_bin_step_mz"] == "1.0221"
    @test getAcquisitionMetadata(plain) === nothing
    @test zt_mode(zt)
    @test !zt_mode(sw)
    @test !zt_mode(plain)
    @test zt_search_k(6) == 3 && zt_search_k(5) == 3 && zt_search_k(7) == 4
end
