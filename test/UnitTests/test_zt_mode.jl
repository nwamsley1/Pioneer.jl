# ZT mode selection: config `acquisition.scanning_quad` wins; otherwise the file's acquisition metadata decides.

using Test, Arrow
using Pioneer
using Pioneer: zt_mode, getAcquisitionMetadata, BasicMassSpecData

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
    none = (nce = 26,)
    @test zt_mode(none, zt) == (true, "file metadata", true)
    @test zt_mode(none, sw) == (false, "file metadata", false)
    @test zt_mode(none, plain) == (false, "file metadata", false)
    @test zt_mode((nce = 26, scanning_quad = false), zt) == (false, "config", true)
    @test zt_mode((nce = 26, scanning_quad = true), plain) == (true, "config", false)
end
