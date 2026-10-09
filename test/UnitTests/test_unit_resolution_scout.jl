# Unit-resolution scout: the converter's `mass_resolution` header key sets the wide-scout fragment window.

using Test, Arrow
using Pioneer
using Pioneer: scout_unit_resolution_tol, BasicMassSpecData

function unit_res_arrow(path; metadata = nothing)
    tbl = (mz_array = [Union{Missing,Float32}[100f0]], intensity_array = [Union{Missing,Float32}[1f0]],
           scanHeader = [""], scanNumber = Int32[1], packetType = Int32[0], retentionTime = Float32[0.1],
           lowMz = Float32[100], highMz = Float32[1500], TIC = Float32[1], centerMz = Union{Missing,Float32}[missing],
           isolationWidthMz = Union{Missing,Float32}[missing], collisionEnergyField = Union{Missing,Float32}[missing],
           collisionEnergyEvField = Float32[0], msOrder = UInt8[1])
    metadata === nothing ? Arrow.write(path, tbl) : Arrow.write(path, tbl; metadata = metadata)
    BasicMassSpecData(path)
end

@testset "scout_unit_resolution_tol" begin
    d = mktempdir()
    # PioneerConverter header for a Stellar (ion-trap MS2) file
    stellar = unit_res_arrow(joinpath(d, "stellar.arrow"); metadata = Dict(
        "instrument_model" => "Stellar", "ms2_mass_analyzer" => "ITMS", "mass_resolution" => "0.5"))
    # Orbitrap/Astral files carry no mass_resolution key
    astral = unit_res_arrow(joinpath(d, "astral.arrow"); metadata = Dict(
        "instrument_model" => "Orbitrap Astral", "ms2_mass_analyzer" => "ASTMS"))
    plain = unit_res_arrow(joinpath(d, "plain.arrow"))
    bad = unit_res_arrow(joinpath(d, "bad.arrow"); metadata = Dict("mass_resolution" => "n/a"))
    zero_res = unit_res_arrow(joinpath(d, "zero.arrow"); metadata = Dict("mass_resolution" => "0"))
    @test scout_unit_resolution_tol(stellar) === 0.5f0
    @test scout_unit_resolution_tol(astral) === nothing
    @test scout_unit_resolution_tol(plain) === nothing
    @test scout_unit_resolution_tol(bad) === nothing
    @test scout_unit_resolution_tol(zero_res) === nothing
end
