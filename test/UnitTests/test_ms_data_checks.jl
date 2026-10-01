# The search refuses, before any search work, to mix timsTOF .tdfs runs with .arrow files and to search .tdfs data
# with a library that has no predicted ion mobility.

using Test
using Pioneer
using Pioneer: check_ms_data_vendors, check_library_ion_mobility

@testset "timsTOF input checks" begin
    root = mktempdir()
    tdfs = mkpath(joinpath(root, "a.tdfs"))
    arrow = joinpath(root, "b.arrow"); write(arrow, "x")
    # vendors
    @test check_ms_data_vendors([tdfs]) === nothing
    @test check_ms_data_vendors([arrow]) === nothing
    err = try check_ms_data_vendors([tdfs, arrow]); nothing catch e; e end
    @test err isa ErrorException && occursin("mixes 1 timsTOF .tdfs run with 1 .arrow file", err.msg)
    # library ion mobility
    @test check_library_ion_mobility([tdfs], true, "lib.poin") === nothing
    @test check_library_ion_mobility([arrow], false, "lib.poin") === nothing   # Thermo data needs none
    err = try check_library_ion_mobility([tdfs], false, "lib.poin"); nothing catch e; e end
    @test err isa ErrorException && occursin("no ion-mobility predictions", err.msg) &&
          occursin("alphapept_ccs", err.msg)
end
