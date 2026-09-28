# SciexWiff (src/vendor/SciexWiff, folded in from SciexWiff.jl). Wrapped in a module so the submodule's exports
# (ConvertParams, bin_to_mz, ...) do not leak into the shared test namespace or clash with TimsSlices'.
module SciexWiffTests
using Test
using Pioneer
using Pioneer.SciexWiff

# Real-data tests run only when a data directory is given (PIONEER_FORMATS_TEST_DATA, or SciexWiff's SCIEXWIFF_DATA).
const DATA_DIR = get(ENV, "PIONEER_FORMATS_TEST_DATA", get(ENV, "SCIEXWIFF_DATA", ""))
# The truncated SCIEX fixture: downloaded from Zenodo into temp/zenodo by .github/actions/precompile-data
# (PIONEER_FORMATS_FIXTURES overrides the directory).
const FIXTURES = get(ENV, "PIONEER_FORMATS_FIXTURES", joinpath(pkgdir(Pioneer), "temp", "zenodo"))
const FIXTURE_WIFF = joinpath(FIXTURES, "sciex_wiff_fixture", "BenchSample_B_nswath4_25ng.wiff")

@testset "SciexWiff" begin
    include("test_cfb.jl")
    include("test_scan.jl")
    include("test_container.jl")
    if isfile(FIXTURE_WIFF)
        include("test_fixture.jl")
    else
        @warn "SCIEX .wiff fixture not found at $FIXTURE_WIFF (Zenodo download); skipping fixture tests"
    end
end

end # module SciexWiffTests
