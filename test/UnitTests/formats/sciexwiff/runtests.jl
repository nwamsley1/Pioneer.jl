# SciexWiff (src/vendor/SciexWiff, folded in from SciexWiff.jl). Wrapped in a module so the submodule's exports
# (ConvertParams, bin_to_mz, ...) do not leak into the shared test namespace or clash with TimsSlices'.
module SciexWiffTests
using Test
using Pioneer
using Pioneer.SciexWiff

# Real-data tests run only when a data directory is given (PIONEER_FORMATS_TEST_DATA, or SciexWiff's SCIEXWIFF_DATA).
const DATA_DIR = get(ENV, "PIONEER_FORMATS_TEST_DATA", get(ENV, "SCIEXWIFF_DATA", ""))

@testset "SciexWiff" begin
    include("test_cfb.jl")
    include("test_scan.jl")
    include("test_container.jl")
end

end # module SciexWiffTests
