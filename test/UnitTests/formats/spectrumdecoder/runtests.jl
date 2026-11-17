# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# SpectrumDecoder (src/vendor/SpectrumDecoder). Wrapped in a module so the submodule's exports do not leak into the
# shared test namespace.
module SpectrumDecoderTests
using Test
using Pioneer
using Pioneer.SpectrumDecoder
import Pioneer.SpectrumDecoder: _fit_peak, _median_positive, _median_positive!

# Real-data tests run only when a data directory is given (PIONEER_FORMATS_TEST_DATA): every `.raw` file in it.
const DATA_DIR = get(ENV, "PIONEER_FORMATS_TEST_DATA", "")
const RAW_FILES = isdir(DATA_DIR) ? filter(p -> endswith(lowercase(p), ".raw"), readdir(DATA_DIR; join = true)) : String[]

@testset "SpectrumDecoder" begin
    include("test_centroid.jl")
    include("test_synthetic.jl")
    include("test_pioneer_input.jl")
    if isempty(RAW_FILES)
        @info "No .raw files in PIONEER_FORMATS_TEST_DATA; skipping SpectrumDecoder real-data tests"
    else
        include("test_files.jl")
    end
end

end # module SpectrumDecoderTests
