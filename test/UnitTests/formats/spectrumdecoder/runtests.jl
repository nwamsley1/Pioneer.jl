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

@testset "SpectrumDecoder" begin
    include("test_centroid.jl")
end
end
