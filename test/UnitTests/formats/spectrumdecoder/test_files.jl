# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Acquisition files in PIONEER_FORMATS_TEST_DATA: every consistency check passes, and the written `.arrow` opens in
# Pioneer with one row per scan.

@testset "acquisition files" begin
    for p in RAW_FILES
        @testset "$(basename(p))" begin
            f = open_acquisition(p)
            R = check_consistency(f)
            failed = [r.name for r in R if !SpectrumDecoder.ok(r)]
            @test isempty(failed)
            isempty(failed) || @info "$(basename(p)): checks not passing" failed
            out = write_pioneer_arrow(f, joinpath(synthetic_dir, splitext(basename(p))[1] * ".arrow"))
            d = Pioneer.loadMassSpecData(out)
            @test length(d) == length(f)
            rm(out; force = true)
        end
    end
end
