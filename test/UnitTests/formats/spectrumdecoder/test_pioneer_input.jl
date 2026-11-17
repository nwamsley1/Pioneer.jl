# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# The written `.arrow` opens through Pioneer's own loader with the expected values.

@testset "Pioneer reads the written file" begin
    p = write_synthetic_acquisition(joinpath(synthetic_dir, "pioneer.raw"); cycles = 5, windows = 10)
    f = open_acquisition(p)
    out = write_pioneer_arrow(f, joinpath(synthetic_dir, "pioneer.arrow"))
    d = Pioneer.loadMassSpecData(out)
    @test d isa Pioneer.BasicNonIonMobilityMassSpecData{Float32}
    @test length(d) == 55
    @test Pioneer.getMsOrder(d, 1) == 0x01 && Pioneer.getMsOrder(d, 2) == 0x02
    @test collect(Pioneer.getMzArray(d, 2)) == peaks(f, 2)[1]
    @test collect(Pioneer.getIntensityArray(d, 2)) == peaks(f, 2)[2]
    @test Pioneer.getCenterMz(d, 2) == 405f0 && ismissing(Pioneer.getCenterMz(d, 1))
    @test Pioneer.getRetentionTime(d, 3) == Float32(scans(f)[3].retention_time)
    @test Pioneer.getCycleIdx(d, 13) == 2
end
