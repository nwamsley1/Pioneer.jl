# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

using Arrow

@testset "acquisition mode from the method summary" begin
    mode(t) = (a = Pioneer.SpectrumDecoder.acquisition_mode_from_method(t); a === nothing ? nothing : a.mode)
    exp(kinds...) = join(["Experiment $(i)\nExperiment Name = x\nScan $(k)\nMSn Level = 1" for (i, k) in enumerate(kinds)], "\n")
    @test mode(exp("MasterScan", "DIAScan")) == "dia"
    @test mode(exp("MasterScan", "tMSn")) == "dia"
    @test mode(exp("MasterScan", "ddMSnScan")) == "dda"
    @test mode(exp("MasterScan", "ddMSnScan") * "\nFilter TargetedInclusionMassList") == "targeted"
    @test mode(exp("MasterScan")) == "not_searchable"
    @test mode(exp("tSIM")) == "not_searchable"
    @test mode(exp("MasterScan", "SomeOtherScan")) === nothing
    @test mode("   Experiments\n\nFULL MS — SIM\nFull MS — SIM\n\nDIA\nGeneral\nDIA\n") == "dia"
    @test mode("   Experiment\n\nFULL MS / DD-MS² (TOPN)\nFull MS\ndd-MS² / dd-SIM\nTopN   20\n") == "dda"
    @test mode("Number of Scan Events: 6\nData Dependent Settings:\n  Charge state rejection enabled\n") == "dda"
    @test mode("") === nothing
    @test mode("Scan Range (m/z) = 350-1500") === nothing
    a = Pioneer.SpectrumDecoder.acquisition_mode_from_method(exp("MasterScan", "DIAScan"))
    @test a.evidence == "method: Scan DIAScan"
end

@testset "acquisition mode label in the .arrow metadata" begin
    mktempdir() do d
        f = open_acquisition(write_synthetic_acquisition(joinpath(d, "s.raw")))
        @test acquisition_mode(f) === nothing          # no embedded method
        out = write_pioneer_arrow(f, joinpath(d, "s.arrow"))
        @test Arrow.getmetadata(Arrow.Table(out)) === nothing
    end
end
