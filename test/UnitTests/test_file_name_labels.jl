# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

@testset "file names: full for tables, distinguishing for plot labels" begin
    paths = ["/d/20240101_lab_rep1.arrow", "/d/20240101_lab_rep2.arrow", "/d/20240101_lab_lowsignal.arrow"]
    full = Pioneer.parseFileNames(paths)
    @test full == ["20240101_lab_rep1", "20240101_lab_rep2", "20240101_lab_lowsignal"]
    @test Pioneer.distinguishingFileNames(full) == ["rep1", "rep2", "lowsignal"]

    # Differing token counts: nothing is stripped.
    @test Pioneer.distinguishingFileNames(["a_b_c", "a_b"]) == ["a_b_c", "a_b"]
    # Nothing distinguishes the runs: the full name stands in.
    @test Pioneer.distinguishingFileNames(["same_name", "same_name"]) == ["same_name", "same_name"]
    @test Pioneer.distinguishingFileNames(["only_one"]) == ["only_one"]

    # convertSciex's `.zt` belongs to the extension, not the run name
    @test Pioneer.parseFileNames(["/d/run1.zt.scxs", "/d/run2.zt.arrow", "/d/run3.scxs", "/d/zt.arrow"]) ==
          ["run1", "run2", "run3", "zt"]
end
