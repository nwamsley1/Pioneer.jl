# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# convertRaw: input handling, output naming, --skip-existing, per-file failures and the CLI exit codes.

using Test
using Pioneer
using Pioneer: convertRaw, main_convertRaw

@testset "convertRaw" begin
    root = mktempdir()
    quiet(f) = redirect_stdout(f, devnull)
    @test_throws ArgumentError convertRaw(joinpath(root, "missing"))
    empty = mkpath(joinpath(root, "empty"))
    mkpath(joinpath(empty, "folder.raw"))                     # a folder named .raw is not a file to convert
    @test_throws ArgumentError convertRaw(empty)

    data = mkpath(joinpath(root, "data"))
    Pioneer.SpectrumDecoder.write_synthetic_acquisition(joinpath(data, "a.raw"))
    Pioneer.SpectrumDecoder.write_synthetic_acquisition(joinpath(data, "b.RAW"); ms1_profile = true)
    write(joinpath(data, "notes.txt"), "not an acquisition")

    # folder input, default output folder
    out = quiet(() -> convertRaw(data))
    @test out == [joinpath(data, "arrow_out", "a.arrow"), joinpath(data, "arrow_out", "b.arrow")]
    @test all(isfile, out)
    @test !any(endswith(".partial"), readdir(joinpath(data, "arrow_out")))
    @test length(Pioneer.loadMassSpecData(out[1])) == 55

    # single-file input, default output beside it; --skip-existing leaves an existing output untouched
    single = quiet(() -> convertRaw(joinpath(data, "a.raw")))
    @test single == [joinpath(data, "arrow_out", "a.arrow")]
    mtime_before = mtime(single[1])
    sleep(1.1)
    quiet(() -> convertRaw(joinpath(data, "a.raw"); skip_existing = true))
    @test mtime(single[1]) == mtime_before

    # an unreadable file fails on its own; the other files are still written
    bad = mkpath(joinpath(root, "bad"))
    cp(joinpath(data, "a.raw"), joinpath(bad, "good.raw"))
    write(joinpath(bad, "broken.raw"), zeros(UInt8, 64))
    @test_throws ErrorException quiet(() -> convertRaw(bad; output_dir = joinpath(root, "bad_out")))
    @test readdir(joinpath(root, "bad_out")) == ["good.arrow"]

    # the CLI returns 1 rather than throwing, and 1 for no arguments (help)
    @test quiet(() -> main_convertRaw(String[])) == 1
    @test quiet(() -> main_convertRaw([empty])) == 1
    @test quiet(() -> main_convertRaw(["--bogus"])) == 1
    @test quiet(() -> main_convertRaw([data, "--threads-per-file", "0"])) == 1
    @test quiet(() -> main_convertRaw(["--help"])) == 0
    @test quiet(() -> main_convertRaw([data, "-o", joinpath(root, "cli"), "--threads-per-file", "2"])) == 0
    @test readdir(joinpath(root, "cli")) == ["a.arrow", "b.arrow"]
end
