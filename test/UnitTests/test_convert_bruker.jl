# convertBruker input handling: what counts as a Bruker bundle, and the errors for inputs that are not. The
# conversion itself is TimsSlices.convert, tested in TimsSlices.jl and end to end by the timsTOF searches.

using Test
using Pioneer
using Pioneer: convertBruker, main_convertBruker

@testset "convertBruker input handling" begin
    root = mktempdir()
    # not a path at all
    @test_throws ArgumentError convertBruker(joinpath(root, "missing"))
    # a folder with no .d bundles in it (a file named .d does not count)
    empty = mkpath(joinpath(root, "empty"))
    write(joinpath(empty, "notes.d"), "a file, not a bundle")
    @test_throws ArgumentError convertBruker(empty)
    # the CLI returns 1 rather than throwing, and 1 for no arguments (help)
    @test redirect_stdout(devnull) do
        main_convertBruker(String[])
    end == 1
    @test redirect_stdout(devnull) do
        main_convertBruker([empty])
    end == 1
    @test redirect_stdout(devnull) do
        main_convertBruker(["--bogus"])
    end == 1
    @test redirect_stdout(devnull) do
        main_convertBruker(["--help"])
    end == 0
end
