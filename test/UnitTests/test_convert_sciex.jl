# convertSciex input handling: what counts as a SCIEX run, and the errors for inputs that are not. The conversion
# itself is SciexWiff.convert, tested in SciexWiff.jl. With PIONEER_TEST_WIFF set to a .wiff (its .wiff.scan beside
# it), also converts it and opens the result as MS data.

using Test
using Pioneer
using Pioneer: convertSciex, main_convertSciex, ask_zt_scan

@testset "convertSciex input handling" begin
    root = mktempdir()
    # not a path at all
    @test_throws ArgumentError convertSciex(joinpath(root, "missing"))
    # a folder with no .wiff files in it (a folder named .wiff and a .wiff2 do not count)
    empty = mkpath(joinpath(root, "empty"))
    mkpath(joinpath(empty, "run.wiff"))
    write(joinpath(empty, "run.wiff2"), "encrypted, unsupported")
    @test_throws ArgumentError convertSciex(empty)
    # the CLI returns 1 rather than throwing, and 1 for no arguments (help)
    @test redirect_stdout(devnull) do
        main_convertSciex(String[])
    end == 1
    @test redirect_stdout(devnull) do
        main_convertSciex([empty])
    end == 1
    @test redirect_stdout(devnull) do
        main_convertSciex(["--bogus"])
    end == 1
    @test redirect_stdout(devnull) do
        main_convertSciex(["--help"])
    end == 0
    @test redirect_stdout(devnull) do
        main_convertSciex(["--zt", "--no-zt", empty])
    end == 1

    # ZT Scan DIA is always asked: y / n / Enter (No); an error with no terminal to ask on
    ask(s) = ask_zt_scan(; input = IOBuffer(s), output = devnull, interactive = true)
    @test ask("\n") == false && ask("n\n") == false && ask("y\n") == true && ask("YES\n") == true
    @test ask("maybe\ny\n") == true
    @test_throws ArgumentError ask_zt_scan(; input = IOBuffer("y\n"), output = devnull, interactive = false)

    wiff = get(ENV, "PIONEER_TEST_WIFF", "")
    if !isempty(wiff)
        out = mktempdir()
        paths = redirect_stdout(devnull) do
            convertSciex(wiff; output_dir = out, zt_scan = false)
        end
        @test length(paths) == 1 && Pioneer.is_scxs_path(only(paths))
        d = Pioneer.loadMassSpecData(only(paths))
        @test d isa Pioneer.ScxsMassSpecData && length(d) > 0
        @test any(==(UInt8(1)), Pioneer.getMsOrders(d)) && any(==(UInt8(2)), Pioneer.getMsOrders(d))
    end
end
