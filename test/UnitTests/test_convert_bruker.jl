# convertBruker input handling: what counts as a Bruker bundle, and the errors for inputs that are not. The
# conversion itself is TimsSlices.convert_run (src/vendor/TimsSlices), tested in test/UnitTests/formats/timsslices and end to end by the timsTOF searches.

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

@testset "TimsSlices locked read (the Windows pread! path)" begin
    # pread! uses this on Windows, where there is no positioned read. Run it on every OS: it is the
    # path that broke the Windows build (lock(f, ::IOStream) has no method; lock(::IOStream) is a no-op).
    mktemp() do path, io
        payload = UInt8.(0:255)
        write(io, payload); flush(io); close(io)
        open(path, "r") do stream
            l = ReentrantLock()
            dst = zeros(UInt8, 8)
            Pioneer.TimsSlices._locked_read!(dst, stream, l, 100, 8)
            @test dst == payload[101:108]
            results = Vector{Vector{UInt8}}(undef, 30)
            Threads.@threads for t in 1:30   # 8t + 16 <= 256 bytes
                buf = zeros(UInt8, 16)
                Pioneer.TimsSlices._locked_read!(buf, stream, l, 8t, 16)
                results[t] = buf
            end
            @test all(results[t] == payload[8t+1:8t+16] for t in 1:30)
        end
    end
end
