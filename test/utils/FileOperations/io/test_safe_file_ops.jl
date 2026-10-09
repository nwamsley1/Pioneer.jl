using Test
using DataFrames, Arrow, Tables

using Pioneer: safeRm, writeArrow, _windows_path

@testset "safeRm path handling" begin
    mktempdir() do temp_dir
        nested_dir = joinpath(temp_dir, "path with spaces", "data")
        mkpath(nested_dir)
        file_path = joinpath(nested_dir, "temporary.arrow")
        write(file_path, "temporary")

        # Git Bash and JSON inputs can supply forward-slash paths on Windows.
        forward_slash_path = replace(file_path, "\\" => "/")
        safeRm(forward_slash_path; force=true)
        @test !isfile(file_path)

        # Missing paths are deliberately a no-op.
        safeRm(forward_slash_path; force=true)
        @test !isfile(file_path)
    end
end

const LONG_PATH_PREFIX = "\\\\?\\"   # \\?\

@testset "Windows path normalization" begin
    relative_path = joinpath("relative", "data", "file with spaces.arrow")
    win_path = _windows_path(relative_path)
    @test !occursin("/", win_path)
    @test endswith(win_path, "relative\\data\\file with spaces.arrow")

    if Sys.iswindows()   # abspath leaves drive and UNC paths alone only on Windows
        long_path = "C:\\" * join(fill("d" ^ 50, 6), "\\") * "\\file.arrow"
        @test _windows_path(long_path) == LONG_PATH_PREFIX * long_path
        long_unc = "\\\\server\\share\\" * join(fill("d" ^ 50, 6), "\\") * "\\file.arrow"
        @test _windows_path(long_unc) == LONG_PATH_PREFIX * "UNC\\" * long_unc[3:end]
        @test _windows_path("C:\\short\\file.arrow") == "C:\\short\\file.arrow"
    end
end

if Sys.iswindows()
    @testset "safeRm on Windows: read-only, memory-mapped and long paths" begin
        mktempdir() do temp_dir
            read_only = joinpath(temp_dir, "read only.arrow")
            write(read_only, "temporary")
            chmod(read_only, 0o444)
            safeRm(read_only)
            @test !isfile(read_only)

            # Julia's rm refuses a file this process has mapped; safeRm must not.
            mapped = joinpath(temp_dir, "mapped.arrow")
            Arrow.write(mapped, (x = collect(1:1000),))
            table = Arrow.Table(mapped)
            @test sum(table.x) == 500500
            safeRm(mapped)
            @test !isfile(mapped)
            @test sum(table.x) == 500500    # the mapping itself stays readable

            deep = joinpath(temp_dir, fill("d" ^ 50, 5)...)
            long_file = joinpath(deep, "file.arrow")
            @test length(long_file) >= 260
            mkpath(LONG_PATH_PREFIX * deep)
            write(LONG_PATH_PREFIX * long_file, "temporary")
            safeRm(long_file)
            @test !isfile(LONG_PATH_PREFIX * long_file)

            # A file held with no sharing cannot be deleted, renamed or collected away:
            # safeRm warns and returns, as the cmd.exe del path it replaced did.
            locked = joinpath(temp_dir, "locked.arrow")
            write(locked, "temporary")
            handle = ccall((:CreateFileW, "kernel32"), stdcall, Ptr{Cvoid},
                (Cwstring, UInt32, UInt32, Ptr{Cvoid}, UInt32, UInt32, Ptr{Cvoid}),
                locked, 0x80000000, 0, C_NULL, 3, 0x80, C_NULL)   # GENERIC_READ, no sharing
            try
                @test (safeRm(locked); true)
                @test isfile(locked)
            finally
                ccall((:CloseHandle, "kernel32"), stdcall, Cint, (Ptr{Cvoid},), handle)
            end
            safeRm(locked)
            @test !isfile(locked)
        end
    end
end

@testset "writeArrow replaces existing files through safeRm" begin
    mktempdir() do temp_dir
        output_path = joinpath(temp_dir, "replace me.arrow")
        writeArrow(output_path, DataFrame(value=Int32[1]))
        writeArrow(output_path, DataFrame(value=Int32[2]))

        result = Arrow.Table(output_path)
        @test collect(Tables.getcolumn(result, :value)) == Int32[2]
    end
end
