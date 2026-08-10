using Test
using DataFrames, Arrow, Tables

using Pioneer: safeRm, writeArrow, _windows_delete_command

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

@testset "Windows delete command normalization" begin
    relative_path = joinpath("relative", "data", "file with spaces.arrow")

    regular_cmd = _windows_delete_command(relative_path)
    force_cmd = _windows_delete_command(relative_path; force=true)

    @test regular_cmd.exec[1:5] == ["cmd.exe", "/d", "/c", "del", "/q"]
    @test "/f" ∉ regular_cmd.exec
    @test "/f" ∈ force_cmd.exec
    @test !occursin("/", regular_cmd.exec[end])
    @test endswith(regular_cmd.exec[end], "relative\\data\\file with spaces.arrow")
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
