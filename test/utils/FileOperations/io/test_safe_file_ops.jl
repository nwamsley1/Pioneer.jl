using Test
using DataFrames, Arrow, Tables

using Pioneer: safeRm, writeArrow, _windows_delete_command

struct FailingArrowColumn <: AbstractVector{Int32} end
Base.size(::FailingArrowColumn) = (1,)
Base.getindex(::FailingArrowColumn, ::Int) = error("injected Arrow serialization failure")

@testset "writeArrow failure preserves existing output" begin
    mktempdir() do dir
        temp_dir = joinpath(dir, "staging")
        mkdir(temp_dir)
        output_path = joinpath(dir, "output.arrow")
        writeArrow(output_path, DataFrame(value=Int32[7]); temp_dir)
        original = read(output_path)
        df = DataFrame(value=FailingArrowColumn(); copycols=false)
        @test_throws Exception writeArrow(output_path, df; temp_dir)
        @test read(output_path) == original
        @test isempty(readdir(temp_dir))
    end
end

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

@testset "writeArrow removes temporary files after replacement failure" begin
    mktempdir() do dir
        temp_dir = joinpath(dir, "staging")
        mkdir(temp_dir)
        output_path = joinpath(dir, "missing", "output.arrow")
        @test_throws Exception writeArrow(output_path, DataFrame(value=Int32[1]); temp_dir)
        @test isempty(readdir(temp_dir))
        @test !isfile(output_path)
    end
end

@testset "Windows delete command normalization" begin
    relative_path = joinpath("relative", "data", "file with spaces.arrow")

    delete_cmd = _windows_delete_command(relative_path)

    @test delete_cmd.exec[1:6] == ["cmd.exe", "/d", "/c", "del", "/f", "/q"]
    @test !occursin("/", delete_cmd.exec[end])
    @test endswith(delete_cmd.exec[end], "relative\\data\\file with spaces.arrow")
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
