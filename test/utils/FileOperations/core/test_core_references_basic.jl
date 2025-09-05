using Test
using DataFrames, Arrow

using Pioneer: PSMFileReference, create_reference, FileSchema, has_column

@testset "FileReferences basic" begin
    tmp = mktempdir()
    missing_path = joinpath(tmp, "missing.arrow")
    ref_missing = PSMFileReference(missing_path)
    @test ref_missing.file_exists == false
    @test ref_missing.row_count == 0

    # Write a small file
    path = joinpath(tmp, "psms.arrow")
    df = DataFrame(id=1:3, score=Float32[0.9,0.8,0.7])
    Arrow.write(path, df)

    ref = create_reference(path, PSMFileReference)
    @test ref.file_exists == true
    @test ref.row_count == 3
    @test has_column(ref.schema, :id)
    @test has_column(ref.schema, :score)
end

