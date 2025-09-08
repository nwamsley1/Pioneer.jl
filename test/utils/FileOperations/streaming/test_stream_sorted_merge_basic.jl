using Test
using DataFrames, Arrow

using Pioneer: PSMFileReference, stream_sorted_merge, load_dataframe, sort_file_by_keys!

@testset "stream_sorted_merge basic" begin
    tmp = mktempdir()
    f1 = joinpath(tmp, "left.arrow")
    f2 = joinpath(tmp, "right.arrow")
    out = joinpath(tmp, "merged.arrow")

    # Prepare two sorted files by :id
    df1 = DataFrame(id = [1,3,5], score = Float32[0.9, 0.7, 0.8])
    df2 = DataFrame(id = [2,4,6], score = Float32[0.95, 0.6, 0.85])
    Arrow.write(f1, df1)
    Arrow.write(f2, df2)

    ref1 = PSMFileReference(f1)
    ref2 = PSMFileReference(f2)
    # Ensure sorted_by metadata is set; sort if needed
    sort_file_by_keys!(ref1, :id)
    sort_file_by_keys!(ref2, :id)

    merged_ref = stream_sorted_merge([ref1, ref2], out, :id)
    dfm = load_dataframe(merged_ref)

    @test dfm.id == collect(1:6)
    @test length(dfm.score) == 6
end
