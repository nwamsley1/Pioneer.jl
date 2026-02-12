using Test
using DataFrames, Arrow

using Pioneer: PSMFileReference, stream_sorted_merge, stream_sorted_merge_chunked,
               load_dataframe, sort_file_by_keys!

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

@testset "hierarchical merge with many files" begin
    tmp = mktempdir()
    n_files = 20  # triggers multiple stages with max_fanin=4
    refs = PSMFileReference[]

    for i in 1:n_files
        path = joinpath(tmp, "file_$i.arrow")
        df = DataFrame(id = Int64[i, i + n_files], score = Float32[rand(), rand()])
        Arrow.write(path, df)
        ref = PSMFileReference(path)
        sort_file_by_keys!(ref, :id)
        push!(refs, ref)
    end

    out = joinpath(tmp, "merged.arrow")
    merged_ref = stream_sorted_merge(refs, out, :id; max_fanin=4)
    dfm = load_dataframe(merged_ref)

    @test dfm.id == collect(1:2*n_files)
    @test length(dfm.score) == 2 * n_files
end

@testset "hierarchical chunked merge with many files" begin
    tmp = mktempdir()
    n_files = 20
    refs = PSMFileReference[]

    # Each file has rows for two groups, sorted by (group, id)
    for i in 1:n_files
        path = joinpath(tmp, "file_$i.arrow")
        df = DataFrame(
            group = String["A", "B"],
            id    = Int64[i, i],
            score = Float32[rand(), rand()]
        )
        Arrow.write(path, df)
        ref = PSMFileReference(path)
        sort_file_by_keys!(ref, :group, :id)
        push!(refs, ref)
    end

    out_dir = joinpath(tmp, "chunks")
    chunk_refs = stream_sorted_merge_chunked(
        refs, out_dir, :group, :group, :id; max_fanin=4
    )

    # Collect all rows across chunks
    all_dfs = [load_dataframe(cr) for cr in chunk_refs]
    combined = vcat(all_dfs...)

    @test nrow(combined) == 2 * n_files
    @test sort(unique(combined.group)) == ["A", "B"]
    # Within each group, ids should be sorted
    for g in ["A", "B"]
        gdf = filter(r -> r.group == g, combined)
        @test gdf.id == sort(gdf.id)
    end
end
