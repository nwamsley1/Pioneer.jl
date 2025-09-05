using Test
using DataFrames, Arrow, Tables

using Pioneer: PSMFileReference,
               sort_file_by_keys!, write_arrow_file, transform_and_write!,
               load_dataframe, column_names, has_columns, create_reference

@testset "ArrowOperations basics" begin
    temp_dir = mktempdir()
    f1 = joinpath(temp_dir, "a.arrow")
    f2 = joinpath(temp_dir, "b.arrow")

    df = DataFrame(id=1:5,
                   score=Float32[0.5, 0.9, 0.7, 0.6, 0.8],
                   target=Bool[true, false, true, false, true])
    Arrow.write(f1, df)

    # Create reference and basic helpers
    ref = PSMFileReference(f1)
    @test has_columns(ref, :id, :score, :target)
    cols = column_names(ref)
    @test :id ∈ cols && :score ∈ cols && :target ∈ cols

    # sort_file_by_keys! single file
    sort_file_by_keys!(ref, :score, :target; reverse=[true, true])
    df_sorted = load_dataframe(ref)
    @test issorted(df_sorted.score; rev=true)
    # ties on score should be sorted by target desc
    if any(diff(df_sorted.score) .== 0)
        t_inds = findall(x->x==df_sorted.score[1], df_sorted.score)
        @test issorted(df_sorted.target[t_inds]; rev=true)
    end

    # write_arrow_file updates metadata
    df2 = DataFrame(id=6:8, score=Float32[0.1,0.2,0.3], target=Bool[true,false,true])
    write_arrow_file(ref, df2)
    df_reload = load_dataframe(ref)
    @test nrow(df_reload) == 3
    @test ref.sorted_by == ()  # reset after write

    # transform_and_write! (in-place)
    transform_and_write!(ref) do d
        select(d, :id)
    end
    df_only_id = load_dataframe(ref)
    @test names(df_only_id) == [:id]

    # transform_and_write! to new path (no change to original)
    new_ref = transform_and_write!(d->d, ref, f2)
    @test load_dataframe(new_ref) == df_only_id
end

