using Test
using DataFrames, Arrow, Tables

using Pioneer: PSMFileReference,
               sort_file_by_keys!, write_arrow_file, transform_and_write!,
               load_dataframe, column_names, has_columns, create_reference,
               _ensure_typed_missing_file_columns!

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
    @test ref.sorted_by == (:score, :target)
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
    @test Symbol.(names(df_only_id)) == [:id]

    # transform_and_write! to new path (no change to original)
    new_ref = transform_and_write!(d->d, ref, f2)
    @test load_dataframe(new_ref) == df_only_id
end

@testset "Arrow record-batch missing Float32 stability" begin
    temp_dir = mktempdir()
    out_path = joinpath(temp_dir, "missing_float32_batches.arrow")

    open(Arrow.Writer, out_path; file=false) do writer
        batch1 = DataFrame(
            id = 1:2,
            run_col = Union{Missing, Float32}[1.0f0, missing]
        )
        Arrow.write(writer, batch1)

        # Match Pioneer write path: add absent run columns as typed Missing Float32.
        batch2 = DataFrame(id = 3:4)
        _ensure_typed_missing_file_columns!(batch2, ["run_col"], Float32)
        Arrow.write(writer, batch2)
    end

    table = Arrow.Table(out_path)
    c = Tables.getcolumn(table, :run_col)

    @test eltype(c) == Union{Missing, Float32}
    @test c[1] == 1.0f0
    @test ismissing(c[2])
    @test ismissing(c[3])
    @test ismissing(c[4])
end

@testset "Arrow record-batch unstack all-missing stability" begin
    temp_dir = mktempdir()
    out_path = joinpath(temp_dir, "unstack_all_missing_batches.arrow")

    long1 = DataFrame(
        precursor_idx = UInt32[1, 2],
        file_name = ["run_col", "run_col"],
        peak_area = Union{Missing, Float32}[1.0f0, missing]
    )
    long2 = DataFrame(
        precursor_idx = UInt32[3, 4],
        file_name = ["run_col", "run_col"],
        peak_area = Union{Missing, Float32}[missing, missing]
    )

    batch1 = unstack(long1, [:precursor_idx], :file_name, :peak_area; combine=sum)
    batch2 = unstack(long2, [:precursor_idx], :file_name, :peak_area; combine=sum)
    @test eltype(batch2.run_col) == Missing

    _ensure_typed_missing_file_columns!(batch1, ["run_col"], Float32)
    _ensure_typed_missing_file_columns!(batch2, ["run_col"], Float32)

    open(Arrow.Writer, out_path; file=false) do writer
        Arrow.write(writer, batch1)
        Arrow.write(writer, batch2)
    end

    table = Arrow.Table(out_path)
    c = Tables.getcolumn(table, :run_col)

    @test eltype(c) == Union{Missing, Float32}
    @test c[1] == 1.0f0
    @test ismissing(c[2])
    @test ismissing(c[3])
    @test ismissing(c[4])
end
