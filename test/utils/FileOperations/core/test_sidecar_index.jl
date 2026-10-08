using Test, Arrow, DataFrames, Pioneer

function _sidecar_reference_signature(ref)
    return (
        ref.file_path, ref.file_exists, ref.row_count, ref.sorted_by,
        ref.schema.columns, [(side.path, side.cols) for side in ref.sidecars],
    )
end

function _check_indexed_sidecar_parity(paths)
    index = Pioneer.index_sidecar_paths(paths)
    for path in paths
        ordinary = Pioneer.PSMFileReference(path)
        indexed = Pioneer.PSMFileReference(path; sidecar_paths=index[path])
        @test _sidecar_reference_signature(indexed) == _sidecar_reference_signature(ordinary)
        if ordinary.file_exists
            @test isequal(Pioneer.load_with_sidecars(indexed), Pioneer.load_with_sidecars(ordinary))
        end
    end
end

@testset "Shared sidecar discovery preserves references and values" begin
    directory = mktempdir()
    try
        paths = String[]
        for (subdir, name, rows) in (
            ("one", "é.測試.arrow", 3),
            ("two", "zero.arrow", 0),
            ("one", "a.arrow", 3),
            ("one", "a.arrow.more.arrow", 3),
            ("two", "a.arrow", 3),
            ("two", "self.sidecar.arrow", 3),
        )
            path = joinpath(mkpath(joinpath(directory, subdir)), name)
            Arrow.write(path, (id=UInt32.(1:rows), base=Float32.(1:rows)))
            Arrow.write(path * ".a.sidecar.arrow",
                (base=fill(99f0, rows), score=fill(1f0, rows)))
            Arrow.write(path * ".b.sidecar.arrow",
                (score=fill(2f0, rows), extra=fill(3f0, rows)))
            Arrow.write(path * ".wrong.sidecar.arrow", (bad=Float32[1, 2, 3, 4],))
            write(path * ".malformed.sidecar.arrow", "not Arrow")
            push!(paths, path)
        end
        plain_path = joinpath(directory, "two", "plain.arrow")
        Arrow.write(plain_path, (id=UInt32[1],))
        push!(paths, plain_path)
        push!(paths, joinpath(directory, "one", "missing.arrow"))
        push!(paths, joinpath(directory, "absent_directory", "missing.arrow"))
        push!(paths, first(paths))
        write(joinpath(directory, "one", "unrelated.txt"), "unrelated")

        _check_indexed_sidecar_parity(paths)
        index = Pioneer.index_sidecar_paths(paths)
        @test index[plain_path] == String[]
        @test index[paths[end - 1]] == String[]
        @test isempty(Pioneer.index_sidecar_paths(String[]))
        @test index[first(paths)] == sort(String[
            first(paths) * ".a.sidecar.arrow", first(paths) * ".b.sidecar.arrow",
            first(paths) * ".wrong.sidecar.arrow", first(paths) * ".malformed.sidecar.arrow",
        ])
        # Main columns win; among sidecars, sorted discovery order wins.
        ref = Pioneer.PSMFileReference(first(paths); sidecar_paths=index[first(paths)])
        values = Pioneer.load_with_sidecars(ref)
        @test values.base == Float32[1, 2, 3]
        @test values.score == fill(1f0, 3)
        @test values.extra == fill(3f0, 3)
        @test !hasproperty(values, :bad)
    finally
        # Existing Arrow constructors may leave mappings pending finalization
        # on Windows. Release those for test cleanup, not inside discovery.
        GC.gc()
        GC.gc()
        rm(directory; recursive=true)
    end
end

@testset "Sidecar snapshots refresh between writes and own no handles" begin
    directory = mktempdir()
    try
        path = joinpath(directory, "snapshot.arrow")
        Arrow.write(path, (id=UInt32[1, 2, 3],))
        paths = String[path]
        original = Pioneer.index_sidecar_paths(paths)
        late = path * ".late.sidecar.arrow"
        Arrow.write(late, (late=Float32[4, 5, 6],))
        @test original[path] == String[]
        @test !Pioneer.has_column_anywhere(
            Pioneer.PSMFileReference(path; sidecar_paths=original[path]), :late)
        _check_indexed_sidecar_parity(paths)
        @test Pioneer.has_column_anywhere(Pioneer.PSMFileReference(path), :late)

        # The index below only stores filenames. Never read the candidate
        # before removal/replacement, so this also works with Windows locks.
        changed = path * ".changed.sidecar.arrow"
        Arrow.write(changed, (old=Float32[1, 2, 3],))
        stale = Pioneer.index_sidecar_paths(paths)
        @test changed in stale[path]
        rm(changed)
        @test !isfile(changed)
        @test !Pioneer.has_column_anywhere(
            Pioneer.PSMFileReference(path; sidecar_paths=stale[path]), :old)
        Arrow.write(changed, (replacement=Float32[7, 8, 9],))
        fresh = Pioneer.index_sidecar_paths(paths)
        @test Pioneer.has_column_anywhere(
            Pioneer.PSMFileReference(path; sidecar_paths=fresh[path]), :replacement)
        _check_indexed_sidecar_parity(paths)
    finally
        GC.gc()
        GC.gc()
        rm(directory; recursive=true)
    end
end
