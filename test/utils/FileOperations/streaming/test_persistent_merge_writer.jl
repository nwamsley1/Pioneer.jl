using Test
using Arrow, DataFrames, Tables
using Pioneer: PSMFileReference, stream_sorted_merge, stream_sorted_merge_chunked,
               mark_sorted!, row_count, file_path, is_sorted_by,
               _with_merge_output, _write_batch_typed, _open_merge_output,
               _close_merge_output!, _publish_merge_output!, _cleanup_merge_output!

# Read bytes instead of mapping the output, so handle tests cannot be confused
# by a reader owned by this test. The merge's input mappings are a separate issue.
_read_merge_result(path) = DataFrame(Arrow.Table(read(path)))

function _with_merge_test_dir(f)
    mktempdir() do dir
        try
            f(dir)
        finally
            # The production merge reads Arrow.Table inputs. After f returns,
            # their unreachable mappings may need collection before Windows can
            # remove the fixture directory. Do not use GC to pass writer tests.
            if Sys.iswindows()
                lock(Pioneer._WINDOWS_DELETE_GC_LOCK) do
                    GC.gc(true)
                end
            end
        end
    end
end

function _merge_fixture_ref(path, df, keys...)
    Arrow.write(path, df)
    ref = PSMFileReference(path; sidecar_paths=String[])
    mark_sorted!(ref, keys...)
    return ref
end

@testset "persistent writer snapshots and immediate handle release" begin
    _with_merge_test_dir() do dir
        path = joinpath(dir, "résult with spaces.arrow")
        batch = DataFrame(id=zeros(Int64, 4), value=zeros(Float32, 4))
        session = Ref{Any}()
        _with_merge_output(path) do output
            session[] = output
            for i in 1:256
                batch.id .= (4i-3):(4i)
                batch.value .= Float32(i)
                _write_batch_typed(output, batch, 4)
                # Exercise buffer reuse immediately, including the last batch.
                batch.id .= -1
                batch.value .= -1
            end
        end
        @test !isopen(session[].io)
        @test session[].published
        @test !ispath(session[].temp_path)
        result = _read_merge_result(path)
        @test result.id == collect(Int64, 1:1024)
        @test result.value == repeat(Float32.(1:256); inner=4)

        # Use ordinary rename/remove, without retries or GC: the writer must
        # release its handle even while the session itself remains reachable.
        renamed = joinpath(dir, "renamed.arrow")
        mv(path, renamed)
        rm(renamed)
        @test !ispath(renamed)
    end
end

@testset "failure closes output and preserves destination" begin
    _with_merge_test_dir() do dir
        path = joinpath(dir, "existing.arrow")
        Arrow.write(path, DataFrame(id=Int64[99]))
        before = read(path)
        session = Ref{Any}()
        failure = ErrorException("injected merge failure")
        caught = try
            _with_merge_output(path) do output
                session[] = output
                _write_batch_typed(output, DataFrame(id=Int64[1, 2]), 2)
                throw(failure)
            end
        catch error
            error
        end
        @test caught === failure
        @test !isopen(session[].io)
        @test !ispath(session[].temp_path)
        @test read(path) == before
        @test readdir(dir) == ["existing.arrow"]
    end
end

@testset "publication failure retains closed recoverable stream" begin
    _with_merge_test_dir() do dir
        # A directory is not a replaceable Arrow destination on either OS.
        destination = joinpath(dir, "blocked.arrow")
        mkdir(destination)
        write(joinpath(destination, "keep.txt"), "keep")
        session = Ref{Any}()
        @test_logs (:warn, r"completed temporary output retained") begin
            @test_throws Exception _with_merge_output(destination) do output
                session[] = output
                _write_batch_typed(output, DataFrame(id=Int64[7]), 1)
            end
        end
        @test !isopen(session[].io)
        @test session[].complete
        @test !session[].published
        @test isfile(joinpath(destination, "keep.txt"))
        @test _read_merge_result(session[].temp_path).id == [7]
        rm(session[].temp_path)
    end
end

@testset "writer finalization failure still closes IO" begin
    _with_merge_test_dir() do dir
        output = _open_merge_output(joinpath(dir, "bad close.arrow"))
        _write_batch_typed(output, DataFrame(id=Int64[1]), 1)
        # Inject an IO failure before Arrow finalizes its stream. Its consumer
        # may also encounter this failure; the underlying file must still close.
        close(output.io)
        @test_throws Exception _close_merge_output!(output)
        @test !isopen(output.io)
        @test !output.complete
        @test_throws Exception _publish_merge_output!(output)
        _cleanup_merge_output!(output)
        @test !ispath(output.temp_path)
    end
end

@testset "consumer IO failure wakes blocked producers" begin
    _with_merge_test_dir() do dir
        output = _open_merge_output(joinpath(dir, "consumer failure.arrow"))
        close(output.io)
        @test_throws Exception begin
            # One message may be accepted before its consumer encounters the
            # closed IO. A subsequent submission must throw rather than block.
            _write_batch_typed(output, DataFrame(id=Int64[1]), 1)
            _write_batch_typed(output, DataFrame(id=Int64[2]), 1)
        end
        _cleanup_merge_output!(output)
        @test !isopen(output.io)
        @test !ispath(output.temp_path)
        @test !ispath(output.path)
    end
end

@testset "consumer failure on the last batch is not published" begin
    _with_merge_test_dir() do dir
        destination = joinpath(dir, "last batch failure.arrow")
        output = _open_merge_output(destination)
        _write_batch_typed(output, DataFrame(id=Int64[1, 2]), 2)
        # The final submission is accepted by the unbuffered channel before the
        # consumer writes it; make that write fail and finalize straight after.
        close(output.io)
        try
            _write_batch_typed(output, DataFrame(id=Int64[3]), 1)
        catch
        end
        timedwait(() -> istaskdone(output.writer.task), 5.0)
        @test istaskfailed(output.writer.task)
        @test_throws Exception _close_merge_output!(output)
        @test !output.complete
        @test_throws Exception _publish_merge_output!(output)
        _cleanup_merge_output!(output)
        @test !ispath(destination)
        @test !ispath(output.temp_path)
    end
end

if Sys.iswindows()
    @testset "external Windows lock leaves old destination intact" begin
        _with_merge_test_dir() do dir
            path = joinpath(dir, "locked output.arrow")
            Arrow.write(path, DataFrame(id=Int64[99]))
            before = read(path)
            wide_path = [transcode(UInt16, abspath(path)); UInt16(0)]
            handle = ccall((:CreateFileW, "kernel32"), Ptr{Cvoid},
                           (Ptr{UInt16}, UInt32, UInt32, Ptr{Cvoid}, UInt32, UInt32, Ptr{Cvoid}),
                           wide_path, 0x80000000, 0, C_NULL, 3, 0x80, C_NULL)
            @assert handle != Ptr{Cvoid}(typemax(UInt)) "Could not create Windows lock fixture"
            session = Ref{Any}()
            try
                @test_logs (:warn, r"completed temporary output retained") begin
                    @test_throws Exception _with_merge_output(path) do output
                        session[] = output
                        _write_batch_typed(output, DataFrame(id=Int64[7]), 1)
                    end
                end
                @test !isopen(session[].io)
                @test !session[].published
                @test session[].complete
            finally
                ccall((:CloseHandle, "kernel32"), Int32, (Ptr{Cvoid},), handle)
            end
            @test read(path) == before
            @test _read_merge_result(session[].temp_path).id == [7]
            rm(session[].temp_path)
        end
    end
end

@testset "persistent merge values, schema, partial and exact batches" begin
    _with_merge_test_dir() do dir
        expected = DataFrame(
            id=collect(Int64, 1:33),
            score=Union{Missing,Float32}[i % 5 == 0 ? missing : Float32(i / 10) for i in 1:33],
            name=Union{Missing,String}[i % 7 == 0 ? missing : "peptide_$i" for i in 1:33],
            target=[isodd(i) for i in 1:33],
            fragments=[Float32[i, i + 0.5] for i in 1:33],
        )
        refs = [_merge_fixture_ref(joinpath(dir, "source_$i.arrow"), expected[i:3:end, :], :id) for i in 1:3]
        for batch_size in (1, 3, 4, 33, 50)
            path = joinpath(dir, "merged_$batch_size.arrow")
            ref = stream_sorted_merge(refs, path, :id; batch_size)
            result = _read_merge_result(path)
            @test isequal(result, expected)
            @test eltype(result.id) == Int64
            @test eltype(result.score) == Union{Missing,Float32}
            @test eltype(result.name) == Union{Missing,String}
            @test row_count(ref) == 33
            @test is_sorted_by(ref, :id)
            @test length(collect(Arrow.Stream(read(path)))) == cld(33, batch_size)
        end
        # Replacing an existing destination cannot append or duplicate rows.
        path = joinpath(dir, "replaced.arrow")
        for _ in 1:3
            stream_sorted_merge(refs, path, :id; batch_size=2)
            @test isequal(_read_merge_result(path), expected)
        end
    end
end

@testset "multi-key direction and ties survive small batches" begin
    _with_merge_test_dir() do dir
        left = DataFrame(group=["A", "A", "B"], score=Int32[4, 2, 3], id=Int64[1, 2, 5])
        right = DataFrame(group=["A", "A", "B"], score=Int32[4, 1, 2], id=Int64[3, 4, 6])
        refs = [_merge_fixture_ref(joinpath(dir, "mixed_$i.arrow"), df, :group, :score)
                for (i, df) in enumerate((left, right))]
        output = joinpath(dir, "mixed.arrow")
        stream_sorted_merge(refs, output, :group, :score; reverse=[false, true], batch_size=1)
        @test _read_merge_result(output).id == [1, 3, 2, 4, 5, 6]

        descending = sort(vcat(left, right), :id; rev=true)
        ref = _merge_fixture_ref(joinpath(dir, "descending.arrow"), descending, :id)
        stream_sorted_merge([ref], output, :id; reverse=true, batch_size=2)
        @test isequal(_read_merge_result(output), descending)
    end
end

@testset "empty sources produce readable schema-preserving outputs" begin
    _with_merge_test_dir() do dir
        empty = DataFrame(id=Int64[], group=String[], score=Union{Missing,Float32}[])
        ref = _merge_fixture_ref(joinpath(dir, "empty.arrow"), empty, :id)
        output = joinpath(dir, "empty merged.arrow")
        result_ref = stream_sorted_merge([ref], output, :id; batch_size=2)
        @test isequal(_read_merge_result(output), empty)
        @test row_count(result_ref) == 0
        chunks = stream_sorted_merge_chunked([ref], joinpath(dir, "empty chunks"), :group, :id; batch_size=2)
        @test length(chunks) == 1
        @test isequal(_read_merge_result(file_path(only(chunks))), empty)

        full = DataFrame(id=Int64[1, 2], group=["A", "B"], score=Union{Missing,Float32}[1, missing])
        full_ref = _merge_fixture_ref(joinpath(dir, "full.arrow"), full, :id)
        stream_sorted_merge([ref, full_ref], output, :id; batch_size=1)
        @test isequal(_read_merge_result(output), full)
    end
end

@testset "chunk rotation, missing and oversized groups, hierarchical staging" begin
    _with_merge_test_dir() do dir
        # Group A is larger than the target and spans many batches. All groups,
        # including missing, must remain whole despite rotating writers.
        expected = DataFrame(
            group=Union{Missing,String}[fill("A", 20); fill("B", 8); fill(missing, 4)],
            id=collect(Int64, 1:32),
            score=Float32.(1:32),
        )
        refs = [_merge_fixture_ref(joinpath(dir, "group_$i.arrow"), expected[i:8:end, :], :group, :id) for i in 1:8]
        output_dir = joinpath(dir, "chunk outputs")
        for _ in 1:2
            chunks = stream_sorted_merge_chunked(refs, output_dir, :group, :group, :id;
                                                batch_size=3, max_chunk_bytes=1, max_fanin=2)
            @test length(chunks) == 3
            results = [_read_merge_result(file_path(ref)) for ref in chunks]
            @test isequal(vcat(results...), expected)
            @test nrow.(results) == [20, 8, 4]
            @test all(ref -> is_sorted_by(ref, :group, :id), chunks)
            @test all(df -> length(unique(df.group)) == 1, results)
            @test all(path -> endswith(path, ".arrow"), readdir(output_dir))
        end
        output = joinpath(dir, "hierarchical.arrow")
        stream_sorted_merge(refs, output, :group, :id; batch_size=2, max_fanin=2)
        @test isequal(_read_merge_result(output), expected)
    end
end

@testset "merge rejects aliases and invalid bounds before writing" begin
    _with_merge_test_dir() do dir
        source = joinpath(dir, "source.arrow")
        ref = _merge_fixture_ref(source, DataFrame(id=Int64[1, 2]), :id)
        before = read(source)
        @test_throws ArgumentError stream_sorted_merge([ref], source, :id; batch_size=1)
        @test_throws ArgumentError stream_sorted_merge([ref], joinpath(dir, ".", "source.arrow"), :id)
        if Sys.iswindows()
            @test_throws ArgumentError stream_sorted_merge([ref], uppercase(source), :id)
        else
            alias = joinpath(dir, "alias.arrow")
            symlink(source, alias)
            @test_throws ArgumentError stream_sorted_merge([ref], alias, :id)
        end
        hard_alias = joinpath(dir, "hard alias.arrow")
        hardlink(source, hard_alias)
        @test_throws ArgumentError stream_sorted_merge([ref], hard_alias, :id)
        @test read(source) == before
        output = joinpath(dir, "invalid.arrow")
        @test_throws ArgumentError stream_sorted_merge([ref], output, :id; batch_size=0)
        @test_throws ArgumentError stream_sorted_merge([ref], output, :id; max_fanin=1)
        @test_throws ArgumentError stream_sorted_merge_chunked([ref], dir, :id, :id; max_chunk_bytes=0)
        @test !ispath(output)

        chunk_path = joinpath(dir, "chunk_0001.arrow")
        chunk_ref = _merge_fixture_ref(chunk_path, DataFrame(id=Int64[1]), :id)
        @test_throws ArgumentError stream_sorted_merge_chunked([chunk_ref], dir, :id, :id)
        @test _read_merge_result(chunk_path).id == [1]
    end
end
