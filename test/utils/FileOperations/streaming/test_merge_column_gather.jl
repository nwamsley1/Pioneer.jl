# The k-way merge gathers each output column from typed source columns (one lookup per column
# per batch). These tests pin the gathered values to an in-memory merge for every column type the
# merge writes, and keep a mismatched input column an error.

using Test, Arrow, DataFrames, Random
using Pioneer: PSMFileReference, sort_file_by_keys!, stream_sorted_merge, stream_sorted_merge_chunked

const GATHER_KEYS = (:group, :target, :entrapment_group_id, :precursor_idx)
const GATHER_REV = [false, true, true, true]

function _gather_frame(rng, n; score_type = Float32)
    DataFrame(
        group = rand(rng, ["P1", "P2;P3", "P4", "P10"], n),                 # String sort and group key
        target = rand(rng, Bool, n),
        entrapment_group_id = rand(rng, UInt8(0):UInt8(2), n),
        precursor_idx = UInt32.(rand(rng, 1:50, n)),
        score = rand(rng, score_type, n),                                     # Real
        qval = allowmissing(rand(rng, Float32, n)),                           # Union{Missing, Real}
        note = Union{Missing, String}[isodd(i) ? "n$i" : missing for i in 1:n],  # Union{Missing, String}
        name = ["row$i" for i in 1:n],                                        # String
        flag = rand(rng, UInt8, n),
        mz_range = [(rand(rng, Float32), rand(rng, Float32)) for _ in 1:n],   # Tuple
        weights = [rand(rng, Float32, rand(rng, 0:3)) for _ in 1:n],          # list column
    )
end

function _write_sorted(dir, frames; dictencode = fill(false, length(frames)))
    map(enumerate(frames)) do (i, df)
        path = joinpath(dir, "in_$(lpad(i, 3, '0')).arrow")
        Arrow.write(path, df; dictencode = dictencode[i])
        ref = PSMFileReference(path)
        sort_file_by_keys!(ref, GATHER_KEYS...; reverse = GATHER_REV)
        ref
    end
end

# Same rows as the in-memory merge (ties may be ordered differently), in merge-key order.
function _check_against_memory(got::DataFrame, frames)
    expected = reduce(vcat, frames)
    @test nrow(got) == nrow(expected)
    @test names(got) == names(expected)
    @test issorted(got, collect(GATHER_KEYS); rev = GATHER_REV)
    canonical(df) = sort(df, names(df); by = x -> x isa AbstractVector ? collect(x) : x)
    @test isequal(canonical(got), canonical(expected))
end

@testset "k-way merge gathers every column type exactly" begin
    rng = MersenneTwister(20261005)
    mktempdir() do dir
        frames = [_gather_frame(rng, n) for n in (0, 1, 7, 2_003, 517)]   # empty, tiny, uneven
        refs = _write_sorted(dir, frames)
        out = joinpath(dir, "merged.arrow")
        stream_sorted_merge(refs, out, GATHER_KEYS...; batch_size = 97, reverse = GATHER_REV)
        _check_against_memory(DataFrame(Arrow.Table(read(out))), frames)
    end
end

@testset "inputs with different Arrow encodings of the same column" begin
    rng = MersenneTwister(7)
    mktempdir() do dir
        frames = [_gather_frame(rng, 400) for _ in 1:3]
        refs = _write_sorted(dir, frames; dictencode = [false, true, false])   # untyped gather path
        out = joinpath(dir, "merged.arrow")
        stream_sorted_merge(refs, out, GATHER_KEYS...; batch_size = 128, reverse = GATHER_REV)
        _check_against_memory(DataFrame(Arrow.Table(read(out))), frames)
    end
end

@testset "nullable input column into a non-nullable batch column" begin
    rng = MersenneTwister(11)
    mktempdir() do dir
        a = _gather_frame(rng, 50)
        b = _gather_frame(rng, 50); b.score = allowmissing(b.score)        # nullable type, no missing values
        refs = _write_sorted(dir, [a, b])
        out = joinpath(dir, "merged.arrow")
        stream_sorted_merge(refs, out, GATHER_KEYS...; batch_size = 16, reverse = GATHER_REV)
        _check_against_memory(DataFrame(Arrow.Table(read(out))), [a, disallowmissing(b, :score)])

        c = _gather_frame(rng, 50); c.score = allowmissing(c.score); c.score[3] = missing
        refs = _write_sorted(joinpath(dir, ""), [a, c])
        @test_throws Exception stream_sorted_merge(refs, joinpath(dir, "bad.arrow"), GATHER_KEYS...;
                                                   batch_size = 16, reverse = GATHER_REV)
    end
end

@testset "mismatched input column type is an error" begin
    rng = MersenneTwister(13)
    mktempdir() do dir
        refs = _write_sorted(dir, [_gather_frame(rng, 40), _gather_frame(rng, 40; score_type = Float64)])
        @test_throws ArgumentError stream_sorted_merge(refs, joinpath(dir, "out.arrow"), GATHER_KEYS...;
                                                       batch_size = 16, reverse = GATHER_REV)
    end
end

@testset "chunked merge keeps groups whole and gathers exactly, through staging" begin
    rng = MersenneTwister(17)
    mktempdir() do dir
        frames = [_gather_frame(rng, rand(rng, 0:60)) for _ in 1:70]       # > 64 inputs: staged
        refs = _write_sorted(dir, frames)
        chunk_dir = joinpath(dir, "chunks"); mkpath(chunk_dir)
        chunks = stream_sorted_merge_chunked(refs, chunk_dir, :group, GATHER_KEYS...;
                                             batch_size = 50, reverse = GATHER_REV, max_chunk_bytes = 2_000)
        parts = [DataFrame(Arrow.Table(read(Pioneer.file_path(c)))) for c in chunks]
        @test length(parts) > 1
        groups = [Set(p.group) for p in parts if nrow(p) > 0]
        @test all(isempty(intersect(groups[i], groups[j])) for i in eachindex(groups) for j in eachindex(groups) if i < j)
        _check_against_memory(reduce(vcat, parts), frames)
    end
end
