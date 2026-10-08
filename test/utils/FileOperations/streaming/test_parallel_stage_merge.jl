# Staging groups of a hierarchical merge are merged concurrently. The result must be identical, row for row and
# including the order of tied rows, to merging the groups one after another (the previous behaviour).

using Test, Arrow, DataFrames, Random
using Pioneer: PSMFileReference, sort_file_by_keys!, stream_sorted_merge, stream_sorted_merge_chunked, file_path

const STAGE_KEYS = (:group, :target, :entrapment_group_id, :precursor_idx)
const STAGE_REV = [false, true, true, true]

# Few distinct key values, so many rows tie across files and tie order is exercised.
_stage_frame(rng, n, file) = DataFrame(
    group = rand(rng, ["P1", "P2", "P3"], n),
    target = rand(rng, Bool, n),
    entrapment_group_id = rand(rng, UInt8(0):UInt8(1), n),
    precursor_idx = UInt32.(rand(rng, 1:4, n)),
    source = fill(UInt16(file), n),          # identifies the input file of each row
    row = collect(UInt32, 1:n),
)

function _stage_inputs(dir, rng; n_files = 9)
    map(1:n_files) do i
        path = joinpath(dir, "in_$(lpad(i, 2, '0')).arrow")
        Arrow.write(path, _stage_frame(rng, rand(rng, 0:300), i))
        ref = PSMFileReference(path)
        sort_file_by_keys!(ref, STAGE_KEYS...; reverse = STAGE_REV)
        ref
    end
end

# The previous algorithm: merge groups of `fanin` one after another, level by level, then merge the rest.
function _sequential_staged(refs, out, fanin)
    level = 0
    while length(refs) > fanin
        level += 1
        refs = [stream_sorted_merge(collect(g), joinpath(dirname(out), "seq_$(level)_$i.arrow"), STAGE_KEYS...;
                                    reverse = STAGE_REV, batch_size = 37, max_fanin = fanin)
                for (i, g) in enumerate(Iterators.partition(refs, fanin))]
    end
    stream_sorted_merge(refs, out, STAGE_KEYS...; reverse = STAGE_REV, batch_size = 37, max_fanin = fanin)
end

@testset "parallel merge staging matches sequential staging exactly (threads=$(Threads.nthreads()))" begin
    rng = MersenneTwister(20261008)
    for fanin in (2, 3)
        mktempdir() do dir
            refs = _stage_inputs(dir, rng)
            expected = DataFrame(Arrow.Table(read(file_path(_sequential_staged(refs, joinpath(dir, "seq.arrow"), fanin)))))

            out = joinpath(dir, "par.arrow")
            stream_sorted_merge(refs, out, STAGE_KEYS...; reverse = STAGE_REV, batch_size = 37, max_fanin = fanin)
            @test isequal(DataFrame(Arrow.Table(read(out))), expected)

            chunks = stream_sorted_merge_chunked(refs, joinpath(dir, "chunks"), :group, STAGE_KEYS...;
                                                 reverse = STAGE_REV, batch_size = 37, max_chunk_bytes = 2_000,
                                                 max_fanin = fanin)
            @test isequal(reduce(vcat, [DataFrame(Arrow.Table(read(file_path(c)))) for c in chunks]), expected)
            @test nrow(expected) > 0
        end
    end
end
