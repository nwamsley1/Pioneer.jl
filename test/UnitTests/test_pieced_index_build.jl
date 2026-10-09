# write_fragment_indexes from a selection: an index over the piece limit is written as pieces whose partitions,
# in order, are the single index's (built from a spilled, memory-mapped selection).

using Test, Random, JSON
using Pioneer

@testset "pieced fragment index build" begin
    rng = Xoshiro(11)
    n = 3000
    prec_mzs = Float32.(400 .+ 600 .* rand(rng, n)); prec_irts = Float32.(100 .* rand(rng, n))
    counts = UInt8.(rand(rng, 0:8, n))
    starts = Vector{Int}(undef, n); frags = Pioneer.SimpleFrag{Float32}[]
    for pid in 1:n
        starts[pid] = length(frags) + 1
        for r in 1:counts[pid]
            push!(frags, Pioneer.SimpleFrag{Float32}(Float32(150 + 1500 * rand(rng)), UInt32(pid), prec_mzs[pid],
                                                     prec_irts[pid], 0x00, UInt8(1) << UInt8(r - 1)))
        end
    end
    sel = Pioneer.IndexFragSelection(frags, starts, counts, prec_mzs)
    single, pieced = mktempdir(), mktempdir()
    kw = (frag_bin_tol_ppm = 10.0f0, frag_bin_tol_mda = 2.0f0, rt_bin_tol = 3.0f0)
    Pioneer.write_fragment_indexes(single, sel, (5.0f0,), "UInt32"; kw...)
    Pioneer.write_fragment_indexes(pieced, sel, (5.0f0,), "UInt32"; kw..., max_piece_bytes = 20_000)
    entry = JSON.parsefile(joinpath(pieced, Pioneer.FRAGMENT_INDEX_DESCRIPTOR))["indexes"][1]
    @test entry["pieced"] == true && entry["main"] == "partitioned_fragment_index_pieces"
    @test !haskey(JSON.parsefile(joinpath(single, Pioneer.FRAGMENT_INDEX_DESCRIPTOR))["indexes"][1], "pieced")
    @test isempty(filter(f -> endswith(f, ".tmp"), readdir(pieced)))          # spilled selection removed
    same_partition(a, b) = all(fieldnames(typeof(a))) do f
        x = getfield(a, f); y = getfield(b, f)
        f === :fragments ? [(g.local_id, g.score) for g in x] == [(g.local_id, g.score) for g in y] :
        x isa Pioneer.SoAFragBins ? all(g -> getfield(x, g) == getfield(y, g), fieldnames(typeof(x))) : x == y
    end
    # the raw index file reads back as the index built in memory
    built = Pioneer.build_partitioned_index_from_selection(sel; partition_width = 5.0f0, kw..., id_type = UInt32)
    stored = Pioneer.load_fragment_index(joinpath(single, "partitioned_fragment_index.bin"))
    @test stored.partition_bounds == built.partition_bounds && stored.n_partitions == built.n_partitions
    @test all(k -> same_partition(stored.partitions[k], built.partitions[k]), 1:built.n_partitions)
    for (file, dir) in (("partitioned_fragment_index.bin", "partitioned_fragment_index_pieces"),
                        ("presearch_partitioned_fragment_index.bin", "presearch_partitioned_fragment_index_pieces"))
        whole = Pioneer.load_fragment_index(joinpath(single, file))
        pfi = Pioneer.load_fragment_index(joinpath(pieced, dir))
        @test pfi isa Pioneer.PiecedFragmentIndex
        @test length(pfi.pieces) > 1
        parts = []; bounds = Tuple{Float32, Float32}[]
        for piece in pfi.pieces
            idx = Pioneer.read_index_piece(joinpath(pfi.dir, piece.file))
            append!(parts, idx.partitions); append!(bounds, idx.partition_bounds)
        end
        @test bounds == whole.partition_bounds
        @test length(parts) == whole.n_partitions
        @test all(1:length(parts)) do k
            a, b = parts[k], whole.partitions[k]
            all(f -> (x = getfield(a, f); y = getfield(b, f);
                      x isa Pioneer.SoAFragBins ? all(g -> getfield(x, g) == getfield(y, g), fieldnames(typeof(x))) : x == y),
                fieldnames(typeof(a)))
        end
    end
end
