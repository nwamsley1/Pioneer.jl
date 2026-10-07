# EXPERIMENT check: a pieced index (build_index_pieces) has exactly the single index's partitions, and pieces
# round-trip through write/read_index_piece. usage: julia test_pieces.jl LIB_DIR OUT_DIR MAX_PIECE_BYTES ID_TYPE
using Pioneer
const P = Pioneer
include(joinpath(@__DIR__, "build_synth.jl"))
lib_dir, out, maxb, idreq = ARGS[1], ARGS[2], parse(Int, ARGS[3]), ARGS[4]
sel = P.select_index_fragments(load_rank_ordered_lib(lib_dir))
I = idreq == "UInt32" ? UInt32 : UInt16
single = P.build_partitioned_index_from_selection(sel; partition_width = 5.0f0, rt_bin_tol = 3.0f0, id_type = I)
t = @timed P.build_index_pieces(sel, out; partition_width = 5.0f0, rt_bin_tol = 3.0f0, id_type_request = idreq,
                                max_piece_bytes = maxb)
pfi = P.load_pieced_index(out)
println("pieces: ", length(pfi.pieces), "  sizes GB: ", [round(p.bytes / 1e9, digits = 3) for p in pfi.pieces],
        "  build ", round(t.time, digits = 1), " s, single index ", round(P.index_bytes(single) / 1e9, digits = 3), " GB")
parts = reduce(vcat, [P.read_index_piece(joinpath(out, p.file)).partitions for p in pfi.pieces])
bounds = reduce(vcat, [P.read_index_piece(joinpath(out, p.file)).partition_bounds for p in pfi.pieces])
same(a, b) = a.fragments == b.fragments && a.local_to_global == b.local_to_global && a.n_local_precs == b.n_local_precs &&
    a.fragment_bins.lows == b.fragment_bins.lows && a.fragment_bins.highs == b.fragment_bins.highs &&
    a.fragment_bins.first_bins == b.fragment_bins.first_bins && a.fragment_bins.last_bins == b.fragment_bins.last_bins &&
    a.skip_hints == b.skip_hints &&
    [(r.lb, r.ub, r.first_bin, r.last_bin) for r in a.rt_bins] == [(r.lb, r.ub, r.first_bin, r.last_bin) for r in b.rt_bins]
ok = length(parts) == single.n_partitions && bounds == single.partition_bounds && all(same(parts[k], single.partitions[k]) for k in eachindex(parts))
println("PIECES IDENTICAL TO SINGLE INDEX: ", ok)
ok || exit(1)
