# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# ── Pieced fragment index ─────────────────────────────────────────────────────
#
# A fragment index too large to hold in memory at once is stored as pieces: each piece is an ordinary partitioned
# index over a consecutive run of the initial precursor-m/z bins, so the pieces together have exactly the
# partitions of the single index over all bins. The search (`searchFragmentIndexPartitionMajorHinted`) loads one
# piece at a time, searches every scan whose isolation window overlaps it, and frees it; a scan's candidates are
# the union over pieces (pieces hold disjoint precursors). Pieces are stored as raw arrays (`write_index_piece`),
# which load at disk speed.

const INDEX_PIECES_MANIFEST = "index_pieces.json"
const INDEX_PIECE_MAGIC = UInt64(0x5049_4f4e_5049_4543)   # "PIONPIEC"
const INDEX_PIECE_VERSION = UInt32(1)

"One piece of a `PiecedFragmentIndex`: its file and the precursor m/z range it covers."
struct IndexPiece
    file::String
    prec_mz_min::Float32
    prec_mz_max::Float32
    n_partitions::Int
    n_precursors::Int
    n_fragments::Int
    bytes::Int
    id_type::String
end

"""
    PiecedFragmentIndex

A fragment index stored as pieces in `dir`, listed in `INDEX_PIECES_MANIFEST`; see the comment at the top of
`pieces.jl`. Searched by `searchFragmentIndexPartitionMajorHinted` like a single index.
"""
struct PiecedFragmentIndex
    dir::String
    pieces::Vector{IndexPiece}
    partition_width::Float32
end

"Bytes held by a partitioned index's arrays (what it occupies in memory once loaded, excluding object headers)."
function index_bytes(pfi::AbstractLocalPartitionedFragmentIndex)
    b = sizeof(pfi.partition_bounds)
    for p in getPartitions(pfi)
        fb = p.fragment_bins
        b += sizeof(fb.lows) + sizeof(fb.highs) + sizeof(fb.first_bins) + sizeof(fb.last_bins) +
             sizeof(p.rt_bins) + sizeof(p.fragments) + sizeof(p.local_to_global) + sizeof(p.skip_hints)
    end
    return b
end

_write_vec(io::IO, v::Vector) = (write(io, Int64(length(v))); write(io, v))
_read_vec(io::IO, ::Type{T}) where {T} = read!(io, Vector{T}(undef, read(io, Int64)))

"""
    write_index_piece(path, pfi)

Write a partitioned index as raw arrays: a header (magic, version, local ID type, partition count), the partition
bounds, then per partition its local precursor count and arrays. Read back with `read_index_piece`.
"""
function write_index_piece(path::AbstractString, pfi::AbstractLocalPartitionedFragmentIndex{Float32})
    open(path, "w") do io
        write(io, INDEX_PIECE_MAGIC, INDEX_PIECE_VERSION, UInt8(sizeof(local_id_type(pfi))), Int64(pfi.n_partitions))
        _write_vec(io, pfi.partition_bounds)
        for p in getPartitions(pfi)
            write(io, UInt32(p.n_local_precs))
            fb = p.fragment_bins
            _write_vec(io, fb.lows); _write_vec(io, fb.highs); _write_vec(io, fb.first_bins); _write_vec(io, fb.last_bins)
            _write_vec(io, p.rt_bins); _write_vec(io, p.fragments); _write_vec(io, p.local_to_global)
            _write_vec(io, p.skip_hints)
        end
    end
    return path
end

"Read a partitioned index written by `write_index_piece`."
function read_index_piece(path::AbstractString)
    open(path, "r") do io
        read(io, UInt64) == INDEX_PIECE_MAGIC || error("$path is not a fragment index piece")
        (v = read(io, UInt32)) == INDEX_PIECE_VERSION || error("$path: unsupported index piece version $v")
        id_bytes = read(io, UInt8)
        id_bytes == 2 ? _read_index_piece(io, UInt16) :
        id_bytes == 4 ? _read_index_piece(io, UInt32) : error("$path: bad local ID width $id_bytes")
    end
end

function _read_index_piece(io::IO, ::Type{I}) where {I<:Unsigned}
    n_partitions = Int(read(io, Int64))
    bounds = _read_vec(io, Tuple{Float32, Float32})
    PT = local_partition_type(I){Float32}
    partitions = Vector{PT}(undef, n_partitions)
    for k in 1:n_partitions
        n_local = I(read(io, UInt32))
        fb = SoAFragBins{Float32}(_read_vec(io, Float32), _read_vec(io, Float32), _read_vec(io, UInt32), _read_vec(io, UInt32))
        rt_bins = _read_vec(io, FragIndexBin{Float32})
        fragments = _read_vec(io, local_fragment_type(I))
        l2g = _read_vec(io, UInt32)
        hints = _read_vec(io, UInt16)
        partitions[k] = PT(fb, rt_bins, fragments, l2g, n_local, hints)
    end
    return local_index_type(I){Float32}(partitions, bounds, n_partitions)
end

"""
    build_index_pieces(sel, dir; partition_width, frag_bin_tol_ppm, frag_bin_tol_mda, rt_bin_tol,
                       id_type_request = "auto", max_piece_bytes = 4_000_000_000) -> PiecedFragmentIndex

Build the partitioned index of `sel` as pieces of at most `max_piece_bytes` (`index_bytes`) each, write them to
`dir` with the manifest `INDEX_PIECES_MANIFEST`, and return the `PiecedFragmentIndex`. Consecutive initial m/z bins
are grouped by an upper-bound size estimate; a piece that still comes out too large is split in half and rebuilt.
Each piece resolves its own local ID type (`resolve_local_id_type`). Only one piece is held in memory at a time.
"""
function build_index_pieces(sel::IndexFragSelection, dir::AbstractString;
                            partition_width::Float32 = 5.0f0,
                            frag_bin_tol_ppm::Float32 = 0.0f0,
                            frag_bin_tol_mda::Float32 = 2.0f0,
                            rt_bin_tol::Float32 = 3.0f0,
                            id_type_request::AbstractString = "auto",
                            max_piece_bytes::Integer = 4_000_000_000)
    mkpath(dir)
    bins = initial_partitions(sel.prec_mzs, partition_width)
    nfrag(pids) = sum(pid -> n_index_frags(sel, pid), pids; init = 0)
    # upper bound: 8-byte fragments, 4-byte local->global ids, one 18-byte m/z bin per two fragments
    est(k) = 17 * nfrag(bins[k]) + 4 * length(bins[k])

    groups = UnitRange{Int}[]
    lo, acc = 1, 0
    for k in eachindex(bins)
        e = est(k)
        if acc > 0 && acc + e > max_piece_bytes
            push!(groups, lo:(k - 1)); lo, acc = k, 0
        end
        acc += e
    end
    push!(groups, lo:length(bins))

    pieces = IndexPiece[]
    function build_group!(g::UnitRange{Int})
        pids = reduce(vcat, bins[g]; init = UInt32[])
        I = first(resolve_local_id_type(id_type_request, view(sel.prec_mzs, pids), partition_width))
        idx = build_partitioned_index_from_selection(sel; partition_width = partition_width,
            frag_bin_tol_ppm = frag_bin_tol_ppm, frag_bin_tol_mda = frag_bin_tol_mda, rt_bin_tol = rt_bin_tol,
            id_type = I, initial_partition_pids = [copy(bins[k]) for k in g])
        b = index_bytes(idx)
        if b > max_piece_bytes && length(g) > 1
            idx = nothing
            mid = first(g) + length(g) ÷ 2 - 1
            build_group!(first(g):mid); build_group!((mid + 1):last(g))
            return
        end
        b > max_piece_bytes && @user_warn "Fragment index piece of one $(partition_width) Da bin is $(round(b / 1e9, digits = 2)) GB, over the $(round(max_piece_bytes / 1e9, digits = 2)) GB limit"
        file = @sprintf("piece_%03d.bin", length(pieces) + 1)
        write_index_piece(joinpath(dir, file), idx)
        nonempty = filter(bd -> bd[1] <= bd[2], idx.partition_bounds)
        push!(pieces, IndexPiece(file,
            isempty(nonempty) ? Inf32 : minimum(first, nonempty), isempty(nonempty) ? -Inf32 : maximum(last, nonempty),
            idx.n_partitions, length(pids), nfrag(pids), b, string(I)))
        @debug_l1 "fragment index piece $(file): bins $(g), $(length(pids)) precursors, $(I), $(round(b / 1e9, digits = 2)) GB"
        return
    end
    for g in groups
        build_group!(g)
        GC.gc()
    end

    open(joinpath(dir, INDEX_PIECES_MANIFEST), "w") do io
        JSON.print(io, Dict{String, Any}("format_version" => 1, "partition_width_da" => partition_width,
            "pieces" => [Dict{String, Any}(string(f) => getfield(p, f) for f in fieldnames(IndexPiece)) for p in pieces]), 2)
    end
    return PiecedFragmentIndex(String(dir), pieces, partition_width)
end

"Open the pieced fragment index in `dir` (reads only the manifest; pieces are loaded when searched)."
function load_pieced_index(dir::AbstractString)
    m = JSON.parsefile(joinpath(dir, INDEX_PIECES_MANIFEST))
    pieces = [IndexPiece(String(p["file"]), Float32(p["prec_mz_min"]), Float32(p["prec_mz_max"]),
                         Int(p["n_partitions"]), Int(p["n_precursors"]), Int(p["n_fragments"]), Int(p["bytes"]),
                         String(p["id_type"])) for p in m["pieces"]]
    return PiecedFragmentIndex(String(dir), pieces, Float32(m["partition_width_da"]))
end

"""
    for_each_index_piece(f, pfi::PiecedFragmentIndex, scan_prec_min, scan_prec_max)

Load each piece whose precursor m/z range overlaps some scan's window, call `f(index)`, and free it before loading
the next, so at most one piece is in memory.
"""
function for_each_index_piece(f, pfi::PiecedFragmentIndex, scan_prec_min::Vector{Float32}, scan_prec_max::Vector{Float32})
    for piece in pfi.pieces
        any(i -> scan_prec_min[i] <= piece.prec_mz_max && scan_prec_max[i] >= piece.prec_mz_min,
            eachindex(scan_prec_min)) || continue
        index = read_index_piece(joinpath(pfi.dir, piece.file))
        f(index)
        index = nothing
        GC.gc(false)
    end
    return nothing
end
