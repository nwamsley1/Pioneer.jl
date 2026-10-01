# Determinism of BuildSpecLib (content-based).
#
# Verifies that, under the deterministic SyntheticKoinaClient:
#   * the same `seed` reproduces a byte-for-content-identical library
#     (precursors/proteins/fragments + the *content* of the fragment index —
#     the index .jls serialization layout is not byte-stable, but its content
#     is), and
#   * a different `seed` produces different shuffled decoy sequences (so the
#     library content changes).
#
# Run standalone:  julia --project=. -e 'using Pioneer; include("test/Routines/BuildSpecLib/test_build_determinism.jl")'

using Test
using Pioneer

const _REPO_ROOT = abspath(joinpath(@__DIR__, "..", "..", ".."))

# Content fingerprint of a built .poin. Byte-hashes the byte-deterministic
# payload files and content-hashes the fragment indexes (whose serialized bytes
# vary by object-encoding but whose content is deterministic). Uses Julia's
# content-based `hash` (not randomized across processes).
function _lib_content_fingerprint(lib::AbstractString)
    h = UInt(0)
    for f in ("precursors_table.arrow", "proteins_table.arrow", "detailed_fragments.jls")
        h = hash(read(joinpath(lib, f)), h)
    end
    for f in ("partitioned_fragment_index.jls", "presearch_partitioned_fragment_index.jls")
        idx = Pioneer.deserialize_from_jls(joinpath(lib, f))
        for p in idx.partitions
            s = p.fragment_bins
            h = hash(s.lows, h); h = hash(s.highs, h)
            h = hash(s.first_bins, h); h = hash(s.last_bins, h)
            h = hash(getfield.(p.fragments, 1), h)   # prec local id
            h = hash(getfield.(p.fragments, 2), h)   # score
            for fi in 1:4                            # FragIndexBin: lb, ub, first_bin, last_bin
                h = hash(getfield.(p.rt_bins, fi), h)
            end
            h = hash(p.local_to_global, h)
            h = hash(p.skip_hints, h)
            h = hash(p.n_local_precs, h)
        end
        h = hash(idx.partition_bounds, h)
        h = hash(idx.n_partitions, h)
    end
    return h
end

# Build keap1 offline with a given seed into <outdir>/<name>.poin; returns the path.
# `library_params_extra` (JSON members, e.g. `"prec_partition_width": 10.0`) is prepended to library_params.
function _build_keap1(outdir::AbstractString, name::AbstractString, seed::Int; library_params_extra::AbstractString = "")
    base = read(joinpath(_REPO_ROOT, "test", "integration", "build_keap1.json"), String)
    fwd(p) = replace(p, "\\" => "/")
    s = base
    s = replace(s, r"\"out_dir\"\s*:\s*\"[^\"]*\""      => "\"out_dir\": \"$(fwd(outdir))\"")
    s = replace(s, r"\"lib_name\"\s*:\s*\"[^\"]*\""     => "\"lib_name\": \"$name\"")
    s = replace(s, r"\"new_lib_name\"\s*:\s*\"[^\"]*\"" => "\"new_lib_name\": \"$name\"")
    s = replace(s, r"\"library_path\"\s*:\s*\"[^\"]*\"" => "\"library_path\": \"$(fwd(outdir))/$name\"")
    # Inject/override the seed (add it before the closing brace if absent).
    if occursin(r"\"seed\"\s*:", s)
        s = replace(s, r"\"seed\"\s*:\s*\d+" => "\"seed\": $seed")
    else
        s = replace(s, r"\}\s*$" => ",\n  \"seed\": $seed\n}")
    end
    if !isempty(library_params_extra)
        s = replace(s, r"\"library_params\"\s*:\s*\{" => "\"library_params\": {" * library_params_extra * ",")
    end
    cfg = joinpath(outdir, "cfg_$name.json")
    open(io -> write(io, s), cfg, "w")
    cd(_REPO_ROOT) do
        Pioneer.with_koina_client(Pioneer.SyntheticKoinaClient()) do
            BuildSpecLib(cfg)
        end
    end
    return joinpath(outdir, "$name.poin")
end

@testset "BuildSpecLib determinism (content-based, seeded)" begin
    tmp = mktempdir()
    try
        libA = _build_keap1(tmp, "det_a", 1844)
        libB = _build_keap1(tmp, "det_b", 1844)   # same seed
        libC = _build_keap1(tmp, "det_c", 9999)   # different seed

        fpA = _lib_content_fingerprint(libA)
        fpB = _lib_content_fingerprint(libB)
        fpC = _lib_content_fingerprint(libC)

        @test fpA == fpB        # same seed → identical library content
        @test fpA != fpC        # different seed → different shuffled decoys
    finally
        # Best-effort cleanup; on Windows the just-built Arrow files may still be
        # mmap-locked in-process, so don't fail the test on a cleanup ENOTEMPTY.
        GC.gc()
        try; rm(tmp; recursive=true, force=true); catch; end
    end
end

# (global precursor ID, rank bitmask) of every fragment in an index, sorted: the index content independent of how
# precursors are grouped into partitions and of the local ID width.
function _index_prec_scores(idx)
    v = Tuple{UInt32, UInt8}[]
    for p in idx.partitions, f in p.fragments
        push!(v, (p.local_to_global[f.local_id], f.score))
    end
    sort!(v)
end

@testset "BuildSpecLib fragment indexes (5 + 10 Da, local ID type, choice at search time)" begin
    tmp = mktempdir()
    try
        lib16 = _build_keap1(tmp, "fi16", 1844)
        lib32 = _build_keap1(tmp, "fi32", 1844;
            library_params_extra = "\"frag_index_local_id_type\": \"UInt32\", \"prec_partition_width\": 10.0")
        # "auto" (the default) on this small library: every 5 Da bin fits, so UInt16; the resolved type and width
        # are recorded in config.json
        cfg16 = Pioneer.JSON.parsefile(joinpath(lib16, "config.json"))["library_params"]
        @test cfg16["frag_index_local_id_type"] == "auto"
        @test cfg16["frag_index_local_id_type_resolved"] == "UInt16"
        @test cfg16["prec_partition_width_resolved"] == 5.0
        cfg32 = Pioneer.JSON.parsefile(joinpath(lib32, "config.json"))["library_params"]
        @test cfg32["frag_index_local_id_type_resolved"] == "UInt32" && cfg32["prec_partition_width_resolved"] == 10.0
        # the library payload does not depend on the index variant
        for f in ("precursors_table.arrow", "proteins_table.arrow", "detailed_fragments.jls")
            @test read(joinpath(lib16, f)) == read(joinpath(lib32, f))
        end
        for f in ("partitioned_fragment_index.jls", "presearch_partitioned_fragment_index.jls")
            i16 = Pioneer.deserialize_from_jls(joinpath(lib16, f))
            i32 = Pioneer.deserialize_from_jls(joinpath(lib32, f))
            @test i16 isa Pioneer.LocalPartitionedFragmentIndex{Float32}
            @test i32 isa Pioneer.LocalPartitionedFragmentIndex32{Float32}
            @test eltype(i32.partitions[1].fragments) === Pioneer.LocalFragment32
            @test _index_prec_scores(i16) == _index_prec_scores(i32)   # same fragments and rank bits per precursor
            # 10 Da partitions: every partition spans < 10 Da of precursor m/z
            @test all(b -> b[2] - b[1] < 10.0f0, i32.partition_bounds)
            # the library loads with either index type
            @test Pioneer.local_id_type(i32) === UInt32
        end

        # A default build has a 5 and a 10 Da index, listed in fragment_indices.json; the 5 Da pair keeps the
        # historical names. Both hold the same fragments and rank bits per precursor.
        desc = Pioneer.JSON.parsefile(joinpath(lib16, Pioneer.FRAGMENT_INDEX_DESCRIPTOR))
        @test desc["format_version"] == 1
        @test [e["partition_width_da"] for e in desc["indexes"]] == [5.0, 10.0]
        @test desc["indexes"][1]["main"] == "partitioned_fragment_index.jls"
        @test desc["indexes"][2]["main"] == "partitioned_fragment_index_w10.jls"
        @test desc["indexes"][2]["presearch"] == "presearch_partitioned_fragment_index_w10.jls"
        for e in desc["indexes"], k in ("main", "presearch")
            @test isfile(joinpath(lib16, e[k]))
        end
        i5 = Pioneer.deserialize_from_jls(joinpath(lib16, "partitioned_fragment_index.jls"))
        i10 = Pioneer.deserialize_from_jls(joinpath(lib16, "partitioned_fragment_index_w10.jls"))
        @test _index_prec_scores(i10) == _index_prec_scores(i5)
        @test all(b -> b[2] - b[1] < 10.0f0, i10.partition_bounds)
        # the hidden prec_partition_width override builds that one width, under the historical names
        desc32 = Pioneer.JSON.parsefile(joinpath(lib32, Pioneer.FRAGMENT_INDEX_DESCRIPTOR))
        @test length(desc32["indexes"]) == 1 && desc32["indexes"][1]["partition_width_da"] == 10.0
        @test desc32["indexes"][1]["main"] == "partitioned_fragment_index.jls"

        # SearchDIA's choice from the data: the median MS2 isolation width of the first file, in size classes
        function _ms_file(path, width)
            n = 6
            Pioneer.Arrow.write(path, (mz_array = [Union{Missing, Float32}[500f0] for _ in 1:n],
                intensity_array = [Union{Missing, Float32}[1f0] for _ in 1:n], scanHeader = fill("", n),
                scanNumber = Int32.(1:n), packetType = zeros(Int32, n), retentionTime = Float32.(1:n),
                lowMz = fill(100f0, n), highMz = fill(1700f0, n), TIC = ones(Float32, n),
                centerMz = Union{Missing, Float32}[missing; fill(612.5f0, n - 1)],
                isolationWidthMz = Union{Missing, Float32}[missing; fill(Float32(width), n - 1)],
                collisionEnergyField = Vector{Union{Missing, Float32}}(fill(30f0, n)),
                collisionEnergyEvField = zeros(Float32, n), msOrder = UInt8[1; fill(0x02, n - 1)],
                cycle_idx = Int32.(1:n)))
            path
        end
        wide = _ms_file(joinpath(tmp, "wide.arrow"), 25.0)      # timsTOF diaPASEF-like
        narrow = _ms_file(joinpath(tmp, "narrow.arrow"), 2.9)   # SCIEX / Astral-like
        c = Pioneer.choose_fragment_index(lib16, [wide])
        @test c.width == 10.0 && c.main == "partitioned_fragment_index_w10.jls" && c.window == 25.0
        c = Pioneer.choose_fragment_index(lib16, [narrow, wide])  # the first file decides
        @test c.width == 5.0 && c.main == "partitioned_fragment_index.jls" && c.window == Float64(2.9f0)
        @test Pioneer.choose_fragment_index(lib16, String[]).width == 5.0       # no data: the first index
        @test Pioneer.choose_fragment_index(lib32, [narrow]).width == 10.0      # a single index is always used
        # timsTOF data always gets the widest class, whatever its window width; the .tdfs is not even read
        tdfs = mkpath(joinpath(tmp, "run.tdfs"))
        c = Pioneer.choose_fragment_index(lib16, [tdfs])
        @test c.width == 10.0 && c.tims && c.window === nothing
        @test Pioneer.target_partition_width(2.9, true) == 10.0
        @test Pioneer.target_partition_width(nothing, true) == 10.0
        @test [Pioneer.target_partition_width(w, false) for w in (2.0, 3.7, 3.75, 7.4, 7.5, 25.0)] == [2.5, 2.5, 5.0, 5.0, 10.0, 10.0]
        @test Pioneer.target_partition_width(nothing, false) === nothing
        legacy = mktempdir()                                                    # no descriptor: an older library
        c = Pioneer.choose_fragment_index(legacy, [wide])
        @test c.width === nothing && c.main == "partitioned_fragment_index.jls"

        # add_fragment_indexes! on a library built before the descriptor existed (only the 5 Da pair): it adds the
        # 10 Da pair, rebuilt from the stored m/z-sorted fragments, identical in content to the one a fresh build
        # writes. (Content, not bytes: serialisation writes the uninitialised padding byte of each LocalFragment.)
        function _same_index(x, y)
            typeof(x) == typeof(y) && x.n_partitions == y.n_partitions && x.partition_bounds == y.partition_bounds ||
                return false
            for (p, q) in zip(x.partitions, y.partitions), f in fieldnames(typeof(p))
                u, v = getfield(p, f), getfield(q, f)
                if f === :fragments
                    [(g.local_id, g.score) for g in u] == [(g.local_id, g.score) for g in v] || return false
                elseif u isa Pioneer.SoAFragBins
                    all(getfield(u, g) == getfield(v, g) for g in fieldnames(typeof(u))) || return false
                else
                    u == v || return false
                end
            end
            return true
        end
        old = joinpath(tmp, "old.poin"); cp(lib16, old)
        for f in (Pioneer.FRAGMENT_INDEX_DESCRIPTOR, "partitioned_fragment_index_w10.jls",
                  "presearch_partitioned_fragment_index_w10.jls")
            rm(joinpath(old, f))
        end
        @test Pioneer.add_fragment_indexes!(old) == Float32[10.0]
        d = Pioneer.JSON.parsefile(joinpath(old, Pioneer.FRAGMENT_INDEX_DESCRIPTOR))
        @test [e["partition_width_da"] for e in d["indexes"]] == [5.0, 10.0]
        @test d["indexes"][1]["main"] == "partitioned_fragment_index.jls"
        for f in ("partitioned_fragment_index_w10.jls", "presearch_partitioned_fragment_index_w10.jls")
            @test _same_index(Pioneer.deserialize_from_jls(joinpath(old, f)), Pioneer.deserialize_from_jls(joinpath(lib16, f)))
        end
        @test isempty(Pioneer.add_fragment_indexes!(old))                       # nothing left to add
    finally
        GC.gc()
        try; rm(tmp; recursive=true, force=true); catch; end
    end
end
