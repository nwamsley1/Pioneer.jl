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

@testset "BuildSpecLib UInt32 fragment index (frag_index_local_id_type)" begin
    tmp = mktempdir()
    try
        lib16 = _build_keap1(tmp, "fi16", 1844)
        lib32 = _build_keap1(tmp, "fi32", 1844;
            library_params_extra = "\"frag_index_local_id_type\": \"UInt32\", \"prec_partition_width\": 10.0")
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
    finally
        GC.gc()
        try; rm(tmp; recursive=true, force=true); catch; end
    end
end
