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
function _build_keap1(outdir::AbstractString, name::AbstractString, seed::Int)
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
