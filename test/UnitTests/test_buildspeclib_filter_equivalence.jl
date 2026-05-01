# Byte-equivalence regression test for BuildSpecLib's keap1 synthetic build.
#
# Originally written to gate §9 commit 3 (move filters from build to
# predict time) — the reference was a snapshot from before that
# migration, and the test pinned that the migration didn't change
# library content. §9 has long since shipped, and several intentional
# changes have landed since (1fa48851: emit Compact frags directly;
# 55e5abc9: m/z-sort detailed_fragments before write). The reference
# at test/integration/keap1_reference/detailed_fragments.jls has been
# refreshed against current code; the test now pins "no further
# unintentional changes in the keap1 synthetic build's
# detailed_fragments.jls".
#
# To refresh after an INTENTIONAL change to BuildSpecLib output:
#
#     julia --project=. -e '
#         using Pioneer
#         Pioneer.with_koina_client(Pioneer.SyntheticKoinaClient()) do
#             BuildSpecLib("test/integration/build_keap1.json")
#         end'
#     cp test/integration/keap1_lib.poin/detailed_fragments.jls \
#        test/integration/keap1_reference/detailed_fragments.jls
#
# Review the diff (deserialize both, structural compare) before
# committing the new reference.

using Test
using Pioneer

@testset "BuildSpecLib early-filter equivalence (keap1, synthetic)" begin
    repo_root = abspath(joinpath(@__DIR__, "..", ".."))
    build_params = joinpath(repo_root, "test", "integration", "build_keap1.json")
    lib_path     = joinpath(repo_root, "test", "integration", "keap1_lib.poin")
    ref_dir      = joinpath(repo_root, "test", "integration", "keap1_reference")

    @test isfile(build_params)
    @test isfile(joinpath(ref_dir, "detailed_fragments.jls"))

    # Clean stale build output so the test runs from scratch.
    isdir(lib_path) && rm(lib_path; recursive=true, force=true)

    cd(repo_root) do
        Pioneer.with_koina_client(Pioneer.SyntheticKoinaClient()) do
            @test BuildSpecLib(build_params) === nothing
        end
    end

    # Byte-for-byte equivalence on detailed_fragments.jls — this is THE
    # invariant. The partitioned-index .jls files have different byte
    # representations on identical content (Serialization.serialize uses
    # non-deterministic object-id backreferences for some immutable structs
    # with vector fields like SoAFragBins), so we don't byte-compare them.
    # Since the indices are deterministic transforms of detailed_fragments
    # (by index-building code that §9 does not touch), pinning
    # detailed_fragments.jls bytes is sufficient to gate the §9 change.
    @test read(joinpath(ref_dir, "detailed_fragments.jls")) ==
          read(joinpath(lib_path, "detailed_fragments.jls"))

    # Cleanup
    isdir(lib_path) && rm(lib_path; recursive=true, force=true)
end
