# Unit-test wrapper for the offline BuildSpecLib smoke test.
#
# Runs BuildSpecLib end-to-end under SyntheticKoinaClient against the
# small keap1.fasta input. Pins the contract that no Koina endpoint is
# hit, that the full pipeline (digest → predict → parse → build →
# index) doesn't crash, and that the produced .poin has the expected
# files. Output dir is cleaned up afterwards (gitignored anyway).

using Test
using Pioneer

@testset "BuildSpecLib offline (SyntheticKoinaClient) — keap1.fasta" begin
    repo_root    = abspath(joinpath(@__DIR__, "..", ".."))
    build_params = joinpath(repo_root, "test", "integration", "build_keap1.json")
    lib_path     = joinpath(repo_root, "test", "integration", "keap1_lib.poin")
    fasta_path   = joinpath(repo_root, "data", "fasta", "keap1.fasta")

    @test isfile(fasta_path)
    @test isfile(build_params)

    # Clean any previous output so the test is deterministic.
    safe_rmdir(lib_path)

    cd(repo_root) do
        Pioneer.with_koina_client(Pioneer.SyntheticKoinaClient()) do
            @test BuildSpecLib(build_params) === nothing
        end
    end

    @test isdir(lib_path)
    for f in ("config.json", "detailed_fragments.jls",
              "partitioned_fragment_index.jls",
              "presearch_partitioned_fragment_index.jls",
              "precursors_table.arrow", "proteins_table.arrow")
        @test isfile(joinpath(lib_path, f))
    end

    # Cleanup
    safe_rmdir(lib_path)
end
