# Pioneer test suite — part 1 (heavy machinery).
#
# Runs in its own Julia subprocess so it has a clean compiler and heap
# state. The integration test and the BuildSpecLib testset both compile
# a lot of code; in a single-process suite this corrupts compiler /
# heap state and downstream files crash. Splitting into subprocesses
# isolates them.
#
# This file is invoked by test/runtests.jl.

include(joinpath(@__DIR__, "test_setup.jl"))

# No outer @testset wrapper here on purpose — each include below has
# its own @testset and prints its own summary. Wrapping them in a
# parent @testset triggers a Julia 1.12 compiler crash at the finish
# print_test_results step (the cumulative compile load in this
# subprocess is enough to corrupt inference state). The driver's
# exit code is what matters; an unbroken sub-summary list is more
# useful than an aggregate that crashes.

println("dir ", @__DIR__)

@testset "process_test" begin
    # Integration test: SearchDIA only (no library build).
    # The library at test/integration/ecoli_lib.poin/ is committed to
    # the repo so CI / GitHub Actions can run this test without needing
    # Koina. To rebuild locally:
    #   julia --project=. -t auto test/integration/run.jl --build-only
    @test SearchDIA(joinpath(@__DIR__, "integration", "search_ecoli.json")) === nothing
end

# FASTA parameter enhancement
include("./Routines/BuildSpecLib/params/test_fasta_params.jl")

# BuildSpecLib end-to-end (uses the offline Koina fixture at
# test/Routines/BuildSpecLib/koina_fixtures.json).
include("./Routines/BuildSpecLib/test_build_spec_lib.jl")

# Coverage data is normally flushed to .cov files by an atexit hook
# that `_exit` skips. Flush it explicitly before terminating, otherwise
# `--code-coverage=...` produces zero output and codecov reports 0%.
Base.JLOptions().code_coverage != 0 &&
    ccall(:jl_write_coverage_data, Cvoid, (Cstring,), C_NULL)

# Bypass Julia's normal shutdown — see runtests_part2_units.jl for why.
ccall(:_exit, Cvoid, (Cint,), 0)
