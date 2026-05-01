# Integration test driver: runnable from a fresh shell, no REPL state.
#
#   julia --project=. -t auto test/integration/run.jl              # search only (uses committed .poin)
#   julia --project=. -t auto test/integration/run.jl --build      # rebuild library locally, then search
#   julia --project=. -t auto test/integration/run.jl --build-only # rebuild library locally, skip search
#
# The library at test/integration/ecoli_lib.poin/ IS committed to the repo
# so CI / GitHub Actions can run the search-only path without Koina (which
# requires internet for fragment prediction). Rebuild only when you want
# to refresh the committed library — then `git add` the changed files.
#
# Exit codes: 0 success, 1 failure.

using Pioneer

const REPO_ROOT     = abspath(joinpath(@__DIR__, "..", ".."))
const BUILD_PARAMS  = joinpath(REPO_ROOT, "test", "integration", "build_ecoli.json")
const SEARCH_PARAMS = joinpath(REPO_ROOT, "test", "integration", "search_ecoli.json")
const LIB_PATH      = joinpath(REPO_ROOT, "test", "integration", "ecoli_lib.poin")
const FASTA_PATH    = joinpath(REPO_ROOT, "data", "precompile", "ecoli.fasta")

# Run from REPO_ROOT so the relative paths inside the JSON configs resolve.
cd(REPO_ROOT)

# CLI flags
do_build = "--build" in ARGS || "--build-only" in ARGS
do_search = !( "--build-only" in ARGS )

# If no .poin exists and the user didn't ask for build, surface the issue
# clearly rather than letting SearchDIA fail with an opaque error.
if do_search && !do_build && !isdir(LIB_PATH)
    @error "No library at $LIB_PATH. Re-run with --build, or build separately with BuildSpecLib."
    exit(1)
end

if do_build
    if !isfile(FASTA_PATH)
        @error "Missing FASTA at $FASTA_PATH"
        exit(1)
    end
    @info "Building library: $BUILD_PARAMS → $LIB_PATH"
    try
        BuildSpecLib(BUILD_PARAMS)
    catch e
        @error "BuildSpecLib failed" exception=(e, catch_backtrace())
        exit(1)
    end
    if !isdir(LIB_PATH)
        @error "BuildSpecLib finished but $LIB_PATH was not produced"
        exit(1)
    end
    @info "Library built successfully"
end

if do_search
    @info "Running SearchDIA: $SEARCH_PARAMS"
    try
        result = SearchDIA(SEARCH_PARAMS)
        if result !== nothing
            @error "SearchDIA returned non-nothing result: $result"
            exit(1)
        end
    catch e
        @error "SearchDIA failed" exception=(e, catch_backtrace())
        exit(1)
    end
    @info "Integration test PASSED"
end

exit(0)
