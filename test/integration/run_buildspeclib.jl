# Offline BuildSpecLib smoke test.
#
#   julia --project=. -t auto test/integration/run_buildspeclib.jl
#
# Builds a synthetic keap1.fasta library using SyntheticKoinaClient — no
# network required. The output library is NOT biologically meaningful
# (intensities/iRTs are hashed from the peptide string), so it must NOT
# be used as input to SearchDIA in production. Its purpose is purely
# end-to-end pipeline coverage: every code path from FASTA digestion
# through fragment_predict → fragment_parse → buildPionLib runs against
# deterministic data.
#
# Exit codes: 0 success, 1 failure.

using Pioneer

const REPO_ROOT     = abspath(joinpath(@__DIR__, "..", ".."))
const BUILD_PARAMS  = joinpath(REPO_ROOT, "test", "integration", "build_keap1.json")
const LIB_PATH      = joinpath(REPO_ROOT, "test", "integration", "keap1_lib.poin")
const FASTA_PATH    = joinpath(REPO_ROOT, "data", "fasta", "keap1.fasta")

cd(REPO_ROOT)

if !isfile(FASTA_PATH)
    @error "Missing FASTA at $FASTA_PATH"
    exit(1)
end

# Clean any previous output so we genuinely smoke-test from scratch.
isdir(LIB_PATH) && rm(LIB_PATH; recursive=true, force=true)

@info "Building synthetic keap1 library: $BUILD_PARAMS → $LIB_PATH"
try
    Pioneer.with_koina_client(Pioneer.SyntheticKoinaClient()) do
        BuildSpecLib(BUILD_PARAMS)
    end
catch e
    @error "BuildSpecLib failed under SyntheticKoinaClient" exception=(e, catch_backtrace())
    exit(1)
end

# Sanity checks on the produced .poin
if !isdir(LIB_PATH)
    @error "BuildSpecLib finished but $LIB_PATH was not produced"
    exit(1)
end

required = ["config.json", "detailed_fragments.jls", "partitioned_fragment_index.jls",
            "presearch_partitioned_fragment_index.jls", "precursors_table.arrow",
            "proteins_table.arrow"]
missing_files = [f for f in required if !isfile(joinpath(LIB_PATH, f))]
if !isempty(missing_files)
    @error "Library is missing expected files: $(join(missing_files, ", "))"
    exit(1)
end

@info "Synthetic BuildSpecLib smoke test PASSED"
exit(0)
