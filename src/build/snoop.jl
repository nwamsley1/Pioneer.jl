using Pioneer

root = joinpath(@__DIR__)
data_dir = joinpath(root, "data")

# Generate library build parameters using bundled test FASTA
GetBuildLibParams(mktempdir(), "test_lib", joinpath(data_dir, "fasta"))

# Generate search parameters using bundled ecoli test data
GetSearchParams(
    joinpath(data_dir, "ecoli_test", "altimeter_ecoli.poin"),
    joinpath(data_dir, "ecoli_test", "raw"),
    mktempdir(),
)

# Parse example empirical library to compile parsing code
ParseSpecLib(joinpath(data_dir, "library_test", "defaultParseEmpiricalLibParams2.json"))


# Build a tiny Prosit library and search it with low memory thresholds
BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_prosit.json"))
#SearchDIA(joinpath(data_dir, "precompile", "search_ecoli_prosit.json"))

# Build a tiny Altimeter library
#BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_prosit.json"))

# Execute a small Altimeter search
SearchDIA(joinpath(data_dir, "ecoli_test", "ecoli_test_params.json"))









#convertMzML("data/mzML/")