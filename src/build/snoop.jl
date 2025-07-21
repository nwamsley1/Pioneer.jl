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

# Execute a small search to precompile search routines
SearchDIA(joinpath(data_dir, "ecoli_test", "ecoli_test_params.json"))



#BuildSpecLib("build/buildspeclib_params.json")
#convertMzML("data/mzML/")