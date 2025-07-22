using Pioneer

root = joinpath(@__DIR__)
data_dir = joinpath(root, "..", "..", "data")

##########################################
# Generate params
##########################################

# search
GetSearchParams(
    joinpath(data_dir, "ecoli_test", "altimeter_ecoli.poin"),
    joinpath(data_dir, "ecoli_test", "raw"),
    mktempdir(),
)
# predict
GetBuildLibParams(mktempdir(), "test_lib", joinpath(data_dir, "fasta"))
# empirical
GetParseSpecLibParams( "test_lib", mktempdir())


##########################################
# Empirical libraries
##########################################
ParseSpecLib(joinpath(data_dir, "precompile", "build_empirical.json"))


##########################################
# Predict libraries
##########################################

# Build a tiny Prosit library and search it with low memory thresholds
BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_prosit.json"))
# TODO
# Build a tiny Altimeter library
#BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_prosit.json"))


##########################################
# Search
##########################################

# prosit
SearchDIA(joinpath(data_dir, "precompile", "search_ecoli_prosit.json"))
# altimeter + MBR
SearchDIA(joinpath(data_dir, "precompile", "search_yeast_altimeter.json"))
# altimeter + MBR + OOM
SearchDIA(joinpath(data_dir, "precompile", "search_yeast_altimeter_OOM.json"))



##########################################
# ConvertMzML
##########################################
#convertMzML("data/mzML/")