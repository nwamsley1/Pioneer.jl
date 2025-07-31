using Pioneer

root = joinpath(@__DIR__)
data_dir = joinpath(root, "..", "..", "data")

cmd = get(ENV, "PIONEER_CMD", nothing)

function maybe_run(name, f)
    if cmd === nothing || cmd == name
        try
            f()
        catch e
            @warn "Error executing $name during precompile of $cmd " exception=(e, catch_backtrace()) 
        end
    end
end

##########################################
# Generate params
##########################################

# search
maybe_run("GetSearchParams") do
    GetSearchParams(
        joinpath(data_dir, "ecoli_test", "altimeter_ecoli.poin"),
        joinpath(data_dir, "ecoli_test", "raw"),
        mktempdir(),
    )
end

# predict
maybe_run("GetBuildLibParams") do
    GetBuildLibParams(mktempdir(), "test_lib", joinpath(data_dir, "fasta"))
end

# empirical
maybe_run("GetParseSpecLibParams") do
    GetParseSpecLibParams("test_lib", mktempdir())
end


##########################################
# Empirical libraries
##########################################
maybe_run("ParseSpecLib") do
    ParseSpecLib(joinpath(data_dir, "precompile", "build_empirical.json"))
end



##########################################
# Predict libraries
##########################################

# Build a tiny Prosit library and search it with low memory thresholds
maybe_run("BuildSpecLib") do
    BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_prosit.json"))
end
# TODO
# Build a tiny Altimeter library
#maybe_run("BuildSpecLib") do
#    BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_altimeter.json"))
#end


##########################################
# Search
##########################################

maybe_run("SearchDIA") do
    SearchDIA(joinpath(data_dir, "precompile", "search_ecoli_prosit.json"))         # prosit
    SearchDIA(joinpath(data_dir, "precompile", "search_yeast_altimeter.json"))      # altimeter + MBR
    SearchDIA(joinpath(data_dir, "precompile", "search_yeast_altimeter_OOM.json"))  # altimeter + MBR + OOM
end


##########################################
# ConvertMzML
##########################################
maybe_run("convertMzML") do
    convertMzML(joinpath(data_dir, "precompile", "convert_example.mzML"))
end
