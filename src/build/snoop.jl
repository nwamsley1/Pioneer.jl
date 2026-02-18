using Pioneer

root = joinpath(@__DIR__)
data_dir = joinpath(root, "..", "..", "data")

cmd = get(ENV, "PIONEER_CMD", nothing)

function maybe_run(f, name)
    if cmd === nothing || cmd == name
        try
            f()
        catch e
            bt = catch_backtrace()
            target_cmd = cmd === nothing ? "<all>" : cmd
            @user_warn "Error executing $name during precompile of $target_cmd: $(sprint(showerror, e))"
            @warn "Precompile exception details for $name" exception=(e, bt)
        end
    end
end

##########################################
# Generate params
##########################################

# search
maybe_run("GetSearchParams") do
    Pioneer.GetSearchParams(
        joinpath(data_dir, "ecoli_test", "altimeter_ecoli.poin"),
        joinpath(data_dir, "ecoli_test", "raw"),
        mktempdir(),
    )
end

# predict
maybe_run("GetBuildLibParams") do
    Pioneer.GetBuildLibParams(mktempdir(), "test_lib", joinpath(data_dir, "fasta"))
end

# empirical
maybe_run("GetParseSpecLibParams") do
    Pioneer.GetParseSpecLibParams("test_lib", mktempdir())
end


##########################################
# Empirical libraries
##########################################
maybe_run("ParseSpecLib") do
    Pioneer.ParseSpecLib(joinpath(data_dir, "precompile", "build_empirical.json"))
end



##########################################
# Predict libraries
##########################################

# Build a tiny Prosit library and search it with low memory thresholds
maybe_run("BuildSpecLib") do
    Pioneer.BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_prosit.json"))
end
# Build a tiny Altimeter library
maybe_run("BuildSpecLib") do
    Pioneer.BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_altimeter.json"))
end


##########################################
# Search
##########################################

maybe_run("SearchDIA") do
    Pioneer.SearchDIA(joinpath(data_dir, "precompile", "search_ecoli_prosit.json"))         # prosit
    Pioneer.SearchDIA(joinpath(data_dir, "precompile", "search_yeast_altimeter.json"))      # altimeter + MBR
    Pioneer.SearchDIA(joinpath(data_dir, "precompile", "search_yeast_altimeter_OOM.json"))  # altimeter + MBR + OOM
end


##########################################
# ConvertMzML
##########################################
maybe_run("convertMzML") do
    Pioneer.convertMzML(joinpath(data_dir, "precompile", "convert_example.mzML"))
end
