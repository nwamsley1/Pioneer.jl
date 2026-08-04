using Pioneer

root = joinpath(@__DIR__)
data_dir = joinpath(root, "..", "..", "data")

cmd = get(ENV, "PIONEER_CMD", nothing)

# A failing target used to be swallowed into a warning, so the build stayed green while shipping a
# binary with that target's code paths uncompiled. That is exactly the failure we cannot detect from
# the outside: `SearchDIA` silently stopped precompiling when 6bab0dadb changed the .poin layout, and
# nothing surfaced it. Failures are now collected and rethrown at the end of the script, which fails
# `create_app` rather than degrading the artifact.
#
# Set PIONEER_SNOOP_ALLOW_FAILURES=1 to restore warn-only behaviour when *debugging* the script
# itself. Never set it in a release build.
const SNOOP_FAILURES = Tuple{String, String}[]
const SNOOP_ALLOW_FAILURES = get(ENV, "PIONEER_SNOOP_ALLOW_FAILURES", "0") == "1"

function maybe_run(f, name)
    if cmd === nothing || cmd == name
        start_time = time()
        @user_info "Starting precompile target $name"
        try
            f()
            elapsed = round(time() - start_time; digits=1)
            @user_info "Finished precompile target $name in $(elapsed)s"
        catch e
            bt = catch_backtrace()
            target_cmd = cmd === nothing ? "<all>" : cmd
            elapsed = round(time() - start_time; digits=1)
            @user_warn "Error executing $name during precompile of $target_cmd: $(sprint(showerror, e))"
            @warn "Precompile exception details for $name after $(elapsed)s" exception=(e, bt)
            push!(SNOOP_FAILURES, (name, sprint(showerror, e)))
        end
    end
end

# Called once at the end of the script, after every target has had its turn, so a single run reports
# *all* broken targets rather than stopping at the first -- which matters when each target takes
# minutes.
function report_snoop_failures()
    isempty(SNOOP_FAILURES) && return nothing
    summary = join(("  $n: $msg" for (n, msg) in SNOOP_FAILURES), "\n")
    if SNOOP_ALLOW_FAILURES
        @user_warn "PIONEER_SNOOP_ALLOW_FAILURES=1; $(length(SNOOP_FAILURES)) precompile target(s) failed:\n$summary"
        return nothing
    end
    error("$(length(SNOOP_FAILURES)) precompile target(s) failed; the app would ship without them " *
          "compiled:\n$summary")
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
# ParseSpecLib is currently disabled in importScripts.jl due EmpiricalLibrary loading issues.
# Keep this precompile target disabled until ParseSpecLib is re-enabled in the module.
# maybe_run("ParseSpecLib") do
#     Pioneer.ParseSpecLib(joinpath(data_dir, "precompile", "build_empirical.json"))
# end



##########################################
# Predict libraries
##########################################

# Build a tiny Prosit library and search it with low memory thresholds
#maybe_run("BuildSpecLib") do
#    Pioneer.BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_prosit.json"))
#end
# Build a tiny Altimeter library
maybe_run("BuildSpecLib") do
    Pioneer.BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_altimeter.json"))
end


##########################################
# Search
##########################################

# Real-world Arrow files are written in several record batches, so a column reads back as a
# `SentinelArrays.ChainedVector`. The checked-in ecoli fixtures are single-batch, so their columns
# read back as bare `Arrow.Primitive`. The per-scan getters in src/structs/MassSpecData/getters.jl
# are annotated with a two-member Union over exactly those two representations
# (`_MSDataCol` / `_MSDataPeakCol`), which means the two layouts specialise every method downstream
# of a getter differently.
#
# Precompiling only the single-batch arm therefore leaves the multi-batch arm -- the one real user
# data actually takes -- uncompiled. Measured against a real 6-file Exploris search, that was a 43%
# gap in Pioneer-owned methods, including `extract_chromatograms` and the deconvolution solver.
#
# So materialise a second copy of one fixture written as two record batches, and search a directory
# holding both layouts. Generated at snoop time rather than checked in, so it cannot drift from the
# fixtures it is derived from.
function prepare_mixed_batch_ms_data(raw_dir::String, out_dir::String)
    rm(out_dir; force = true, recursive = true)
    mkpath(out_dir)
    arrows = sort(filter(f -> endswith(f, ".arrow"), readdir(raw_dir)))
    for f in arrows
        cp(joinpath(raw_dir, f), joinpath(out_dir, f))
    end
    isempty(arrows) && return out_dir

    src = joinpath(raw_dir, first(arrows))
    tbl = Pioneer.Arrow.Table(src)
    cols = Pioneer.Tables.columnnames(tbl)
    n = length(Pioneer.Tables.getcolumn(tbl, first(cols)))
    n < 2 && return out_dir
    mid = n ÷ 2
    chunk(r) = (; (k => Pioneer.Tables.getcolumn(tbl, k)[r] for k in cols)...)
    Pioneer.Arrow.write(
        joinpath(out_dir, "ecoli_multibatch.arrow"),
        Pioneer.Tables.partitioner([chunk(1:mid), chunk(mid+1:n)]),
    )
    return out_dir
end

# Searches the library the BuildSpecLib target above just produced, NOT the Zenodo yeast.poin.
#
# The Zenodo artifact (record 16289168) was built before 6bab0dadb changed the .poin layout to
# `partitioned_fragment_index.jls`, so searching it now throws
#     SystemError: opening file ".../yeast.poin/partitioned_fragment_index.jls"
# and `maybe_run` swallows that into a warning -- the build still succeeds, but ships with none of
# SearchDIA precompiled, which is the single biggest target here. Sourcing the library from the
# preceding BuildSpecLib target keeps this self-consistent: whatever code builds the library is the
# same code that searches it, so a future format change cannot silently disable this target again.
#
# To go back to yeast, regenerate the Zenodo artifact with the current BuildSpecLib and restore the
# search_yeast_altimeter*.json configs -- the yeast data is larger and exercises more scans per file.
maybe_run("SearchDIA") do
    prepare_mixed_batch_ms_data(
        joinpath(data_dir, "ecoli_test", "raw"),
        joinpath(data_dir, "precompile", "ecoli_raw_mixed"),
    )
    Pioneer.SearchDIA(joinpath(data_dir, "precompile", "search_ecoli_altimeter.json"))      # altimeter + MBR
    Pioneer.SearchDIA(joinpath(data_dir, "precompile", "search_ecoli_altimeter_OOM.json"))  # altimeter + MBR + OOM
end


##########################################
# ConvertMzML
##########################################
maybe_run("convertMzML") do
    Pioneer.convertMzML(joinpath(data_dir, "precompile", "convert_example.mzML"))
end


##########################################
# Fail the build if any target did not run
##########################################
report_snoop_failures()
