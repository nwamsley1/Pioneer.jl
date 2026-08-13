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

# Build a tiny Altimeter library
maybe_run("BuildSpecLib") do
    Pioneer.BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_altimeter.json"))
end

# Build the same library through the Prosit (InstrumentAgnosticModel) path.
#
# This is not redundant with the Altimeter target. Prosit returns scalar intensities
# rather than spline coefficients, so buildPionLib produces a StandardFragmentLookup
# (ConstantType) where Altimeter produces a SplineFragmentLookup, and no spline knots
# are written. Different concrete type, so every method downstream of a fragment
# getter specialises differently -- the same trap documented for the
# ChainedVector/Primitive split below, which measured a 43% gap in Pioneer-owned
# methods when only one arm was precompiled.
#
# Same FASTA and m/z window as the Altimeter fixture, so the two are directly
# comparable. Measured locally: 117s to build (about half the Altimeter target,
# since there are no spline coefficients to fit), and searching it recovers 1,707
# precursors / 855 protein groups against Altimeter's 1,822 / 902 -- 89% Jaccard
# overlap on precursors, which is the expected level of agreement between two
# independent fragment predictors rather than a degenerate subset.
maybe_run("BuildSpecLib_prosit") do
    Pioneer.BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_prosit.json"))
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
    # Two snoop targets share this directory: SearchDIA builds it, then
    # SearchDIA_prosit asks for it again. By then the first target's searches
    # have memory-mapped the Arrow files inside it, and Windows refuses to
    # delete a mapped file -- the rm fails with ENOTEMPTY and takes the whole
    # app build down with it. POSIX unlinks an open file happily, which is why
    # this only ever failed on Windows.
    #
    # The contents are deterministic, so if they are already there, they are
    # already correct: reuse them rather than deleting and rewriting.
    expected = sort(filter(f -> endswith(f, ".arrow"), readdir(raw_dir)))
    if isdir(out_dir)
        have = Set(readdir(out_dir))
        if all(in(have), expected) && ("ecoli_multibatch.arrow" in have)
            return out_dir
        end
        # Present but incomplete: a previous build died partway. Rebuilding is
        # worth attempting, but a mapped leftover must not fail the build.
        try
            rm(out_dir; force = true, recursive = true)
        catch err
            @warn "could not clear $out_dir; reusing what is there" exception = err
        end
    end
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
maybe_run("SearchDIA") do
    prepare_mixed_batch_ms_data(
        joinpath(data_dir, "ecoli_test", "raw"),
        joinpath(data_dir, "precompile", "ecoli_raw_mixed"),
    )
    Pioneer.SearchDIA(joinpath(data_dir, "precompile", "search_ecoli_altimeter.json"))      # altimeter + MBR
    Pioneer.SearchDIA(joinpath(data_dir, "precompile", "search_ecoli_altimeter_OOM.json"))  # altimeter + MBR + OOM
end

# Search the Prosit library built above. Building it is not enough: the search side
# reads a StandardFragmentLookup (ConstantType) where every other target reads a
# SplineFragmentLookup, so without this the whole search pipeline stays uncompiled
# for anyone using a Prosit library and re-JITs on first use.
#
# Regenerates the mixed-batch fixture rather than assuming the SearchDIA target ran,
# so this target still works when invoked alone via the `cmd` filter. The function
# rm -rf's its output first, so calling it twice is safe.
maybe_run("SearchDIA_prosit") do
    prepare_mixed_batch_ms_data(
        joinpath(data_dir, "ecoli_test", "raw"),
        joinpath(data_dir, "precompile", "ecoli_raw_mixed"),
    )
    Pioneer.SearchDIA(joinpath(data_dir, "precompile", "search_ecoli_prosit.json"))
end


##########################################
# Search -- yeast
##########################################

# The ecoli fixtures cannot reach two things, and both matter.
#
# 1. MS1. `ecoli_filtered_01.arrow` is 19,275 scans, every one msOrder=2. MainSearch computes a
#    whole MS1 feature block (features.jl: ms1_m0_mass_err_ppm, ms1_m1_to_m0_ratio,
#    ms1_isotope_dotp_m0_m1_m2, ...) via add_ms1_lookup_features!. Julia compiles whole method
#    bodies, so those methods do compile on ecoli even though they early-return -- but the
#    MS1-present path instantiates types the MS1-absent path never does.
#
# 2. Scale and shape. Measured against the snoop trace, a real yeast search compiles 301 statements
#    (177 naming Pioneer) that the ecoli fixtures do not, and 203 of those are not covered by a real
#    Olsen search either. One dataset is not a substitute for another: a binary precompiled from one
#    and run on the other re-JITs ~150 Pioneer methods.
#
# The fixture is deliberately small: all MS1 scans kept, MS2 restricted to a 500-600 m/z band (the
# same convention build_ecoli_altimeter.json uses), full RT range -- squeezing RT would make the RT
# and quad spline fits degenerate, and those Interpolations paths are part of what we want compiled.
#
# Three files, and the split matters. rep1/rep2 are ALTERNATING CYCLES of a single run, so both span
# the full 0-4.98 min gradient and both converge in parameter tuning (4,049 and 4,171 PSMs against
# min_samples 3500) and quad tuning. lowsignal is 1 MB of void volume that yields zero PSMs; it is
# there on purpose, because it is the only file that compiles the fallback paths -- parameter-tuning
# fallback, SquareQuadModel, RT-recalibration skip -- which is what real users with a bad run hit.
# Together: 3 files, 2 converged, 1 fallback; 17,444 precursors / 5,001 protein groups.
#
# Produced by split_yeast.jl; see .github/actions/precompile-data for the Zenodo artifact.
#
# This target reads a PRE-BUILT library rather than building one, because Koina is slow and
# occasionally down, and a build-time network dependency is a poor thing to put on the critical path
# of a release. That is the arrangement that produced bug #4 -- a stale Zenodo .poin whose failure
# `maybe_run` swallowed into a warning, shipping a binary with SearchDIA uncompiled. It is safe now
# only because report_snoop_failures() rethrows: a stale library fails the build instead of
# degrading the artifact. Bump the precompile-data cache-key whenever the artifact is regenerated.
maybe_run("SearchDIA_yeast") do
    Pioneer.SearchDIA(joinpath(data_dir, "precompile", "search_yeast_altimeter.json"))
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
