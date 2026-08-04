# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# Entry point for PackageCompiler
function main_SearchDIA(argv=ARGS)::Cint
    
    settings = ArgParseSettings(; autofix_names = true)
    @add_arg_table! settings begin
        "params_path"
            help = "Path to search parameters JSON file"
            arg_type = String
    end
    parsed_args = parse_args(argv, settings; as_symbols = true)
    
    try
        SearchDIA(parsed_args[:params_path])
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

"""
    asset_path(parts...)
    Return the path to a bundled asset. The function first checks the
    compile-time `assets/` directory and falls back to a path relative to the
    installed executable.
"""
function asset_path(parts...)
    compile_dir = joinpath(@__DIR__, "..", "..", "assets", parts...)
    if ispath(compile_dir)
        return compile_dir
    end
    exe = PROGRAM_FILE
    if !isabspath(exe)
        exe_full = Sys.which(exe)
        exe = exe_full !== nothing ? exe_full : exe
    end
    exe_dir = abspath(dirname(realpath(exe)))
    return joinpath(exe_dir, "..", "data", parts...)


    
end

"""
    Locate the isotope spline XML file bundled with the application.
"""
isotope_spline_path() = asset_path("IsotopeSplines_10kDa_21isotopes.xml")


"""
    SearchDIA(params_path::String)

Main entry point for the DIA (Data-Independent Acquisition) search workflow.
Executes a series of `SearchMethod`s and generates performance metrics.

Parameters:
- params_path: Path to JSON configuration file containing search parameters

Output:
- Generates a log file in the results directory
- Long and wide-formatted tables (.arrow and .csv) for protein-group and precursor level id's and quantitation.
- Reports timing and memory usage statistics

Example:
```julia
julia> SearchDIA("/path/to/config.json")
==========================================================================================
Sarting SearchDIA
==========================================================================================

Starting search at: 2024-12-30T14:01:01.510
Output directory: ./../data/ecoli_test/ecoli_test_results
[ Info: Loading Parameters...
[ Info: Loading Spectral Library...
 .
 .
 .
```
If it does not already exist, SearchDIA creates the user-specified results_dir and generates quality control plots, data tables, and logs.
```
results_dir/
├── pioneer_search_log.txt
├── qc_plots/
│   ├── collision_energy_alignment/
│   │   └── nce_alignment_plots.pdf
│   ├── quad_transmission_model/
│   │   ├── quad_data
│   │   │   └── quad_data_plots.pdf
│   │   └── quad_models
│   │       └── quad_model_plots.pdf
│   ├── rt_alignment_plots/
│   │   └── rt_alignment_plots.pdf
│   ├── mass_error_plots/
│   │   └── mass_error_plots.pdf
│   └── QC_PLOTS.pdf
├── precursors_long.arrow
├── precursors_long.tsv
├── precursors_wide.arrow
├── precurosrs_wide.tsv
├── protein_groups_long.arrow
├── protein_groups_long.tsv
├── protein_groups_wide.arrow
└── protein_groups_wide.tsv
```
"""
function SearchDIA(params_path::String)
    # Clean up any old file handlers in case the program crashed
    GC.gc()
     #params_path = normpath("C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\OlsenAstralMS1SpeedTest\\OlsenAstralMS1SpeedTest.json.json")
   
    # === Initialize logging ===
    checkParams(params_path)
    params_string = read(params_path, String)

    params = parse_pioneer_parameters(params_path)
    mkpath(params.paths[:results])  # Ensure directory exists

    debug_level = haskey(params.logging, :debug_console_level) ?
                  params.logging.debug_console_level : 0
    # Independent of the console: the debug log is useful by default (see DEBUG_FILE_LEVEL).
    debug_file_level = haskey(params.logging, :debug_file_level) ?
                       params.logging.debug_file_level : 1
    max_bytes = haskey(params.logging, :max_message_bytes) ?
                try Int(params.logging.max_message_bytes) catch; 4096 end : 4096
    warnings_full_path = init_pioneer_logging(
        params.paths[:results],
        "Pioneer Search Log";
        debug_console_level = debug_level,
        debug_file_level    = debug_file_level,
        max_message_bytes   = max_bytes,
    )

    try
        @user_print "\n" * repeat("=", 90)
        @user_print "Starting SearchDIA"
        @user_print repeat("=", 90)
        @user_info "\nStarting search at: $(Dates.now())"
        @user_info "Output directory: $(params.paths[:results])"
        

        # === Initialize performance tracking ===
        timings = Dict{String, Any}()
        
        # === Load and validate data paths ===
        @user_info "Loading Parameters..."
        params_timing = @timed begin
            # Setup MS data directory
            MS_DATA_DIR = params.paths[:ms_data]
            if !isabspath(MS_DATA_DIR)
                MS_DATA_DIR = joinpath(@__DIR__, "../../", MS_DATA_DIR)
            end

            # Setup spectral library directory
            SPEC_LIB_DIR = params.paths[:library]
            if !isabspath(SPEC_LIB_DIR)
                SPEC_LIB_DIR = joinpath(@__DIR__, "../../", SPEC_LIB_DIR)
            end

            if !isdir(MS_DATA_DIR)
                @user_error "ms_data directory does not exist: " * MS_DATA_DIR
                return
            end

            # Find all Arrow files in MS data directory
            MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) 
                            for file in readdir(MS_DATA_DIR)
                            if isfile(joinpath(MS_DATA_DIR, file)) && 
                               match(r"\.arrow$", file) != nothing]

            if length(MS_TABLE_PATHS) <= 0
                @user_error "No .arrow files found in ms_data directory: " * MS_DATA_DIR
                return
            end

            nothing
        end
        timings["Parameter Loading"] = params_timing
        @user_info "Searching $(length(MS_TABLE_PATHS)) MS file$(length(MS_TABLE_PATHS) == 1 ? "" : "s") from $MS_DATA_DIR"

        # === Initialize spectral library and search context ===
        @user_info "Loading Spectral Library..."
        lib_timing = @timed begin
            SPEC_LIB = loadSpectralLibrary(SPEC_LIB_DIR, params)
            nothing
        end

        timings["Spectral Library Loading"] = lib_timing

        # Initialize Search Context
        @user_info "Initializing Search Context..."
        context_timing = @timed begin
            # Load isotope splines and initialize search context
            SEARCH_CONTEXT = initSearchContext(
                SPEC_LIB,
                parseIsoXML(isotope_spline_path()),
                ArrowTableReference(MS_TABLE_PATHS),
                Threads.nthreads(),
                250000 # Default temp array batch size
            )
            setDataOutDir!(SEARCH_CONTEXT, params.paths[:results])

            # Ensure temporary files are written to the results directory
            ENV["TMPDIR"] = params.paths[:results]

            # Write complete merged config (user + defaults) to output directory
            merged_config = params_to_dict(params)
            merged_json = JSON.json(merged_config, 4)  # Pretty print with 4-space indent
            write(joinpath(normpath(params.paths[:results]), "config.json"), merged_json)
            nothing
        end
        [rm(fpath) for fpath in readdir(getDataOutDir(SEARCH_CONTEXT), join=true) if endswith(fpath, ".tsv")]
        [rm(fpath) for fpath in readdir(getDataOutDir(SEARCH_CONTEXT), join=true) if endswith(fpath, ".arrow")]
        timings["Search Context Initialization"] = context_timing

        # === Execute search pipeline ===
        # Define search phases in order of execution
        searches = [
            ("Parameter Tuning", ParameterTuningSearch()),
            ("Quadrupole Tuning", QuadTuningSearch()),
            ("BitVec Calibration", BitVecCalibrationSearch()),
            ("Main Search", MainSearch()),
            ("Precursor Scoring", PrecursorScoringSearch()),
        ]
        if get(params.optimization.chromatogram_integration, :deconvolution_solver, "huber") == "huber"
            push!(searches, ("Huber Calibration", HuberTuningSearch()))
        end
        append!(searches, [
            ("Chromatogram Integration", IntegrateChromatogramSearch()),
            ("Protein Inference", ProteinInferenceSearch()),
            ("Protein Scoring", ProteinScoringSearch()),
            ("Quantification & Output", MaxLFQSearch())
        ])

        # Execute each search phase and record timing + peak RSS delta
        rss_deltas = Dict{String, Float64}()
        for (name, search) in searches
            @user_info "Executing $name..."
            rss_before = peak_rss()
            search_timing = @timed execute_search(search, SEARCH_CONTEXT, params)
            rss_after = peak_rss()
            timings[name] = search_timing
            rss_deltas[name] = (rss_after - rss_before) / (1024^3)
        end

        # === Generate performance report ===
        print_performance_report(timings, MS_TABLE_PATHS, SEARCH_CONTEXT, rss_deltas)
        
    catch e
        error_msg = try
            "$(typeof(e)): $(e.msg)"
        catch
            "$(typeof(e))"
        end
        @user_error "Search failed with error: $error_msg"
        @user_error "Stacktrace: $(stacktrace(catch_backtrace()))"
        rethrow(e)
    finally
        finalize_pioneer_logging(warnings_full_path; banner_title = "Search completed")
        # End-of-search runtime cleanup. Forces a full GC pass + flushes any
        # buffered C-side stdio (libomp/LightGBM emit warnings via printf).
        # Cheap (~tens of ms) and gives the runtime a clean state before the
        # next SearchDIA call in the same REPL — mitigates the libomp
        # re-entrance assertion that has historically aborted multi-call runs.
        GC.gc(true)
        Libc.flush_cstdio()
    end

    return nothing
end


"""
Helper function to print formatted performance metrics
"""
function print_performance_report(timings, ms_table_paths, search_context, rss_deltas=Dict{String,Float64}())
    # Header
    @user_print "\n" * repeat("=", 102)
    @user_print "DIA Search Performance Report"
    @user_print repeat("=", 102)

    # Detailed analysis
    @user_print "\nDetailed Step Analysis (in execution order):"
    @user_print repeat("-", 102)
    @user_print rpad("Step", 30) * " " *
            rpad("Time (s)", 12) * " " *
            rpad("Mem Alloc (GB)", 16) * " " *
            rpad("GC Time (s)", 12) * " " *
            rpad("GC %", 8) * " " *
            rpad("RSS Δ (GB)", 12)
    @user_print repeat("-", 102)

    # Calculate totals
    peak_memory = peak_rss()
    total_time = 0.0
    total_memory = 0
    total_gc = 0.0

    # Print step-by-step metrics in execution order
    # (insertion-ordered by SearchDIA driver; falls back to sorted keys for any extras).
    execution_order = [
        "Parameter Loading", "Spectral Library Loading", "Search Context Initialization",
        "Parameter Tuning", "Quadrupole Tuning", "BitVec Calibration",
        "Main Search", "Precursor Scoring", "Huber Calibration", "Chromatogram Integration",
        "Protein Inference", "Protein Scoring", "Quantification & Output",
    ]
    seen = Set{String}()
    sorted_steps = String[]
    for s in execution_order
        if haskey(timings, s)
            push!(sorted_steps, s)
            push!(seen, s)
        end
    end
    for s in sort(collect(keys(timings)))
        s in seen && continue
        push!(sorted_steps, s)
    end
    for step in sorted_steps
        timing = timings[step]
        time_s = timing.time
        mem_gb = timing.bytes / (1024^3)
        gc_s = timing.gctime
        gc_pct = (gc_s / time_s) * 100

        total_time += time_s
        total_memory += timing.bytes
        total_gc += gc_s

        rss_delta_str = haskey(rss_deltas, step) ? @sprintf("%+.2f", rss_deltas[step]) : "-"
        @user_print rpad(step, 30) * " " *
            rpad(@sprintf("%.2f", time_s), 12) * " " *
            rpad(@sprintf("%.2f", mem_gb), 16) * " " *
            rpad(@sprintf("%.2f", gc_s), 12) * " " *
            rpad(@sprintf("%.1f", gc_pct), 8) * " " *
            rpad(rss_delta_str, 12)
    end

    # Print summary statistics
    print_summary_statistics(
        total_time, total_memory, peak_memory, total_gc,
        length(timings), length(ms_table_paths),
        search_context, rss_deltas
    )
end

"""
Helper function to print summary statistics
"""
function print_summary_statistics(total_time, total_memory, peak_memory, total_gc,
                                n_steps, n_files, search_context, rss_deltas=Dict{String,Float64}())
    # Print totals
    @user_print repeat("-", 102)
    @user_print rpad("TOTAL", 30) * " " *
            rpad(@sprintf("%.2f", total_time), 12) * " " *
            rpad(@sprintf("%.2f", total_memory/(1024^3)), 16) * " " *
            rpad(@sprintf("%.2f", total_gc), 12) * " " *
            rpad(@sprintf("%.1f",(total_gc/total_time)*100), 8)

    # Memory summary
    # - "Allocated" sums GC bytes across stages (cumulative; can exceed system RAM
    #   because the same buffers get re-allocated and freed).
    # - "Working set (RSS)" is the OS-level high-water mark from Sys.maxrss().
    @user_print "\nMemory Usage Summary:"
    @user_print repeat("-", 102)
    current_mem = Sys.total_memory() / 1024^3
    @user_print "Total alloc (cumulative):   $(round(total_memory/1024^3, digits=2)) GB"
    @user_print "Working set (peak RSS):     $(round(peak_memory/1024^3, digits=2)) GB"
    @user_print "System RAM:                 $(round(current_mem, digits=2)) GB"
    if !isempty(rss_deltas)
        peak_step = argmax(rss_deltas)
        @user_print "Largest RSS jump:           +$(@sprintf("%.2f", rss_deltas[peak_step])) GB during $peak_step"
    end

    # Runtime summary
    @user_print "\nRuntime Summary:"
    @user_print repeat("-", 102)
    @user_print "Total Runtime: $(round(total_time/60, digits=2)) min ($(round(total_time, digits=2)) s)"
    @user_print "Per Raw File:  $(round(total_time/n_files, digits=2)) s"

    # === Final results block — count rows in the headline output files ===
    out_dir = getDataOutDir(search_context)
    prec_long = joinpath(out_dir, "precursors_long.arrow")
    pg_long   = joinpath(out_dir, "protein_groups_long.arrow")
    n_prec_rows = isfile(prec_long) ? length(Arrow.Table(prec_long)[1]) : nothing
    n_pg_rows   = isfile(pg_long)   ? length(Arrow.Table(pg_long)[1])   : nothing
    if n_prec_rows !== nothing || n_pg_rows !== nothing
        @user_print "\nResults:"
        @user_print repeat("-", 102)
        n_prec_rows !== nothing && @user_print rpad("Precursor rows (long):", 28) * "$n_prec_rows"
        n_pg_rows   !== nothing && @user_print rpad("Protein-group rows (long):", 28) * "$n_pg_rows"
        @user_print rpad("Output directory:", 28) * "$out_dir"
    end
    @user_print "\n" * repeat("=", 102)
end
