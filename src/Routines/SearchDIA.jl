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
    
    # Set debug console level from parameters (default to 0)
    if haskey(params.logging, :debug_console_level)
        DEBUG_CONSOLE_LEVEL[] = params.logging.debug_console_level
    else
        DEBUG_CONSOLE_LEVEL[] = 0  # Default: no debug on console
    end
    
    # Configure max log message bytes (clamped with ENV override)
    local max_bytes_val::Int = 4096
    if haskey(params.logging, :max_message_bytes)
        try
            max_bytes_val = Int(params.logging.max_message_bytes)
        catch
            max_bytes_val = 4096
        end
    end
    if haskey(ENV, "PIONEER_MAX_LOG_MSG_BYTES")
        env_val = tryparse(Int, ENV["PIONEER_MAX_LOG_MSG_BYTES"])
        if env_val !== nothing
            max_bytes_val = env_val
        end
    end
    max_bytes_val = clamp(max_bytes_val, 1024, 1048576)
    Pioneer.MAX_LOG_MSG_BYTES[] = max_bytes_val
    
    # Open FOUR log files
    essential_path = joinpath(params.paths[:results], "pioneer_search_report.txt")
    console_path = joinpath(params.paths[:results], "pioneer_search_log.log")
    debug_path = joinpath(params.paths[:results], "pioneer_search_debug.log")
    warnings_path = joinpath(params.paths[:results], "pioneer_warnings.log")
    
    Pioneer.ESSENTIAL_FILE[] = open(essential_path, "w")
    Pioneer.CONSOLE_FILE[] = open(console_path, "w")
    Pioneer.DEBUG_FILE[] = open(debug_path, "w")
    Pioneer.WARNINGS_FILE[] = open(warnings_path, "w")
    
    # Get Pioneer version from Project.toml
    project_toml = joinpath(pkgdir(Pioneer), "Project.toml")
    pioneer_version = "unknown"
    if isfile(project_toml)
        content = read(project_toml, String)
        version_match = match(r"version\s*=\s*\"([^\"]+)\"", content)
        if version_match !== nothing
            pioneer_version = version_match[1]
        end
    end
    
    # Essential file - clean header (like dual_println)
    essential_header = [
        "=" ^ 90,
        "Pioneer Search Log",
        "Version: $pioneer_version",
        "Started: $(Dates.now())",
        "Output Directory: $(params.paths[:results])",
        "=" ^ 90,
        ""
    ]
    for line in essential_header
        println(Pioneer.ESSENTIAL_FILE[], line)
    end
    
    # Console file - same header
    for line in essential_header
        println(Pioneer.CONSOLE_FILE[], line)
    end
    
    # Debug file - detailed header
    debug_header = [
        "=" ^ 90,
        "Pioneer Debug Log (Full Trace)",
        "Version: $pioneer_version",
        "Started: $(Dates.now())",
        "Output Directory: $(params.paths[:results])",
        "Julia Version: $(VERSION)",
        "Threads: $(Threads.nthreads())",
        "=" ^ 90,
        ""
    ]
    for line in debug_header
        println(Pioneer.DEBUG_FILE[], line)
    end
    
    # Warnings file - simple header
    warnings_header = [
        "=" ^ 90,
        "Pioneer Warnings Log",
        "Started: $(Dates.now())",
        "=" ^ 90,
        ""
    ]
    for line in warnings_header
        println(Pioneer.WARNINGS_FILE[], line)
    end
    
    # Log initialization message
    @user_info "Pioneer logging system initialized"
    
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

            # Disable match-between-runs if only a single run is provided
            if length(MS_TABLE_PATHS) == 1 && params.global_settings.match_between_runs
                @user_warn "Only one run detected; disabling match_between_runs (MBR)."

                params_dict = JSON.parse(params_string)
                params_dict["global"]["match_between_runs"] = false
                params_string = JSON.json(params_dict, 4)

                updated_global = (; params.global_settings..., match_between_runs=false)
                params = PioneerParameters(
                    updated_global,
                    params.parameter_tuning,
                    params.first_search,
                    params.quant_search,
                    params.acquisition,
                    params.rt_alignment,
                    params.optimization,
                    params.protein_inference,
                    params.maxLFQ,
                    params.output,
                    params.logging,
                    params.paths,
                )
            end
            
            nothing
        end
        timings["Parameter Loading"] = params_timing

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
            ("NCE Tuning", NceTuningSearch()),
            ("Quadrupole Tuning", QuadTuningSearch()),
            ("First Pass Search", FirstPassSearch()),
            ("Huber Tuning", HuberTuningSearch()),
            ("Second Pass Search", SecondPassSearch()),
            ("Scoring", ScoringSearch()),
            ("Chromatogram Integration", IntegrateChromatogramSearch()),
            ("MaxLFQ", MaxLFQSearch())
        ]

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
        # Count warnings if any exist
        warning_count = 0
        warnings_full_path = ""  # Store full path for display
        if Pioneer.WARNINGS_FILE[] !== nothing
            # Close and get full path
            close(Pioneer.WARNINGS_FILE[])
            warnings_full_path = abspath(warnings_path)  # Get absolute path
            if isfile(warnings_path)
                warning_count = max(0, countlines(warnings_path) - 5)  # Subtract header lines
            end
            Pioneer.WARNINGS_FILE[] = nothing
        end
        
        # Close all four files with appropriate footers
        if Pioneer.ESSENTIAL_FILE[] !== nothing
            footer = ["", "=" ^ 90]
            if warning_count > 0
                push!(footer, "⚠️  $warning_count warnings were generated during search")
            end
            push!(footer, "Search completed at: $(Dates.now())")
            push!(footer, "=" ^ 90)
            for line in footer
                println(Pioneer.ESSENTIAL_FILE[], line)
            end
            close(Pioneer.ESSENTIAL_FILE[])
            Pioneer.ESSENTIAL_FILE[] = nothing
        end
        
        if Pioneer.CONSOLE_FILE[] !== nothing
            footer = ["", "=" ^ 90]
            if warning_count > 0
                push!(footer, "⚠️  $warning_count warnings were generated during search")
            end
            push!(footer, "Search completed at: $(Dates.now())")
            push!(footer, "=" ^ 90)
            for line in footer
                println(Pioneer.CONSOLE_FILE[], line)
            end
            close(Pioneer.CONSOLE_FILE[])
            Pioneer.CONSOLE_FILE[] = nothing
        end
        
        if Pioneer.DEBUG_FILE[] !== nothing
            footer = ["", "=" ^ 90, "Search completed at: $(Dates.now())", "=" ^ 90]
            for line in footer
                println(Pioneer.DEBUG_FILE[], line)
            end
            close(Pioneer.DEBUG_FILE[])
            Pioneer.DEBUG_FILE[] = nothing
        end
        
        # Print warning count to console if there were any
        if warning_count > 0
            printstyled("┌ ", color=:yellow)
            printstyled("Warning:", bold=true, color=:yellow)
            println(" $warning_count warnings were generated during search - see $warnings_full_path")
        end
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
    @user_print "\nDetailed Step Analysis:"
    @user_print repeat("-", 102)
    @user_print rpad("Step", 30) * " " *
            rpad("Time (s)", 12) * " " *
            rpad("Memory (GB)", 12) * " " *
            rpad("GC Time (s)", 12) * " " *
            rpad("GC %", 12) * " " *
            rpad("Peak RSS +", 12)
    @user_print repeat("-", 102)

    # Calculate totals
    peak_memory = peak_rss()
    total_time = 0.0
    total_memory = 0
    total_gc = 0.0

    # Print step-by-step metrics
    sorted_steps = sort(collect(keys(timings)), by=x->timings[x][:time])
    for step in sorted_steps
        timing = timings[step]
        time_s = timing.time
        mem_gb = timing.bytes / (1024^3)
        gc_s = timing.gctime
        gc_pct = (gc_s / time_s) * 100

        total_time += time_s
        total_memory += timing.bytes
        total_gc += gc_s

        rss_delta_str = haskey(rss_deltas, step) ? @sprintf("%.2f GB", rss_deltas[step]) : "-"
        @user_print rpad(step, 30) * " " *
            rpad(@sprintf("%.2f", time_s), 12) * " " *
            rpad(@sprintf("%.2f", mem_gb), 12) * " " *
            rpad(@sprintf("%.2f", gc_s), 12) * " " *
            rpad(@sprintf("%.1f", gc_pct), 12) * " " *
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
            rpad(@sprintf("%.2f", total_memory/(1024^3)), 12) * " " *
            rpad(@sprintf("%.2f", total_gc), 12) * " " *
            rpad(@sprintf("%.1f",(total_gc/total_time)*100), 12)

    # Memory summary
    @user_print "\nMemory Usage Summary:"
    @user_print repeat("-", 102)
    current_mem = Sys.total_memory() / 1024^3
    @user_print "Total Memory Allocated: $(round(total_memory/1024^3, digits=2)) GB"
    @user_print "Peak  Memory Allocated: $(round(peak_memory/1024^3, digits=2)) GB"
    @user_print "Total Available Memory: $(round(current_mem, digits=2)) GB"
    if !isempty(rss_deltas)
        peak_step = argmax(rss_deltas)
        @user_print "Peak RSS Step: $peak_step (+$(@sprintf("%.2f", rss_deltas[peak_step])) GB)"
    end

    # Runtime summary
    @user_print "\nRuntime Summary:"
    @user_print repeat("-", 102)
    @user_print "Total Runtime: $(round(total_time/60, digits=2)) minutes"
    @user_print "Average Runtime per Step: $(round(total_time/n_steps, digits=2)) seconds"
    @user_print "Average Runtime per Raw File: $(round(total_time/n_files, digits=2)) seconds"
    @user_print "\n" * repeat("=", 102)
    @user_print "Huber δ: " * string(getHuberDelta(search_context))
end
