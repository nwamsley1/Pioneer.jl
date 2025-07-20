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
function main_SearchDIA()::Cint
    try
        SearchDIA(ARGS[1])
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

"""
    asset_path(parts...)
    Return the path to a bundled asset. The function first checks the
    compile-time `data/` directory and falls back to a path relative to the
    installed executable.
"""
function asset_path(parts...)
    compile_dir = joinpath(@__DIR__, "..", "..", "data", parts...)
    if ispath(compile_dir)
        return compile_dir
    end
    exe_dir = abspath(dirname(PROGRAM_FILE))
    println("ISO PATH:", joinpath(exe_dir, "..", "data", parts...))
    return joinpath(exe_dir, "..", "data", parts...)
end

"""
    Locate the isotope spline XML file bundled with the application.
"""
isotope_spline_path() = asset_path("IsotopeSplines", "IsotopeSplines_10kDa_21isotopes-1.xml")


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
    log_path = joinpath(params.paths[:results], "pioneer_search_log.txt")
    mkpath(dirname(log_path))  # Ensure directory exists


    open(log_path, "w") do log_file
        # Utility functions for dual console/file output
        function dual_println(args...)
            println(stdout, args...)  # Write to console
            println(log_file, args...)  # Write to file
        end
        
        function dual_print(args...)
            print(stdout, args...)  # Write to console
            print(log_file, args...)  # Write to file
        end

        dual_println("\n", repeat("=", 90))
        dual_println("Sarting SearchDIA")
        dual_println(repeat("=", 90))
        dual_println("\nStarting search at: ", Dates.now())
        dual_println("Output directory: ", params.paths[:results])
        

        # === Initialize performance tracking ===
        timings = Dict{String, Any}()
        
        # === Load and validate data paths ===
        @info "Loading Parameters..."
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
                @error "ms_data directory does not exist: " * MS_DATA_DIR
                return
            end

            # Find all Arrow files in MS data directory
            MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) 
                            for file in readdir(MS_DATA_DIR)
                            if isfile(joinpath(MS_DATA_DIR, file)) && 
                               match(r"\.arrow$", file) != nothing]

            if length(MS_TABLE_PATHS) <= 0
                @error "No .arrow files found in ms_data directory: " * MS_DATA_DIR
                return
            end

            # Disable match-between-runs if only a single run is provided
            if length(MS_TABLE_PATHS) == 1 && params.global_settings.match_between_runs
                @warn "Only one run detected; disabling match_between_runs (MBR)."

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
                    params.paths,
                )
            end
            
            nothing
        end
        timings["Parameter Loading"] = params_timing

        # === Initialize spectral library and search context ===
        @info "Loading Spectral Library..."
        lib_timing = @timed begin
            SPEC_LIB = loadSpectralLibrary(SPEC_LIB_DIR, params)
            nothing
        end

        timings["Spectral Library Loading"] = lib_timing

        # Initialize Search Context
        @info "Initializing Search Context..."
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
           
            write( joinpath(normpath(params.paths[:results]), "config.json"), params_string)
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

        # Execute each search phase and record timing
        for (name, search) in searches
            @info "Executing $name..."
            search_timing = @timed execute_search(search, SEARCH_CONTEXT, params)
            timings[name] = search_timing
        end

        # === Generate performance report ===
        print_performance_report(timings, MS_TABLE_PATHS, SEARCH_CONTEXT, dual_println)
    end
    return nothing
end


"""
Helper function to print formatted performance metrics
"""
function print_performance_report(timings, ms_table_paths, search_context, println_func)
    # Header
    println_func("\n" * repeat("=", 90))
    println_func("DIA Search Performance Report")
    println_func(repeat("=", 90))

    # Detailed analysis
    println_func("\nDetailed Step Analysis:")
    println_func(repeat("-", 90))
    println_func(rpad("Step", 30), " ", 
            rpad("Time (s)", 12), " ",
            rpad("Memory (GB)", 12), " ", 
            rpad("GC Time (s)", 12), " ",
            rpad("GC %", 12))
    println_func(repeat("-", 90))

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

        println_func(rpad(step, 30), " ",
            rpad(@sprintf("%.2f", time_s), 12), " ",
            rpad(@sprintf("%.2f", mem_gb), 12), " ",
            rpad(@sprintf("%.2f", gc_s), 12), " ",
            rpad(@sprintf("%.1f", gc_pct), 12))
    end

    # Print summary statistics
    print_summary_statistics(
        total_time, total_memory, peak_memory, total_gc,
        length(timings), length(ms_table_paths),
        println_func, search_context
    )
end

"""
Helper function to print summary statistics
"""
function print_summary_statistics(total_time, total_memory, peak_memory, total_gc, 
                                n_steps, n_files, println_func, search_context)
    # Print totals
    println_func(repeat("-", 90))
    println_func(rpad("TOTAL", 30), " ",
            rpad(@sprintf("%.2f", total_time), 12), " ",
            rpad(@sprintf("%.2f", total_memory/(1024^3)), 12), " ",
            rpad(@sprintf("%.2f", total_gc), 12), " ",
            rpad(@sprintf("%.1f",(total_gc/total_time)*100), 12))

    # Memory summary
    println_func("\nMemory Usage Summary:")
    println_func(repeat("-", 90))
    current_mem = Sys.total_memory() / 1024^3
    println_func("Total Memory Allocated: $(round(total_memory/1024^3, digits=2)) GB")
    println_func("Peak  Memory Allocated: $(round(peak_memory/1024^3, digits=2)) GB")
    println_func("Total Available Memory: $(round(current_mem, digits=2)) GB")
    
    # Runtime summary
    println_func("\nRuntime Summary:")
    println_func(repeat("-", 90))
    println_func("Total Runtime: $(round(total_time/60, digits=2)) minutes")
    println_func("Average Runtime per Step: $(round(total_time/n_steps, digits=2)) seconds")
    println_func("Average Runtime per Raw File: $(round(total_time/n_files, digits=2)) seconds")
    println_func("\n" * repeat("=", 90))
    println_func("Huber δ: ", getHuberDelta(search_context))
end
