"""
    SearchDIA(params_path::String)

Main entry point for DIA (Data-Independent Acquisition) search workflow.
Executes a series of `SearchMethod`s and generates performance metrics.

The workflow consists of:
1. Parameter loading and validation
2. Spectral library initialization
3. Multiple search phases (parameter tuning, NCE tuning, etc.)
4. Performance reporting

Parameters:
- params_path: Path to JSON configuration file containing search parameters

Output:
- Generates a log file in the results directory
- Executes all search phases
- Reports timing and memory usage statistics
"""
function SearchDIA(params_path::String)
    # === Initialize logging ===
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

        # === Initialize performance tracking ===
        timings = Dict{String, NamedTuple{(:value, :time, :bytes, :gctime, :gcstats),Tuple{Nothing, Float64,Int64,Float64,Base.GC_Diff}}}()
        
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

            # Find all Arrow files in MS data directory
            MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) 
                            for file in readdir(MS_DATA_DIR)
                            if isfile(joinpath(MS_DATA_DIR, file)) && 
                               match(r"\.arrow$", file) != nothing]
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
                parseIsoXML(joinpath(@__DIR__,"../../data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml")),
                ArrowTableReference(MS_TABLE_PATHS),
                Threads.nthreads(),
                250000 # Default temp array batch size 
            )
            setDataOutDir!(SEARCH_CONTEXT, params.paths[:results])
            nothing
        end
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
        print_performance_report(timings, MS_TABLE_PATHS, dual_println)
    end
    return nothing
end


"""
Helper function to print formatted performance metrics
"""
function print_performance_report(timings, ms_table_paths, println_func)
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
        total_time, total_memory, total_gc,
        length(timings), length(ms_table_paths),
        println_func
    )
end

"""
Helper function to print summary statistics
"""
function print_summary_statistics(total_time, total_memory, total_gc, 
                                n_steps, n_files, println_func)
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
    println_func("Total Available Memory: $(round(current_mem, digits=2)) GB")
    
    # Runtime summary
    println_func("\nRuntime Summary:")
    println_func(repeat("-", 90))
    println_func("Total Runtime: $(round(total_time/60, digits=2)) minutes")
    println_func("Average Runtime per Step: $(round(total_time/n_steps, digits=2)) seconds")
    println_func("Average Runtime per Raw File: $(round(total_time/n_files, digits=2)) seconds")

    println_func("\n" * repeat("=", 90))
end
