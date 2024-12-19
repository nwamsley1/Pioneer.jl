function SearchDIA(params_path::String)
    # Create log file in the results directory
    params = parse_pioneer_parameters(params_path)
    log_path = joinpath(params.paths[:results], "pioneer_search_log.txt")
    mkpath(dirname(log_path))  # Ensure directory exists
    
    # Open log file and redirect output
    open(log_path, "w") do log_file
        # Create a function to write to both console and file
        function dual_println(args...)
            println(stdout, args...)  # Write to console
            println(log_file, args...)  # Write to file
        end
        
        function dual_print(args...)
            print(stdout, args...)  # Write to console
            print(log_file, args...)  # Write to file
        end

        timings = Dict{String, NamedTuple{(:value, :time, :bytes, :gctime, :gcstats),Tuple{Nothing, Float64,Int64,Float64,Base.GC_Diff}}}()
        
        # Load Parameters
        @info "Loading Parameters..."
        params_timing = @timed begin
            if !isabspath(params_path)
                params_path = joinpath(@__DIR__, "../../", params_path)
            end
            MS_DATA_DIR = params.paths[:ms_data]
            SPEC_LIB_DIR = params.paths[:library]
            if !isabspath(SPEC_LIB_DIR)
                SPEC_LIB_DIR = joinpath(@__DIR__, "../../", SPEC_LIB_DIR)
            end

            if isabspath(MS_DATA_DIR)
                MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))]
            else
                MS_DATA_DIR = joinpath(@__DIR__, "../../", MS_DATA_DIR)
                MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))]
            end
            nothing
        end
        timings["Parameter Loading"] = params_timing

        # Load Spectral Library
        @info "Loading Spectral Library..."
        lib_timing = @timed begin
            SPEC_LIB = loadSpectralLibrary(SPEC_LIB_DIR, params)
            dual_println("c ", typeof(SPEC_LIB))
            nothing
        end
        timings["Spectral Library Loading"] = lib_timing

        # Initialize Search Context
        @info "Initializing Search Context..."
        context_timing = @timed begin
            SEARCH_CONTEXT = initSearchContext(
                SPEC_LIB,
                parseIsoXML(joinpath(@__DIR__,"../../data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml")),
                ArrowTableReference(MS_TABLE_PATHS),
                Threads.nthreads(),
                250000
            )
            setDataOutDir!(SEARCH_CONTEXT, params.paths[:results])
            nothing
        end
        timings["Search Context Initialization"] = context_timing

        # Define and execute searches
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

        for (name, search) in searches
            @info "Executing $name..."
            search_timing = @timed execute_search(search, SEARCH_CONTEXT, params)
            timings[name] = search_timing
        end

        # Print performance report
        dual_println("\n" * repeat("=", 90))
        dual_println("DIA Search Performance Report")
        dual_println(repeat("=", 90))

        # Step-by-step analysis
        dual_println("\nDetailed Step Analysis:")
        dual_println(repeat("-", 90))
        dual_println(rpad("Step", 30), " ", 
                rpad("Time (s)", 12), " ",
                rpad("Memory (GB)", 12), " ", 
                rpad("GC Time (s)", 12), " ",
                rpad("GC %", 12))
        dual_println(repeat("-", 90))

        total_time = 0.0
        total_memory = 0
        total_gc = 0.0

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

            dual_println(rpad(step, 30), " ",
                rpad(@sprintf("%.2f", time_s), 12), " ",
                rpad(@sprintf("%.2f", mem_gb), 12), " ",
                rpad(@sprintf("%.2f", gc_s), 12), " ",
                rpad(@sprintf("%.1f", gc_pct), 12))
        end

        dual_println(repeat("-", 90))
        dual_println(rpad("TOTAL", 30), " ",
                rpad(@sprintf("%.2f", total_time), 12), " ",
                rpad(@sprintf("%.2f", total_memory/(1024^3)), 12), " ",
                rpad(@sprintf("%.2f", total_gc), 12), " ",
                rpad(@sprintf("%.1f",(total_gc/total_time)*100), 12))

        # Memory Usage Summary
        dual_println("\nMemory Usage Summary:")
        dual_println(repeat("-", 90))
        current_mem = Sys.total_memory() / 1024^3
        dual_println("Total Memory Allocated: $(round(total_memory/1024^3, digits=2)) GB")
        dual_println("Total Available Memory: $(round(current_mem, digits=2)) GB")
        
        # Runtime Summary
        dual_println("\nRuntime Summary:")
        dual_println(repeat("-", 90))
        dual_println("Total Runtime: $(round(total_time/60, digits=2)) minutes")
        dual_println("Average Runtime per Step: $(round(total_time/length(timings), digits=2)) seconds")
        dual_println("Average Runtime per Raw File: $(round(total_time/length(MS_TABLE_PATHS), digits=2)) seconds")

        dual_println("\n" * repeat("=", 90))
    end

    return nothing
end