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

# src/BuildSpecLib.jl

# Entry point for PackageCompiler
function main_BuildSpecLib(argv=ARGS)::Cint
    
    settings = ArgParseSettings(; autofix_names = true)
    @add_arg_table! settings begin
        "params_path"
            help = "Path to BuildSpecLib parameters JSON file"
            arg_type = String
    end
    parsed_args = parse_args(argv, settings; as_symbols = true)
    
    try
        BuildSpecLib(parsed_args[:params_path])
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

"""
    BuildSpecLib(params_path::String)

Main function to build a spectral library from parameters. Executes a series of steps:
1. Parameter validation and directory setup
2. Fragment bound detection
3. Retention time prediction (optional)
4. Fragment prediction (optional)
5. Library index building

Parameters:
- params_path: Path to JSON configuration file containing library building parameters

Output:
- Generates a spectral library in the specified output directory
- Creates a detailed log file with timing and performance metrics
- Returns nothing
"""
function BuildSpecLib(params_path::String)
    # Clean up any old file handlers in case the program crashed
    GC.gc()
    Random.seed!(1844)
    # Initialize timing dictionary for performance tracking
    #params_path = "/Users/n.t.wamsley/RIS_temp/koina_testing/config.json"
    timings = Dict{String, Any}()
    
    # Read and validate parameters
    params_timing = @timed begin
        params_string = read(params_path, String)
        params = check_params_bsp(params_string)

        # Get library directory (already has .poin extension added in check_params)
        lib_dir = params["_lib_dir"]
        mkpath(lib_dir)

        # Setup logging
        log_path = joinpath(lib_dir, "build_log.txt")
        params_out_path = joinpath(lib_dir, "config.json")

        # Write complete merged parameters to config.json (not just user input)
        params_json = JSON.json(params, 2)  # Pretty-print with 2-space indent
        write(params_out_path, params_json)
        nothing
    end
    timings["Parameter Loading"] = params_timing

    # Open log file for writing
    open(log_path, "w") do log_file

        # Utility functions for dual console/file output
        function dual_println(args...)
            println(stdout, args...)
            println(log_file, args...)
        end
        
        function dual_print(args...)
            print(stdout, args...)
            print(log_file, args...)
        end

    
        dual_println("\n", repeat("=", 90))
        dual_println("Spectral Library Building Process")
        dual_println(repeat("=", 90))
        dual_println("\nStarting library build at: ", Dates.now())
        dual_println("Output directory: ", lib_dir)
        
        # Create temporary chronologer directory
        setup_timing = @timed begin
            chronologer_dir = joinpath(lib_dir, "chronologer_temp")
            mkpath(chronologer_dir)
            chronologer_in_path = joinpath(chronologer_dir, "precursors_for_chronologer.arrow")
            chronologer_out_path = joinpath(chronologer_dir, "precursors_for_chronologer_rt.arrow")
            
            # Format parameters into named tuple
            _params = (
                fasta_digest_params = Dict{String, Any}(k => v for (k, v) in params["fasta_digest_params"]),
                nce_params = Dict{String, Any}(k => v for (k, v) in params["nce_params"]),
                library_params = Dict{String, Any}(k => v for (k, v) in params["library_params"])
            )
            nothing
        end
        timings["Directory Setup"] = setup_timing

        # Get fragment bounds
        dual_println("\nDetecting fragment bounds...")
        bounds_timing = @timed begin
            # calibration_raw_file is optional - use empty string if not provided
            calibration_file = get(_params, "calibration_raw_file", "")
            frag_bounds, prec_mz_min, prec_mz_max = get_fragment_bounds(
                _params.library_params["auto_detect_frag_bounds"],
                calibration_file,
                (Float32(_params.library_params["frag_mz_min"]), Float32(_params.library_params["frag_mz_max"])),
                (Float32(_params.library_params["prec_mz_min"]), Float32(_params.library_params["prec_mz_max"]))
            )
            dual_println("Fragment m/z range: ", frag_bounds)
            dual_println("Precursor m/z range: ", (prec_mz_min, prec_mz_max))
            nothing
        end
        timings["Fragment Bound Detection"] = bounds_timing

        N_PRECURSORS = 0
        N_FRAGMENTS = 0
        N_TARGETS = 0
        N_DECOYS = 0
        if params["predict_fragments"]
            # Fragment prediction workflow
            dual_println("\nStarting fragment prediction workflow...")
            
            # Validate prediction model
            model_timing = @timed begin
                prediction_model = params["library_params"]["prediction_model"]
                PREDICTION_MODELS = collect(keys(KOINA_URLS))
                if !(prediction_model ∈ PREDICTION_MODELS)
                    error("Invalid prediction model: $prediction_model. Valid options: $(join(PREDICTION_MODELS, ", "))")
                end
                dual_println("Using prediction model: $prediction_model")

                # Get annotation type and model type
                frag_annotation_type = MODEL_CONFIGS[prediction_model].annotation_type
                koina_model_type = MODEL_CONFIGS[prediction_model].model_type
                nothing
            end
            timings["Model Validation"] = model_timing

            # Collision energy interpolation
            ce_timing = @timed begin
                mz_to_ev_interp = missing
                if occursin("unispec", prediction_model)
                    try
                        mz_to_ev_interp = get_mz_to_ev_interp(
                            _params.calibration_raw_file,
                            lib_dir
                        )
                        dual_println("Successfully created collision energy interpolator")
                    catch
                        @user_warn "Could not estimate mz to ev conversion. Using default NCE."
                        dual_println("Warning: Using default NCE (collision energy interpolation failed)")
                    end
                end

                instrument_type = ismissing(mz_to_ev_interp) ? 
                                _params.library_params["instrument_type"] : 
                                "NONE"
                nothing
            end
            timings["Collision Energy Setup"] = ce_timing

            # Chronologer prediction workflow
            dual_println("\nPreparing chronologer workflow...")
            chrono_prep_timing = @timed begin
                prepare_chronologer_input(params,
                                        mz_to_ev_interp,
                                        prec_mz_min,
                                        prec_mz_max,
                                        chronologer_in_path,
                                        joinpath(lib_dir, "proteins_table.arrow"))
                nothing
            end
            timings["Chronologer Preparation"] = chrono_prep_timing

            dual_println("Predicting retention times...")
            rt_timing = @timed begin
                predict_retention_times(chronologer_in_path, chronologer_out_path)
                nothing
            end
            timings["Retention Time Prediction"] = rt_timing
            # Parse results and prepare for fragment prediction
            parse_timing = @timed begin
                iso_mod_to_mass = Dict{String, Float32}()
                precursors_arrow_path = parse_chronologer_output(
                    chronologer_out_path,
                    lib_dir,
                    Dict{String,Int8}(),
                    iso_mod_to_mass,
                    params["isotope_mod_groups"],
                    Float32(_params.library_params["rt_bin_tol"])
                )
                #println("precursors_arrow_path $precursors_arrow_path")
                # Cleanup temporary files
                GC.gc()
                rm(chronologer_in_path, force=true)
                rm(chronologer_out_path, force=true)
                dir, filename = splitdir(precursors_arrow_path)
                raw_fragments_arrow_path = joinpath(dir, "raw_fragments.arrow")
                rm(raw_fragments_arrow_path, force=true)
                nothing
            end
            timings["Chronologer Output Processing"] = parse_timing

            # Fragment prediction
            dual_println("\nPredicting fragment ion intensities...")
            frag_predict_timing = @timed begin
                predict_fragments(
                    precursors_arrow_path,
                    raw_fragments_arrow_path,
                    koina_model_type,
                    instrument_type,
                    params["max_koina_requests"],
                    params["max_koina_batch"],
                    prediction_model
                )
                nothing
            end
            timings["Fragment Prediction"] = frag_predict_timing

            # Process predictions
            process_timing = @timed begin
                # Load tables
                precursors_table = Arrow.Table(precursors_arrow_path)
                fragments_table = Arrow.Table(raw_fragments_arrow_path)
                #Record the spline knots 
                try
                    spl_knots = copy(fragments_table[:knot_vector][1])
                    jldsave(
                        joinpath(lib_dir, "spline_knots.jld2");
                        spl_knots
                    )
                    ion_dictionary = get_altimeter_ion_dict(asset_path("ion_dictionary.txt"))

                    parse_altimeter_fragments(
                        precursors_table,
                        fragments_table,
                        frag_annotation_type,
                        ion_dictionary,
                        10000,
                        asset_path("immonium.txt"),
                        lib_dir,
                        Dict{String, Int8}(),
                        iso_mod_to_mass,
                        koina_model_type
                    )

                catch
                    #println("No spline knots. static library")

                    # Process ion annotations
                    ion_annotation_set = get_ion_annotation_set(fragments_table[:annotation])
                    frag_name_to_idx = Dict(ion => UInt16(i) for (i, ion) in enumerate(ion_annotation_set))

                    ion_annotation_dict = parse_koina_fragments(
                        precursors_table,
                        fragments_table,
                        frag_annotation_type,
                        ion_annotation_set,
                        frag_name_to_idx,
                        10000,
                        asset_path("immonium.txt"),
                        lib_dir,
                        Dict{String, Int8}(),
                        iso_mod_to_mass,
                        koina_model_type
                    )
                end
                

                # Process precursor table
                N_FRAGMENTS = length(fragments_table[:mz])
                N_PRECURSORS = length(precursors_table[:mz])
                N_DECOYS  = sum(precursors_table[:decoy])
                N_TARGETS = N_PRECURSORS - N_DECOYS
                fragments_table = nothing
                precursors_table = DataFrame(precursors_table)
                rename!(precursors_table, [
                    :accession_number => :accession_numbers,
                    :precursor_charge => :prec_charge,
                    :decoy => :is_decoy,
                    :mods => :structural_mods
                ])

                # Convert types
                precursors_table[!, :missed_cleavages] = UInt8.(precursors_table[!, :missed_cleavages])
                precursors_table[!, :prec_charge] = UInt8.(precursors_table[!, :prec_charge])
                precursors_table[!, :mz] = Float32.(precursors_table[!, :mz])
                precursors_table[!, :irt] = Float32.(precursors_table[!, :irt])
                precursors_table[!, :start_idx] = UInt32.(precursors_table[!, :start_idx])

                # Save processed precursor table
                println("   Before add_pair_indices!: $(nrow(precursors_table)) precursors")
                println("   Unique pair_ids available: $(length(unique(precursors_table.pair_id)))")
                add_pair_indices!(precursors_table)  # Add partner indices AFTER all sorting is complete
                partner_col = precursors_table.partner_precursor_idx
                max_partner = all(ismissing, partner_col) ? 0 : Int64(maximum(skipmissing(partner_col)))
                println("   After add_pair_indices!: max partner_idx = $max_partner (table size: $(nrow(precursors_table)))")
                println("   Valid indices: $(max_partner <= nrow(precursors_table) ? "✅ YES" : "❌ NO")")
                
                # Add entrapment target indices if entrapment_pair_id column exists AND entrapment is enabled
                entrapment_r = get(params["fasta_digest_params"], "entrapment_r", 0)
                if hasproperty(precursors_table, :entrapment_pair_id) && entrapment_r > 0
                    println("   Adding entrapment target indices...")
                    add_entrapment_indices!(precursors_table)
                    ent_col = precursors_table.entrapment_target_idx
                    max_entrap_target = all(ismissing, ent_col) ? 0 : Int64(maximum(skipmissing(ent_col)))
                    println("   After add_entrapment_indices!: max entrapment_target_idx = $max_entrap_target")
                    println("   Entrapment indices valid: $(max_entrap_target <= nrow(precursors_table) ? "✅ YES" : "❌ NO")")
                end
                Arrow.write(
                    joinpath(lib_dir, "precursors_table.arrow"),
                    precursors_table
                )

                GC.gc()
                nothing
            end
            timings["Prediction Processing"] = process_timing
        end

        # Verify required files
        verify_timing = @timed begin
            required_files = ["fragments_table.arrow", "prec_to_frag.arrow", "precursors_table.arrow", "proteins_table.arrow"]
            if !all(isfile.(joinpath.(lib_dir, required_files)))
                error("Missing required files in $lib_dir. Try running with predict_fragments=true")
            end
            nothing
        end
        timings["File Verification"] = verify_timing

        # Build final indices
        dual_println("\nBuilding final library indices...")
        index_timing = @timed begin
            buildPionLib(
                lib_dir,
                UInt8(_params.library_params["y_start_index"]),
                UInt8(_params.library_params["y_start"]),
                UInt8(_params.library_params["b_start_index"]),
                UInt8(_params.library_params["b_start"]),
                _params.library_params["include_p_index"],
                _params.library_params["include_p"],
                _params.library_params["include_isotope"],
                _params.library_params["include_immonium"],
                _params.library_params["include_internal"],
                _params.library_params["include_neutral_diff"],
                UInt8(_params.library_params["max_frag_charge"]),
                UInt8(_params.library_params["max_frag_rank"]),
                Float32(_params.library_params["length_to_frag_count_multiple"]),
                Float32(_params.library_params["min_frag_intensity"]),
                UInt8.(_params.library_params["rank_to_score"]),
                frag_bounds,
                Float32(_params.library_params["frag_bin_tol_ppm"]),
                Float32(_params.library_params["rt_bin_tol"]),
                koina_model_type
            )          

            nothing
        end
        timings["Index Building"] = index_timing

        # Print performance report
        print_performance_report(timings, dual_println,
        N_PRECURSORS = N_PRECURSORS,
        N_FRAGMENTS = N_FRAGMENTS,
        N_TARGETS = N_TARGETS,
        N_DECOYS = N_DECOYS)
        
        dual_println("\nLibrary building completed at: ", Dates.now())
        dual_println("\nLibrary location: ", lib_dir)
    end

    # cleanup temp files
    GC.gc()
    rm(joinpath(lib_dir, "raw_fragments.arrow"), force=true)
    rm(joinpath(lib_dir,"fragments_table.arrow"), force=true);
    rm(joinpath(lib_dir,"prec_to_frag.arrow"), force=true);
    rm(joinpath(lib_dir,"precursors.arrow"), force=true);
    GC.gc()

    return nothing
end

"""
Helper function to print formatted performance metrics
"""
function print_performance_report(timings, println_func; kwargs...)
    # Header
    println_func("\n", repeat("=", 90))
    println_func("Library Building Performance Report")
    println_func(repeat("=", 90))

    println_func("\nDetailed Step Analysis:")

    println_func(repeat("-", 90))
    println_func(rpad("# Precursors", 12), " ",
                rpad("# Fragments", 12), " ",
                rpad("# Targets", 12), " ",
                rpad("# Decoys", 12), " "
    )
    println_func(repeat("-", 90))
    println_func(lpad(@sprintf("%.2f", kwargs[:N_PRECURSORS]), 12), " ",
                 lpad(@sprintf("%.2f", kwargs[:N_FRAGMENTS]), 12), " ",
                 lpad(@sprintf("%.2f", kwargs[:N_TARGETS]), 12), " ",
                 lpad(@sprintf("%.2f", kwargs[:N_DECOYS]), 12), " ",
    )
    # Print detailed analysis
    println_func("\nDetailed Step Analysis:")
    println_func(repeat("-", 90))
    println_func(rpad("Step", 40), " ",
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
    sorted_steps = sort(collect(keys(timings)), by=x->timings[x].time)
    for step in sorted_steps
        timing = timings[step]
        time_s = timing.time
        mem_gb = timing.bytes / (1024^3)
        gc_s = timing.gctime
        gc_pct = (gc_s / time_s) * 100

        total_time += time_s
        total_memory += timing.bytes
        total_gc += gc_s

        println_func(rpad(step, 40), " ",
                    rpad(@sprintf("%.2f", time_s), 12), " ",
                    rpad(@sprintf("%.2f", mem_gb), 12), " ",
                    rpad(@sprintf("%.2f", gc_s), 12), " ",
                    rpad(@sprintf("%.1f", gc_pct), 12))
    end

    # Print totals
    println_func(repeat("-", 90))
    println_func(rpad("TOTAL", 40), " ",
                rpad(@sprintf("%.2f", total_time), 12), " ",
                rpad(@sprintf("%.2f", total_memory/(1024^3)), 12), " ",
                rpad(@sprintf("%.2f", total_gc), 12), " ",
                rpad(@sprintf("%.1f", (total_gc/total_time)*100), 12))

    # Print summary statistics
    println_func("\nPerformance Summary:")
    println_func(repeat("-", 90))
    println_func("Total Runtime: $(round(total_time/60, digits=2)) minutes")
    println_func("Total Garbage Collection Time: $(round(total_gc, digits=2)) seconds")
    println_func("Number of Steps: $(length(timings))")
    println_func("Average Step Runtime: $(round(total_time/length(timings), digits=2)) seconds")
    
    # Memory statistics
    println_func("\nMemory Usage Summary:")
    println_func(repeat("-", 90))
    current_mem = Sys.total_memory() / 1024^3
    println_func("Total Memory Allocated: $(round(total_memory/1024^3, digits=2)) GB")
    println_func("Peak  Memory Usage: $(round(peak_memory/1024^3, digits=2)) GB")
    println_func("Total Available System Memory: $(round(current_mem, digits=2)) GB")
    
    # Process statistics
    println_func("\nProcess Information:")
    println_func(repeat("-", 90))
    println_func("Number of Threads: $(Threads.nthreads())")
    println_func("Julia Version: $(VERSION)")
    println_func("System: $(Sys.MACHINE)")
    
    println_func("\n", repeat("=", 90))
end
