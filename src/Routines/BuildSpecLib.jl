# src/BuildSpecLib.jl

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

    # Initialize timing dictionary for performance tracking
    #params_path = "/Users/n.t.wamsley/RIS_temp/koina_testing/config.json"
    timings = Dict{String, NamedTuple{(:value, :time, :bytes, :gctime, :gcstats),
                                    Tuple{Nothing, Float64, Int64, Float64, Base.GC_Diff}}}()
    
    # Read and validate parameters
    params_timing = @timed begin
        params_string = read(params_path, String)
        params = check_params_bsp(params_string)
        
        # Create output directories
        lib_out_dir = params["out_dir"]
        mkpath(lib_out_dir)
        
        # Library directory (.poin extension)
        lib_dir = joinpath(lib_out_dir, params["lib_name"] * ".poin")
        mkpath(lib_dir)
        
        # Setup logging
        log_path = joinpath(lib_dir, "build_log.txt")
        params_out_path = joinpath(lib_dir, "config.json")
        
        write(params_out_path, params_string)
        #dual_println("Saved parameter file to: ", params_out_path)
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
            chronologer_out_path = joinpath(chronologer_dir, "precursors_for_chronologer.arrow")
            
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
            frag_bounds, prec_mz_min, prec_mz_max = get_fragment_bounds(
                _params.library_params["auto_detect_frag_bounds"],
                _params.library_params["calibration_raw_file"],
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
                            _params.library_params["calibration_raw_file"],
                            lib_dir
                        )
                        dual_println("Successfully created collision energy interpolator")
                    catch
                        @warn "Could not estimate mz to ev conversion. Using default NCE."
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
                                        chronologer_out_path)
                nothing
            end
            timings["Chronologer Preparation"] = chrono_prep_timing

            dual_println("Predicting retention times...")
            rt_timing = @timed begin
                predict_retention_times(chronologer_out_path)
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
                
                # Cleanup temporary files
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
                catch
                    println("No spline knots. static library")
                end
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
                    joinpath(@__DIR__, "../../data/immonium.txt"),
                    lib_dir,
                    Dict{String, Int8}(),
                    iso_mod_to_mass,
                    koina_model_type
                )

                # Process precursor table
                N_FRAGMENTS = length(fragments_table[:mz])
                N_PRECURSORS = length(precursors_table[:mz])
                N_TARGETS = sum(precursors_table[:decoy])
                N_DECOYS = N_PRECURSORS - N_TARGETS
                fragments_table = nothing
                rm(raw_fragments_arrow_path)
                precursors_table = DataFrame(precursors_table)
                rename!(precursors_table, [
                    :accession_number => :accession_numbers,
                    :precursor_charge => :prec_charge,
                    :decoy => :is_decoy,
                    :mods => :structural_mods,
                    :isotope_mods => :isotopic_mods
                ])

                # Convert types
                precursors_table[!, :missed_cleavages] = UInt8.(precursors_table[!, :missed_cleavages])
                precursors_table[!, :prec_charge] = UInt8.(precursors_table[!, :prec_charge])
                precursors_table[!, :mz] = Float32.(precursors_table[!, :mz])
                precursors_table[!, :irt] = Float32.(precursors_table[!, :irt])

                # Save processed precursor table
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
            required_files = ["fragments_table.arrow", "prec_to_frag.arrow", "precursors_table.arrow"]
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
    println_func("Peak Memory Usage: $(round(total_memory/1024^3, digits=2)) GB")
    println_func("Total Garbage Collection Time: $(round(total_gc, digits=2)) seconds")
    println_func("Number of Steps: $(length(timings))")
    println_func("Average Step Runtime: $(round(total_time/length(timings), digits=2)) seconds")
    
    # Memory statistics
    println_func("\nMemory Usage Summary:")
    println_func(repeat("-", 90))
    current_mem = Sys.total_memory() / 1024^3
    println_func("Total Memory Allocated: $(round(total_memory/1024^3, digits=2)) GB")
    println_func("Total Available System Memory: $(round(current_mem, digits=2)) GB")
    println_func("Peak Memory Usage: $(round(maximum([t.bytes for t in values(timings)])/1024^3, digits=2)) GB")
    
    # Process statistics
    println_func("\nProcess Information:")
    println_func(repeat("-", 90))
    println_func("Number of Threads: $(Threads.nthreads())")
    println_func("Julia Version: $(VERSION)")
    println_func("System: $(Sys.MACHINE)")
    
    println_func("\n", repeat("=", 90))
end


#=
using Test 
old_prec = Arrow.Table("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Altimeter111124_MixedSpecies_OlsenAstral_NoEntrapment_SplineTest.poin/precursors_table.arrow")
new_prec = Arrow.Table("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Altimeter4M_ThreeProteome_NoNL_122624.poin/precursors.arrow")

new_seqs_sorted = sort(new_prec[:sequence]);
old_seqs_sorted = sort(old_prec[:sequence]);
new_seqs_set = Set(new_seqs_sorted);
old_seqs_set = Set(old_seqs_sorted);
length(setdiff(old_seqs_set, new_seqs_set))/length(new_seqs_set)
length(setdiff(old_seqs_set, new_seqs_set))/length(old_seqs_set)
length((old_seqs_set ∩ new_seqs_set))/length(old_seqs_set)

seqs_mismatch = (new_seqs_sorted.==old_seqs_sorted).==false;

findall(seqs_mismatch)

setdiff(Set(new_seqs_sorted), Set(old_seqs_sorted))

setdiff(Set(old_seqs_sorted), Set(new_seqs_sorted))

any((sort(old_prec[:sequence]).==sort(new_prec[:sequence])).==false)




precursors_table = DataFrame(Arrow.Table("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Keap1Altimeter.poin/precursors_table.arrow"))
findall((precursors_table[!,:sequence].=="NTGDFGGVAYLLRNLVAVGVGIR").&(precursors_table[!,:prec_charge].==3))
pid_to_fid = load("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Keap1Altimeter.poin/precursor_to_fragment_indices.jld2")["pid_to_fid"]
frags = load("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Keap1Altimeter.poin/detailed_fragments.jld2")["data"]
knots = load("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Keap1Altimeter.poin/spline_knots.jld2")["spl_knots"]
function getIntensity(
                    pf::Pioneer.SplineDetailedFrag{N, T},
                    intensity_type::SplineType{M, T}
                    ) where {M, N, T<:AbstractFloat}
   return splevl(getNCE(intensity_type), getKnots(intensity_type), pf.intensity, getDegree(intensity_type))::T
end
spline_features = SplineType(Tuple(knots), 25.0f0, 3)
frag_name_to_idx = load("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Keap1Altimeter.poin/frag_name_to_idx.jld2")["frag_name_to_idx"];
frag_idx_to_name = Dict(zip(values(frag_name_to_idx), keys(frag_name_to_idx)))
intensities = []
frag_name = []
mz = []
for i in range(pid_to_fid[474],pid_to_fid[475]-1)
    push!(intensities, getIntensity(frags[i], spline_features))
    push!(frag_name, frag_idx_to_name[frags[i].ion_type])
    push!(mz, frags[i].mz)
end
test_data = DataFrame(Dict(:intensity=>intensities,:annotation=>frag_name,:mz=>mz))
CSV.write("/Users/n.t.wamsley/Desktop/test_NTGDFGGVAYLLRNLVAVGVGIR_3_13m.csv", test_data)


x = test_data[!,:mz]
y = test_data[!,:intensity]
p = plot()
for i in 1:length(x)
    plot!([x[i], x[i]], [0, y[i]], 
          color=:black, 
          label=false,
          linewidth=2)
end

# Add labels and title
plot!(
    xlabel="m/z",
    ylabel="intensity",
    title="NTGDFGGVAYLLRNLVAVGVGIR_3"
)



#===========================================
old model
===========================================#


precursors_table = DataFrame(Arrow.Table("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Altimeter111124_MixedSpecies_OlsenAstral_NoEntrapment_SplineTest.poin/precursors_table.arrow"))
prec_idx = findall((precursors_table[!,:sequence].=="NTGDFGGVAYLLRNLVAVGVGIR").&(precursors_table[!,:prec_charge].==3))[1]
pid_to_fid = load("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Altimeter111124_MixedSpecies_OlsenAstral_NoEntrapment_SplineTest.poin/precursor_to_fragment_indices.jld2")["pid_to_fid"]
frags = load("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Altimeter111124_MixedSpecies_OlsenAstral_NoEntrapment_SplineTest.poin/detailed_fragments.jld2")["data"]
knots = load("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Altimeter111124_MixedSpecies_OlsenAstral_NoEntrapment_SplineTest.poin/spline_knots.jld2")["spl_knots"]
function getIntensity(
                    pf::Pioneer.SplineDetailedFrag{N, T},
                    intensity_type::SplineType{M, T}
                    ) where {M, N, T<:AbstractFloat}
   return splevl(getNCE(intensity_type), getKnots(intensity_type), pf.intensity, getDegree(intensity_type))::T
end
spline_features = SplineType(Tuple(knots), 25.0f0, 3)
#frag_name_to_idx = load("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Keap1Altimeter.poin/frag_name_to_idx.jld2")["frag_name_to_idx"];
#frag_idx_to_name = Dict(zip(values(frag_name_to_idx), keys(frag_name_to_idx)))
intensities = []
#frag_name = []
mz = []
for i in range(pid_to_fid[prec_idx],pid_to_fid[prec_idx + 1]-1)
    push!(intensities, getIntensity(frags[i], spline_features))
    #push!(frag_name, frag_idx_to_name[frags[i].ion_type])
    push!(mz, frags[i].mz)
end
test_data = DataFrame(Dict(:intensity=>intensities,:mz=>mz))
CSV.write("/Users/n.t.wamsley/Desktop/test_NTGDFGGVAYLLRNLVAVGVGIR_3_old_model.csv", test_data)


x = test_data[!,:mz]
y = test_data[!,:intensity]
p = plot()
for i in 1:length(x)
    plot!([x[i], x[i]], [0, y[i]], 
          color=:black, 
          label=false,
          linewidth=2)
end

# Add labels and title
plot!(
    xlabel="m/z",
    ylabel="intensity",
    title="NTGDFGGVAYLLRNLVAVGVGIR_3"
)
=#