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
    timings = Dict{String, Any}()

    # Read and validate parameters
    params_timing = @timed begin
        params_string = read(params_path, String)
        params = check_params_bsp(params_string)

        # Deterministic RNG for decoy/entrapment shuffling. Configurable via the
        # `seed` build param (default 1844, backward-compatible): a fixed seed
        # makes the library fully reproducible; a different seed produces
        # different shuffled decoy/entrapment sequences.
        build_seed = Int(get(params, "seed", 1844))
        Random.seed!(build_seed)

        # Get library directory (already has .poin extension added in check_params)
        lib_dir = params["_lib_dir"]
        mkpath(lib_dir)

        # Write complete merged parameters to config.json (not just user input)
        params_out_path = joinpath(lib_dir, "config.json")
        params_json = JSON.json(params, 2)  # Pretty-print with 2-space indent
        write(params_out_path, params_json)
        nothing
    end
    timings["Parameter Loading"] = params_timing

    # Initialize the shared Pioneer logging system: same four files
    # (pioneer_search_report.txt / pioneer_search_log.log /
    #  pioneer_search_debug.log / pioneer_warnings.log) and the same
    # @user_info / @user_warn macros that SearchDIA uses.
    debug_level = haskey(params, "logging") && haskey(params["logging"], "debug_console_level") ?
                  Int(params["logging"]["debug_console_level"]) : 0
    # Same knob as SearchDIA: the debug log is populated independently of the console.
    debug_file_lvl = haskey(params, "logging") && haskey(params["logging"], "debug_file_level") ?
                     Int(params["logging"]["debug_file_level"]) : 1
    max_bytes = haskey(params, "logging") && haskey(params["logging"], "max_message_bytes") ?
                Int(params["logging"]["max_message_bytes"]) : 4096
    warnings_full_path = init_pioneer_logging(
        lib_dir,
        "Pioneer Library Build Log";
        debug_console_level = debug_level,
        debug_file_level    = debug_file_lvl,
        max_message_bytes   = max_bytes,
    )

    try
        @user_print "\n" * repeat("=", 90)
        @user_print "Spectral Library Building Process"
        @user_print repeat("=", 90)
        @user_info "\nStarting library build at: $(Dates.now())"
        @user_info "Output directory: $lib_dir"
        @user_info "RNG seed: $build_seed"

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
        @user_info "Detecting fragment bounds..."
        bounds_timing = @timed begin
            # calibration_raw_file is optional - use empty string if not provided
            calibration_file = get(params, "calibration_raw_file", "")
            frag_bounds, prec_mz_min, prec_mz_max = get_fragment_bounds(
                _params.library_params["auto_detect_frag_bounds"],
                calibration_file,
                (Float32(_params.library_params["frag_mz_min"]), Float32(_params.library_params["frag_mz_max"])),
                (Float32(_params.library_params["prec_mz_min"]), Float32(_params.library_params["prec_mz_max"]))
            )
            @user_info "Fragment m/z range: $frag_bounds"
            @user_info "Precursor m/z range: ($prec_mz_min, $prec_mz_max)"
            nothing
        end
        timings["Fragment Bound Detection"] = bounds_timing

        N_PRECURSORS = 0
        N_FRAGMENTS = 0
        N_TARGETS = 0
        N_DECOYS = 0
        # In-memory handoff for the fused raw->detailed_frags path (set in the
        # prediction block below). When predict_fragments is false (resume from
        # existing files) they stay nothing and buildPionLib uses the disk path.
        detailed_frags = nothing
        pid_to_fid = nothing
        if params["predict_fragments"]
            # Fragment prediction workflow
            @user_info "Starting fragment prediction workflow..."

            # Altimeter (SplineCoefficientModel) is the only supported fragment-prediction model.
            model_timing = @timed begin
                prediction_model = "altimeter"
                @user_info "Using prediction model: $prediction_model"
                frag_annotation_type = MODEL_CONFIGS[prediction_model].annotation_type
                koina_model_type = MODEL_CONFIGS[prediction_model].model_type
                nothing
            end
            timings["Model Validation"] = model_timing

            mz_to_ev_interp = missing

            # Chronologer prediction workflow
            @user_info "Preparing chronologer workflow..."
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

            @user_info "Predicting retention times..."
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
                    3.0f0  # rt_bin_tol for precursor sorting (matches fragment index)
                )
                #println("precursors_arrow_path $precursors_arrow_path")
                # Cleanup temporary files. Use safeRm (GC + retry + rename
                # fallback): on Windows the just-read Arrow files are still
                # mmap-locked, so a plain rm() throws EACCES (GC.gc() alone is
                # not enough to release the mapping).
                GC.gc()
                safeRm(chronologer_in_path; force=true)
                safeRm(chronologer_out_path; force=true)
                dir, filename = splitdir(precursors_arrow_path)
                raw_fragments_arrow_path = joinpath(dir, "raw_fragments.arrow")
                safeRm(raw_fragments_arrow_path; force=true)
                nothing
            end
            timings["Chronologer Output Processing"] = parse_timing

            # Fragment-filter knobs. These are hardcoded for the Altimeter
            # spline path (the JSON's `library_params.max_frag_rank` etc. are
            # NOT read here — see the buildPionLib call below where the same
            # constants are passed). They live as locals so the early filter
            # in `predict_fragments` and the late path in `buildPionLib`
            # share a single source of truth.
            const_y_start                       = UInt8(3)
            const_b_start                       = UInt8(2)
            const_include_p                     = false
            const_include_isotope               = false
            const_include_immonium              = false
            const_include_internal              = false
            const_include_neutral_diff          = true
            const_max_frag_charge               = UInt8(3)
            const_max_frag_rank                 = UInt8(10)
            const_length_to_frag_count_multiple = Float32(1000)
            const_min_frag_intensity            = 0.0f0

            # Fragment prediction
            @user_info "Predicting fragment ion intensities..."
            frag_predict_timing = @timed begin
                # Build the per-batch filter context that lets predict_fragments
                # apply the same metadata + rank filters that getDetailedFrags
                # otherwise applies after raw_fragments.arrow is on disk.
                precursors_df_for_ctx = DataFrame(Arrow.Table(precursors_arrow_path))
                filter_ctx = build_spline_frag_filter_ctx(
                    frag_bounds,
                    precursors_df_for_ctx,
                    asset_path("ion_dictionary.txt"),
                    asset_path("immonium.txt"),
                    frag_annotation_type;
                    y_start = const_y_start,
                    b_start = const_b_start,
                    include_p = const_include_p,
                    include_isotope = const_include_isotope,
                    include_immonium = const_include_immonium,
                    include_internal = const_include_internal,
                    include_neutral_diff = const_include_neutral_diff,
                    max_frag_charge = const_max_frag_charge,
                    max_frag_rank = const_max_frag_rank,
                    length_to_frag_count_multiple = const_length_to_frag_count_multiple,
                    min_frag_intensity = const_min_frag_intensity,
                )
                predict_fragments(
                    precursors_arrow_path,
                    raw_fragments_arrow_path,
                    koina_model_type,
                    filter_ctx,
                    24,     # max_koina_requests
                    1000,   # max_koina_batch
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
                # Record the spline knots (stored as Arrow metadata, not per-row column).
                # Altimeter is the only supported fragment-prediction backend, so the knot
                # metadata is always present.
                meta = Arrow.getmetadata(fragments_table)
                spl_knots = Vector{Float32}(JSON.parse(meta["knot_vector"]))
                serialize_to_jls(
                    joinpath(lib_dir, "spline_knots.jls"),
                    spl_knots
                )
                ion_dictionary = get_altimeter_ion_dict(asset_path("ion_dictionary.txt"))

                # Fused: decode raw_fragments directly into the final
                # Vector{SplineCompactFrag} + CSR, in memory, skipping the
                # fragments_table.arrow re-encode/re-read. Handed to buildPionLib
                # below without touching disk.
                detailed_frags, pid_to_fid = build_detailed_frags_from_raw(
                    precursors_table,
                    fragments_table,
                    frag_annotation_type,
                    ion_dictionary,
                    asset_path("immonium.txt"),
                    lib_dir,
                    Dict{String, Int8}(),
                    iso_mod_to_mass,
                    koina_model_type
                )



                # Process precursor table
                N_FRAGMENTS = length(fragments_table[:mz])
                N_PRECURSORS = length(precursors_table[:mz])
                N_DECOYS  = sum(precursors_table[:decoy])
                N_TARGETS = N_PRECURSORS - N_DECOYS
                fragments_table = nothing
                # raw_fragments.arrow is fully consumed by parse_altimeter_fragments
                # (process_spline_batch! decode + rebuild_prec_to_frag_index!); its
                # mmap handle (fragments_table) is released just above. Delete it now —
                # before the precursor DataFrame work + precursors_table.arrow write —
                # so raw no longer coexists with fragments_table.arrow at the peak.
                GC.gc()
                safeRm(raw_fragments_arrow_path; force=true)
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
                @debug_l1 "  Before add_pair_indices!: $(nrow(precursors_table)) precursors"
                @debug_l1 "  Unique pair_ids available: $(length(unique(precursors_table.pair_id)))"
                add_pair_indices!(precursors_table)  # Add partner indices AFTER all sorting is complete
                partner_col = precursors_table.partner_precursor_idx
                max_partner = all(ismissing, partner_col) ? 0 : Int64(maximum(skipmissing(partner_col)))
                @debug_l1 "  After add_pair_indices!: max partner_idx = $max_partner (table size: $(nrow(precursors_table)))"
                @debug_l1 "  Valid indices: $(max_partner <= nrow(precursors_table) ? "✅ YES" : "❌ NO")"

                # Add entrapment target indices if entrapment_pair_id column exists AND entrapment is enabled
                entrapment_r = get(params["fasta_digest_params"], "entrapment_r", 0)
                if hasproperty(precursors_table, :entrapment_pair_id) && entrapment_r > 0
                    @debug_l1 "  Adding entrapment target indices..."
                    add_entrapment_indices!(precursors_table)
                    ent_col = precursors_table.entrapment_target_idx
                    max_entrap_target = all(ismissing, ent_col) ? 0 : Int64(maximum(skipmissing(ent_col)))
                    @debug_l1 "  After add_entrapment_indices!: max entrapment_target_idx = $max_entrap_target"
                    @debug_l1 "  Entrapment indices valid: $(max_entrap_target <= nrow(precursors_table) ? "✅ YES" : "❌ NO")"
                end
                Arrow.write(
                    joinpath(lib_dir, "precursors_table.arrow"),
                    precursors_table
                )

                GC.gc()
                # precursors.arrow is now fully materialized into precursors_table
                # (DataFrame copy + Arrow.write above), so delete it before buildPionLib
                # + File Verification. (raw_fragments.arrow was deleted earlier, right
                # after parse_altimeter_fragments.) safeRm = GC + retry for the Windows
                # mmap lock.
                safeRm(precursors_arrow_path; force=true)
                nothing
            end
            timings["Prediction Processing"] = process_timing
        end

        # Verify required files
        verify_timing = @timed begin
            # The fused path hands detailed_frags to buildPionLib in memory, so it
            # needs only the precursor/protein tables. The disk path (resume with
            # predict_fragments=false) still requires the fragment intermediates.
            required_files = detailed_frags === nothing ?
                ["fragments_table.arrow", "prec_to_frag.arrow", "precursors_table.arrow", "proteins_table.arrow"] :
                ["precursors_table.arrow", "proteins_table.arrow"]
            if !all(isfile.(joinpath.(lib_dir, required_files)))
                error("Missing required files in $lib_dir. Try running with predict_fragments=true")
            end
            nothing
        end
        timings["File Verification"] = verify_timing

        # Build final indices
        @user_info "Building final library indices..."
        index_timing = @timed begin
            buildPionLib(
                lib_dir,
                UInt8(4),       # y_start_index
                const_y_start,
                UInt8(3),       # b_start_index
                const_b_start,
                false,          # include_p_index
                const_include_p,
                const_include_isotope,
                const_include_immonium,
                const_include_internal,
                const_include_neutral_diff,
                const_max_frag_charge,
                const_max_frag_rank,
                const_length_to_frag_count_multiple,
                const_min_frag_intensity,
                frag_bounds,
                Float32(get(_params.library_params, "frag_bin_tol_ppm", 0.0)),  # frag_bin_tol_ppm (0 = use mDa mode)
                3.0f0,          # rt_bin_tol
                koina_model_type;
                frag_bin_tol_mda = Float32(get(_params.library_params, "frag_bin_tol_mda", 2.0)),
                detailed_frags = detailed_frags,
                pid_to_fid = pid_to_fid
            )

            nothing
        end
        timings["Index Building"] = index_timing

        # Apply DIA-NN-style decoy conversion if requested
        decoy_method = get(params["fasta_digest_params"], "decoy_method", "shuffle")
        if decoy_method == "diann_mutation"
            @user_info "Applying DIA-NN-style decoy conversion..."
            diann_timing = @timed begin
                apply_diann_decoy_style!(lib_dir)
                nothing
            end
            timings["DIA-NN Decoy Conversion"] = diann_timing
        end

        # Print performance report
        print_performance_report(timings;
            N_PRECURSORS = N_PRECURSORS,
            N_FRAGMENTS  = N_FRAGMENTS,
            N_TARGETS    = N_TARGETS,
            N_DECOYS     = N_DECOYS,
        )

        @user_info "\nLibrary building completed at: $(Dates.now())"
        @user_info "Library location: $lib_dir"
    catch e
        error_msg = try
            "$(typeof(e)): $(e.msg)"
        catch
            "$(typeof(e))"
        end
        @user_error "BuildSpecLib failed: $error_msg"
        @user_error "Stacktrace: $(stacktrace(catch_backtrace()))"
        rethrow(e)
    finally
        finalize_pioneer_logging(warnings_full_path; banner_title = "Library build completed")
    end

    # Clean up temp files after build; on Windows these may still be transiently locked.
    cleanUpLibrary(lib_dir)
    GC.gc()

    return nothing
end

"""
    print_performance_report(timings; N_PRECURSORS, N_FRAGMENTS, N_TARGETS, N_DECOYS)

BuildSpecLib-specific performance report. Writes through `@user_print`,
so it lands on console, the essential log, the console mirror, and the
debug log via the Pioneer logging system.
"""
function print_performance_report(timings; kwargs...)
    @user_print "\n" * repeat("=", 90)
    @user_print "Library Building Performance Report"
    @user_print repeat("=", 90)

    @user_print "\nLibrary Composition:"
    @user_print repeat("-", 90)
    @user_print rpad("# Precursors", 14) * " " *
                rpad("# Fragments", 14) * " " *
                rpad("# Targets", 14) * " " *
                rpad("# Decoys", 14)
    @user_print repeat("-", 90)
    @user_print lpad(string(kwargs[:N_PRECURSORS]), 14) * " " *
                lpad(string(kwargs[:N_FRAGMENTS]),  14) * " " *
                lpad(string(kwargs[:N_TARGETS]),    14) * " " *
                lpad(string(kwargs[:N_DECOYS]),     14)

    @user_print "\nDetailed Step Analysis:"
    @user_print repeat("-", 90)
    @user_print rpad("Step", 40) * " " *
                rpad("Time (s)", 12) * " " *
                rpad("Mem Alloc (GB)", 16) * " " *
                rpad("GC Time (s)", 12) * " " *
                rpad("GC %", 8)
    @user_print repeat("-", 90)

    peak_memory = peak_rss()
    total_time = 0.0
    total_memory = 0
    total_gc = 0.0

    sorted_steps = sort(collect(keys(timings)), by = x -> timings[x].time)
    for step in sorted_steps
        timing = timings[step]
        time_s = timing.time
        mem_gb = timing.bytes / (1024^3)
        gc_s   = timing.gctime
        gc_pct = time_s > 0 ? (gc_s / time_s) * 100 : 0.0

        total_time += time_s
        total_memory += timing.bytes
        total_gc += gc_s

        @user_print rpad(step, 40) * " " *
                    rpad(@sprintf("%.2f", time_s), 12) * " " *
                    rpad(@sprintf("%.2f", mem_gb), 16) * " " *
                    rpad(@sprintf("%.2f", gc_s), 12) * " " *
                    rpad(@sprintf("%.1f", gc_pct), 8)
    end

    @user_print repeat("-", 90)
    @user_print rpad("TOTAL", 40) * " " *
                rpad(@sprintf("%.2f", total_time), 12) * " " *
                rpad(@sprintf("%.2f", total_memory/(1024^3)), 16) * " " *
                rpad(@sprintf("%.2f", total_gc), 12) * " " *
                rpad(@sprintf("%.1f", total_time > 0 ? (total_gc/total_time)*100 : 0.0), 8)

    @user_print "\nMemory Usage Summary:"
    @user_print repeat("-", 90)
    current_mem = Sys.total_memory() / 1024^3
    @user_print "Total alloc (cumulative):   $(round(total_memory/1024^3, digits=2)) GB"
    @user_print "Working set (peak RSS):     $(round(peak_memory/1024^3, digits=2)) GB"
    @user_print "System RAM:                 $(round(current_mem, digits=2)) GB"

    @user_print "\nRuntime Summary:"
    @user_print repeat("-", 90)
    @user_print "Total Runtime: $(round(total_time/60, digits=2)) min ($(round(total_time, digits=2)) s)"
    @user_print "Threads: $(Threads.nthreads())"
    @user_print "\n" * repeat("=", 90)
end
