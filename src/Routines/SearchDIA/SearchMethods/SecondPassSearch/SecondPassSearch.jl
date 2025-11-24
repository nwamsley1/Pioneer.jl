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

#==========================================================
MS1 Feature Schema Definition
==========================================================#
"""
Complete MS1 feature schema with column names (without _ms1 suffix), types,
and sentinel values. These are the column names as they appear in ms1_psms
before joining. The leftjoin operation will automatically add the _ms1 suffix.

This ensures consistent Arrow schemas across all files, regardless of whether
MS1 PSMs were found during search.

Sentinel values indicate missing MS1 data:
- Float types: -1.0 (no MS1 measurement)
- UInt types: 0 (no MS1 data)
- Bool: false (feature not present)

Schema derived from Ms1ScoredPSM struct + parseMs1Psms additions.
"""
const MS1_BASE_SCHEMA = [
    # From Ms1ScoredPSM struct (core MS1 scoring features)
    (:m0, Bool, false),
    (:n_iso, UInt8, UInt8(0)),
    (:big_iso, UInt8, UInt8(0)),
    (:m0_error, Float16, Float16(-1)),
    (:error, Float16, Float16(-1)),
    (:spectral_contrast, Float16, Float16(-1)),
    (:fitted_spectral_contrast, Float16, Float16(-1)),
    (:gof, Float16, Float16(-1)),
    (:max_matched_residual, Float16, Float16(-1)),
    (:max_unmatched_residual, Float16, Float16(-1)),
    (:fitted_manhattan_distance, Float16, Float16(-1)),
    (:matched_ratio, Float16, Float16(-1)),
    (:weight, Float32, Float32(-1)),
    (:ms_file_idx, UInt32, UInt32(0)),
    (:scan_idx, UInt32, UInt32(0)),
    # From parseMs1Psms (RT-related features)
    (:rt, Float32, Float32(-1)),
    (:rt_max_intensity, Float32, Float32(-1)),
    (:rt_diff_max_intensity, Float32, Float32(-1)),
    # From SecondPassSearch join (precursor pairing)
    (:pair_idx, UInt32, UInt32(0)),
]

"""
    get_expected_column_order(psms::DataFrame) -> Vector{Symbol}

Returns the expected column order for PSMs DataFrame with MS1 features.
Ensures consistent Arrow schema across all files by enforcing deterministic column positions.

# Column Ordering Strategy
1. All non-MS1-related columns in their current order
2. MS1 feature columns in MS1_BASE_SCHEMA order (with _ms1 suffix)
3. MS1-derived computed columns (ms1_ms2_rt_diff, ms1_features_missing)

# Rationale
Arrow.jl requires not just matching column names and types, but also matching column **positions**
across files. The leftjoin operation preserves the column order from the right DataFrame (ms1_psms),
which can vary depending on what parseMs1Psms() returns. This function enforces a consistent order.

# Performance
Uses select!() for in-place reordering - only reorders column metadata pointers, does not copy data.
Overhead: <1ms per file.
"""
function get_expected_column_order(psms::DataFrame)
    all_cols = propertynames(psms)

    # MS1 columns in schema order (with _ms1 suffix added)
    ms1_cols_ordered = [Symbol(String(col_name) * "_ms1") for (col_name, _, _) in MS1_BASE_SCHEMA]

    # MS1-derived computed columns
    ms1_computed = [:ms1_ms2_rt_diff, :ms1_features_missing]

    # All MS1-related columns
    ms1_related = Set(vcat(ms1_cols_ordered, ms1_computed))

    # Non-MS1 columns in their current order
    non_ms1_cols = [c for c in all_cols if c ∉ ms1_related]

    # Final order: non-MS1 columns, then MS1 columns in schema order, then computed columns
    return vcat(
        non_ms1_cols,
        ms1_cols_ordered,
        ms1_computed
    )
end

"""
    SecondPassSearch

Second pass search method using optimized parameters from initial searches.

This search:
1. Uses optimized parameters from first pass search
2. Performs PSM identification with previously calculated Huber delta
3. Tracks retention time windows for efficient searching
4. Records chromatogram information for later analysis

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :second_pass_params => Dict(
        "min_y_count" => 1,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "n_frag_isotopes" => 2,
        "max_frag_rank" => 10
    )
)

# Execute search
results = execute_search(SecondPassSearch(), search_context, params)
```
"""
struct SecondPassSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for second pass search.
"""
struct SecondPassSearchResults <: SearchResults
    psms::Base.Ref{DataFrame}          # PSMs for each file
    ms1_psms::Base.Ref{DataFrame}
end

"""
Parameters for second pass search.
"""
struct SecondPassSearchParameters{P<:PrecEstimation, I<:IsotopeTraceType} <: FragmentIndexSearchParameters
    # Core parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_fraction_transmitted::Float32
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    spec_order::Set{Int64}
    match_between_runs::Bool

    # Deconvolution parameters (MS2)
    lambda::Float32
    reg_type::RegularizationType
    max_iter_newton::Int64
    max_iter_bisection::Int64
    max_iter_outer::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    max_diff::Float32

    # MS1 deconvolution parameters
    ms1_lambda::Float32
    ms1_reg_type::RegularizationType
    ms1_huber_delta::Float32

    # PSM filtering
    min_y_count::Int64
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::Int64
    
    # Precursor estimation strategy
    isotope_tracetype::I
    prec_estimation::P

    # Collect MS1 data?
    ms1_scoring::Bool

    function SecondPassSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        global_params = params.global_settings
        quant_params = params.quant_search
        frag_params = quant_params.fragment_settings
        deconv_params = params.optimization.deconvolution
        
        # Determine isotope trace type based on global settings
        isotope_trace_type = if haskey(global_params.isotope_settings, :combine_traces) && 
                               global_params.isotope_settings.combine_traces
            SeperateTraces() #CombineTraces(0.0f0)  # Default min_fraction_transmitted
        else
            SeperateTraces()
        end

        isotope_bounds = global_params.isotope_settings.err_bounds_quant_search
        min_fraction_transmitted = global_params.isotope_settings.min_fraction_transmitted
        prec_estimation = global_params.isotope_settings.partial_capture ? PartialPrecCapture() : FullPrecCapture()

        # Parse MS2 regularization type
        reg_type = deconv_params.ms2.reg_type
        if reg_type == "none"
            reg_type = NoNorm()
        elseif reg_type == "l1"
            reg_type = L1Norm()
        elseif reg_type == "l2"
            reg_type = L2Norm()
        else
            reg_type = NoNorm()
            @user_warn "Warning. MS2 reg type `$reg_type` not recognized. Using NoNorm. Accepted types are `none`, `l1`, `l2`"
        end

        # Parse MS1 regularization type
        ms1_reg_type = deconv_params.ms1.reg_type
        if ms1_reg_type == "none"
            ms1_reg_type = NoNorm()
        elseif ms1_reg_type == "l1"
            ms1_reg_type = L1Norm()
        elseif ms1_reg_type == "l2"
            ms1_reg_type = L2Norm()
        else
            ms1_reg_type = NoNorm()
            @user_warn "Warning. MS1 reg type `$ms1_reg_type` not recognized. Using NoNorm. Accepted types are `none`, `l1`, `l2`"
        end

        ms1_scoring = Bool(global_params.ms1_scoring)
        new{typeof(prec_estimation), typeof(isotope_trace_type)}(
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            Float32(min_fraction_transmitted),
            Int64(frag_params.n_isotopes),
            UInt8(frag_params.max_rank),
            Set{Int64}([2]),
            Bool(global_params.match_between_runs),
            
            Float32(deconv_params.ms2.lambda),
            reg_type,
            Int64(deconv_params.newton_iters),
            Int64(deconv_params.bisection_iters),
            Int64(deconv_params.outer_iters),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.max_diff),

            Float32(deconv_params.ms1.lambda),
            ms1_reg_type,
            Float32(deconv_params.ms1.huber_delta),

            Int64(frag_params.min_y_count),
            Int64(frag_params.min_count),
            Float32(frag_params.min_spectral_contrast),
            Float32(frag_params.min_log2_ratio),
            (Int64(first(frag_params.min_top_n)), Int64(last(frag_params.min_top_n))),
            Int64(frag_params.max_rank),
            
            isotope_trace_type,
            prec_estimation,

            ms1_scoring
        )
    end
end

getIsotopeTraceType(p::SecondPassSearchParameters) = p.isotope_tracetype
#==========================================================
Interface Implementation
==========================================================#

get_parameters(::SecondPassSearch, params::Any) = SecondPassSearchParameters(params)

function init_search_results(::P, search_context::SearchContext) where {P<:SecondPassSearchParameters}
    second_pass_psms = joinpath(getDataOutDir(search_context), "temp_data", "second_pass_psms")
    !isdir(second_pass_psms) && mkdir(second_pass_psms)
    return SecondPassSearchResults(
        DataFrame(),
        DataFrame()
    )
end


#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single file for second pass search.
"""
function process_file!(
    results::SecondPassSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:SecondPassSearchParameters}

    # Check if file should be skipped due to previous failure
    if check_and_skip_failed_file(search_context, ms_file_idx, "SecondPassSearch")
        return results  # Return early with unchanged results
    end
    
    try
        # Get RT index
        rt_index = buildRtIndex(
            DataFrame(Arrow.Table(getRtIndex(getMSData(search_context), ms_file_idx))),
            bin_rt_size = 0.1)

        # Perform second pass search
        psms = perform_second_pass_search(
            spectra,
            rt_index,
            search_context,
            params,
            ms_file_idx,
            MS2CHROM()
        )
        if params.ms1_scoring
            precursors_passing = unique(psms[!,:precursor_idx])

            # Check if we have any passing precursors for MS1 scoring
            if !isempty(precursors_passing)
                precursors = getPrecursors(getSpecLib(search_context));
                seqs = [getSequence(precursors)[pid] for pid in precursors_passing]
                pids = [pid for pid in precursors_passing]
                pcharge = [getCharge(precursors)[pid] for pid in precursors_passing]
                pmz = [getMz(precursors)[pid] for pid in precursors_passing]
                isotopes_dict = getIsotopes(seqs, pmz, pids, pcharge, QRoots(5), 5)
                precursors_passing = Set(precursors_passing)
                # Perform MS1 search (diagnostic timing)
                ms1_psms = perform_second_pass_search(
                    spectra,
                    rt_index,
                    search_context,
                    params,
                    ms_file_idx,
                    precursors_passing,
                    isotopes_dict,
                    MS1CHROM()
                )
                pair_idx = getPairIdx(precursors);
                is_decoy = getIsDecoy(precursors);
                partner_idx = getPartnerPrecursorIdx(precursors);
                ms1_psms[!,:pair_idx] = [pair_idx[pid] for pid in ms1_psms[!,:precursor_idx]]
                ms1_psms_partner = copy(ms1_psms)
                for i in range(1, size(ms1_psms_partner, 1))
                    p = partner_idx[ms1_psms_partner[i,:precursor_idx]]
                    if !ismissing(p)
                        ms1_psms_partner[i,:precursor_idx] = p
                    else
                        ms1_psms_partner[i,:precursor_idx] = 0
                    end
                end
                filter!(x->!iszero(x.precursor_idx), ms1_psms_partner);
                ms1_psms = vcat([ms1_psms, ms1_psms_partner]...)
            else
                # No passing precursors for MS1 scoring
                @debug_l2 "No passing precursors found for MS1 scoring in file $ms_file_idx. Skipping MS1 isotope calculations."
                ms1_psms = DataFrame()
            end
            #rt_irt_model = getRtIrtModel(search_context, ms_file_idx)
            #ms1_psms[!,:irt] = zeros(Float32, size(ms1_psms, 1))
            #ms1_psms[!,:pred_irt] = zeros(Float32, size(ms1_psms, 1))
            #for i in range(1, size(ms1_psms, 1))
            #    ms1_psms[i,:irt] = rt_irt_model(getRetentionTime(spectra, ms1_psms[i,:scan_idx]))
            #    ms1_psms[i,:pred_irt] = getIrt(precursors)[ms1_psms[i,:precursor_idx]]
            #end

            #ms1_psms[!,:is_decoy_] = [is_decoy[pid] for pid in ms1_psms[!,:precursor_idx]]
            #end
                #pprof();
            #end
            #println("times ")
            #println("time.time ", times.time)
            #println("times.bytes ", times.bytes)
            #println("times.gctime ", times.gctime)
            #println("times.gcstats", times.gcstats)
            #println("loc_conflicts ", times.lock_conflicts)
            #println("compile_time ", times.compile_time)
            #println("recompile_time ", times.recompile_time)
            #println("\n")
        else
            ms1_psms = DataFrame()
        end

        results.psms[] = psms
        results.ms1_psms[] = ms1_psms

    catch e
        # Handle failures gracefully using helper function (logs full stacktrace)
        handle_search_error!(search_context, ms_file_idx, "SecondPassSearch", e, createFallbackResults!, results)
    end

    return results
end

"""
Create empty PSM results for a failed file in SecondPassSearch.
"""
function createFallbackResults!(results::SecondPassSearchResults, ms_file_idx::Int64)
    # Create empty PSM DataFrame with proper schema for regular PSMs
    empty_psms = DataFrame(
        ms_file_idx = UInt32[],
        scan_idx = UInt32[], 
        precursor_idx = UInt32[],
        rt = Float32[],
        q_value = Float32[],
        score = Float32[], 
        prob = Float32[]
    )
    
    # Create empty MS1 PSM DataFrame with proper schema
    empty_ms1_psms = DataFrame(
        ms_file_idx = UInt32[],
        scan_idx = UInt32[], 
        precursor_idx = UInt32[],
        rt = Float32[],
        q_value = Float32[],
        score = Float32[], 
        prob = Float32[]
    )
    
    # Set empty results (don't append since this file failed)
    results.psms[] = empty_psms
    results.ms1_psms[] = empty_ms1_psms
end

function process_search_results!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext, 
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:SecondPassSearchParameters}

    # Check if file should be skipped due to previous failure
    if check_and_skip_failed_file(search_context, ms_file_idx, "SecondPassSearch results processing")
        return nothing  # Return early 
    end

    try
        # Get PSMs from results container
        psms = results.psms[]
        ms1_psms = results.ms1_psms[]
        # Add basic search columns (RT, charge, target/decoy status)
        add_second_search_columns!(psms, 
            getRetentionTimes(spectra),
            getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge], 
            getIsDecoy(getPrecursors(getSpecLib(search_context))),#[:is_decoy],
            getPrecursors(getSpecLib(search_context))
            );

        # Determine which precursor isotopes are captured in each scan's isolation window
        get_isotopes_captured!(
            psms,
            getIsotopeTraceType(params),#.isotope_tracetype,
            getQuadTransmissionModel(search_context, ms_file_idx),
            getSearchData(search_context),
            psms[!, :scan_idx],
            getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge],
            getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
            getSulfurCount(getPrecursors(getSpecLib(search_context))),
            getCenterMzs(spectra),
            getIsolationWidthMzs(spectra)
        )

        # Remove PSMs where only M2+ isotopes are captured (expect poor quantification)
        #excluded_isotopes = (Int8(-1), Int8(-1))
        #filter!(row -> row.isotopes_captured != excluded_isotopes, psms)
        filter!(row -> row.precursor_fraction_transmitted >= params.min_fraction_transmitted, psms)
        #filter!(row -> first(row.isotopes_captured) > 2, psms)

        # Initialize columns for best scan selection and summary statistics
        psms[!,:best_scan] = zeros(Bool, size(psms, 1));
        init_summary_columns!(psms);

        # Calculate summary scores for each PSM group
        for (key, gpsms) in pairs(groupby(psms, getPsmGroupbyCols(getIsotopeTraceType(params))))
            get_summary_scores!(
                gpsms, 
                gpsms[!,:weight],
                gpsms[!,:gof],
                gpsms[!,:matched_ratio],
                gpsms[!,:fitted_manhattan_distance],
                gpsms[!,:fitted_spectral_contrast],
                gpsms[!,:scribe],
                gpsms[!,:y_count],
                getRtToRefinedIrtModel(search_context, ms_file_idx)
            );
        end
        # Keep only apex scans for each PSM group
        filter!(x->x.best_scan, psms);

        # Build MS2 RT lookup for efficient MS1 alignment
        ms2_rt_lookup = Dict{UInt32, Float32}(
            row.precursor_idx => row.rt for row in eachrow(psms)
        )

        # Apply hybrid MS1 selection (RT proximity + max intensity features)
        ms1_psms = parseMs1Psms(ms1_psms, spectra, ms2_rt_lookup)

        # Standardize MS1 schema to ensure consistent Arrow output across all files
        if size(ms1_psms, 1) > 0
            # Ensure ms1_psms has correct types for all columns in schema
            for (col_name, col_type, sentinel_value) in MS1_BASE_SCHEMA
                if hasproperty(ms1_psms, col_name)
                    # Convert to correct type
                    ms1_psms[!, col_name] = col_type.(ms1_psms[!, col_name])
                else
                    # Add missing column with sentinel values
                    ms1_psms[!, col_name] = fill(sentinel_value, nrow(ms1_psms))
                end
            end

            # Join - leftjoin automatically adds _ms1 suffix via renamecols
            psms = leftjoin(
                psms,
                ms1_psms,
                on = :precursor_idx,
                makeunique = true,
                renamecols = "" => "_ms1"
            )

            # Fill missing values (precursors without MS1 match) with sentinels
            for (col_name, col_type, sentinel_value) in MS1_BASE_SCHEMA
                ms1_col = Symbol(String(col_name) * "_ms1")
                psms[!, ms1_col] = coalesce.(psms[!, ms1_col], sentinel_value)
                disallowmissing!(psms, ms1_col)
            end

            # Track which PSMs have missing MS1 data
            miss_mask = ismissing.(psms[!, Symbol(String(first(MS1_BASE_SCHEMA)[1]) * "_ms1")])
        else
            # No MS1 data - create all columns with _ms1 suffix using sentinels
            for (col_name, col_type, sentinel_value) in MS1_BASE_SCHEMA
                ms1_col = Symbol(String(col_name) * "_ms1")
                psms[!, ms1_col] = fill(sentinel_value, nrow(psms))
            end
            miss_mask = trues(nrow(psms))
        end

        # Calculate MS1-MS2 RT difference in refined iRT space with explicit Float32 conversion
        rt_to_refined_irt_model = getRtToRefinedIrtModel(search_context, ms_file_idx)
        psms[!,:ms1_ms2_rt_diff] = Float32.(ifelse.(psms[!,:rt_ms1] .== Float32(-1),
                          Float32(-1),
                          abs.(rt_to_refined_irt_model.(psms[!,:rt]) .- rt_to_refined_irt_model.(psms[!,:rt_ms1]))))

        psms[!, :ms1_features_missing] = miss_mask

        # Critical: Enforce consistent MS1 column ordering for Arrow schema compatibility
        # Arrow requires matching column positions across files, not just names/types.
        # The leftjoin operation preserves input order from ms1_psms, which can vary.
        expected_order = get_expected_column_order(psms)
        select!(psms, expected_order)

        #Add additional features for final analysis
        add_features!(
            psms,
            search_context,
            getTICs(spectra),
            getMzArrays(spectra),
            ms_file_idx,
            getRtIrtModel(search_context, ms_file_idx),
            getRtToRefinedIrtModel(search_context, ms_file_idx),
            getPrecursorDict(search_context)
        )

        # Initialize probability scores (will be calculated later)
        initialize_prob_group_features!(psms, params.match_between_runs)

        # Only save results if we have actual PSMs
        if nrow(psms) > 0
            # Save processed results
            temp_path = joinpath(
                getDataOutDir(search_context), "temp_data",
                "second_pass_psms",
                getParsedFileName(search_context, ms_file_idx) * ".arrow"
            )
            writeArrow(temp_path, psms)
            setSecondPassPsms!(getMSData(search_context), ms_file_idx, temp_path)
        else
            # No PSMs found - mark as empty but don't fail
            @debug_l2 "No PSMs found for file $ms_file_idx in SecondPassSearch, setting empty path"
            setSecondPassPsms!(getMSData(search_context), ms_file_idx, "")
        end
    catch e
        throw(e)
        # Mark file as failed and handle gracefully
        file_name = try
            getMassSpecData(search_context).file_id_to_name[ms_file_idx]
        catch
            "file_$ms_file_idx"
        end

        reason = "SecondPassSearch failed: $(typeof(e))"
        markFileFailed!(search_context, ms_file_idx, reason)
        # Also mark in ArrowTableReference for downstream methods
        setFailedIndicator!(getMSData(search_context), ms_file_idx, true)

        # Log detailed error information with immediate flush
        @user_error "Second pass search failed for MS data file: $file_name"
        flush(stdout); flush(stderr)
        @user_error "Actual error: $(typeof(e)): $e"
        flush(stdout); flush(stderr)
        bt = catch_backtrace()
        @user_error sprint(showerror, e, bt)
        flush(stdout); flush(stderr)

        # Print to console directly as backup
        println(stderr, "\n=== DETAILED ERROR ===")
        println(stderr, "File: $file_name")
        println(stderr, "Error type: $(typeof(e))")
        println(stderr, "Error message: $e")
        println(stderr, "\nStacktrace:")
        showerror(stderr, e, bt)
        println(stderr, "\n=== END ERROR ===\n")
        flush(stderr)

        @user_warn "Creating empty results to continue pipeline."

        # Set empty path for failed file
        setSecondPassPsms!(getMSData(search_context), ms_file_idx, "")

        # TEMPORARY: Rethrow to expose error during debugging
        rethrow(e)

        # Don't rethrow - continue with next file
    end

    return nothing
end


"""
Reset results containers.
"""
function reset_results!(results::SecondPassSearchResults)
    results.psms[] = DataFrame()
end

function summarize_results!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:SecondPassSearchParameters}

    return nothing
end
