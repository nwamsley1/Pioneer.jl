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

"""
    IntegrateChromatogramSearch

Search method for analyzing chromatograms to get quantitative information.

This search:
1. Uses precursor and trace information from previous searches
2. Builds chromatograms for each precursor
3. Integrates areas for quantification
4. Incorporates isotope pattern information
"""
struct IntegrateChromatogramSearch <: SearchMethod end

#==========================================================
Type Definitions  
==========================================================#

"""
Results container for chromatogram integration search.
"""
struct IntegrateChromatogramSearchResults <: SearchResults
    psms::Base.Ref{DataFrame}  # Chromatogram data per file
end

"""
Parameters for chromatogram integration search.
"""
struct IntegrateChromatogramSearchParameters{P<:PrecEstimation, I<:IsotopeTraceType} <: FragmentIndexSearchParameters
    # Core parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_fraction_transmitted::Float32
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    spec_order::Set{Int64}
    ms1_quant::Bool

    # Chromatogram parameters
    wh_smoothing_strength::Float32
    n_pad::Int64
    max_apex_offset::Int64
    write_decoys::Bool
    
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

    # Analysis strategies
    isotope_tracetype::I
    prec_estimation::P

    function IntegrateChromatogramSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        global_params = params.global_settings
        quant_params = params.quant_search
        frag_params = quant_params.fragment_settings
        chrom_params = quant_params.chromatogram
        deconv_params = params.optimization.deconvolution
        output_params = params.output
        
        # Determine isotope trace type
        isotope_trace_type = if haskey(global_params.isotope_settings, :combine_traces) && 
                               global_params.isotope_settings.combine_traces
            CombineTraces(0.0f0)  # Default min_fraction_transmitted
        else
            SeperateTraces()
        end

        isotope_bounds = global_params.isotope_settings.err_bounds_quant_search
        min_fraction_transmitted = global_params.isotope_settings.min_fraction_transmitted
        # Always use partial precursor capture for integrate chromatogram
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
        new{typeof(prec_estimation), typeof(isotope_trace_type)}(
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            Float32(min_fraction_transmitted),
            Int64(frag_params.n_isotopes),
            UInt8(frag_params.max_rank),
            Set{Int64}([2]),
            global_params.ms1_quant,

            Float32(chrom_params.smoothing_strength),
            Int64(chrom_params.padding),
            Int64(chrom_params.max_apex_offset),
            Bool(output_params.write_decoys),

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

            isotope_trace_type,
            prec_estimation
        )
    end
end


#==========================================================
Interface Implementation
==========================================================#

get_parameters(::IntegrateChromatogramSearch, params::Any) = IntegrateChromatogramSearchParameters(params)

function init_search_results(::IntegrateChromatogramSearchParameters, search_context::SearchContext)
    return IntegrateChromatogramSearchResults(
        Ref(DataFrame())
    )
end

"""
Process a single file for chromatogram integration.
"""
function process_file!(
    results::IntegrateChromatogramSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData) where {P<:IntegrateChromatogramSearchParameters}

    # Check if file should be skipped due to previous failure
    if check_and_skip_failed_file(search_context, ms_file_idx, "IntegrateChromatogramSearch")
        return results  # Return early with unchanged results
    end

    try
        # Check if required files exist (not empty paths from failed files)
        rt_index_path = getRtIndex(getMSData(search_context), ms_file_idx)
        passing_psms_path = getPassingPsms(getMSData(search_context), ms_file_idx)

        if isempty(rt_index_path) || isempty(passing_psms_path)
            file_name = try
                getMassSpecData(search_context).file_id_to_name[ms_file_idx]
            catch
                "file_$ms_file_idx"
            end
            @debug_l2 "Skipping IntegrateChromatogramSearch for file $file_name - missing required files from previous steps"
            return results
        end

        # Build retention time index for efficient precursor lookup
        # Groups RTs into bins of 0.1 minute width
        rt_index = buildRtIndex(
            DataFrame(Arrow.Table(rt_index_path)),
            bin_rt_size = 0.1)

        # Load PSMs that passed previous filtering steps
        # Convert to DataFrame for processing
        passing_psms = DataFrame(Tables.columntable(Arrow.Table(passing_psms_path)))#load_passing_psms(search_context, parsed_fname)
        
        # Keep only target (non-decoy) PSMs
        if !params.write_decoys
            filter!(row -> row.target, passing_psms)
        end

        # Initialize columns to store integration results
        # peak_area: Integrated area of chromatographic peak
        # new_best_scan: Updated apex scan after refinement
        passing_psms[!, :peak_area] = zeros(Float32, nrow(passing_psms))
        passing_psms[!, :new_best_scan] = zeros(UInt32, nrow(passing_psms))
        passing_psms[!, :points_integrated] = zeros(UInt32, nrow(passing_psms))
        passing_psms[!, :precursor_fraction_transmitted_traces] = fill("", nrow(passing_psms))
        passing_psms[!, :isotopes_captured_traces] = fill("", nrow(passing_psms))
        if params.ms1_quant==true
            passing_psms[!, :new_best_scan] = zeros(UInt32, nrow(passing_psms))
            passing_psms[!, :peak_area_ms1] = zeros(Float32, nrow(passing_psms)) 
            passing_psms[!, :ms1_best_scan] = zeros(UInt32, nrow(passing_psms))
            passing_psms[!, :ms1_points_integrated] = zeros(UInt32, nrow(passing_psms))
        end
        # Extract chromatograms for all passing PSMs
        # Builds chromatograms using parallel processing across scan ranges
        chromatograms = extract_chromatograms(
            spectra,
            passing_psms,
            rt_index,
            search_context,
            params,
            ms_file_idx,
            MS2CHROM(),
        )
        #out_dir = getDataOutDir(search_context)
        #Arrow.write(joinpath(out_dir, "test_chroms_ms2.arrow"), chromatograms)
        #jldsave("/Users/nathanwamsley/Desktop/rt_index.jld2"; rt_index)
        if params.ms1_quant==true
            ms1_chromatograms = extract_chromatograms(
                spectra,
                passing_psms,
                rt_index,
                search_context,
                params,
                ms_file_idx,
                MS1CHROM(),
            )
            # MS1 always uses CombineTraces, so sort by [:precursor_idx, :rt]
            sort!(ms1_chromatograms, [:precursor_idx, :rt])
            ms1_chromatograms[!,:precursor_fraction_transmitted] = ones(Float32, size(ms1_chromatograms, 1))
        end
        #Arrow.write(joinpath(out_dir, "test_chroms_ms1.arrow"), ms1_chromatograms)
        #jldsave("/Users/nathanwamsley/Desktop/test_chroms_ms1.jld2"; ms1_chromatograms)
        # Determine which isotopes are captured in each isolation window
        # Uses quadrupole transmission model to check isotope coverage
        get_isotopes_captured!(
            chromatograms,
            params.isotope_tracetype,
            getQuadTransmissionModel(search_context, ms_file_idx),
            getSearchData(search_context),
            chromatograms[!, :scan_idx],
            getCharge(getActivePrecursors(search_context)),
            getMz(getActivePrecursors(search_context)),
            getSulfurCount(getActivePrecursors(search_context)),
            getCenterMzs(spectra),
            getIsolationWidthMzs(spectra)
        )

        # Pre-sort chromatograms to avoid concurrent sorting in threads
        # Groups will inherit this ordering, eliminating need for per-group sorting
        if seperateTraces(params.isotope_tracetype)
            sort!(chromatograms, [:precursor_idx, :isotopes_captured, :rt])
        else
            sort!(chromatograms, [:precursor_idx, :rt])
        end

        # Integrate chromatographic peaks for each precursor
        # Updates peak_area and new_best_scan in passing_psms
        integrate_precursors(
            chromatograms,
            params.isotope_tracetype,
            params.min_fraction_transmitted,
            passing_psms[!, :precursor_idx],
            passing_psms[!, :isotopes_captured],
            passing_psms[!, :scan_idx],
            passing_psms[!, :peak_area],
            passing_psms[!, :new_best_scan],
            passing_psms[!, :points_integrated],
            passing_psms[!, :precursor_fraction_transmitted_traces],
            passing_psms[!, :isotopes_captured_traces],
            ms_file_idx,
            λ = params.wh_smoothing_strength,
            n_pad = params.n_pad,
            max_apex_offset = params.max_apex_offset
        )
        if params.ms1_quant==true
            integrate_precursors(
                ms1_chromatograms,
                CombineTraces(zero(Float32)),
                params.min_fraction_transmitted,
                passing_psms[!, :precursor_idx],
                passing_psms[!, :isotopes_captured],
                passing_psms[!, :scan_idx],
                passing_psms[!, :peak_area_ms1],
                passing_psms[!, :ms1_best_scan],
                passing_psms[!, :ms1_points_integrated],
                missing, # precursor_fraction_transmitted_traces
                missing, # isotopes_captured_traces
                ms_file_idx,
                λ = params.wh_smoothing_strength,
                n_pad = params.n_pad,
                max_apex_offset = 5,#typemax(Int64),
                test_print = true
            )
        end
        # Clear chromatograms to free memory
        chromatograms = nothing

        # Store processed PSMs in results
        results.psms[] = passing_psms
    catch e
        # Handle failures gracefully using helper function
        handle_search_error!(search_context, ms_file_idx, "IntegrateChromatogramSearch", e, createFallbackResults!, results)
    end

    return results
end

"""
Create empty chromatogram results for a failed file in IntegrateChromatogramSearch.
"""
function createFallbackResults!(results::IntegrateChromatogramSearchResults, ms_file_idx::Int64)
    # Create empty PSM DataFrame with proper schema
    empty_psms = DataFrame(
        ms_file_idx = UInt32[],
        scan_idx = UInt32[], 
        precursor_idx = UInt32[],
        rt = Float32[],
        q_value = Float32[],
        score = Float32[], 
        prob = Float32[],
        peak_area = Float32[],
        peak_height = Float32[],
        fwhm = Float32[]
    )
    
    # Set empty results (don't append since this file failed)
    results.psms[] = empty_psms
end

function process_search_results!(
    results::IntegrateChromatogramSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:IntegrateChromatogramSearchParameters}

    # Check if file should be skipped due to previous failure
    if check_and_skip_failed_file(search_context, ms_file_idx, "IntegrateChromatogramSearch results processing")
        return nothing  # Return early 
    end

    try
        passing_psms = results.psms[]

        # Skip processing if no PSMs (empty DataFrame from failed search)
        if nrow(passing_psms) == 0 || ncol(passing_psms) == 0
            @debug_l2 "No PSMs to process for file $ms_file_idx in IntegrateChromatogramSearch results"
            return nothing
        end

        parsed_fname = getFileIdToName(getMSData(search_context), ms_file_idx)
        # Process final PSMs
        process_final_psms!(
            passing_psms,
            search_context,
            parsed_fname,
            ms_file_idx
        )
        # Save results
        writeArrow(getPassingPsms(getMSData(search_context))[ms_file_idx], passing_psms)
    catch e
        @user_warn "Chromatogram processing failed" ms_file_idx exception=e
        rethrow(e)
    end
    return nothing
end

function reset_results!(results::IntegrateChromatogramSearchResults)
    results.psms[] = DataFrame()
    GC.gc()
    return nothing
end

function summarize_results!(
    ::IntegrateChromatogramSearchResults,
    ::P,
    ::SearchContext
) where {P<:IntegrateChromatogramSearchParameters}
    return nothing
end
