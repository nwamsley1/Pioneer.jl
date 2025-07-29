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
    sample_rate::Float32
    spec_order::Set{Int64}
    match_between_runs::Bool

    # Deconvolution parameters
    lambda::Float32
    reg_type::RegularizationType
    max_iter_newton::Int64
    max_iter_bisection::Int64
    max_iter_outer::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    max_diff::Float32

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

        reg_type = deconv_params.reg_type
        if reg_type == "none"
            reg_type = NoNorm()
        elseif reg_type == "l1"
            reg_type = L1Norm()
        elseif reg_type == "l2"
            reg_type = L2Norm()
        else
            reg_type = NoNorm()
            @warn "Warning. Reg type `$reg_type` not recognized. Using NoNorm. Accepted types are `none`, `l1`, `l2`"
        end

        ms1_scoring = Bool(global_params.ms1_scoring)
        new{typeof(prec_estimation), typeof(isotope_trace_type)}(
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            Float32(min_fraction_transmitted),
            Int64(frag_params.n_isotopes),
            UInt8(frag_params.max_rank),
            1.0f0,  # Full sampling rate
            Set{Int64}([2]),
            Bool(global_params.match_between_runs),
            
            Float32(deconv_params.lambda),
            reg_type,
            Int64(deconv_params.newton_iters),
            Int64(deconv_params.bisection_iters),
            Int64(deconv_params.outer_iters),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.max_diff),
            
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

    try
        setNceModel!(
            getFragmentLookupTable(getSpecLib(search_context)), 
            getNceModelModel(search_context, ms_file_idx)
        )
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
            precursors = getPrecursors(getSpecLib(search_context));
            seqs = [getSequence(precursors)[pid] for pid in precursors_passing]
            pids = [pid for pid in precursors_passing]
            pcharge = [getCharge(precursors)[pid] for pid in precursors_passing]
            pmz = [getMz(precursors)[pid] for pid in precursors_passing]
            isotopes_dict = getIsotopes(seqs, pmz, pids, pcharge, QRoots(5), 5)
            precursors_passing = Set(precursors_passing)
            # Perform MS1 search
            #times = @timed begin
            #Profile.clear()
            #@profile begin
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
        @warn "Second pass search failed" ms_file_idx exception=e
        rethrow(e)
    end

    return results
end

function process_search_results!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext, 
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:SecondPassSearchParameters}

    try
        # Get PSMs from results container
        psms = results.psms[]
        ms1_psms = results.ms1_psms[]
        ms1_psms = parseMs1Psms( #Reduce to max intensity scan per precursor_idx
            ms1_psms,
            spectra
        )
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
                getRtIrtModel(search_context, ms_file_idx)
            );
        end
        # Keep only apex scans for each PSM group
        filter!(x->x.best_scan, psms);
        #Need to remove inf gof_ms1?
        #Join MS1 PSMs to MS2 PSMs
        if size(ms1_psms, 1) > 0
            psms = leftjoin(
                psms,
                ms1_psms,
                on = :precursor_idx,
                makeunique = true,
                renamecols = "" => "_ms1",
            )
            ms1_cols = filter(col -> endswith(String(col), "_ms1"), names(psms))
            miss_mask = ismissing.(psms[!, ms1_cols[1]])
        else
            ms1_cols = [
                :weight_ms1, :gof_ms1, :max_matched_residual_ms1,
                :max_unmatched_residual_ms1, :fitted_spectral_contrast_ms1,
                :error_ms1, :m0_error_ms1, :n_iso_ms1, :big_iso_ms1
            ]
            miss_mask = trues(size(psms, 1))
            for col in ms1_cols
                psms[!, col] = -1*ones(Float32, size(psms, 1))
            end
        end

        for col in ms1_cols
            psms[!, col] = coalesce.(psms[!, col], zero(nonmissingtype(eltype(psms[!, col]))))
            disallowmissing!(psms, col)
        end
        psms[!,:rt_diff] = abs.(psms[!,:rt] .- psms[!,:rt_ms1])
        psms[!, :ms1_features_missing] = miss_mask
        #Add additional features for final analysis
        add_features!(
            psms,
            search_context,
            getTICs(spectra),
            getMzArrays(spectra),
            ms_file_idx,
            getRtIrtModel(search_context, ms_file_idx),
            getPrecursorDict(search_context)
        )

        # Initialize probability scores (will be calculated later)
        initialize_prob_group_features!(psms, params.match_between_runs)
        
        # Save processed results
        temp_path = joinpath(
            getDataOutDir(search_context), "temp_data",
            "second_pass_psms",
            getParsedFileName(search_context, ms_file_idx) * ".arrow"
        )
        writeArrow(temp_path, psms)
        setSecondPassPsms!(getMSData(search_context), ms_file_idx, temp_path)
    catch e
        @warn "Failed to process search results" ms_file_idx exception=e
        rethrow(e)
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
