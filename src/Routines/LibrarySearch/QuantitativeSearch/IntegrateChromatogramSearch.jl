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
    psms_paths::Dict{Int64, String}  # Maps file index to PSMs path
    chromatograms::Vector{ChromObject}  # Chromatogram data per file
end

"""
Parameters for chromatogram integration search.
"""
struct IntegrateChromatogramSearchParameters{P<:PrecEstimation, I<:IsotopeTraceType} <: FragmentIndexSearchParameters
    # Core parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    sample_rate::Float32
    spec_order::Set{Int64}

    # Chromatogram parameters
    wh_smoothing_strength::Float32
    n_pad::Int64
    max_apex_offset::Float32
    
    # Deconvolution parameters
    lambda::Float32
    max_iter_newton::Int64
    max_iter_bisection::Int64
    max_iter_outer::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    max_diff::Float32

    # Analysis strategies
    isotope_tracetype::I
    prec_estimation::P

    function IntegrateChromatogramSearchParameters(params::Any)
        qp = params[:quant_search_params]
        dp = params[:deconvolution_params]
        
        # Determine isotope trace type
        isotope_trace_type = if qp["combine_isotope_traces"]
            CombineTraces(Float32(qp["min_fraction_transmitted"]))
        else
            SeperateTraces()
        end

        prec_estimation = PartialPrecCapture()

        new{typeof(prec_estimation), typeof(isotope_trace_type)}(
            (UInt8(3), UInt8(0)),
            Int64(qp["n_frag_isotopes"]),
            UInt8(qp["max_frag_rank"]),
            1.0f0,
            Set(2),
            Float32(qp["WH_smoothing_strength"]),
            Int64(qp["n_pad"]),
            Float32(qp["max_apex_offset"]),
            Float32(dp["lambda"]),
            Int64(dp["max_iter_newton"]),
            Int64(dp["max_iter_bisection"]),
            Int64(dp["max_iter_outer"]),
            Float32(dp["accuracy_newton"]),
            Float32(dp["accuracy_bisection"]),
            Float32(dp["max_diff"]),
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
        Dict{Int64, String}(),
        Vector{ChromObject}(),
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
    spectra::Arrow.Table
) where {P<:IntegrateChromatogramSearchParameters}

    try
        # Get RT index
        #parsed_fname = getParsedFileName(search_context, ms_file_idx)
        rt_df = DataFrame(Arrow.Table(getRtIndexPaths(search_context)[ms_file_idx]))
        rt_index = buildRtIndex(rt_df, bin_rt_size = 0.1)

        # Load passing PSMs
        passing_psms = load_passing_psms(search_context, parsed_fname)
        filter!(row -> row.target, passing_psms)

        # Initialize peak area columns
        passing_psms[!, :peak_area] = zeros(Float32, nrow(passing_psms))
        passing_psms[!, :new_best_scan] = zeros(UInt32, nrow(passing_psms))

        # Get chromatograms
        chromatograms = DataFrame(extract_chromatograms(
            spectra,
            passing_psms,
            rt_index,
            search_context,
            params,
            ms_file_idx
        ))

        # Process isotopes
        getIsotopesCaptured!(
            chromatograms,
            params.isotope_tracetype,
            getQuadTransmissionModel(search_context, ms_file_idx),
            chromatograms[!, :scan_idx],
            getPrecursors(getSpecLib(search_context))[:prec_charge],
            getPrecursors(getSpecLib(search_context))[:mz],
            spectra[:centerMz],
            spectra[:isolationWidthMz]
        )

        # Filter and integrate
        filter!(row -> first(row.isotopes_captured) < 2, chromatograms)
        
        integratePrecursors(
            chromatograms,
            params.isotope_tracetype,
            passing_psms[!, :precursor_idx],
            passing_psms[!, :isotopes_captured],
            passing_psms[!, :scan_idx],
            passing_psms[!, :peak_area],
            passing_psms[!, :new_best_scan],
            λ = params.wh_smoothing_strength,
            n_pad = params.n_pad,
            max_apex_offset = params.max_apex_offset
        )

        # Process final PSMs
        process_final_psms!(
            passing_psms,
            search_context,
            parsed_fname,
            ms_file_idx
        )

        # Save results
        temp_path = joinpath(
            getDataOutDir(search_context),
            "chromatogram_psms",
            parsed_fname * ".arrow"
        )

        Arrow.write(temp_path, passing_psms)
        results.psms_paths[ms_file_idx] = temp_path
        results.chromatograms[ms_file_idx] = chromatograms

    catch e
        @warn "Chromatogram integration failed" ms_file_idx exception=e
        rethrow(e)
    end

    return results
end

function process_search_results!(
    ::IntegrateChromatogramSearchResults,
    ::P, 
    ::SearchContext,
    ::Int64,
    ::Arrow.Table
) where {P<:IntegrateChromatogramSearchParameters}
    return nothing
end

function reset_results!(::IntegrateChromatogramSearchResults)
    return nothing
end

function summarize_results!(
    ::IntegrateChromatogramSearchResults,
    ::P,
    ::SearchContext
) where {P<:IntegrateChromatogramSearchParameters}
    return nothing
end

#==========================================================
Helper Methods
==========================================================#

"""
Extract chromatograms using parallel processing.
"""
function extract_chromatograms(
    spectra::Arrow.Table,
    passing_psms::DataFrame,
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64
)
    thread_tasks = partition_scans(spectra, Threads.nthreads())

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]

            return build_chromatograms(
                spectra,
                last(thread_task),
                Set(passing_psms[!, :precursor_idx]),
                rt_index,
                search_context,
                search_data,
                params,
                ms_file_idx
            )
        end
    end

    return vcat(fetch.(tasks)...)
end

"""
Process final PSMs with additional columns.
"""
function process_final_psms!(
    psms::DataFrame,
    search_context::SearchContext,
    parsed_fname::String,
    ms_file_idx::Int64
)
    # Remove invalid peak areas
    filter!(row -> !isnan(row.peak_area::Float32), psms)
    filter!(row -> row.peak_area::Float32 > 0.0, psms)

    # Add columns
    precursors = getPrecursors(getSpecLib(search_context))
    psms[!, :accession_numbers] = getindex.(Ref(precursors[:accession_numbers]), psms[!, :precursor_idx])
    psms[!, :protein_idx] = getindex.(Ref(getAccessionNumberToId(search_context)), psms[!, :accession_numbers])
    psms[!, :ms_file_idx] = UInt32.(psms[!, :ms_file_idx])
    
    sort!(psms, [:protein_idx, :precursor_idx, :ms_file_idx])
    
    psms[!, :species] = getindex.(Ref(precursors[:proteome_identifiers]), psms[!, :precursor_idx])
    psms[!, :peak_area] = allowmissing(psms[!, :peak_area])
    psms[!, :peak_area_normalized] = allowmissing(zeros(Float32, nrow(psms)))
    
    # Add modification and sequence information
    psms[!, :structural_mods] = allowmissing(getindex.(Ref(precursors[:structural_mods]), psms[!, :precursor_idx]))
    psms[!, :isotopic_mods] = allowmissing(getindex.(Ref(precursors[:isotopic_mods]), psms[!, :precursor_idx]))
    psms[!, :charge] = allowmissing(getindex.(Ref(precursors[:prec_charge]), psms[!, :precursor_idx]))
    psms[!, :sequence] = allowmissing(getindex.(Ref(precursors[:sequence]), psms[!, :precursor_idx]))
    
    psms[!, :file_name] .= parsed_fname
end

"""
Build chromatograms for a set of scans with RT bin caching.
"""
function build_chromatograms(
    spectra::Arrow.Table,
    scan_range::Vector{Int64},
    precursors_passing::Set{UInt32},
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64
)
    # Initialize working arrays
    Hs = getHs(search_data)
    weights = getTempWeights(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    chromatograms = Vector{ChromObject}(undef, 500000)  # Initial size

    # RT bin tracking state
    irt_start, irt_stop = 1, 1
    prec_mz_string = ""
    ion_idx = 0
    rt_idx = 0
    precs_temp = getPrecIds(search_data)  # Use search_data's prec_ids
    prec_temp_size = 0
    
    irt_tol = getIrtErrors(search_context)[ms_file_idx]

    for scan_idx in scan_range
        scan_idx > length(spectra[:mz_array]) && continue
        
        # Process MS2 scans
        msn = spectra[:msOrder][scan_idx]
        msn ∉ params.spec_order && continue

        # Calculate RT window
        irt = getRtIrtModel(search_context, ms_file_idx)(spectra[:retentionTime][scan_idx])
        irt_start_new = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
        irt_stop_new = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

        # Check for m/z change
        prec_mz_string_new = string(spectra[:centerMz][scan_idx])
        prec_mz_string_new = prec_mz_string_new[1:min(length(prec_mz_string_new), 6)]

        # Update transitions if window changed
        if (irt_start_new != irt_start) || (irt_stop_new != irt_stop) || (prec_mz_string_new != prec_mz_string)
            irt_start = irt_start_new
            irt_stop = irt_stop_new
            prec_mz_string = prec_mz_string_new

            quad_func = getQuadTransmissionFunction(
                getQuadTransmissionModel(search_context, ms_file_idx),
                spectra[:centerMz][scan_idx],
                spectra[:isolationWidthMz][scan_idx]
            )

            ion_idx = selectTransitions!(
                getIonTemplates(search_data),
                RTIndexedTransitionSelection(),
                params.prec_estimation,
                getFragmentLookupTable(getSpecLib(search_context)),
                precs_temp,
                getPrecursors(getSpecLib(search_context))[:mz],
                getPrecursors(getSpecLib(search_context))[:prec_charge],
                getPrecursors(getSpecLib(search_context))[:sulfur_count],
                getIsoSplines(search_data),
                quad_func,
                getPrecursorTransmission(search_data),
                getIsotopes(search_data),
                params.n_frag_isotopes,
                params.max_frag_rank,
                rt_index,
                irt_start,
                irt_stop,
                (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
                precursors_passing = precursors_passing,
                isotope_err_bounds = params.isotope_err_bounds,
                block_size = 10000
            )
        end

        # Match peaks
        nmatches, nmisses = matchPeaks!(
            getIonMatches(search_data),
            getIonMisses(search_data),
            getIonTemplates(search_data),
            ion_idx,
            spectra[:mz_array][scan_idx],
            spectra[:intensity_array][scan_idx],
            getMassErrorModel(search_context, ms_file_idx),
            spectra[:highMz][scan_idx],
            UInt32(scan_idx),
            UInt32(ms_file_idx)
        )

        sort!(@view(getIonMatches(search_data)[1:nmatches]), by = x->(x.peak_ind, x.prec_id), alg=QuickSort)

        # Process matches
        if nmatches > 2
            # Build design matrix
            buildDesignMatrix!(
                Hs,
                getIonMatches(search_data),
                getIonMisses(search_data),
                nmatches,
                nmisses,
                getIdToCol(search_data)
            )

            # Handle array resizing
            if getIdToCol(search_data).size > length(weights)
                new_entries = getIdToCol(search_data).size - length(weights) + 1000
                resize!(weights, length(weights) + new_entries)
                resize!(getSpectralScores(search_data), length(getSpectralScores(search_data)) + new_entries)
                append!(getUnscoredPsms(search_data), [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])
            end

            # Initialize weights
            for i in 1:getIdToCol(search_data).size
                weights[getIdToCol(search_data)[getIdToCol(search_data).keys[i]]] = 
                    precursor_weights[getIdToCol(search_data).keys[i]]
            end

            # Solve deconvolution
            initResiduals!(residuals, Hs, weights)
            solveHuber!(
                Hs,
                residuals,
                weights,
                getHuberDelta(search_context),
                params.lambda,
                params.max_iter_newton,
                params.max_iter_bisection,
                params.max_iter_outer,
                params.accuracy_newton,
                params.accuracy_bisection,
                10.0f0,
                params.max_diff
            )

            # Record chromatogram points with weights
            for j in 1:prec_temp_size
                rt_idx += 1
                if rt_idx + 1 > length(chromatograms)
                    append!(chromatograms, Vector{ChromObject}(undef, 500000))
                end

                if !iszero(getIdToCol(search_data)[precs_temp[j]])
                    chromatograms[rt_idx] = ChromObject(
                        Float16(spectra[:retentionTime][scan_idx]),
                        weights[getIdToCol(search_data)[precs_temp[j]]],
                        scan_idx,
                        precs_temp[j]
                    )
                else
                    chromatograms[rt_idx] = ChromObject(
                        Float16(spectra[:retentionTime][scan_idx]),
                        zero(Float32),
                        scan_idx,
                        precs_temp[j]
                    )
                end
            end

            # Update precursor weights
            for i in 1:getIdToCol(search_data).size
                precursor_weights[getIdToCol(search_data).keys[i]] = 
                    weights[getIdToCol(search_data)[getIdToCol(search_data).keys[i]]]
            end

        else
            # Record zero intensity points for non-matching precursors
            for j in 1:prec_temp_size
                rt_idx += 1
                if rt_idx + 1 > length(chromatograms)
                    append!(chromatograms, Vector{ChromObject}(undef, 500000))
                end

                chromatograms[rt_idx] = ChromObject(
                    Float16(spectra[:retentionTime][scan_idx]),
                    zero(Float32),
                    scan_idx,
                    precs_temp[j]
                )
            end
        end

        # Reset arrays
        for i in 1:Hs.n
            getUnscoredPsms(search_data)[i] = eltype(getUnscoredPsms(search_data))()
        end
        reset!(getIdToCol(search_data))
        reset!(Hs)
    end

    return chromatograms[1:rt_idx]#DataFrame(@view(chromatograms[1:rt_idx]))
end
