"""
    MaxLFQSearch

Search method for performing MaxLFQ normalization and protein quantification.

This search:
1. Normalizes quantitative values across runs
2. Performs MaxLFQ protein quantification
3. Generates long and wide format results
4. Creates QC plots
"""
struct MaxLFQSearch <: SearchMethod end

#==========================================================
Type Definitions 
==========================================================#

"""
Results container for MaxLFQ search.
"""
struct MaxLFQSearchResults <: SearchResults
    precursors_long_path::String
    precursors_wide_path::String
    proteins_long_path::String
    proteins_wide_path::String
    normalized_quant::Dict{Int64, Float32}  # Normalization factors
end

"""
Parameters for MaxLFQ search.
"""

struct MaxLFQSearchParameters <: SearchParameters
    # Normalization parameters
    n_rt_bins::Int64
    spline_n_knots::Int64
    run_to_run_normalization::Bool
    
    # LFQ parameters
    q_value_threshold::Float32
    batch_size::Int64
    min_peptides::Int64

    # Output parameters
    write_csv::Bool
    delete_temp::Bool
    params::Any  # Store full parameters for reference

    function MaxLFQSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        norm_params = params.global_settings.normalization
        output_params = params.output
        global_params = params.global_settings
        maxLFQ_params = params.maxLFQ
        protein_inference_params = params.protein_inference
        
        new(
            Int64(norm_params.n_rt_bins),
            Int64(norm_params.spline_n_knots),
            Bool(maxLFQ_params.run_to_run_normalization),
            Float32(global_params.scoring.q_value_threshold),
            Int64(100000),  # Default batch size
            Int64(protein_inference_params.min_peptides),
            Bool(output_params.write_csv),
            Bool(output_params.delete_temp),
            params  # Store full parameters
        )
    end
end


#==========================================================
Interface Implementation
==========================================================#

get_parameters(::MaxLFQSearch, params::Any) = MaxLFQSearchParameters(params)

function init_search_results(::MaxLFQSearchParameters, search_context::SearchContext)
    return MaxLFQSearchResults(
        joinpath(getDataOutDir(search_context), "precursors_long.arrow"),
        joinpath(getDataOutDir(search_context), "precursors_wide.arrow"),
        joinpath(getDataOutDir(search_context), "protein_groups_long.arrow"),
        joinpath(getDataOutDir(search_context), "protein_groups_wide.arrow"),
        Dict{Int64, Float32}()
    )
end

"""
Process a single file for MaxLFQ analysis.
"""
function process_file!(
    ::MaxLFQSearchResults,
    ::MaxLFQSearchParameters,
    ::SearchContext,
    ::Int64,
    ::MassSpecData
) 
    # No per-file processing needed
    return nothing
end

"""
No per-file results processing needed.
"""
function process_search_results!(
    ::MaxLFQSearchResults,
    ::MaxLFQSearchParameters,
    ::SearchContext,
    ::Int64,
    ::MassSpecData
)
    return nothing
end

function reset_results!(::MaxLFQSearchResults)
    return nothing
end

"""
Perform MaxLFQ analysis across all files.
"""
function summarize_results!(
    results::MaxLFQSearchResults,
    params::MaxLFQSearchParameters,
    search_context::SearchContext
)
    try
        # Get paths
        temp_folder = joinpath(getDataOutDir(search_context), "temp_data")
        passing_psms_folder = joinpath(temp_folder, "passing_psms")
        qc_plot_folder = joinpath(getDataOutDir(search_context), "qc_plots")
        precursors_long_path = joinpath(getDataOutDir(search_context), "precursors_long.arrow")
        protein_long_path = joinpath(getDataOutDir(search_context), "protein_groups_long.arrow")
        @info "Performing intensity normalization..."
        # Normalize quantitative values
        normalizeQuant(
            passing_psms_folder,
            :peak_area,
            N = params.n_rt_bins,
            spline_n_knots = params.spline_n_knots
        )

        @info "Merging quantification tables..."
        # Merge quantification tables
        mergeSortedArrowTables(
            passing_psms_folder,
            precursors_long_path,
            (:inferred_protein_group, :precursor_idx),
            N = 1000000
        )

        @info "Writing precursor results..."
        # Create wide format precursor table
        precursors_wide_path = writePrecursorCSV(
            precursors_long_path,
            sort(collect(getParsedFileNames(getMSData(search_context)))),
            params.run_to_run_normalization,
            write_csv = params.write_csv
        )

        @info "Performing MaxLFQ..."
        # Perform MaxLFQ protein quantification
        precursor_quant_col = params.run_to_run_normalization ? :peak_area_normalized : :peak_area

        LFQ(
            DataFrame(Arrow.Table(precursors_long_path)),
            protein_long_path,
            precursor_quant_col,
            collect(getFileIdToName(getMSData(search_context))),
            params.q_value_threshold,
            search_context.global_pg_score_to_qval[],
            search_context.pg_score_to_qval[],#getPGQValueInterp(search_context),
            batch_size = params.batch_size,
            min_peptides = params.min_peptides
        )

        @info "Writing protein group results..."
        # Create wide format protein table
        precursors = getPrecursors(getSpecLib(search_context))
        proteins_wide_path = writeProteinGroupsCSV(
            results.proteins_long_path,
            getSequence(precursors),
            getIsotopicMods(precursors),
            getStructuralMods(precursors),
            getCharge(precursors),
            sort(collect(values(getFileIdToName(getMSData(search_context))))),
            write_csv = params.write_csv
        )

        @info "Creating QC plots..."
        # Create QC plots
        qc_plot_path = joinpath(qc_plot_folder, "QC_PLOTS.pdf")
        isfile(qc_plot_path) && rm(qc_plot_path)
        create_qc_plots(
            precursors_wide_path,
            precursors_long_path,
            proteins_wide_path,
            search_context,
            precursors,
            params
        )

        if params.delete_temp
            @info "Removing temporary data..."
            temp_path = joinpath(getDataOutDir(search_context), "temp_data")
            isdir(temp_path) && rm(temp_path; recursive=true, force=true)
        end

    catch e
        @error "MaxLFQ analysis failed" exception=e
        rethrow(e)
    end

    return nothing
end
