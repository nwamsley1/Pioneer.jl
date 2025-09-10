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
    MaxLFQSearch

Search method for performing MaxLFQ normalization and protein quantification.

This search:
1. Normalizes quantitative values across runs
2. Performs MaxLFQ protein quantification
3. Generates long and wide format results
4. Creates QC plots
"""
struct MaxLFQSearch <: SearchMethod end

# Note: FileReferences, SearchResultReferences, and FileOperations are already
# included by importScripts.jl - no need to include them here

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
    results::MaxLFQSearchResults,
    params::MaxLFQSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) 
    # Check if file should be skipped due to previous failure
    if check_and_skip_failed_file(search_context, ms_file_idx, "MaxLFQSearch")
        return results  # Return early with unchanged results
    end
    
    # No per-file processing needed
    return results
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
        # Normalize quantitative values
        normalizeQuant(
            passing_psms_folder,
            :peak_area,
            N = params.n_rt_bins,
            spline_n_knots = params.spline_n_knots
        )

        # Get PSM paths from MSData and create references
        passing_psm_paths = getPassingPsms(getMSData(search_context))
        psm_refs = [PSMFileReference(path) for path in passing_psm_paths]
        
        # Ensure all PSM files are sorted correctly for MaxLFQ
        sort_keys = (:inferred_protein_group, :target, :entrapment_group_id, :precursor_idx)
        #for (i, psm_ref) in ProgressBar(enumerate(psm_refs))
        #    if !is_sorted_by(psm_ref, sort_keys...)
        #        sort_file_by_keys!(psm_ref, sort_keys...)
        #    end
        #end
        sort_file_by_keys!(psm_refs, :inferred_protein_group, :target, :entrapment_group_id, :precursor_idx;
                           reverse=[true, true, true, true], parallel=true )
        
        # Use reference-based merge with 4 sort keys (all descending)
        @time merged_psm_ref = stream_sorted_merge(psm_refs, precursors_long_path, sort_keys...;
                                           batch_size=1000000, reverse=true)

        # Verify the merged file is sorted (reference-based merge guarantees this)

        # Add FileReference validation for MaxLFQ input
        validate_maxlfq_input(merged_psm_ref)
        
        # Validate MaxLFQ parameters
        validate_maxlfq_parameters(Dict(
            :q_value_threshold => params.q_value_threshold,
            :batch_size => params.batch_size,
            :min_peptides => params.min_peptides
        ))
        # Create wide format precursor table
        precursors_wide_path = writePrecursorCSV(
            precursors_long_path,
            sort(collect(getParsedFileNames(getMSData(search_context)))),
            params.run_to_run_normalization,
            getProteins(getSpecLib(search_context)),
            write_csv = params.write_csv
        )

        @user_info "Performing MaxLFQ..."
        # Perform MaxLFQ protein quantification
        precursor_quant_col = params.run_to_run_normalization ? :peak_area_normalized : :peak_area

        # Use FileReference-based LFQ with TransformPipeline preprocessing
        LFQ(
            merged_psm_ref,  # Use FileReference instead of DataFrame
            protein_long_path,
            precursor_quant_col,
            collect(getFileIdToName(getMSData(search_context))),
            params.q_value_threshold,
            batch_size = params.batch_size
        )
        
        # Create FileReference for output metadata tracking
        protein_ref = ProteinQuantFileReference(protein_long_path)
        @user_info "MaxLFQ completed - output_file: $protein_long_path, n_protein_groups: $(n_protein_groups(protein_ref)), n_experiments: $(n_experiments(protein_ref))"

        @user_info "Writing protein group results..."
        # Create wide format protein table
        precursors = getPrecursors(getSpecLib(search_context))
        proteins_wide_path = writeProteinGroupsCSV(
            results.proteins_long_path,
            getSequence(precursors),
            getIsotopicMods(precursors),
            getStructuralMods(precursors),
            getCharge(precursors),
            sort(collect(values(getFileIdToName(getMSData(search_context))))),
            getProteins(getSpecLib(search_context)),
            write_csv = params.write_csv
        )

        @user_info "Creating QC plots..."
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
            @user_info "Removing temporary data..."
            temp_path = joinpath(getDataOutDir(search_context), "temp_data")
            GC.gc()
            isdir(temp_path) && rm(temp_path; recursive=true, force=true)
        end

    catch e
        @error "MaxLFQ analysis failed" exception=e
        rethrow(e)
    end

    return nothing
end
