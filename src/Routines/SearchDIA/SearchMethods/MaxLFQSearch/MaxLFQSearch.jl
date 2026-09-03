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
    # Run-to-run normalization on/off (the median-spline normalizer's tuning
    # constants — n_rt_bins, spline_n_knots — are hardcoded below since they
    # have no shipping override).
    run_to_run_normalization::Bool

    # LFQ parameters
    q_value_threshold::Float32
    batch_size::Int64
    min_peptides::Int64

    # Chunked merge parameters
    max_chunk_size_mb::Float64

    # Output parameters
    write_csv::Bool
    delete_temp::Bool
    params::Any  # Store full parameters for reference

    function MaxLFQSearchParameters(params::PioneerParameters)
        output_params = params.output
        global_params = params.global_settings
        maxLFQ_params = params.maxLFQ
        protein_scoring_params = params.protein_scoring

        new(
            Bool(maxLFQ_params.run_to_run_normalization),
            _resolve_q_value_threshold(global_params),
            Int64(100000),  # Default batch size
            Int64(protein_scoring_params.min_peptides),
            Float64(get(maxLFQ_params, :max_chunk_size_mb, 1024)),
            Bool(output_params.write_csv),
            Bool(output_params.delete_temp),
            params  # Store full parameters
        )
    end
end

# Hardcoded run-to-run normalization tuning constants (formerly
# global.normalization.n_rt_bins / spline_n_knots). Values match what every
# shipping config used.
const MAXLFQ_NORM_N_RT_BINS = 100
const MAXLFQ_NORM_SPLINE_N_KNOTS = 7


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
    # Get paths
    temp_folder = joinpath(getDataOutDir(search_context), "temp_data")
    passing_psms_folder = joinpath(temp_folder, "passing_psms")
    qc_plot_folder = joinpath(getDataOutDir(search_context), "qc_plots")
    precursors_long_path = joinpath(getDataOutDir(search_context), "precursors_long.arrow")
    protein_long_path = joinpath(getDataOutDir(search_context), "protein_groups_long.arrow")
    spec_lib = getSpecLib(search_context)
    proteins = getProteins(spec_lib)
    precursors = getPrecursors(spec_lib)
    output_schema_policy = getOutputSchemaPolicy(spec_lib)
    # Collect every populated passing PSM path. Skip 0-row files (e.g. files
    # that yielded no PSMs after upstream filtering — they lack the columns
    # MaxLFQ expects).
    indexed_paths = get_all_indexed_paths(getPassingPsms, search_context)
    existing_passing_psm_paths = String[]
    existing_passing_psm_run_ids = UInt32[]
    for (file_idx, path) in indexed_paths
        ref = PSMFileReference(path)
        if row_count(ref) > 0
            push!(existing_passing_psm_paths, path)
            push!(existing_passing_psm_run_ids, UInt32(file_idx))
        end
    end

    # Get all file names (needed for LFQ indexing by experiment ID)
    all_file_names = collect(getMSData(search_context).file_id_to_name)

    if isempty(existing_passing_psm_paths)
        @user_warn "No PSM files found for MaxLFQ analysis"
        return nothing
    end

    psm_refs = [PSMFileReference(path) for path in existing_passing_psm_paths]

    if !params.params.output.write_decoys
        decoy_filter_pipeline = TransformPipeline() |>
            filter_rows(row -> row.target; desc = "keep_target_rows")
        apply_pipeline!(psm_refs, decoy_filter_pipeline)
    end

    # Normalize the same integrated precursor tables that will feed MaxLFQ.
    # Only run when the result will actually be consumed: MaxLFQ's quant
    # column selector + WriteOutputs both gate on run_to_run_normalization.
    # IntegrateChromatograms already populates :peak_area_normalized with
    # placeholder zeros, so a missing column never breaks downstream.
    if params.run_to_run_normalization
        precursor_scoring_results = get_results(
            search_context,
            PrecursorScoringSearch,
        )
        run_similarity_atlas =
            precursor_scoring_results isa PrecursorScoringSearchResults ?
            precursor_scoring_results.run_similarity[] : nothing
        normalizeQuant(
            [file_path(ref) for ref in psm_refs],
            :peak_area,
            N = MAXLFQ_NORM_N_RT_BINS,
            spline_n_knots = MAXLFQ_NORM_SPLINE_N_KNOTS,
            run_ids = existing_passing_psm_run_ids,
            run_similarity_atlas = run_similarity_atlas,
        )
    end

    # Keep inferred protein groups contiguous for MaxLFQ. Ascending group order
    # places `missing` (non-inferred peptides) after all inferred groups; the
    # remaining tie-breakers retain their historical descending order.
    sort_keys = (:inferred_protein_group, :target, :entrapment_group_id, :precursor_idx)
    sort_file_by_keys!(psm_refs, :inferred_protein_group, :target, :entrapment_group_id, :precursor_idx;
                       reverse=[false, true, true, true], parallel=true )

    # Chunked merge: split into protein-group-aligned chunks bounded by max_chunk_size_mb
    chunk_dir = joinpath(temp_folder, "merge_chunks")
    max_chunk_bytes = round(Int, params.max_chunk_size_mb * 1_000_000)
    chunk_refs = stream_sorted_merge_chunked(
        psm_refs, chunk_dir, :inferred_protein_group, sort_keys...;
        batch_size=1_000_000, reverse=[false, true, true, true],
        max_chunk_bytes=max_chunk_bytes
    )

    # Validate first chunk has required columns for MaxLFQ
    validate_maxlfq_input(first(chunk_refs))

    # Validate MaxLFQ parameters
    validate_maxlfq_parameters(Dict(
        :q_value_threshold => params.q_value_threshold,
        :batch_size => params.batch_size,
        :min_peptides => params.min_peptides
    ))

    # Chunked precursor CSV writing (bounded memory per chunk)
    @user_info "Writing precursor tables..."
    precursors_wide_path = writePrecursorCSV_chunked(
        chunk_refs,
        getDataOutDir(search_context),
        all_file_names,
        params.run_to_run_normalization,
        proteins,
        output_schema_policy = output_schema_policy,
        write_csv = params.write_csv
    )

    @user_info "Performing MaxLFQ..."
    # Chunked MaxLFQ protein quantification (bounded memory per chunk)
    precursor_quant_col = params.run_to_run_normalization ? :peak_area_normalized : :peak_area
    LFQ_chunked(
        chunk_refs,
        protein_long_path,
        precursor_quant_col,
        all_file_names,
        getSequence(precursors),
        getStructuralMods(precursors),
        getIsotopicMods(precursors),
        params.q_value_threshold,
        build_accession_to_species(precursors),
        output_schema_policy = output_schema_policy,
        batch_size = params.batch_size
    )
    chunk_paths = [file_path(ref) for ref in chunk_refs]

    # Create FileReference for output metadata tracking
    protein_ref = ProteinQuantFileReference(protein_long_path)
    @user_info "MaxLFQ: $(n_protein_groups(protein_ref)) protein groups across $(n_experiments(protein_ref)) experiments"

    @user_info "Writing protein group results..."
    # Create wide format protein table (protein groups table is small, no chunking needed)
    proteins_wide_path = writeProteinGroupsCSV(
        results.proteins_long_path,
        getSequence(precursors),
        getIsotopicMods(precursors),
        getStructuralMods(precursors),
        getCharge(precursors),
        all_file_names,
        proteins,
        output_schema_policy = output_schema_policy,
        write_csv = params.write_csv
    )

    # Concatenate chunks into a final Arrow file-format export for QC plots.
    @debug_l1 "Concatenating chunks to precursors_long.arrow..."
    isfile(precursors_long_path) && rm(precursors_long_path)
    open(Arrow.Writer, precursors_long_path; file=true) do arrow_writer
        for chunk_ref in chunk_refs
            let tbl = Arrow.Table(file_path(chunk_ref))
                # Dictionary-encode the repeated string columns here, at the final write only --
                # see OUTPUT_DICT_ENCODED_COLUMNS. -14.3% on this file, values unchanged.
                Arrow.write(arrow_writer, dict_encode_output_columns(
                    enabled_output_table(output_schema_policy, :precursors, tbl)))
            end
        end
    end
    chunk_refs = nothing
    GC.gc()

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
        params,
        all_file_names
    )

    # Cleanup chunk files after dropping Arrow.Table references.
    if isdir(chunk_dir)
        GC.gc()
        try
            for chunk_path in chunk_paths
                isfile(chunk_path) && safeRm(chunk_path; force=true)
            end
            rm(chunk_dir; recursive=true, force=true)
        catch e
            # On Windows, a chunk file may still be briefly locked by Arrow mmap state.
            @debug_l1 "merge_chunks cleanup deferred: $e"
        end
    end
    chunk_paths = nothing

    if params.delete_temp
        @user_info "Removing temporary data..."
        temp_path = joinpath(getDataOutDir(search_context), "temp_data")
        GC.gc()
        try
            isdir(temp_path) && rm(temp_path; recursive=true, force=true)
        catch e
            @debug_l1 "Could not fully remove temp_data (likely lingering file handles): $e"
        end
    end

    return nothing
end
