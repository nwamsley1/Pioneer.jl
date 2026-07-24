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
    ProteinScoringSearch

Search method for building protein-group tables, protein rescoring, and
protein-level q-value calculation after chromatogram integration and protein
inference have completed.
"""
struct ProteinScoringSearch <: SearchMethod end

struct ProteinScoringSearchResults <: SearchResults end

# PEP threshold for mining negatives during semi-supervised iterations.
# Hardcoded — develop surfaced this as `min_PEP_neg_threshold_itr` but every
# shipping config used 0.90.
const PROTEIN_SCORING_MIN_PEP_NEG_THRESHOLD_ITR = 0.90f0

struct ProteinScoringSearchParameters <: SearchParameters
    max_in_memory_table_mb::Float64
    min_peptides::Int64
    q_value_threshold::Float32
    min_pep_neg_threshold_itr::Float32
    q_value_interpolation_points_per_bin::Int64
    write_qc_plots::Bool

    function ProteinScoringSearchParameters(params::PioneerParameters)
        ml_params = params.optimization.machine_learning
        global_params = params.global_settings
        protein_scoring_params = params.protein_scoring

        new(
            Float64(ml_params.max_psm_memory_mb),
            Int64(protein_scoring_params.min_peptides),
            _resolve_q_value_threshold(global_params),
            PROTEIN_SCORING_MIN_PEP_NEG_THRESHOLD_ITR,
            Int64(ml_params.pep_bin_size),
            Bool(protein_scoring_params.write_qc_plots)
        )
    end
end

get_parameters(::ProteinScoringSearch, params::Any) = ProteinScoringSearchParameters(params)

function init_search_results(::ProteinScoringSearchParameters, search_context::SearchContext)
    return ProteinScoringSearchResults()
end

function process_file!(
    results::ProteinScoringSearchResults,
    params::ProteinScoringSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    return results
end

function process_search_results!(
    results::ProteinScoringSearchResults,
    params::ProteinScoringSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    return nothing
end

function reset_results!(results::ProteinScoringSearchResults)
    return nothing
end


function summarize_results!(
    ::ProteinScoringSearchResults,
    params::ProteinScoringSearchParameters,
    search_context::SearchContext
)
    indexed_paths = get_all_indexed_paths(getPassingPsms, search_context)
    isempty(indexed_paths) && return nothing

    # Files with zero rows weren't annotated by ProteinInferenceSearch (the
    # pipeline short-circuits on empty inputs), so they lack the
    # `inferred_protein_group` column. Drop them here.
    passing_refs = PSMFileReference[]
    for (_, path) in indexed_paths
        ref = PSMFileReference(path)
        row_count(ref) > 0 && push!(passing_refs, ref)
    end
    isempty(passing_refs) && return nothing

    protein_inference_results = get_results(search_context, ProteinInferenceSearch)
    protein_ambiguity_candidates = if protein_inference_results isa ProteinInferenceSearchResults
        protein_inference_results.protein_ambiguity_candidates
    else
        Dict{UInt32, Vector{ProteinKey}}()
    end
    protein_peptide_opportunities = if protein_inference_results isa ProteinInferenceSearchResults
        protein_inference_results.protein_peptide_opportunities
    else
        Dict{ProteinKey, ProteinPeptideOpportunityCounts}()
    end

    run_protein_scoring!(
        search_context;
        passing_refs = passing_refs,
        protein_ambiguity_candidates = protein_ambiguity_candidates,
        protein_peptide_opportunities = protein_peptide_opportunities,
        max_in_memory_table_mb = params.max_in_memory_table_mb,
        q_value_threshold = params.q_value_threshold,
        min_peptides = params.min_peptides,
        write_qc_plots = params.write_qc_plots,
        min_pep_neg_threshold_itr = params.min_pep_neg_threshold_itr,
        q_value_interpolation_points_per_bin = params.q_value_interpolation_points_per_bin
    )
    empty!(protein_ambiguity_candidates)
    empty!(protein_peptide_opportunities)

    return nothing
end
