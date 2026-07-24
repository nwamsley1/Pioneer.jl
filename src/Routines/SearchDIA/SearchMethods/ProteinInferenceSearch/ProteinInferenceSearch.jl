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
    ProteinInferenceSearch

Search method for annotating integrated passing precursor tables with inferred
protein groups and protein-quant eligibility flags.
"""
struct ProteinInferenceSearch <: SearchMethod end

mutable struct ProteinInferenceSearchResults <: SearchResults
    protein_ambiguity_candidates::Dict{UInt32, Vector{ProteinKey}}
    protein_peptide_opportunities::Dict{
        ProteinKey,
        ProteinPeptideOpportunityCounts
    }
end

struct ProteinInferenceSearchParameters <: SearchParameters
    global_inference::Bool
    function ProteinInferenceSearchParameters(params::PioneerParameters)
        new(Bool(params.protein_scoring.global_protein_inference))
    end
end

get_parameters(::ProteinInferenceSearch, params::Any) = ProteinInferenceSearchParameters(params)

function init_search_results(::ProteinInferenceSearchParameters, search_context::SearchContext)
    return ProteinInferenceSearchResults(
        Dict{UInt32, Vector{ProteinKey}}(),
        Dict{ProteinKey, ProteinPeptideOpportunityCounts}()
    )
end

function process_file!(
    results::ProteinInferenceSearchResults,
    params::ProteinInferenceSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    return results
end

function process_search_results!(
    results::ProteinInferenceSearchResults,
    params::ProteinInferenceSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    return nothing
end

function reset_results!(results::ProteinInferenceSearchResults)
    return nothing
end

function summarize_results!(
    results::ProteinInferenceSearchResults,
    params::ProteinInferenceSearchParameters,
    search_context::SearchContext
)
    indexed_paths = get_all_indexed_paths(getPassingPsms, search_context)
    if isempty(indexed_paths)
        empty!(results.protein_ambiguity_candidates)
        empty!(results.protein_peptide_opportunities)
        store_results!(search_context, ProteinInferenceSearch, results)
        return nothing
    end

    passing_refs = [PSMFileReference(path) for (_, path) in indexed_paths]
    inference_summary = run_protein_inference!(search_context;
        passing_refs = passing_refs,
        global_inference = params.global_inference)
    results.protein_ambiguity_candidates =
        inference_summary.protein_ambiguity_candidates
    results.protein_peptide_opportunities =
        inference_summary.protein_peptide_opportunities
    store_results!(search_context, ProteinInferenceSearch, results)

    return nothing
end
