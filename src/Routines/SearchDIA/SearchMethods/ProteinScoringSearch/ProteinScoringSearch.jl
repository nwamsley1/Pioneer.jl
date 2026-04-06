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

struct ProteinScoringSearchParameters <: SearchParameters
    max_in_memory_table_mb::Float64
    min_peptides::Int64
    q_value_threshold::Float32
    min_pep_neg_threshold_itr::Float32
    q_value_interpolation_points_per_bin::Int64
    write_qc_plots::Bool
    log_feature_importance::Bool

    function ProteinScoringSearchParameters(params::PioneerParameters)
        ml_params = params.optimization.machine_learning
        global_params = params.global_settings
        protein_scoring_params = params.protein_scoring

        new(
            Float64(ml_params.max_in_memory_table_mb),
            Int64(protein_scoring_params.min_peptides),
            Float32(global_params.scoring.q_value_threshold),
            Float32(ml_params.min_PEP_neg_threshold_itr),
            Int64(ml_params.q_value_interpolation_points_per_bin),
            Bool(protein_scoring_params.write_qc_plots),
            Bool(protein_scoring_params.log_feature_importance)
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
    valid_file_data = get_valid_file_paths(search_context, getPassingPsms)
    isempty(valid_file_data) && return nothing

    passing_refs = [PSMFileReference(path) for (_, path) in valid_file_data]

    run_protein_scoring!(
        search_context;
        passing_refs = passing_refs,
        max_in_memory_table_mb = params.max_in_memory_table_mb,
        q_value_threshold = params.q_value_threshold,
        min_peptides = params.min_peptides,
        write_qc_plots = params.write_qc_plots,
        log_feature_importance = params.log_feature_importance,
        min_pep_neg_threshold_itr = params.min_pep_neg_threshold_itr,
        q_value_interpolation_points_per_bin = params.q_value_interpolation_points_per_bin
    )

    return nothing
end
