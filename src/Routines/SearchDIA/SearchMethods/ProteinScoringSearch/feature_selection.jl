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
    remove_zero_variance_columns!(feature_names::Vector{Symbol}, df::AbstractDataFrame)

Remove columns with zero variance from `feature_names` to prevent numerical
issues during model fitting. This includes constant columns, columns with all
missing values, and columns containing Inf/NaN values.

# Arguments
- `feature_names`: Vector of feature symbols to filter
- `df`: DataFrame containing the feature data

# Returns
- `feature_names`: Filtered vector with problematic features removed
"""
function remove_zero_variance_columns!(feature_names::Vector{Symbol}, df::AbstractDataFrame)
    removed_features = Symbol[]

    # Filter out problematic columns
    filter!(feature_names) do feature
        if !hasproperty(df, feature)
            push!(removed_features, feature)
            return false
        end

        col_data = df[!, feature]

        # Check for columns with all missing values
        if all(ismissing, col_data)
            push!(removed_features, feature)
            return false
        end

        # Check for Inf or NaN values
        if any(x -> !ismissing(x) && (isinf(x) || isnan(x)), col_data)
            push!(removed_features, feature)
            return false
        end

        # Check for zero variance (constant columns)
        non_missing_data = collect(skipmissing(col_data))
        if isempty(non_missing_data) || length(unique(non_missing_data)) <= 1 || var(non_missing_data) ≈ 0.0
            push!(removed_features, feature)
            return false
        end

        return true
    end

    # Log removed features if any
    if !isempty(removed_features)
        @user_warn "Removed $(length(removed_features)) problematic features to prevent singular matrix" removed_features = removed_features
    end

    return feature_names
end

"""
    run_level_protein_feature_names()

Return the fixed feature set used for run-level protein rescoring.
"""
function run_level_protein_feature_names()
    return Symbol[
        :pg_score,
        :ambiguous_pg_score,
        :shared_peptide_coverage_logit,
        :shared_coverage_log_ratio,
        :peptide_coverage_logit,
        :all_peptide_coverage_logit,
        :any_common_peps,
        :coverage_log_ratio,
        :precursor_consensus_prefix_shape,
        :shared_precursor_consensus_prefix_shape,
        :single_non_mbr_peptide,
        :single_non_mbr_prefix_shape,
        :mbr_recovered_peptides,
        :mbr_only_protein
    ]
end

"""
    smoothed_coverage_logit(n_peptides, n_possible_peptides)

Compute a smoothed logit transform of peptide coverage to spread out the low end.
"""
function smoothed_coverage_logit(n_peptides::Integer, n_possible_peptides::Integer)::Float32
    n_possible_peptides <= 0 && return 0.0f0
    observed = Float32(n_peptides)
    possible = Float32(n_possible_peptides)
    return log((observed + 0.5f0) / (possible - observed + 0.5f0))
end

"""
    shared_coverage_log_ratio(
        n_shared_peptides,
        n_possible_shared_peptides,
        n_unique_peptides,
        n_possible_unique_peptides
    )

Compare shared- and unique-peptide detection rates within a protein group.
Each rate uses a 0.5 pseudocount before their ratio is log transformed.
Return zero when either peptide class has no theoretical opportunities.
"""
function shared_coverage_log_ratio(
    n_shared_peptides::Integer,
    n_possible_shared_peptides::Integer,
    n_unique_peptides::Integer,
    n_possible_unique_peptides::Integer
)::Float32
    if n_possible_shared_peptides <= 0 || n_possible_unique_peptides <= 0
        return 0.0f0
    end

    shared_rate =
        (Float32(n_shared_peptides) + 0.5f0) /
        (Float32(n_possible_shared_peptides) + 1.0f0)
    unique_rate =
        (Float32(n_unique_peptides) + 0.5f0) /
        (Float32(n_possible_unique_peptides) + 1.0f0)
    return log(shared_rate / unique_rate)
end
