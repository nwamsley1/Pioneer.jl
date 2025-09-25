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
Create QC plots showing quantification metrics.
"""
function create_qc_plots(
    precursors_path::String,
    precursors_long_path::String,
    proteins_path::String,
    search_context::SearchContext,
    precursors::LibraryPrecursors,
    params::Any,
    successful_file_names::Vector{String}
)
    # Create plots showing:
    # - Normalization factors
    # - Missing value patterns
    # - CV distributions
    # - Dynamic range
    # Implementation depends on plotting library
    @user_info "Generating final QC plots"
    # Get only successful files for QC plots
    all_file_paths = collect(getFilePaths(getMSData(search_context)))
    valid_file_indices = get_valid_file_indices(search_context)

    # Filter to get only successful files
    # Note: successful_file_names parameter is actually all_file_names from MaxLFQSearch
    filtered_file_names = [successful_file_names[i] for i in valid_file_indices]
    filtered_file_paths = [all_file_paths[i] for i in valid_file_indices]

    qcPlots(
        precursors_path,
        precursors_long_path,
        proteins_path,
        params.params,
        precursors,
        filtered_file_names,  # Only successful files
        joinpath(getDataOutDir(search_context), "qc_plots"),
        filtered_file_paths,  # Only successful files
        getIrtRtMap(search_context),
        search_context.mass_error_model,
        valid_file_indices  # Pass actual file indices to correctly lookup in dicts
    )
end

"""
Get protein group q-value interpolation function.
"""
function getPGQValueInterp(search_context::SearchContext)
    # Implementation to get protein group q-value interpolation
end

