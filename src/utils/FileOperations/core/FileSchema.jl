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
Schema management for file references.

Provides immutable schema representation and validation utilities
for tracking and validating column information in Arrow files.
"""

using DataFrames

#==========================================================
Schema Management
==========================================================#

"""
    FileSchema

Immutable representation of a file's column schema.
Provides safe access to columns with default support.
"""
struct FileSchema
    columns::Set{Symbol}  # Immutable set of columns
    
    # Constructor validates columns
    FileSchema(cols::Vector{Symbol}) = new(Set(cols))
end

"""
    has_column(schema::FileSchema, col::Symbol) -> Bool

Check if a column exists in the schema.
"""
has_column(schema::FileSchema, col::Symbol) = col in schema.columns

"""
    get_column_or_default(df::DataFrame, schema::FileSchema, col::Symbol, default)

Get a column from DataFrame or return default values if column doesn't exist.
"""
function get_column_or_default(df::DataFrame, schema::FileSchema, col::Symbol, default)
    has_column(schema, col) ? df[!, col] : fill(default, nrow(df))
end

"""
    validate_required_columns(schema::FileSchema, required::Set{Symbol})

Validate that all required columns exist in the schema.
Throws an error with missing columns if validation fails.
"""
function validate_required_columns(schema::FileSchema, required::Set{Symbol})
    missing_cols = setdiff(required, schema.columns)
    isempty(missing_cols) || error("Missing required columns: $missing_cols")
    return nothing
end

# Export schema functions
export FileSchema, has_column, get_column_or_default, validate_required_columns