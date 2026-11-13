# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

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
    validate_required_columns(schema::FileSchema, required::Set{Symbol})

Validate that all required columns exist in the schema.
Throws an error with missing columns if validation fails.
"""
function validate_required_columns(schema::FileSchema, required::Set{Symbol})
    missing_cols = setdiff(required, schema.columns)
    isempty(missing_cols) || error("Missing required columns: $missing_cols")
    return nothing
end
