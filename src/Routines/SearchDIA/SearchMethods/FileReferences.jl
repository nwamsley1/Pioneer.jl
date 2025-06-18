"""
File reference types for managing Arrow files with metadata tracking.

These types provide lightweight references to files on disk with:
- Schema tracking (immutable column information)
- Sort state tracking (which keys the file is sorted by)
- Validation methods to ensure data integrity
"""

using Arrow, DataFrames, Tables

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

#==========================================================
File Reference Types
==========================================================#

"""
    PSMFileReference

Reference to a PSM (Peptide-Spectrum Match) Arrow file with metadata.
Tracks schema, sort state, and basic statistics.
"""
mutable struct PSMFileReference
    file_path::String
    schema::FileSchema            # Immutable schema
    sorted_by::Tuple{Vararg{Symbol}}  # Which keys it's sorted by
    row_count::Int64
    file_exists::Bool
    
    # Constructor that validates file and extracts metadata
    function PSMFileReference(file_path::String)
        if !isfile(file_path)
            return new(file_path, FileSchema(Symbol[]), (), 0, false)
        end
        
        # Read schema from file
        tbl = Arrow.Table(file_path)
        schema = FileSchema(collect(Symbol.(Tables.columnnames(tbl))))
        # Get row count from first column
        row_count = length(Tables.getcolumn(tbl, 1))
        
        new(file_path, schema, (), row_count, true)
    end
end

"""
    ProteinGroupFileReference

Reference to a protein group Arrow file with metadata.
Tracks schema, sort state, and basic statistics.
"""
mutable struct ProteinGroupFileReference  
    file_path::String
    schema::FileSchema
    sorted_by::Tuple{Vararg{Symbol}}
    row_count::Int64
    file_exists::Bool
    
    # Constructor that validates file and extracts metadata
    function ProteinGroupFileReference(file_path::String)
        if !isfile(file_path)
            return new(file_path, FileSchema(Symbol[]), (), 0, false)
        end
        
        # Read schema from file
        tbl = Arrow.Table(file_path)
        schema = FileSchema(collect(Symbol.(Tables.columnnames(tbl))))
        # Get row count from first column
        row_count = length(Tables.getcolumn(tbl, 1))
        
        new(file_path, schema, (), row_count, true)
    end
end

#==========================================================
Paired File References
==========================================================#

"""
    PairedSearchFiles

Enforces 1:1 pairing between PSM and protein group files.
Ensures consistency between related files.
"""
struct PairedSearchFiles
    psm_ref::PSMFileReference
    protein_ref::ProteinGroupFileReference
    ms_file_idx::Int64
    
    # Constructor validates pairing
    function PairedSearchFiles(psm_path::String, protein_path::String, idx::Int64)
        psm_ref = PSMFileReference(psm_path)
        protein_ref = ProteinGroupFileReference(protein_path)
        
        # Validate both files exist or both don't
        if psm_ref.file_exists != protein_ref.file_exists
            error("Inconsistent file existence: PSM exists=$(psm_ref.file_exists), " *
                  "Protein exists=$(protein_ref.file_exists)")
        end
        
        new(psm_ref, protein_ref, idx)
    end
    
    # Alternative constructor from existing references
    function PairedSearchFiles(psm_ref::PSMFileReference, 
                              protein_ref::ProteinGroupFileReference, 
                              idx::Int64)
        # Validate consistency
        if psm_ref.file_exists != protein_ref.file_exists
            error("Inconsistent file existence: PSM exists=$(psm_ref.file_exists), " *
                  "Protein exists=$(protein_ref.file_exists)")
        end
        
        new(psm_ref, protein_ref, idx)
    end
end

#==========================================================
Sort State Management
==========================================================#

"""
    mark_sorted!(ref::Union{PSMFileReference,ProteinGroupFileReference}, keys::Symbol...)

Mark a file reference as sorted by the specified keys.
"""
function mark_sorted!(ref::Union{PSMFileReference,ProteinGroupFileReference}, keys::Symbol...)
    ref.sorted_by = keys
    return ref
end

"""
    is_sorted_by(ref::Union{PSMFileReference,ProteinGroupFileReference}, keys::Symbol...) -> Bool

Check if a file is sorted by the specified keys in the exact order.
"""
function is_sorted_by(ref::Union{PSMFileReference,ProteinGroupFileReference}, keys::Symbol...)
    return ref.sorted_by == keys
end

"""
    ensure_sorted!(ref::PSMFileReference, keys::Symbol...)

Ensure a PSM file is sorted by the specified keys, sorting if necessary.
"""
function ensure_sorted!(ref::PSMFileReference, keys::Symbol...)
    if !is_sorted_by(ref, keys...)
        sort_file_by_keys!(ref, keys...)
    end
    return ref
end

"""
    ensure_sorted!(ref::ProteinGroupFileReference, keys::Symbol...)

Ensure a protein file is sorted by the specified keys, sorting if necessary.
"""
function ensure_sorted!(ref::ProteinGroupFileReference, keys::Symbol...)
    if !is_sorted_by(ref, keys...)
        sort_file_by_keys!(ref, keys...)
    end
    return ref
end

#==========================================================
File Operations
==========================================================#

"""
    validate_exists(ref::Union{PSMFileReference,ProteinGroupFileReference})

Validate that a referenced file exists and update the file_exists flag.
"""
function validate_exists(ref::Union{PSMFileReference,ProteinGroupFileReference})
    ref.file_exists = isfile(ref.file_path)
    ref.file_exists || error("File does not exist: $(ref.file_path)")
    return true
end

"""
    sort_file_by_keys!(ref::PSMFileReference, keys::Symbol...)

Sort a PSM file by the specified keys and update the reference.
"""
function sort_file_by_keys!(ref::PSMFileReference, keys::Symbol...)
    validate_exists(ref)
    
    # Validate all sort keys exist in schema
    for key in keys
        has_column(ref.schema, key) || error("Sort key '$key' not in schema")
    end
    
    # Read into DataFrame (creates a copy), sort, and write back
    df = DataFrame(Tables.columntable(Arrow.Table(ref.file_path)))
    sort!(df, collect(keys))
    
    # Write sorted data
    Arrow.write(ref.file_path, df)
    
    # Update metadata
    mark_sorted!(ref, keys...)
    return ref
end

"""
    sort_file_by_keys!(ref::ProteinGroupFileReference, keys::Symbol...)

Sort a protein group file by the specified keys and update the reference.
"""
function sort_file_by_keys!(ref::ProteinGroupFileReference, keys::Symbol...)
    validate_exists(ref)
    
    # Validate all sort keys exist in schema
    for key in keys
        has_column(ref.schema, key) || error("Sort key '$key' not in schema")
    end
    
    # Read into DataFrame (creates a copy), sort, and write back
    df = DataFrame(Tables.columntable(Arrow.Table(ref.file_path)))
    sort!(df, collect(keys))
    
    # Write sorted data
    Arrow.write(ref.file_path, df)
    
    # Update metadata
    mark_sorted!(ref, keys...)
    return ref
end

"""
    create_psm_reference(file_path::String) -> PSMFileReference

Create a PSM file reference, extracting metadata from the file.
"""
create_psm_reference(file_path::String) = PSMFileReference(file_path)

"""
    create_protein_reference(file_path::String) -> ProteinGroupFileReference

Create a protein group file reference, extracting metadata from the file.
"""
create_protein_reference(file_path::String) = ProteinGroupFileReference(file_path)

"""
    describe_reference(ref::Union{PSMFileReference,ProteinGroupFileReference})

Print a human-readable description of a file reference.
"""
function describe_reference(ref::Union{PSMFileReference,ProteinGroupFileReference})
    println("File Reference:")
    println("  Path: $(ref.file_path)")
    println("  Exists: $(ref.file_exists)")
    println("  Sorted by: $(isempty(ref.sorted_by) ? "none" : join(ref.sorted_by, ", "))")
    println("  Row count: $(ref.row_count)")
    println("  Columns: $(join(sort(collect(ref.schema.columns)), ", "))")
end

# Export all public types and functions
export FileSchema, PSMFileReference, ProteinGroupFileReference, PairedSearchFiles,
       has_column, get_column_or_default, validate_required_columns,
       mark_sorted!, is_sorted_by, ensure_sorted!, validate_exists,
       sort_file_by_keys!, create_psm_reference, create_protein_reference,
       describe_reference