"""
File reference types for managing Arrow files with metadata tracking.

These types provide lightweight references to files on disk with:
- Schema tracking (immutable column information)
- Sort state tracking (which keys the file is sorted by)
- Validation methods to ensure data integrity
"""

using Arrow, DataFrames, Tables

#==========================================================
Abstract Types
==========================================================#

"""
    FileReference

Abstract base type for all file references.
Provides common interface for file metadata and operations.
"""
abstract type FileReference end

# Common interface that all FileReferences must implement
file_path(ref::FileReference) = ref.file_path
schema(ref::FileReference) = ref.schema
sorted_by(ref::FileReference) = ref.sorted_by
row_count(ref::FileReference) = ref.row_count
exists(ref::FileReference) = ref.file_exists

# Specialized accessors for ProteinQuantFileReference
n_protein_groups(ref::ProteinQuantFileReference) = ref.n_protein_groups
n_experiments(ref::ProteinQuantFileReference) = ref.n_experiments

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
    PSMFileReference <: FileReference

Reference to a PSM (Peptide-Spectrum Match) Arrow file with metadata.
Tracks schema, sort state, and basic statistics.
"""
mutable struct PSMFileReference <: FileReference
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
    ProteinGroupFileReference <: FileReference

Reference to a protein group Arrow file with metadata.
Tracks schema, sort state, and basic statistics.
"""
mutable struct ProteinGroupFileReference <: FileReference
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

"""
    ProteinQuantFileReference <: FileReference

Reference to a protein quantification Arrow file with metadata.
Specialized for MaxLFQ quantification results with additional metadata.
"""
mutable struct ProteinQuantFileReference <: FileReference
    file_path::String
    schema::FileSchema
    sorted_by::Tuple{Vararg{Symbol}}
    row_count::Int64
    file_exists::Bool
    n_protein_groups::Int64
    n_experiments::Int64
    
    # Constructor that validates file and extracts metadata
    function ProteinQuantFileReference(file_path::String)
        if !isfile(file_path)
            return new(file_path, FileSchema(Symbol[]), (), 0, false, 0, 0)
        end
        
        # Read schema from file
        tbl = Arrow.Table(file_path)
        schema = FileSchema(collect(Symbol.(Tables.columnnames(tbl))))
        # Get row count from first column
        row_count = length(Tables.getcolumn(tbl, 1))
        
        # Estimate protein groups and experiments for MaxLFQ context
        n_protein_groups = 0
        n_experiments = 0
        
        # Try to extract meaningful counts if the file has appropriate columns
        if has_column(schema, :protein)
            protein_col = Tables.getcolumn(tbl, :protein)
            n_protein_groups = length(unique(skipmissing(protein_col)))
        end
        
        if has_column(schema, :experiments) || has_column(schema, :file_name)
            exp_col = has_column(schema, :experiments) ? 
                     Tables.getcolumn(tbl, :experiments) : 
                     Tables.getcolumn(tbl, :file_name)
            n_experiments = length(unique(skipmissing(exp_col)))
        end
        
        new(file_path, schema, (), row_count, true, n_protein_groups, n_experiments)
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
    mark_sorted!(ref::FileReference, keys::Symbol...)

Mark a file reference as sorted by the specified keys.
"""
function mark_sorted!(ref::FileReference, keys::Symbol...)
    ref.sorted_by = keys
    return ref
end

"""
    is_sorted_by(ref::FileReference, keys::Symbol...) -> Bool

Check if a file is sorted by the specified keys in the exact order.
"""
function is_sorted_by(ref::FileReference, keys::Symbol...)
    return sorted_by(ref) == keys
end

"""
    ensure_sorted!(ref::FileReference, keys::Symbol...)

Ensure a file is sorted by the specified keys, sorting if necessary.
"""
function ensure_sorted!(ref::FileReference, keys::Symbol...)
    if !is_sorted_by(ref, keys...)
        sort_file_by_keys!(ref, keys...)
    end
    return ref
end

"""
    Base.sort!(ref::FileReference, cols::Vector{Symbol}; kwargs...)

Override Base.sort! to automatically track sort state for FileReferences.
This ensures that any sorting operation updates the FileReference metadata.
"""
function Base.sort!(ref::FileReference, cols::Vector{Symbol}; 
                   rev::Union{Bool, Vector{Bool}}=false, kwargs...)
    # Note: This is a placeholder that delegates to sort_file_by_keys!
    # The actual implementation is in FileOperations.jl
    # We just ensure the interface is consistent with Base.sort!
    
    # Convert single bool to appropriate format for sort_file_by_keys!
    reverse = if rev isa Bool
        rev  # sort_file_by_keys! expects a single bool
    else
        # If mixed directions, we can't use sort_file_by_keys! directly
        error("Mixed sort directions not yet supported. Use all ascending or all descending.")
    end
    
    # Delegate to the existing function
    sort_file_by_keys!(ref, cols...; reverse=reverse)
    
    return ref
end

#==========================================================
File Operations
==========================================================#

"""
    validate_exists(ref::FileReference)

Validate that a referenced file exists and update the file_exists flag.
"""
function validate_exists(ref::FileReference)
    ref.file_exists = isfile(file_path(ref))
    exists(ref) || error("File does not exist: $(file_path(ref))")
    return true
end

"""
    sort_file_by_keys!(ref::FileReference, keys::Symbol...)

Sort a file by the specified keys and update the reference.
"""
function sort_file_by_keys!(ref::FileReference, keys::Symbol...)
    validate_exists(ref)
    
    # Validate all sort keys exist in schema
    for key in keys
        has_column(schema(ref), key) || error("Sort key '$key' not in schema")
    end
    
    # Read into DataFrame (creates a copy), sort, and write back
    df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
    sort!(df, collect(keys))
    
    # Write sorted data
    Arrow.write(file_path(ref), df)
    
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
    describe_reference(ref::FileReference)

Print a human-readable description of a file reference.
"""
function describe_reference(ref::FileReference)
    println("File Reference ($(typeof(ref))):")
    println("  Path: $(file_path(ref))")
    println("  Exists: $(exists(ref))")
    println("  Sorted by: $(isempty(sorted_by(ref)) ? "none" : join(sorted_by(ref), ", "))")
    println("  Row count: $(row_count(ref))")
    println("  Columns: $(join(sort(collect(schema(ref).columns)), ", "))")
end

#==========================================================
Additional Helper Functions
==========================================================#

"""
    validate_schema(ref::FileReference, required_cols::Set{Symbol})

Validate that a file reference has all required columns.
"""
function validate_schema(ref::FileReference, required_cols::Set{Symbol})
    validate_required_columns(schema(ref), required_cols)
end

"""
    create_reference(file_path::String, ::Type{T}) where T <: FileReference

Create a file reference of the specified type.
"""
function create_reference(file_path::String, ::Type{PSMFileReference})
    return PSMFileReference(file_path)
end

function create_reference(file_path::String, ::Type{ProteinGroupFileReference})
    return ProteinGroupFileReference(file_path)
end

function create_reference(file_path::String, ::Type{ProteinQuantFileReference})
    return ProteinQuantFileReference(file_path)
end

"""
    create_protein_quant_reference(file_path::String) -> ProteinQuantFileReference

Create a protein quantification file reference, extracting metadata from the file.
"""
create_protein_quant_reference(file_path::String) = ProteinQuantFileReference(file_path)

# Export all public types and functions
export FileReference, FileSchema, PSMFileReference, ProteinGroupFileReference, ProteinQuantFileReference, PairedSearchFiles,
       file_path, schema, sorted_by, row_count, exists, n_protein_groups, n_experiments,
       has_column, get_column_or_default, validate_required_columns, validate_schema,
       mark_sorted!, is_sorted_by, ensure_sorted!, validate_exists,
       sort_file_by_keys!, create_psm_reference, create_protein_reference, create_protein_quant_reference, create_reference,
       describe_reference