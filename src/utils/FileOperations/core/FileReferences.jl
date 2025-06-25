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
        # Get row count from first column if columns exist
        col_names = Tables.columnnames(tbl)
        row_count = isempty(col_names) ? 0 : length(Tables.getcolumn(tbl, 1))
        
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
        # Get row count from first column if columns exist
        col_names = Tables.columnnames(tbl)
        row_count = isempty(col_names) ? 0 : length(Tables.getcolumn(tbl, 1))
        
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
        # Get row count from first column if columns exist
        col_names = Tables.columnnames(tbl)
        row_count = isempty(col_names) ? 0 : length(Tables.getcolumn(tbl, 1))
        
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

# Specialized accessors for ProteinQuantFileReference
n_protein_groups(ref::ProteinQuantFileReference) = ref.n_protein_groups
n_experiments(ref::ProteinQuantFileReference) = ref.n_experiments

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
    create_protein_quant_reference(file_path::String) -> ProteinQuantFileReference

Create a protein quantification file reference, extracting metadata from the file.
"""
create_protein_quant_reference(file_path::String) = ProteinQuantFileReference(file_path)

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

"""
    validate_schema(ref::FileReference, required_cols::Set{Symbol})

Validate that a file reference has all required columns.
"""
function validate_schema(ref::FileReference, required_cols::Set{Symbol})
    validate_required_columns(schema(ref), required_cols)
end

# Export all public types and functions
export FileReference, PSMFileReference, ProteinGroupFileReference, ProteinQuantFileReference, PairedSearchFiles,
       file_path, schema, sorted_by, row_count, exists, n_protein_groups, n_experiments,
       validate_exists, validate_schema,
       create_psm_reference, create_protein_reference, create_protein_quant_reference, create_reference,
       describe_reference