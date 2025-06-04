"""
Test and debug protein group ML scoring
"""

using DataFrames
using Arrow
using Statistics

"""
    test_protein_scoring_pipeline(
        psm_path::String,
        protein_group_path::String,
        precursors::LibraryPrecursors
    )

Test the protein scoring pipeline with detailed debugging output.
"""
function test_protein_scoring_pipeline(
    psm_path::String,
    protein_group_path::String,
    precursors::LibraryPrecursors;
    n_top_precursors::Int = 5
)
    @info "="^60
    @info "Testing Protein Scoring Pipeline"
    @info "="^60
    
    # Step 1: Load and examine PSMs
    @info "\nStep 1: Loading PSMs from $psm_path"
    psms = DataFrame(Arrow.Table(psm_path))
    @info "PSM shape: $(size(psms))"
    @info "PSM columns: $(names(psms))"
    
    # Check for required columns
    required_psm_cols = [:precursor_idx, :prob, :target, :entrapment_group_id]
    missing_psm_cols = setdiff(required_psm_cols, Symbol.(names(psms)))
    if !isempty(missing_psm_cols)
        @error "Missing PSM columns: $missing_psm_cols"
        return
    end
    
    # Step 2: Add necessary columns from precursors
    @info "\nStep 2: Adding precursor information"
    psms[!, :sequence] = [getSequence(precursors)[pid] for pid in psms.precursor_idx]
    psms[!, :charge] = [getCharge(precursors)[pid] for pid in psms.precursor_idx]
    psms[!, :proteins] = [getAccessionNumbers(precursors)[pid] for pid in psms.precursor_idx]
    psms[!, :cv_fold] = [getCvFold(precursors, pid) for pid in psms.precursor_idx]
    
    # Step 3: Examine protein groups
    @info "\nStep 3: Loading protein groups from $protein_group_path"
    pg_table = Arrow.Table(protein_group_path)
    pg_df = DataFrame(pg_table)
    @info "Protein group shape: $(size(pg_df))"
    @info "Protein group columns: $(names(pg_df))"
    
    # Step 4: Build protein groups dictionary
    @info "\nStep 4: Building protein groups dictionary"
    protein_groups = Dictionary{
        @NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8}, 
        @NamedTuple{pg_score::Float32, peptides::Set{String}}
    }()
    
    # Get peptides for each protein
    protein_peptides = Dict{String, Set{String}}()
    for (idx, row) in enumerate(eachrow(psms))
        protein = row.proteins
        peptide = row.sequence
        
        if !haskey(protein_peptides, protein)
            protein_peptides[protein] = Set{String}()
        end
        push!(protein_peptides[protein], peptide)
    end
    
    # Build protein groups from pg_df
    unique_proteins_found = 0
    for row in eachrow(pg_df)
        protein_key = (
            protein_name = row.protein_name,
            target = row.target,
            entrap_id = row.entrap_id
        )
        
        peptides = get(protein_peptides, row.protein_name, Set{String}())
        if !isempty(peptides)
            protein_groups[protein_key] = (
                pg_score = row.pg_score,
                peptides = peptides
            )
            unique_proteins_found += 1
        end
    end
    
    @info "Built protein groups for $unique_proteins_found proteins"
    @info "Total peptides in protein_peptides: $(sum(length(v) for v in values(protein_peptides)))"
    
    # Step 5: Test feature extraction
    @info "\nStep 5: Testing feature extraction"
    try
        protein_features = prepare_protein_group_features(
            psms, protein_groups, precursors, n_top_precursors
        )
        
        @info "Successfully extracted features for $(nrow(protein_features)) protein groups"
        @info "Feature columns: $(names(protein_features))"
        
        # Show sample of features
        if nrow(protein_features) > 0
            @info "\nSample protein features (first 3):"
            for i in 1:min(3, nrow(protein_features))
                row = protein_features[i, :]
                @info "  Protein $(i): $(row.protein_group.protein_name)"
                @info "    - is_target: $(row.is_target)"
                @info "    - n_peptides: $(row.n_peptides)"
                @info "    - n_precursors: $(row.n_precursors)"
                @info "    - current_score: $(row.current_score)"
                @info "    - cv_fold: $(row.cv_fold)"
                @info "    - top_scores: $([row[Symbol("top_score_", j)] for j in 1:n_top_precursors])"
            end
        end
        
        # Check target/decoy balance
        target_dist = countmap(protein_features.is_target)
        @info "\nTarget/Decoy distribution: $target_dist"
        
        # Check CV fold distribution
        fold_dist = countmap(protein_features.cv_fold)
        @info "CV fold distribution: $fold_dist"
        
        return protein_features
        
    catch e
        @error "Failed to extract protein features: $e"
        @error "Stacktrace:" exception=(e, catch_backtrace())
    end
end

"""
    diagnose_protein_scoring_issue(
        passing_psms_paths::Vector{String},
        passing_pg_paths::Vector{String},
        precursors::LibraryPrecursors
    )

Diagnose issues with protein scoring across multiple files.
"""
function diagnose_protein_scoring_issue(
    passing_psms_paths::Vector{String},
    passing_pg_paths::Vector{String},
    precursors::LibraryPrecursors
)
    @info "="^60
    @info "Diagnosing Protein Scoring Issues"
    @info "="^60
    
    @info "Number of PSM files: $(length(passing_psms_paths))"
    @info "Number of PG files: $(length(passing_pg_paths))"
    
    # Check file existence
    missing_psm_files = [p for p in passing_psms_paths if !isfile(p)]
    missing_pg_files = [p for p in passing_pg_paths if !isfile(p)]
    
    if !isempty(missing_psm_files)
        @error "Missing PSM files: $missing_psm_files"
    end
    
    if !isempty(missing_pg_files)
        @error "Missing PG files: $missing_pg_files"
    end
    
    # Test first file pair
    if length(passing_psms_paths) > 0 && length(passing_pg_paths) > 0
        first_psm = passing_psms_paths[1]
        first_pg = passing_pg_paths[1]
        
        if isfile(first_psm) && isfile(first_pg)
            @info "\nTesting first file pair:"
            @info "  PSM: $(basename(first_psm))"
            @info "  PG: $(basename(first_pg))"
            
            test_protein_scoring_pipeline(first_psm, first_pg, precursors)
        end
    end
    
    # Check overall statistics
    total_psms = 0
    total_pgs = 0
    
    for psm_path in passing_psms_paths
        if isfile(psm_path)
            psm_table = Arrow.Table(psm_path)
            total_psms += length(psm_table[1])
        end
    end
    
    for pg_path in passing_pg_paths
        if isfile(pg_path)
            pg_table = Arrow.Table(pg_path)
            total_pgs += length(pg_table[1])
        end
    end
    
    @info "\nOverall statistics:"
    @info "  Total PSMs: $total_psms"
    @info "  Total protein groups: $total_pgs"
end