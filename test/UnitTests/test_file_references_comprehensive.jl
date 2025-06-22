using Test
using Arrow, DataFrames, Tables
using Random: shuffle

# Include the necessary files
cd(@__DIR__)  # Change to test directory
package_root = dirname(dirname(@__DIR__))
include(joinpath(package_root, "src", "utils", "FileOperations", "FileOperations.jl"))

@testset "FileReferences Comprehensive Tests" begin
    
    @testset "FileSchema Advanced Tests" begin
        # Test schema creation with various column types
        complex_schema = FileSchema([:int_col, :float_col, :string_col, :bool_col, :missing_col])
        
        # Test immutability and set operations
        @test complex_schema.columns == Set([:int_col, :float_col, :string_col, :bool_col, :missing_col])
        @test has_column(complex_schema, :int_col)
        @test !has_column(complex_schema, :nonexistent_col)
        
        # Test column validation with complex requirements
        required_subset = Set([:int_col, :float_col])
        @test validate_required_columns(complex_schema, required_subset) === nothing
        
        # Test validation failure with detailed error
        missing_subset = Set([:int_col, :fake_col1, :fake_col2])
        @test_throws ErrorException validate_required_columns(complex_schema, missing_subset)
        
        # Test get_column_or_default with various DataFrame types
        test_df = DataFrame(
            int_col = [1, 2, 3],
            float_col = [1.5, 2.5, 3.5],
            string_col = ["a", "b", "c"]
        )
        # Create schema that matches the actual DataFrame columns
        test_schema = FileSchema([:int_col, :float_col, :string_col])
        
        # Existing column in both schema and DataFrame
        @test get_column_or_default(test_df, test_schema, :int_col, 999) == [1, 2, 3]
        # Column not in schema - should return default
        limited_schema = FileSchema([:int_col, :float_col])
        @test get_column_or_default(test_df, limited_schema, :string_col, "missing") == ["missing", "missing", "missing"]
        # Column not in schema with different default types
        @test get_column_or_default(test_df, limited_schema, :nonexistent_col, true) == [true, true, true]
    end
    
    @testset "PSMFileReference Advanced Tests" begin
        temp_dir = mktempdir()
        
        # Create complex PSM test data with realistic schema
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3, 4, 5],
            prob = Float32[0.9, 0.85, 0.8, 0.75, 0.7],
            inferred_protein_group = ["P1", "P1", "P2", "P2", "P3"],
            target = Bool[true, true, true, false, true],
            entrapment_group_id = UInt8[0, 0, 0, 0, 1],
            pg_qval = Float32[0.01, 0.02, 0.03, 0.04, 0.05],
            global_qval_pg = Float32[0.01, 0.02, 0.03, 0.04, 0.05],
            use_for_protein_quant = Bool[true, true, true, true, true],
            peak_area = Float32[1000.0, 2000.0, 1500.0, 800.0, 1200.0],
            ms_file_idx = UInt16[1, 1, 2, 2, 3]
        )
        
        psm_file = joinpath(temp_dir, "complex_psms.arrow")
        Arrow.write(psm_file, psm_df)
        
        # Test PSMFileReference creation and metadata
        psm_ref = PSMFileReference(psm_file)
        @test psm_ref.file_exists
        @test psm_ref.row_count == 5
        @test psm_ref.sorted_by == ()  # Not sorted initially
        
        # Test schema contains all expected columns
        expected_cols = Set([:precursor_idx, :prob, :inferred_protein_group, :target, 
                           :entrapment_group_id, :pg_qval, :global_qval_pg, 
                           :use_for_protein_quant, :peak_area, :ms_file_idx])
        @test psm_ref.schema.columns == expected_cols
        
        # Test accessor functions
        @test file_path(psm_ref) == psm_file
        @test exists(psm_ref) == true
        @test row_count(psm_ref) == 5
        @test sorted_by(psm_ref) == ()
        
        # Test with empty file
        empty_file = joinpath(temp_dir, "empty_psms.arrow")
        empty_df = DataFrame()
        Arrow.write(empty_file, empty_df)
        empty_ref = PSMFileReference(empty_file)
        @test empty_ref.file_exists
        @test empty_ref.row_count == 0
        @test isempty(empty_ref.schema.columns)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "ProteinQuantFileReference Specialized Tests" begin
        temp_dir = mktempdir()
        
        # Create protein quantification test data
        protein_quant_df = DataFrame(
            file_name = ["file1.raw", "file2.raw", "file1.raw", "file2.raw", "file3.raw"],
            protein = ["P1", "P1", "P2", "P2", "P3"],
            target = Bool[true, true, true, false, true],
            entrap_id = UInt8[0, 0, 0, 0, 1],
            species = ["human", "human", "human", "human", "mouse"],
            abundance = Float32[5000.0, 4500.0, 3000.0, 2800.0, 2000.0],
            n_peptides = Int64[3, 3, 2, 2, 1],
            experiments = UInt32[1, 2, 1, 2, 3]
        )
        
        quant_file = joinpath(temp_dir, "protein_quant.arrow")
        Arrow.write(quant_file, protein_quant_df)
        
        # Test ProteinQuantFileReference creation
        quant_ref = ProteinQuantFileReference(quant_file)
        @test quant_ref.file_exists
        @test quant_ref.row_count == 5
        
        # Test specialized accessors
        @test n_protein_groups(quant_ref) == 3  # P1, P2, P3
        @test n_experiments(quant_ref) == 3     # experiments 1, 2, 3
        
        # Test with file_name column instead of experiments
        protein_quant_df2 = select(protein_quant_df, Not(:experiments))
        quant_file2 = joinpath(temp_dir, "protein_quant2.arrow")
        Arrow.write(quant_file2, protein_quant_df2)
        
        quant_ref2 = ProteinQuantFileReference(quant_file2)
        @test n_experiments(quant_ref2) == 3  # Should count unique file_name values
        
        # Test with missing protein column
        no_protein_df = select(protein_quant_df, Not(:protein))
        no_protein_file = joinpath(temp_dir, "no_protein.arrow")
        Arrow.write(no_protein_file, no_protein_df)
        
        no_protein_ref = ProteinQuantFileReference(no_protein_file)
        @test n_protein_groups(no_protein_ref) == 0
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Sort State Management Advanced Tests" begin
        temp_dir = mktempdir()
        
        # Create test data with multiple sortable columns
        test_df = DataFrame(
            precursor_idx = UInt32[5, 2, 1, 4, 3],
            prob = Float32[0.6, 0.8, 0.9, 0.7, 0.75],
            target = Bool[false, true, true, true, false],
            entrap_id = UInt8[1, 0, 0, 0, 1]
        )
        test_file = joinpath(temp_dir, "sortable.arrow")
        Arrow.write(test_file, test_df)
        
        ref = PSMFileReference(test_file)
        
        # Test complex sort state tracking
        @test !is_sorted_by(ref, :precursor_idx)
        @test !is_sorted_by(ref, :target, :prob)
        
        # Test single key sorting
        mark_sorted!(ref, :precursor_idx)
        @test is_sorted_by(ref, :precursor_idx)
        @test !is_sorted_by(ref, :prob)
        
        # Test multi-key sorting (order matters)
        mark_sorted!(ref, :target, :prob)
        @test is_sorted_by(ref, :target, :prob)
        @test !is_sorted_by(ref, :prob, :target)  # Different order
        @test !is_sorted_by(ref, :target)         # Subset doesn't match
        
        # Test ensure_sorted! functionality
        mark_sorted!(ref, :entrap_id)
        ensure_sorted!(ref, :entrap_id)  # Should not change state
        @test is_sorted_by(ref, :entrap_id)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "File Sorting Operations Extended Tests" begin
        temp_dir = mktempdir()
        
        # Create larger test dataset for realistic sorting
        n_rows = 100
        test_df = DataFrame(
            precursor_idx = UInt32.(shuffle(1:n_rows)),
            prob = rand(Float32, n_rows),
            target = rand(Bool, n_rows),
            score = randn(Float32, n_rows)
        )
        test_file = joinpath(temp_dir, "large_sortable.arrow")
        Arrow.write(test_file, test_df)
        
        ref = PSMFileReference(test_file)
        
        # Test single key sorting with verification
        sort_file_by_keys!(ref, :precursor_idx)
        @test is_sorted_by(ref, :precursor_idx)
        
        # Verify actual file sorting
        sorted_df = DataFrame(Arrow.Table(test_file))
        @test issorted(sorted_df.precursor_idx)
        
        # Test multi-key sorting
        sort_file_by_keys!(ref, :target, :prob)
        @test is_sorted_by(ref, :target, :prob)
        
        # Verify multi-key sort
        sorted_df = DataFrame(Arrow.Table(test_file))
        @test issorted(sorted_df, [:target, :prob])
        
        # Test error handling for invalid columns
        @test_throws ErrorException sort_file_by_keys!(ref, :nonexistent_column)
        
        # Test Base.sort! override (single direction)
        sort!(ref, [:score])
        @test is_sorted_by(ref, :score)
        
        # Test mixed sort directions (should error)
        @test_throws ErrorException sort!(ref, [:target, :prob], rev=[true, false])
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "PairedSearchFiles Advanced Tests" begin
        temp_dir = mktempdir()
        
        # Create matching PSM and protein files
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            inferred_protein_group = ["P1", "P1", "P2"],
            prob = Float32[0.9, 0.8, 0.7]
        )
        protein_df = DataFrame(
            protein_name = ["P1", "P2"],
            n_peptides = Int64[2, 1],
            score = Float32[1.8, 0.7]
        )
        
        psm_file = joinpath(temp_dir, "paired_psms.arrow")
        protein_file = joinpath(temp_dir, "paired_proteins.arrow")
        Arrow.write(psm_file, psm_df)
        Arrow.write(protein_file, protein_df)
        
        # Test successful pairing
        paired = PairedSearchFiles(psm_file, protein_file, 42)
        @test paired.ms_file_idx == 42
        @test paired.psm_ref.file_exists
        @test paired.protein_ref.file_exists
        @test paired.psm_ref.row_count == 3
        @test paired.protein_ref.row_count == 2
        
        # Test pairing with reference objects
        psm_ref = PSMFileReference(psm_file)
        protein_ref = ProteinGroupFileReference(protein_file)
        paired2 = PairedSearchFiles(psm_ref, protein_ref, 24)
        @test paired2.ms_file_idx == 24
        
        # Test validation failure - remove one file
        rm(protein_file)
        @test_throws ErrorException PairedSearchFiles(psm_file, protein_file, 1)
        
        # Test validation failure with existing references
        missing_protein_ref = ProteinGroupFileReference(protein_file)  # Points to deleted file
        @test_throws ErrorException PairedSearchFiles(psm_ref, missing_protein_ref, 1)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Utility Functions Comprehensive Tests" begin
        temp_dir = mktempdir()
        
        # Create various test files
        files_data = [
            ("psm_test.arrow", DataFrame(col1 = [1, 2], col2 = ["a", "b"])),
            ("protein_test.arrow", DataFrame(protein = ["P1"], score = [1.5])),
            ("quant_test.arrow", DataFrame(protein = ["P1", "P2"], 
                                          experiments = UInt32[1, 1], 
                                          abundance = [100.0, 200.0]))
        ]
        
        file_paths = String[]
        for (filename, df) in files_data
            filepath = joinpath(temp_dir, filename)
            Arrow.write(filepath, df)
            push!(file_paths, filepath)
        end
        
        # Test factory methods
        psm_ref = create_psm_reference(file_paths[1])
        @test psm_ref isa PSMFileReference
        @test psm_ref.file_exists
        
        protein_ref = create_protein_reference(file_paths[2])
        @test protein_ref isa ProteinGroupFileReference
        @test protein_ref.file_exists
        
        quant_ref = create_protein_quant_reference(file_paths[3])
        @test quant_ref isa ProteinQuantFileReference
        @test quant_ref.file_exists
        
        # Test generic create_reference with type dispatch
        psm_ref2 = create_reference(file_paths[1], PSMFileReference)
        @test psm_ref2 isa PSMFileReference
        
        protein_ref2 = create_reference(file_paths[2], ProteinGroupFileReference)
        @test protein_ref2 isa ProteinGroupFileReference
        
        quant_ref2 = create_reference(file_paths[3], ProteinQuantFileReference)
        @test quant_ref2 isa ProteinQuantFileReference
        
        # Test describe_reference (should not error and return nothing)
        @test describe_reference(psm_ref) === nothing
        @test describe_reference(protein_ref) === nothing
        @test describe_reference(quant_ref) === nothing
        
        # Test validate_exists
        @test validate_exists(psm_ref) == true
        
        # Test with non-existent file
        fake_ref = PSMFileReference("nonexistent.arrow")
        @test_throws ErrorException validate_exists(fake_ref)
        
        # Test validate_schema
        required_cols = Set([:col1, :col2])
        @test validate_schema(psm_ref, required_cols) === nothing
        
        invalid_cols = Set([:col1, :missing_col])
        @test_throws ErrorException validate_schema(psm_ref, invalid_cols)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Edge Cases and Error Handling" begin
        temp_dir = mktempdir()
        
        # Test with empty DataFrame
        empty_df = DataFrame()
        empty_file = joinpath(temp_dir, "empty.arrow")
        Arrow.write(empty_file, empty_df)
        
        empty_ref = PSMFileReference(empty_file)
        @test empty_ref.file_exists
        @test empty_ref.row_count == 0
        @test isempty(empty_ref.schema.columns)
        
        # Test with single-row DataFrame
        single_df = DataFrame(only_col = [42])
        single_file = joinpath(temp_dir, "single.arrow")
        Arrow.write(single_file, single_df)
        
        single_ref = PSMFileReference(single_file)
        @test single_ref.row_count == 1
        @test has_column(single_ref.schema, :only_col)
        
        # Test with DataFrame containing missing values
        missing_df = DataFrame(
            id = [1, 2, 3],
            value = [1.0, missing, 3.0],
            text = ["a", "b", missing]
        )
        missing_file = joinpath(temp_dir, "missing.arrow")
        Arrow.write(missing_file, missing_df)
        
        missing_ref = PSMFileReference(missing_file)
        @test missing_ref.row_count == 3
        @test has_column(missing_ref.schema, :value)
        @test has_column(missing_ref.schema, :text)
        
        # Test ProteinQuantFileReference with all missing proteins
        all_missing_df = DataFrame(
            protein = [missing, missing],
            experiments = UInt32[1, 2]
        )
        all_missing_file = joinpath(temp_dir, "all_missing.arrow")
        Arrow.write(all_missing_file, all_missing_df)
        
        all_missing_ref = ProteinQuantFileReference(all_missing_file)
        @test n_protein_groups(all_missing_ref) == 0  # Should handle missing gracefully
        @test n_experiments(all_missing_ref) == 2     # Should still count experiments
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end