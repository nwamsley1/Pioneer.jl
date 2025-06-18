using Test
using Arrow, DataFrames

# Include the necessary files
cd(@__DIR__)  # Change to test directory
package_root = dirname(dirname(@__DIR__))
include(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods", "FileReferences.jl"))

@testset "FileReferences Tests" begin
    
    @testset "FileSchema Tests" begin
        # Test schema creation
        schema = FileSchema([:col1, :col2, :col3])
        @test has_column(schema, :col1)
        @test has_column(schema, :col2)
        @test !has_column(schema, :col4)
        
        # Test immutability - columns are stored in a Set
        @test schema.columns == Set([:col1, :col2, :col3])
        
        # Test column validation
        @test_throws ErrorException validate_required_columns(schema, Set([:col1, :col4]))
        @test validate_required_columns(schema, Set([:col1, :col2])) === nothing
        
        # Test get_column_or_default
        df = DataFrame(col1 = [1, 2, 3], col2 = ["a", "b", "c"])
        @test get_column_or_default(df, schema, :col1, 0) == [1, 2, 3]
        @test get_column_or_default(df, schema, :col4, "missing") == ["missing", "missing", "missing"]
    end
    
    @testset "File Reference Creation" begin
        # Create a temporary test file
        temp_dir = mktempdir()
        test_file = joinpath(temp_dir, "test_psms.arrow")
        
        # Create test data
        test_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            prob = Float32[0.9, 0.8, 0.7],
            target = Bool[true, false, true],
            entrapment_group_id = UInt8[0, 0, 1]
        )
        Arrow.write(test_file, test_df)
        
        # Test PSMFileReference creation
        psm_ref = PSMFileReference(test_file)
        @test psm_ref.file_exists
        @test psm_ref.row_count == 3
        @test psm_ref.schema.columns == Set([:precursor_idx, :prob, :target, :entrapment_group_id])
        @test psm_ref.sorted_by == ()  # Not sorted initially
        
        # Test non-existent file
        fake_ref = PSMFileReference("fake_file.arrow")
        @test !fake_ref.file_exists
        @test fake_ref.row_count == 0
        @test isempty(fake_ref.schema.columns)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Sort State Management" begin
        # Create test file
        temp_dir = mktempdir()
        test_file = joinpath(temp_dir, "test_psms.arrow")
        
        test_df = DataFrame(
            precursor_idx = UInt32[3, 1, 2],
            prob = Float32[0.7, 0.9, 0.8],
            target = Bool[true, true, false]
        )
        Arrow.write(test_file, test_df)
        
        # Create reference
        psm_ref = PSMFileReference(test_file)
        
        # Test sort state tracking
        @test !is_sorted_by(psm_ref, :precursor_idx)
        mark_sorted!(psm_ref, :precursor_idx)
        @test is_sorted_by(psm_ref, :precursor_idx)
        @test !is_sorted_by(psm_ref, :prob)
        
        # Test multi-key sorting
        mark_sorted!(psm_ref, :target, :precursor_idx)
        @test is_sorted_by(psm_ref, :target, :precursor_idx)
        @test !is_sorted_by(psm_ref, :precursor_idx, :target)  # Order matters
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "File Sorting Operations" begin
        # Create test file
        temp_dir = mktempdir()
        test_file = joinpath(temp_dir, "test_psms.arrow")
        
        test_df = DataFrame(
            precursor_idx = UInt32[3, 1, 2, 4],
            prob = Float32[0.7, 0.9, 0.8, 0.6],
            target = Bool[true, true, false, true]
        )
        Arrow.write(test_file, test_df)
        
        # Create reference
        psm_ref = PSMFileReference(test_file)
        
        # Sort by single key
        sort_file_by_keys!(psm_ref, :precursor_idx)
        @test is_sorted_by(psm_ref, :precursor_idx)
        
        # Verify file is actually sorted
        sorted_df = DataFrame(Arrow.Table(test_file))
        @test sorted_df.precursor_idx == UInt32[1, 2, 3, 4]
        
        # Sort by multiple keys
        sort_file_by_keys!(psm_ref, :target, :prob)
        @test is_sorted_by(psm_ref, :target, :prob)
        
        # Test ensure_sorted! (should not re-sort if already sorted)
        ensure_sorted!(psm_ref, :target, :prob)
        @test is_sorted_by(psm_ref, :target, :prob)
        
        # Test sorting with non-existent column
        @test_throws ErrorException sort_file_by_keys!(psm_ref, :fake_column)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Paired Files" begin
        # Create test files
        temp_dir = mktempdir()
        psm_file = joinpath(temp_dir, "test_psms.arrow")
        protein_file = joinpath(temp_dir, "test_proteins.arrow")
        
        # Create PSM data
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            inferred_protein_group = ["P1", "P1", "P2"]
        )
        Arrow.write(psm_file, psm_df)
        
        # Create protein data
        protein_df = DataFrame(
            protein_name = ["P1", "P2"],
            n_peptides = Int64[2, 1]
        )
        Arrow.write(protein_file, protein_df)
        
        # Test paired file creation
        paired = PairedSearchFiles(psm_file, protein_file, 1)
        @test paired.ms_file_idx == 1
        @test paired.psm_ref.file_exists
        @test paired.protein_ref.file_exists
        
        # Test pairing validation - one file missing
        rm(protein_file)
        @test_throws ErrorException PairedSearchFiles(psm_file, protein_file, 1)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Utility Functions" begin
        # Create test file
        temp_dir = mktempdir()
        test_file = joinpath(temp_dir, "test_psms.arrow")
        
        test_df = DataFrame(
            col1 = [1, 2, 3],
            col2 = ["a", "b", "c"]
        )
        Arrow.write(test_file, test_df)
        
        # Test reference creation helpers
        psm_ref = create_psm_reference(test_file)
        @test psm_ref isa PSMFileReference
        @test psm_ref.file_exists
        
        # Test describe_reference (just check it doesn't error)
        @test describe_reference(psm_ref) === nothing
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end