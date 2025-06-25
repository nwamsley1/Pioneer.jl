
using Random

@testset "Sort State Management Tests" begin
    
    @testset "Basic Sort State Tracking" begin
        temp_dir = mktempdir()
        
        # Create test data
        test_df = DataFrame(
            precursor_idx = UInt32[5, 2, 1, 4, 3],
            prob = Float32[0.6, 0.8, 0.9, 0.7, 0.75],
            target = Bool[false, true, true, true, false]
        )
        test_file = joinpath(temp_dir, "sortable.arrow")
        Arrow.write(test_file, test_df)
        
        ref = PSMFileReference(test_file)
        
        # Test initial state
        @test !is_sorted_by(ref, :precursor_idx)
        @test !is_sorted_by(ref, :target, :prob)
        @test sorted_by(ref) == ()
        
        # Test single key sorting
        mark_sorted!(ref, :precursor_idx)
        @test is_sorted_by(ref, :precursor_idx)
        @test !is_sorted_by(ref, :prob)
        @test sorted_by(ref) == (:precursor_idx,)
        
        # Test multi-key sorting (order matters)
        mark_sorted!(ref, :target, :prob)
        @test is_sorted_by(ref, :target, :prob)
        @test !is_sorted_by(ref, :prob, :target)  # Different order
        @test !is_sorted_by(ref, :target)         # Subset doesn't match
        @test sorted_by(ref) == (:target, :prob)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Sort State Validation" begin
        temp_dir = mktempdir()
        
        test_df = DataFrame(
            a = [1, 2, 3],
            b = [3, 2, 1],
            c = [1, 1, 2]
        )
        test_file = joinpath(temp_dir, "validation.arrow")
        Arrow.write(test_file, test_df)
        
        ref = PSMFileReference(test_file)
        
        # Test ensure_sorted! functionality
        mark_sorted!(ref, :a)
        ensure_sorted!(ref, :a)  # Should not change state
        @test is_sorted_by(ref, :a)
        
        # Test ensure_sorted! with different column - should sort successfully
        mark_sorted!(ref, :b)  # Mark as sorted by :b
        ensure_sorted!(ref, :c)  # Should sort the file by :c
        @test is_sorted_by(ref, :c)  # File should now be sorted by :c
        
        # Test error condition with non-existent column
        @test_throws ErrorException ensure_sorted!(ref, :nonexistent_column)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Multi-Key Sort Operations" begin
        temp_dir = mktempdir()
        
        # Create larger test dataset
        n_rows = 50
        test_df = DataFrame(
            group = rand(["A", "B", "C"], n_rows),
            score = randn(Float32, n_rows),
            id = UInt32.(shuffle(1:n_rows))
        )
        test_file = joinpath(temp_dir, "multi_sort.arrow")
        Arrow.write(test_file, test_df)
        
        ref = PSMFileReference(test_file)
        
        # Test complex multi-key sort state
        @test !is_sorted_by(ref, :group, :score)
        @test !is_sorted_by(ref, :score, :group)
        
        # Mark as sorted by group then score
        mark_sorted!(ref, :group, :score)
        @test is_sorted_by(ref, :group, :score)
        @test !is_sorted_by(ref, :score, :group)  # Order matters
        @test !is_sorted_by(ref, :group)          # Partial match fails
        
        # Test three-key sorting
        mark_sorted!(ref, :group, :score, :id)
        @test is_sorted_by(ref, :group, :score, :id)
        @test !is_sorted_by(ref, :group, :score)  # Subset fails
        @test !is_sorted_by(ref, :group, :id, :score)  # Wrong order
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Sort State Reset" begin
        temp_dir = mktempdir()
        
        test_df = DataFrame(x = [1, 2, 3], y = [3, 2, 1])
        test_file = joinpath(temp_dir, "reset.arrow")
        Arrow.write(test_file, test_df)
        
        ref = PSMFileReference(test_file)
        
        # Set initial sort state
        mark_sorted!(ref, :x, :y)
        @test is_sorted_by(ref, :x, :y)
        
        # Reset by marking with different keys
        mark_sorted!(ref, :y)
        @test !is_sorted_by(ref, :x, :y)
        @test is_sorted_by(ref, :y)
        
        # Reset to unsorted
        mark_sorted!(ref)  # No arguments = reset to unsorted
        @test !is_sorted_by(ref, :y)
        @test sorted_by(ref) == ()
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end