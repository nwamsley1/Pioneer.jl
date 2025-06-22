using Test
using Arrow, DataFrames

# Include the algorithms module
package_root = dirname(dirname(dirname(dirname(@__DIR__))))
include(joinpath(package_root, "src", "utils", "FileOperations", "FileOperations.jl"))

@testset "PSM Update Algorithm Tests" begin
    
    @testset "Basic PSM Score Updates" begin
        temp_dir = mktempdir()
        
        # Create PSM data
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3, 4],
            prob = Float32[0.8, 0.7, 0.9, 0.6],
            target = Bool[true, true, false, true],
            ms_file_idx = UInt16[1, 1, 2, 2]
        )
        
        # Create score updates
        score_updates_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3, 4],
            global_prob = Float32[0.85, 0.75, 0.88, 0.65],
            new_qval = Float32[0.01, 0.02, 0.05, 0.03]
        )
        
        psm_file = joinpath(temp_dir, "psms.arrow")
        scores_file = joinpath(temp_dir, "scores.arrow")
        Arrow.write(psm_file, psm_df)
        Arrow.write(scores_file, score_updates_df)
        
        psm_ref = PSMFileReference(psm_file)
        scores_ref = PSMFileReference(scores_file)
        
        # Apply score updates
        output_file = joinpath(temp_dir, "updated_psms.arrow")
        updated_ref = update_psms_with_scores(psm_ref, scores_ref, output_file)
        
        @test exists(updated_ref)
        result_df = DataFrame(Arrow.Table(file_path(updated_ref)))
        
        # Should have original columns plus new score columns
        @test hasproperty(result_df, :precursor_idx)
        @test hasproperty(result_df, :prob)
        @test hasproperty(result_df, :global_prob)
        @test hasproperty(result_df, :new_qval)
        @test nrow(result_df) == 4
        
        # Verify score updates were applied
        @test result_df.global_prob == Float32[0.85, 0.75, 0.88, 0.65]
        @test result_df.new_qval == Float32[0.01, 0.02, 0.05, 0.03]
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "PSM Updates with Partial Matches" begin
        temp_dir = mktempdir()
        
        # Create PSM data with more rows than score updates
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3, 4, 5],
            prob = Float32[0.8, 0.7, 0.9, 0.6, 0.85],
            target = Bool[true, true, false, true, true]
        )
        
        # Create score updates for only some PSMs
        score_updates_df = DataFrame(
            precursor_idx = UInt32[1, 3, 5],  # Missing 2 and 4
            global_prob = Float32[0.85, 0.88, 0.87],
            new_feature = ["A", "B", "C"]
        )
        
        psm_file = joinpath(temp_dir, "partial_psms.arrow")
        scores_file = joinpath(temp_dir, "partial_scores.arrow")
        Arrow.write(psm_file, psm_df)
        Arrow.write(scores_file, score_updates_df)
        
        psm_ref = PSMFileReference(psm_file)
        scores_ref = PSMFileReference(scores_file)
        
        # Apply partial updates
        output_file = joinpath(temp_dir, "partial_updated.arrow")
        updated_ref = update_psms_with_scores(psm_ref, scores_ref, output_file)
        
        result_df = DataFrame(Arrow.Table(file_path(updated_ref)))
        
        # Should have all original rows
        @test nrow(result_df) == 5
        
        # Updated rows should have new values, others should have missing
        @test result_df.global_prob[1] == 0.85  # Updated
        @test ismissing(result_df.global_prob[2])  # Not updated
        @test result_df.global_prob[3] == 0.88  # Updated
        @test ismissing(result_df.global_prob[4])  # Not updated
        @test result_df.global_prob[5] == 0.87  # Updated
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "PSM Updates with Conflicting Columns" begin
        temp_dir = mktempdir()
        
        # Create PSM data
        psm_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            prob = Float32[0.8, 0.7, 0.9],
            shared_col = ["original1", "original2", "original3"]
        )
        
        # Create score updates with overlapping column name
        score_updates_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            shared_col = ["updated1", "updated2", "updated3"],  # Same column name
            new_col = Float32[1.0, 2.0, 3.0]
        )
        
        psm_file = joinpath(temp_dir, "conflict_psms.arrow")
        scores_file = joinpath(temp_dir, "conflict_scores.arrow")
        Arrow.write(psm_file, psm_df)
        Arrow.write(scores_file, score_updates_df)
        
        psm_ref = PSMFileReference(psm_file)
        scores_ref = PSMFileReference(scores_file)
        
        # Apply updates - should prefer score file values for conflicts
        output_file = joinpath(temp_dir, "conflict_updated.arrow")
        updated_ref = update_psms_with_scores(psm_ref, scores_ref, output_file)
        
        result_df = DataFrame(Arrow.Table(file_path(updated_ref)))
        
        # Should have updated values from scores file
        @test result_df.shared_col == ["updated1", "updated2", "updated3"]
        @test result_df.new_col == Float32[1.0, 2.0, 3.0]
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "PSM Updates with Different Sort Orders" begin
        temp_dir = mktempdir()
        
        # Create PSM data in one order
        psm_df = DataFrame(
            precursor_idx = UInt32[3, 1, 4, 2],
            prob = Float32[0.9, 0.8, 0.6, 0.7],
            target = Bool[false, true, true, true]
        )
        
        # Create score updates in different order
        score_updates_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3, 4],  # Different order
            global_prob = Float32[0.85, 0.75, 0.88, 0.65]
        )
        
        psm_file = joinpath(temp_dir, "unsorted_psms.arrow")
        scores_file = joinpath(temp_dir, "sorted_scores.arrow")
        Arrow.write(psm_file, psm_df)
        Arrow.write(scores_file, score_updates_df)
        
        psm_ref = PSMFileReference(psm_file)
        scores_ref = PSMFileReference(scores_file)
        
        # Apply updates - should handle different orders
        output_file = joinpath(temp_dir, "reordered_updated.arrow")
        updated_ref = update_psms_with_scores(psm_ref, scores_ref, output_file)
        
        result_df = DataFrame(Arrow.Table(file_path(updated_ref)))
        
        # Should maintain original PSM order but with correct score updates
        @test result_df.precursor_idx == UInt32[3, 1, 4, 2]
        @test result_df.global_prob == Float32[0.88, 0.85, 0.65, 0.75]
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "PSM Updates Error Handling" begin
        temp_dir = mktempdir()
        
        # Create PSM data without precursor_idx
        bad_psm_df = DataFrame(
            some_id = [1, 2, 3],
            prob = Float32[0.8, 0.7, 0.9]
        )
        
        # Create valid score updates
        score_updates_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            global_prob = Float32[0.85, 0.75, 0.88]
        )
        
        bad_psm_file = joinpath(temp_dir, "bad_psms.arrow")
        scores_file = joinpath(temp_dir, "valid_scores.arrow")
        Arrow.write(bad_psm_file, bad_psm_df)
        Arrow.write(scores_file, score_updates_df)
        
        bad_psm_ref = PSMFileReference(bad_psm_file)
        scores_ref = PSMFileReference(scores_file)
        
        output_file = joinpath(temp_dir, "error_output.arrow")
        
        # Should error due to missing precursor_idx column
        @test_throws KeyError update_psms_with_scores(bad_psm_ref, scores_ref, output_file)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Large Scale PSM Updates" begin
        temp_dir = mktempdir()
        
        # Create larger datasets
        n_psms = 10000
        psm_df = DataFrame(
            precursor_idx = UInt32.(1:n_psms),
            prob = rand(Float32, n_psms),
            target = rand(Bool, n_psms),
            ms_file_idx = rand(UInt16(1):UInt16(10), n_psms)
        )
        
        # Create score updates for subset
        n_updates = 7500
        score_updates_df = DataFrame(
            precursor_idx = UInt32.(1:n_updates),
            global_prob = rand(Float32, n_updates),
            new_score = randn(Float32, n_updates)
        )
        
        psm_file = joinpath(temp_dir, "large_psms.arrow")
        scores_file = joinpath(temp_dir, "large_scores.arrow")
        Arrow.write(psm_file, psm_df)
        Arrow.write(scores_file, score_updates_df)
        
        psm_ref = PSMFileReference(psm_file)
        scores_ref = PSMFileReference(scores_file)
        
        # Test performance of large update
        output_file = joinpath(temp_dir, "large_updated.arrow")
        @time updated_ref = update_psms_with_scores(psm_ref, scores_ref, output_file)
        
        @test exists(updated_ref)
        result_df = DataFrame(Arrow.Table(file_path(updated_ref)))
        @test nrow(result_df) == n_psms
        
        # Check that first n_updates rows have scores, rest are missing
        @test !ismissing(result_df.global_prob[1])
        @test !ismissing(result_df.global_prob[n_updates])
        @test ismissing(result_df.global_prob[n_updates + 1])
        @test ismissing(result_df.global_prob[end])
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end