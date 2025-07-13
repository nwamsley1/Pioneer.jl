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
Comprehensive tests for filtering functions in PipelineOperations.jl

This test file covers:
- filter_by_threshold with all comparison operators
- filter_by_multiple_thresholds with AND/OR logic
- Type-stable helper functions for mask updates
- Edge cases and error handling
- Performance with large datasets
"""

@testset "Pipeline Filtering Operations Tests" begin
    
    @testset "filter_by_threshold Tests" begin
        temp_dir = mktempdir()
        
        # Create test data with various numeric values
        test_df = DataFrame(
            id = [1, 2, 3, 4, 5, 6],
            score = Float32[0.9, 0.7, 0.8, 0.6, 0.95, 0.75],
            qval = Float32[0.01, 0.03, 0.02, 0.06, 0.005, 0.04],
            count = Float32[10, 20, 15, 5, 25, 18]
        )
        test_file = joinpath(temp_dir, "threshold_test.arrow")
        Arrow.write(test_file, test_df)
        
        @testset "Comparison operator: $op" for (op, op_sym, expected_ids) in [
            (:(<=), :<=, [1, 2, 3, 5]),      # qval <= 0.03
            (:(<), :<, [1, 3, 5]),           # qval < 0.03
            (:(>=), :>=, [2, 4, 6]),      # qval >= 0.03
            (:(>), :>, [4, 6]),           # qval > 0.03
            (:(==), :(==), [2]),          # qval == 0.03
            (:(!=), :(!=), [1, 3, 4, 5, 6]) # qval != 0.03
        ]
            # Create operation
            filter_op = filter_by_threshold(:qval, 0.03f0; comparison=op)
            desc, op_func = filter_op
            @test contains(desc, "filter_by_threshold")
            @test contains(desc, string(op))
            
            # Apply to DataFrame
            df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
            result_df = op_func(df_copy)
            
            # Verify results
            @test sort(result_df.id) == expected_ids
        end
        
        @testset "Type handling" begin
            # Test with Float32 (standard case)
            filter_op = filter_by_threshold(:score, 0.8f0)
            desc, op_func = filter_op
            
            df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
            result_df = op_func(df_copy)
            @test all(result_df.score .<= 0.8f0)
            @test nrow(result_df) == 4  # scores: 0.7, 0.8, 0.6, 0.75
        end
        
        @testset "Edge cases" begin
            # Empty DataFrame
            empty_df = DataFrame(score = Float32[], qval = Float32[])
            empty_file = joinpath(temp_dir, "empty.arrow")
            Arrow.write(empty_file, empty_df)
            
            filter_op = filter_by_threshold(:score, 0.5f0)
            desc, op_func = filter_op
            
            df_copy = DataFrame(Tables.columntable(Arrow.Table(empty_file)))
            result_df = op_func(df_copy)
            @test nrow(result_df) == 0
            
            # Single row DataFrame
            single_df = DataFrame(score = Float32[0.9], qval = Float32[0.01])
            single_file = joinpath(temp_dir, "single.arrow")
            Arrow.write(single_file, single_df)
            
            df_copy = DataFrame(Tables.columntable(Arrow.Table(single_file)))
            result_df = op_func(df_copy)
            @test nrow(result_df) == 0  # 0.9 > 0.5, so filtered out with <=
            
            # All rows pass filter
            filter_op = filter_by_threshold(:score, 1.0f0)
            desc, op_func = filter_op
            df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
            result_df = op_func(df_copy)
            @test nrow(result_df) == nrow(test_df)
            
            # No rows pass filter
            filter_op = filter_by_threshold(:score, 0.0f0)
            desc, op_func = filter_op
            df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
            result_df = op_func(df_copy)
            @test nrow(result_df) == 0
        end
        
        @testset "Error handling" begin
            # Unsupported comparison operator
            @test_throws Exception filter_by_threshold(:score, 0.5f0; comparison=:unsupported)
        end
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "filter_by_multiple_thresholds Tests" begin
        temp_dir = mktempdir()
        
        # Create test data with mixed types
        test_df = DataFrame(
            id = [1, 2, 3, 4, 5, 6],
            score1 = Float32[0.9, 0.7, 0.8, 0.6, 0.95, 0.75],
            score2 = Float32[0.85, 0.65, 0.9, 0.7, 0.88, 0.8],
            qval = Float32[0.01, 0.03, 0.02, 0.06, 0.005, 0.04],
            category = ["protein", "peptide", "protein", "peptide", "protein", "peptide"],
            status = ["valid", "valid", "invalid", "valid", "valid", "invalid"]
        )
        test_file = joinpath(temp_dir, "multi_threshold_test.arrow")
        Arrow.write(test_file, test_df)
        
        @testset "AND logic" begin
            # Test numeric AND conditions
            filter_op = filter_by_multiple_thresholds([
                (:score1, 0.8f0),
                (:qval, 0.02f0)
            ]; logic=:and)
            desc, op_func = filter_op
            @test contains(desc, "filter_by_multiple_thresholds")
            @test contains(desc, "and")
            
            df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
            result_df = op_func(df_copy)
            @test all(result_df.score1 .<= 0.8f0)
            @test all(result_df.qval .<= 0.02f0)
            @test sort(result_df.id) == [3]  # Only rows 1 and 5 meet both conditions
            
            # Test with mixed numeric and string conditions
            filter_op = filter_by_multiple_thresholds([
                (:qval, 0.03f0),
                (:category, "protein")
            ]; comparison = :(==), logic=:and)
            
            df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
            result_df = op_func(df_copy)
            @test all(result_df.qval .<= 0.03f0)
            @test all(result_df.category .== "protein")
            @test nrow(result_df) == 1  
            
            # Test string equality (case-insensitive)
            filter_op = filter_by_multiple_thresholds([
                (:category, "PROTEIN"),
                (:status, "valid")
            ]; comparison = :(==), logic=:and)
            desc, op_func = filter_op
            
            df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
            result_df = op_func(df_copy)
            @test all(lowercase.(result_df.category) .== "protein")
            @test all(result_df.status .== "valid")
            @test sort(result_df.id) == [1, 5]
        end
        
        @testset "OR logic" begin
            # Test numeric OR conditions
            filter_op = filter_by_multiple_thresholds([
                (:score1, 0.95f0),
                (:score2, 0.9f0)
            ]; logic=:or)
            desc, op_func = filter_op
            @test contains(desc, "or")
            
            df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
            result_df = op_func(df_copy)
            # Should include rows where either score1 <= 0.95 OR score2 <= 0.9
            @test all(row -> row.score1 <= 0.95f0 || row.score2 <= 0.9f0, eachrow(result_df))
            
            # Test OR with no matches
            filter_op = filter_by_multiple_thresholds([
                (:score1, 0.5f0),
                (:score2, 0.6f0)
            ]; logic=:or)
            desc, op_func = filter_op
            
            df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
            result_df = op_func(df_copy)
            @test nrow(result_df) == 0  # No scores that low
        end
        
        @testset "Different comparison operators" begin
            # Test with >= operator
            filter_op = filter_by_multiple_thresholds([
                (:score1, 0.8f0),
                (:score2, 0.85f0)
            ]; comparison = :(>=), logic=:and)
            desc, op_func = filter_op
            
            df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
            result_df = op_func(df_copy)
            @test all(result_df.score1 .>= 0.8f0)
            @test all(result_df.score2 .>= 0.85f0)
            
            # Test with != operator for strings
            filter_op = filter_by_multiple_thresholds([
                (:category, "peptide"),
                (:status, "invalid")
            ]; comparison = :(!=), logic=:and)
            desc, op_func = filter_op
            
            df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
            result_df = op_func(df_copy)
            @test all(lowercase.(result_df.category) .!= "peptide")
            @test all(lowercase.(result_df.status) .!= "invalid")
        end
        
        @testset "Edge cases" begin
            # Single condition (should work like filter_by_threshold)
            filter_op = filter_by_multiple_thresholds([(:qval, 0.02f0)])
            desc, op_func = filter_op
            
            df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
            result_df = op_func(df_copy)
            @test all(result_df.qval .<= 0.02f0)
            
            # All conditions on same column
            filter_op = filter_by_multiple_thresholds([
                (:qval, 0.01f0),
                (:qval, 0.02f0)
            ]; comparison = :(>=), logic=:or)
            desc, op_func = filter_op
            
            df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
            result_df = op_func(df_copy)
            @test all(result_df.qval .>= 0.01f0)  # All values >= 0.01
        end
        
        @testset "Error handling" begin
            # Unsupported logic
            @test_throws Exception begin
                filter_op = filter_by_multiple_thresholds([(:qval, 0.01f0)]; logic=:xor)
                desc, op_func = filter_op
                df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
                op_func(df_copy)
            end
            
            # Unsupported comparison for strings
            @test_throws Exception begin
                filter_op = filter_by_multiple_thresholds([
                    (:category, "protein")
                ]; comparison = :(<))
                desc, op_func = filter_op
                df_copy = DataFrame(Tables.columntable(Arrow.Table(test_file)))
                op_func(df_copy)
            end
        end
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Helper Functions Tests" begin
        # Test update_mask_and! and update_mask_or! functions
        
        @testset "Numeric type helpers" begin
            # Test data
            col = Float32[0.1, 0.5, 0.8, 0.3, 0.9]
            threshold = 0.5f0
            
            # Test AND logic
            mask = trues(5)
            Pioneer.update_mask_and!(mask, col, threshold, <=)
            @test mask == [true, true, false, true, false]
            
            # Test OR logic
            mask = falses(5)
            Pioneer.update_mask_or!(mask, col, threshold, >)
            @test mask == [false, false, true, false, true]
        end
        
        @testset "Union{Missing,T} type helpers" begin
            # Test data with missing values
            col = Union{Missing,Float32}[0.1, missing, 0.8, 0.3, missing]
            threshold = 0.5f0
            
            # Test AND logic - missings should become false
            mask = trues(5)
            Pioneer.update_mask_and!(mask, col, threshold, <=)
            @test mask == [true, false, false, true, false]
            
            # Test OR logic - missings stay false
            mask = falses(5)
            Pioneer.update_mask_or!(mask, col, threshold, >)
            @test mask == [false, false, true, false, false]
        end
        
        @testset "String type helpers" begin
            # Test data
            col = ["Apple", "banana", "APPLE", "orange", "Banana"]
            
            # Test AND logic with == (case-insensitive)
            mask = trues(5)
            Pioneer.update_mask_and!(mask, col, "apple", ==)
            @test mask == [true, false, true, false, false]
            
            # Test OR logic with != (case-insensitive)
            mask = falses(5)
            Pioneer.update_mask_or!(mask, col, "BANANA", !=)
            @test mask == [true, false, true, true, false]
            
            # Test error for unsupported string comparison
            mask = trues(5)
            @test_throws Exception Pioneer.update_mask_and!(mask, col, "apple", <)
        end
    end
    
    @testset "Integration with Pipeline API" begin
        temp_dir = mktempdir()
        
        # Create complex test data
        test_df = DataFrame(
            id = 1:100,
            score = Float32.(0.5 .+ 0.5 .* rand(100)),
            qval = Float32.(0.001 .+ 0.099 .* rand(100)),
            category = rand(["A", "B", "C"], 100),
            flag = rand(Bool, 100)
        )
        test_file = joinpath(temp_dir, "integration_test.arrow")
        Arrow.write(test_file, test_df)
        test_ref = PSMFileReference(test_file)
        
        # Build complex pipeline with multiple filters
        pipeline = TransformPipeline() |>
            filter_by_threshold(:score, 0.7f0; comparison = :(>=)) |>
            filter_by_multiple_thresholds([
                (:qval, 0.05f0),
                (:flag, true)
            ]; comparison = :(==), logic = :and) |>
            add_column(:score_rank, row -> row.score > 0.9 ? "high" : "medium") |>
            filter_by_multiple_thresholds([
                (:category, "A"),
                (:score_rank, "high")
            ]; comparison = :(==), logic = :or) |>
            sort_by([:score]; rev=[true])
        
        # Apply pipeline
        apply_pipeline!(test_ref, pipeline)
        
        result_df = DataFrame(Tables.columntable(Arrow.Table(test_file)))
        
        # Verify all filters were applied correctly
        @test all(result_df.score .>= 0.7f0)
        @test all(result_df.qval .<= 0.05f0)
        @test all(result_df.flag)
        @test all(row -> lowercase(row.category) == "a" || row.score_rank == "high", eachrow(result_df))
        @test issorted(result_df.score, rev=true)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Performance with Large Datasets" begin
        temp_dir = mktempdir()
        
        # Create large dataset
        n_rows = 100_000
        large_df = DataFrame(
            id = 1:n_rows,
            score1 = Float32.(rand(n_rows)),
            score2 = Float32.(rand(n_rows)), 
            score3 = Float32.(rand(n_rows)),
            category = rand(["cat1", "cat2", "cat3", "cat4"], n_rows),
            status = rand(["active", "inactive"], n_rows)
        )
        large_file = joinpath(temp_dir, "large_dataset.arrow")
        Arrow.write(large_file, large_df)
        
        # Test filter_by_threshold performance
        filter_op = filter_by_threshold(:score1, 0.5f0)
        desc, op_func = filter_op
        
        df_copy = DataFrame(Tables.columntable(Arrow.Table(large_file)))
        elapsed = @elapsed result_df = op_func(df_copy)
        @test elapsed < 1.0  # Should complete in under 1 second
        @test all(result_df.score1 .<= 0.5f0)
        
        # Test filter_by_multiple_thresholds performance
        filter_op = filter_by_multiple_thresholds([
            (:score1, 0.7f0),
            (:score2, 0.6f0),
            (:score3, 0.8f0),
            (:category, "cat1"),
            (:status, "active")
        ]; comparison = :(==), logic = :and)
        desc, op_func = filter_op
        
        df_copy = DataFrame(Tables.columntable(Arrow.Table(large_file)))
        elapsed = @elapsed result_df = op_func(df_copy)
        @test elapsed < 2.0  # Should complete in under 2 seconds for complex filter
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end