using Test
using Arrow, DataFrames, Tables
using DataStructures  # For heap operations

# Include the necessary files
cd(@__DIR__)  # Change to test directory
package_root = dirname(dirname(@__DIR__))
include(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods", "FileReferences.jl"))
include(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods", "FileOperations.jl"))

@testset "FileOperations Comprehensive Tests" begin
    
    @testset "TransformPipeline Core Functionality" begin
        # Test pipeline creation and composition
        pipeline = TransformPipeline()
        @test length(pipeline.operations) == 0
        
        # Test operation addition with |> operator
        pipeline = TransformPipeline() |>
            add_column(:computed, row -> row.a + row.b) |>
            filter_rows(row -> row.computed > 5)
        
        @test length(pipeline.operations) == 2
        # Check that operation descriptions exist
        descriptions = [first(op) for op in pipeline.operations]
        @test "add_column_computed" in descriptions
        @test "filter" in descriptions
        
        # Test pipeline execution on simple DataFrame
        test_df = DataFrame(a = [1, 2, 3, 4], b = [2, 3, 4, 5])
        
        # Apply pipeline operations sequentially
        result_df = test_df
        for op_pair in pipeline.operations
            desc, op = op_pair.first, op_pair.second
            result_df = op(result_df)
        end
        
        # Should have computed column and filtered rows where computed > 5
        @test hasproperty(result_df, :computed)
        @test all(result_df.computed .> 5)
        @test nrow(result_df) == 2  # Only rows with a+b > 5
    end
    
    @testset "Column Operations" begin
        temp_dir = mktempdir()
        
        # Create test data
        test_df = DataFrame(
            id = [1, 2, 3, 4],
            value_a = [10, 20, 30, 40],
            value_b = [1.5, 2.5, 3.5, 4.5],
            category = ["A", "B", "A", "B"]
        )
        test_file = joinpath(temp_dir, "column_ops.arrow")
        Arrow.write(test_file, test_df)
        test_ref = PSMFileReference(test_file)
        
        # Test add_column operation
        add_op = add_column(:sum_values, row -> row.value_a + row.value_b)
        desc, op_func = add_op
        @test startswith(desc, "add_column")
        
        result_df = op_func(DataFrame(Arrow.Table(test_file)))
        @test hasproperty(result_df, :sum_values)
        @test result_df.sum_values ≈ [11.5, 22.5, 33.5, 44.5]
        
        # Test rename_column operation
        rename_op = rename_column(:value_a, :renamed_value)
        desc, op_func = rename_op
        @test startswith(desc, "rename")
        
        result_df = op_func(DataFrame(Arrow.Table(test_file)))
        @test hasproperty(result_df, :renamed_value)
        @test !hasproperty(result_df, :value_a)
        
        # Test select_columns operation
        select_op = select_columns([:id, :value_a])
        desc, op_func = select_op
        @test startswith(desc, "select")
        
        result_df = op_func(DataFrame(Arrow.Table(test_file)))
        @test ncol(result_df) == 2
        @test hasproperty(result_df, :id)
        @test hasproperty(result_df, :value_a)
        @test !hasproperty(result_df, :value_b)
        
        # Test remove_columns operation
        remove_op = remove_columns(:value_b, :category)
        desc, op_func = remove_op
        @test startswith(desc, "remove")
        
        result_df = op_func(DataFrame(Arrow.Table(test_file)))
        @test !hasproperty(result_df, :value_b)
        @test !hasproperty(result_df, :category)
        @test hasproperty(result_df, :id)
        @test hasproperty(result_df, :value_a)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Row Operations" begin
        temp_dir = mktempdir()
        
        # Create test data with filtering scenarios
        test_df = DataFrame(
            id = [1, 2, 3, 4, 5],
            score = [0.9, 0.7, 0.8, 0.6, 0.95],
            qval = [0.01, 0.03, 0.02, 0.06, 0.005],
            target = [true, true, false, true, true],
            use_for_quant = [true, true, true, false, true]
        )
        test_file = joinpath(temp_dir, "row_ops.arrow")
        Arrow.write(test_file, test_df)
        
        # Test basic filter_rows operation
        filter_op = filter_rows(row -> row.score > 0.75; desc="high_score_filter")
        desc, op_func = filter_op
        @test desc == "high_score_filter"
        
        result_df = op_func(DataFrame(Arrow.Table(test_file)))
        @test nrow(result_df) == 3  # Scores: 0.9, 0.8, 0.95
        @test all(result_df.score .> 0.75)
        
        # Test filter_by_multiple_thresholds
        multi_filter_op = filter_by_multiple_thresholds([
            (:score, 0.75),
            (:qval, 0.03)
        ])
        desc, op_func = multi_filter_op
        @test startswith(desc, "filter")
        
        result_df = op_func(DataFrame(Arrow.Table(test_file)))
        @test all(result_df.score .>= 0.75)
        @test all(result_df.qval .<= 0.03)
        @test nrow(result_df) == 2  # Only rows 1 and 3 meet both criteria
        
        # Test complex filter with multiple conditions
        complex_filter_op = filter_rows(row -> row.target && row.use_for_quant && row.qval < 0.04)
        desc, op_func = complex_filter_op
        
        result_df = op_func(DataFrame(Arrow.Table(test_file)))
        @test all(result_df.target)
        @test all(result_df.use_for_quant)
        @test all(result_df.qval .< 0.04)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Sort Operations with State Tracking" begin
        temp_dir = mktempdir()
        
        # Create unsorted test data
        test_df = DataFrame(
            primary = [3, 1, 4, 2],
            secondary = [0.8, 0.9, 0.6, 0.7],
            tertiary = ["D", "A", "C", "B"]
        )
        test_file = joinpath(temp_dir, "sort_ops.arrow")
        Arrow.write(test_file, test_df)
        test_ref = PSMFileReference(test_file)
        
        # Test single column sort
        sort_op = sort_by([:primary])
        desc, op_func = sort_op
        @test startswith(desc, "sort")
        
        result_df = op_func(DataFrame(Arrow.Table(test_file)))
        @test issorted(result_df.primary)
        
        # Test multi-column sort
        multi_sort_op = sort_by([:secondary, :primary]; rev=true)
        desc, op_func = multi_sort_op
        
        result_df = op_func(DataFrame(Arrow.Table(test_file)))
        @test issorted(result_df, [:secondary, :primary], rev=true)
        
        # Test sort with mixed directions (should work in DataFrame)
        mixed_sort_op = sort_by([:secondary, :primary]; rev=[false, true])
        desc, op_func = mixed_sort_op
        
        result_df = op_func(DataFrame(Arrow.Table(test_file)))
        @test issorted(result_df, [:secondary, :primary], rev=[false, true])
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Pipeline Execution and apply_pipeline!" begin
        temp_dir = mktempdir()
        
        # Create test data for comprehensive pipeline
        test_df = DataFrame(
            id = [1, 2, 3, 4, 5, 6],
            score_a = [0.9, 0.7, 0.8, 0.6, 0.95, 0.85],
            score_b = [0.8, 0.9, 0.6, 0.7, 0.9, 0.8],
            qval = [0.01, 0.05, 0.02, 0.08, 0.003, 0.025],
            target = [true, true, false, true, true, true],
            old_name = ["A", "B", "C", "D", "E", "F"]
        )
        test_file = joinpath(temp_dir, "pipeline_test.arrow")
        Arrow.write(test_file, test_df)
        test_ref = PSMFileReference(test_file)
        
        # Create comprehensive pipeline
        pipeline = TransformPipeline() |>
            add_column(:combined_score, row -> (row.score_a + row.score_b) / 2) |>
            filter_by_multiple_thresholds([(:qval, 0.03), (:combined_score, 0.75)]) |>
            filter_rows(row -> row.target) |>
            rename_column(:old_name, :new_name) |>
            select_columns([:id, :combined_score, :new_name]) |>
            sort_by([:combined_score]; rev=true)
        
        # Apply pipeline using apply_pipeline!
        apply_pipeline!(test_ref, pipeline)
        
        # Verify results
        result_df = DataFrame(Arrow.Table(test_file))
        
        # Should have new computed column
        @test hasproperty(result_df, :combined_score)
        @test !hasproperty(result_df, :score_a)  # Removed by select_columns
        @test !hasproperty(result_df, :score_b)  # Removed by select_columns
        
        # Should have renamed column
        @test hasproperty(result_df, :new_name)
        @test !hasproperty(result_df, :old_name)
        
        # Should be filtered and sorted
        @test all(result_df.combined_score .>= 0.75)
        @test issorted(result_df.combined_score, rev=true)
        
        # Should only have 3 columns as specified in select_columns
        @test ncol(result_df) == 3
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Streaming Sort Operations" begin
        temp_dir = mktempdir()
        
        # Create test data that needs sorting
        test_df = DataFrame(
            precursor_idx = UInt32[5, 2, 8, 1, 3],
            prob = Float32[0.6, 0.9, 0.5, 0.95, 0.8],
            target = Bool[true, false, true, true, false],
            protein = ["P3", "P1", "P5", "P1", "P2"]
        )
        test_file = joinpath(temp_dir, "sort_test.arrow")
        Arrow.write(test_file, test_df)
        test_ref = PSMFileReference(test_file)
        
        # Test sort_file_by_keys! with single key
        sort_file_by_keys!(test_ref, :precursor_idx)
        @test is_sorted_by(test_ref, :precursor_idx)
        
        # Verify actual sorting in file
        sorted_df = DataFrame(Arrow.Table(test_file))
        @test issorted(sorted_df.precursor_idx)
        @test sorted_df.precursor_idx == UInt32[1, 2, 3, 5, 8]
        
        # Test multi-key sorting
        sort_file_by_keys!(test_ref, :target, :prob)
        @test is_sorted_by(test_ref, :target, :prob)
        
        # Verify multi-key sort
        sorted_df = DataFrame(Arrow.Table(test_file))
        @test issorted(sorted_df, [:target, :prob])
        
        # Test reverse sorting
        sort_file_by_keys!(test_ref, :prob; reverse=true)
        @test is_sorted_by(test_ref, :prob)
        
        sorted_df = DataFrame(Arrow.Table(test_file))
        @test issorted(sorted_df.prob, rev=true)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Stream Sorted Merge Operations" begin
        temp_dir = mktempdir()
        
        # Create multiple sorted files to merge
        files_data = [
            DataFrame(
                id = UInt32[1, 3, 5],
                score = Float32[0.9, 0.7, 0.5],
                group = ["A", "B", "C"]
            ),
            DataFrame(
                id = UInt32[2, 4, 6],
                score = Float32[0.8, 0.6, 0.4],
                group = ["X", "Y", "Z"]
            ),
            DataFrame(
                id = UInt32[7, 8],
                score = Float32[0.3, 0.2],
                group = ["M", "N"]
            )
        ]
        
        # Write files and create references
        file_refs = PSMFileReference[]
        for (i, df) in enumerate(files_data)
            filepath = joinpath(temp_dir, "merge_file_$i.arrow")
            Arrow.write(filepath, df)
            ref = PSMFileReference(filepath)
            mark_sorted!(ref, :id)  # Mark as sorted by id
            push!(file_refs, ref)
        end
        
        # Test merge operation
        output_file = joinpath(temp_dir, "merged_output.arrow")
        merged_ref = stream_sorted_merge(file_refs, output_file, :id; batch_size=2)
        
        # Verify merge results
        @test exists(merged_ref)
        @test is_sorted_by(merged_ref, :id)
        
        merged_df = DataFrame(Arrow.Table(output_file))
        @test nrow(merged_df) == 8  # Total rows from all files
        @test issorted(merged_df.id)
        @test merged_df.id == UInt32[1, 2, 3, 4, 5, 6, 7, 8]
        
        # Test merge with multiple sort keys
        # First re-sort files by score (descending)
        for ref in file_refs
            sort_file_by_keys!(ref, :score; reverse=true)
        end
        
        output_file2 = joinpath(temp_dir, "merged_output2.arrow")
        merged_ref2 = stream_sorted_merge(file_refs, output_file2, :score; 
                                        batch_size=3, reverse=true)
        
        merged_df2 = DataFrame(Arrow.Table(output_file2))
        @test issorted(merged_df2.score, rev=true)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "MaxLFQ Validation Functions" begin
        temp_dir = mktempdir()
        
        # Create valid MaxLFQ input data
        valid_maxlfq_df = DataFrame(
            inferred_protein_group = ["P1", "P1", "P2", "P2"],
            target = Bool[true, true, false, true],
            entrapment_group_id = UInt8[0, 0, 0, 1],
            precursor_idx = UInt32[1, 2, 3, 4],
            pg_qval = Float32[0.01, 0.02, 0.05, 0.03],
            global_qval_pg = Float32[0.01, 0.02, 0.05, 0.03],
            use_for_protein_quant = Bool[true, true, true, false],
            peak_area = Float32[1000.0, 2000.0, 1500.0, 800.0]
        )
        valid_file = joinpath(temp_dir, "valid_maxlfq.arrow")
        Arrow.write(valid_file, valid_maxlfq_df)
        valid_ref = PSMFileReference(valid_file)
        
        # Test valid input passes validation
        @test validate_maxlfq_input(valid_ref) === nothing
        
        # Test missing required columns
        invalid_df = select(valid_maxlfq_df, Not(:peak_area))
        invalid_file = joinpath(temp_dir, "invalid_maxlfq.arrow")
        Arrow.write(invalid_file, invalid_df)
        invalid_ref = PSMFileReference(invalid_file)
        
        @test_throws ErrorException validate_maxlfq_input(invalid_ref)
        
        # Test parameter validation
        valid_params = Dict(
            :q_value_threshold => 0.05,
            :batch_size => 100000,
            :min_peptides => 1
        )
        @test validate_maxlfq_parameters(valid_params) === nothing
        
        # Test invalid parameters
        invalid_params = Dict(
            :q_value_threshold => 1.5,  # > 1.0
            :batch_size => -100,        # negative
            :min_peptides => 0          # < 1
        )
        @test_throws ErrorException validate_maxlfq_parameters(invalid_params)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Memory-Efficient Operations" begin
        temp_dir = mktempdir()
        
        # Create larger dataset for memory testing
        n_rows = 1000
        large_df = DataFrame(
            id = 1:n_rows,
            value = rand(Float32, n_rows),
            category = rand(["A", "B", "C", "D"], n_rows),
            score = rand(Float32, n_rows),
            flag = rand(Bool, n_rows)
        )
        large_file = joinpath(temp_dir, "large_dataset.arrow")
        Arrow.write(large_file, large_df)
        large_ref = PSMFileReference(large_file)
        
        # Test that operations work with larger datasets
        pipeline = TransformPipeline() |>
            filter_rows(row -> row.value > 0.5) |>
            add_column(:score_category, row -> string(row.category, "_", round(row.score, digits=2))) |>
            sort_by([:score]; rev=true)
        
        # This should complete without memory issues
        apply_pipeline!(large_ref, pipeline)
        
        result_df = DataFrame(Arrow.Table(large_file))
        @test hasproperty(result_df, :score_category)
        @test issorted(result_df.score, rev=true)
        @test all(result_df.value .> 0.5)
        
        # Test merge with multiple larger files
        # Create additional files
        additional_files = PSMFileReference[]
        for i in 1:3
            additional_df = DataFrame(
                id = (i*1000):(i*1000+500),
                value = rand(Float32, 501),
                category = rand(["E", "F"], 501),
                score = rand(Float32, 501),
                flag = rand(Bool, 501)
            )
            additional_file = joinpath(temp_dir, "additional_$i.arrow")
            Arrow.write(additional_file, additional_df)
            additional_ref = PSMFileReference(additional_file)
            sort_file_by_keys!(additional_ref, :id)
            push!(additional_files, additional_ref)
        end
        
        # Sort original file and add to merge list
        sort_file_by_keys!(large_ref, :id)
        all_refs = [large_ref; additional_files]
        
        # Test merge of multiple large files
        merged_file = joinpath(temp_dir, "large_merged.arrow")
        merged_ref = stream_sorted_merge(all_refs, merged_file, :id; batch_size=500)
        
        merged_df = DataFrame(Arrow.Table(merged_file))
        @test issorted(merged_df.id)
        @test nrow(merged_df) > 2000  # Should have all rows
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Error Handling and Edge Cases" begin
        temp_dir = mktempdir()
        
        # Test with empty DataFrame
        empty_df = DataFrame()
        empty_file = joinpath(temp_dir, "empty.arrow")
        Arrow.write(empty_file, empty_df)
        empty_ref = PSMFileReference(empty_file)
        
        # Pipeline operations should handle empty DataFrames gracefully
        pipeline = TransformPipeline() |>
            filter_rows(row -> true) |>
            sort_by(Symbol[])  # Empty sort should be handled
        
        # Should not error on empty DataFrame
        apply_pipeline!(empty_ref, pipeline)
        result_df = DataFrame(Arrow.Table(empty_file))
        @test nrow(result_df) == 0
        
        # Test pipeline with operations that might fail
        test_df = DataFrame(
            a = [1, 2, 3],
            b = [0, 1, 2]  # Note: b[1] = 0 for division test
        )
        test_file = joinpath(temp_dir, "error_test.arrow")
        Arrow.write(test_file, test_df)
        test_ref = PSMFileReference(test_file)
        
        # Test operation that might cause division by zero
        error_pipeline = TransformPipeline() |>
            add_column(:ratio, row -> row.a / row.b)  # Will cause Inf for first row
        
        # Should handle mathematical edge cases
        apply_pipeline!(test_ref, error_pipeline)
        result_df = DataFrame(Arrow.Table(test_file))
        @test hasproperty(result_df, :ratio)
        @test isinf(result_df.ratio[1])  # a[1]/b[1] = 1/0 = Inf
        
        # Test invalid sort keys
        invalid_pipeline = TransformPipeline() |>
            sort_by([:nonexistent_column])
        
        @test_throws Exception apply_pipeline!(test_ref, invalid_pipeline)
        
        # Test merge with incompatible schemas
        df1 = DataFrame(a = [1, 2], b = ["x", "y"])
        df2 = DataFrame(a = [3, 4], c = [1.0, 2.0])  # Different column 'c' instead of 'b'
        
        file1 = joinpath(temp_dir, "schema1.arrow")
        file2 = joinpath(temp_dir, "schema2.arrow")
        Arrow.write(file1, df1)
        Arrow.write(file2, df2)
        
        ref1 = PSMFileReference(file1)
        ref2 = PSMFileReference(file2)
        mark_sorted!(ref1, :a)
        mark_sorted!(ref2, :a)
        
        output_file = joinpath(temp_dir, "merge_error.arrow")
        
        # Should handle schema mismatches gracefully
        @test_throws Exception stream_sorted_merge([ref1, ref2], output_file, :a)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Performance and Scalability" begin
        temp_dir = mktempdir()
        
        # Test with varying batch sizes
        test_df = DataFrame(
            id = 1:100,
            value = rand(100)
        )
        test_file = joinpath(temp_dir, "batch_test.arrow")
        Arrow.write(test_file, test_df)
        test_ref = PSMFileReference(test_file)
        
        # Test different batch sizes for merge operations
        # Create multiple files for merging
        merge_refs = PSMFileReference[]
        for i in 1:5
            df = DataFrame(
                id = (i*100):((i+1)*100-1),
                value = rand(100)
            )
            filepath = joinpath(temp_dir, "batch_file_$i.arrow")
            Arrow.write(filepath, df)
            ref = PSMFileReference(filepath)
            sort_file_by_keys!(ref, :id)
            push!(merge_refs, ref)
        end
        
        # Test with small batch size
        output1 = joinpath(temp_dir, "merged_small_batch.arrow")
        merged1 = stream_sorted_merge(merge_refs, output1, :id; batch_size=10)
        
        # Test with large batch size
        output2 = joinpath(temp_dir, "merged_large_batch.arrow")
        merged2 = stream_sorted_merge(merge_refs, output2, :id; batch_size=200)
        
        # Both should produce identical results
        df1 = DataFrame(Arrow.Table(output1))
        df2 = DataFrame(Arrow.Table(output2))
        @test df1 == df2
        @test issorted(df1.id)
        @test issorted(df2.id)
        
        # Test pipeline performance with complex operations
        complex_pipeline = TransformPipeline() |>
            add_column(:computed1, row -> sin(row.value)) |>
            add_column(:computed2, row -> cos(row.value)) |>
            add_column(:computed3, row -> row.computed1 + row.computed2) |>
            filter_rows(row -> row.computed3 > -1.0) |>
            sort_by([:computed3]; rev=true)
        
        # Should complete in reasonable time
        start_time = time()
        apply_pipeline!(test_ref, complex_pipeline)
        elapsed = time() - start_time
        
        @test elapsed < 5.0  # Should complete within 5 seconds
        
        result_df = DataFrame(Arrow.Table(test_file))
        @test hasproperty(result_df, :computed3)
        @test issorted(result_df.computed3, rev=true)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end