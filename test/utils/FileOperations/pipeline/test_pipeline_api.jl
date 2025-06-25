
@testset "Pipeline API Tests" begin
    
    @testset "TransformPipeline Core Functionality" begin
        # Test pipeline creation and composition
        pipeline = TransformPipeline()
        @test length(pipeline.operations) == 0
        @test length(pipeline.post_actions) == 0
        
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
    
    @testset "Column Operations in Pipeline" begin
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
        
        result_df = op_func(DataFrame(Tables.columntable(Arrow.Table(test_file))))
        @test hasproperty(result_df, :sum_values)
        @test result_df.sum_values ≈ [11.5, 22.5, 33.5, 44.5]
        
        # Test rename_column operation
        rename_op = rename_column(:value_a, :renamed_value)
        desc, op_func = rename_op
        @test startswith(desc, "rename")
        
        result_df = op_func(DataFrame(Tables.columntable(Arrow.Table(test_file))))
        @test hasproperty(result_df, :renamed_value)
        @test !hasproperty(result_df, :value_a)
        
        # Test select_columns operation
        select_op = select_columns([:id, :value_a])
        desc, op_func = select_op
        @test startswith(desc, "select")
        
        result_df = op_func(DataFrame(Tables.columntable(Arrow.Table(test_file))))
        @test ncol(result_df) == 2
        @test hasproperty(result_df, :id)
        @test hasproperty(result_df, :value_a)
        @test !hasproperty(result_df, :value_b)
        
        # Test remove_columns operation
        remove_op = remove_columns(:value_b, :category)
        desc, op_func = remove_op
        @test startswith(desc, "remove")
        
        result_df = op_func(DataFrame(Tables.columntable(Arrow.Table(test_file))))
        @test !hasproperty(result_df, :value_b)
        @test !hasproperty(result_df, :category)
        @test hasproperty(result_df, :id)
        @test hasproperty(result_df, :value_a)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Row Operations in Pipeline" begin
        temp_dir = mktempdir()
        
        # Create test data with filtering scenarios
        test_df = DataFrame(
            id = [1, 2, 3, 4, 5],
            score = Float32[0.9, 0.7, 0.8, 0.6, 0.95],
            qval = Float32[0.01, 0.03, 0.02, 0.06, 0.005],
            target = [true, true, false, true, true],
            use_for_quant = [true, true, true, false, true]
        )
        test_file = joinpath(temp_dir, "row_ops.arrow")
        Arrow.write(test_file, test_df)
        
        # Test basic filter_rows operation
        filter_op = filter_rows(row -> row.score > 0.75; desc="high_score_filter")
        desc, op_func = filter_op
        @test desc == "high_score_filter"
        
        result_df = op_func(DataFrame(Tables.columntable(Arrow.Table(test_file))))
        @test nrow(result_df) == 3  # Scores: 0.9, 0.8, 0.95
        @test all(result_df.score .> 0.75)
        
        # Test combined filtering with pipeline composition
        # Use separate filters for mixed comparison operators (score >= 0.75 AND qval <= 0.03)
        combined_pipeline = TransformPipeline() |>
            filter_rows(row -> row.score >= 0.75) |>
            filter_rows(row -> row.qval <= 0.03)
        
        apply_pipeline!(PSMFileReference(test_file), combined_pipeline)
        result_df = DataFrame(Tables.columntable(Arrow.Table(test_file)))
        @test all(result_df.score .>= 0.75)
        @test all(result_df.qval .<= 0.03)
        @test nrow(result_df) == 3  # Rows 1, 3, and 5 meet both criteria (0.9,0.01), (0.8,0.02), (0.95,0.005)
        
        # Test complex filter with multiple conditions
        complex_filter_op = filter_rows(row -> row.target && row.use_for_quant && row.qval < 0.04)
        desc, op_func = complex_filter_op
        
        result_df = op_func(DataFrame(Tables.columntable(Arrow.Table(test_file))))
        @test all(result_df.target)
        @test all(result_df.use_for_quant)
        @test all(result_df.qval .< 0.04)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Sort Operations in Pipeline" begin
        temp_dir = mktempdir()
        
        # Create unsorted test data
        test_df = DataFrame(
            id = [3, 1, 4, 2],
            score = [0.7, 0.9, 0.6, 0.8],
            category = ["B", "A", "C", "A"]
        )
        test_file = joinpath(temp_dir, "sort_ops.arrow")
        Arrow.write(test_file, test_df)
        test_ref = PSMFileReference(test_file)
        
        # Test single key sort
        sort_op = sort_by([:id])
        desc, op_func = sort_op.operation
        @test startswith(desc, "sort")
        
        result_df = op_func(DataFrame(Tables.columntable(Arrow.Table(test_file))))
        @test issorted(result_df.id)
        @test result_df.id == [1, 2, 3, 4]
        
        # Test multi-key sort with mixed directions
        multi_sort_op = sort_by([:category, :score]; rev=[false, true])
        desc, op_func = multi_sort_op.operation
        
        result_df = op_func(DataFrame(Tables.columntable(Arrow.Table(test_file))))
        @test issorted(result_df, [:category, :score], rev=[false, true])
        
        # Test reverse sort
        rev_sort_op = sort_by([:score]; rev=[true])
        desc, op_func = rev_sort_op.operation
        
        result_df = op_func(DataFrame(Tables.columntable(Arrow.Table(test_file))))
        @test issorted(result_df.score, rev=true)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Complex Pipeline Composition" begin
        temp_dir = mktempdir()
        
        # Create complex test dataset
        test_df = DataFrame(
            id = [1, 2, 3, 4, 5, 6],
            score_a = [0.8, 0.6, 0.9, 0.7, 0.5, 0.85],
            score_b = [0.7, 0.8, 0.6, 0.9, 0.7, 0.8],
            qval = [0.01, 0.05, 0.02, 0.03, 0.08, 0.01],
            target = [true, true, false, true, true, false],
            category = ["A", "B", "A", "B", "A", "B"]
        )
        test_file = joinpath(temp_dir, "complex.arrow")
        Arrow.write(test_file, test_df)
        test_ref = PSMFileReference(test_file)
        
        # Build complex pipeline
        pipeline = TransformPipeline() |>
            add_column(:avg_score, row -> (row.score_a + row.score_b) / 2) |>
            filter_rows(row -> row.target && row.qval < 0.05) |>
            add_column(:score_rank, row -> row.avg_score > 0.75 ? "high" : "medium") |>
            select_columns([:id, :avg_score, :score_rank, :category, :target, :qval]) |>
            sort_by([:avg_score]; rev=[true])
        
        # Apply pipeline
        apply_pipeline!(test_ref, pipeline)
        
        result_df = DataFrame(Tables.columntable(Arrow.Table(test_file)))
        
        # Verify pipeline effects
        @test hasproperty(result_df, :avg_score)
        @test hasproperty(result_df, :score_rank)
        @test !hasproperty(result_df, :score_a)  # Should be removed by select
        @test !hasproperty(result_df, :score_b)
        @test all(result_df.target)  # Only target=true should remain
        @test all(result_df.qval .< 0.05)  # Only qval < 0.05 should remain
        @test issorted(result_df.avg_score, rev=true)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Pipeline Error Handling" begin
        temp_dir = mktempdir()
        
        test_df = DataFrame(a = [1, 2], b = [3, 4])
        test_file = joinpath(temp_dir, "error_test.arrow")
        Arrow.write(test_file, test_df)
        test_ref = PSMFileReference(test_file)
        
        # Test pipeline with invalid column reference
        bad_pipeline = TransformPipeline() |>
            add_column(:result, row -> row.nonexistent_col + 1)
        
        @test_throws Exception apply_pipeline!(test_ref, bad_pipeline)
        
        # Test pipeline with invalid sort column
        bad_sort_pipeline = TransformPipeline() |>
            sort_by([:nonexistent_col])
        
        @test_throws ArgumentError apply_pipeline!(test_ref, bad_sort_pipeline)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Pipeline Performance with Large Data" begin
        temp_dir = mktempdir()
        
        # Create larger dataset
        n_rows = 10000
        large_df = DataFrame(
            id = 1:n_rows,
            value = rand(n_rows),
            category = rand(["A", "B", "C", "D"], n_rows),
            score = rand(n_rows),
            flag = rand(Bool, n_rows)
        )
        large_file = joinpath(temp_dir, "large_dataset.arrow")
        Arrow.write(large_file, large_df)
        large_ref = PSMFileReference(large_file)
        
        # Test performance of complex pipeline
        pipeline = TransformPipeline() |>
            filter_rows(row -> row.value > 0.5) |>
            add_column(:score_category, row -> string(row.category, "_", round(row.score, digits=2))) |>
            sort_by([:score]; rev=[true])
        
        # This should complete efficiently
        @time apply_pipeline!(large_ref, pipeline)
        
        result_df = DataFrame(Tables.columntable(Arrow.Table(large_file)))
        @test hasproperty(result_df, :score_category)
        @test issorted(result_df.score, rev=true)
        @test all(result_df.value .> 0.5)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end