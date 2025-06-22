using Test
using Arrow, DataFrames

# Include the streaming operations module
package_root = dirname(dirname(dirname(dirname(@__DIR__))))
include(joinpath(package_root, "src", "utils", "FileOperations", "FileOperations.jl"))

@testset "Stream Operations Tests" begin
    
    @testset "Basic Streaming Operations" begin
        temp_dir = mktempdir()
        
        # Create test data
        test_df = DataFrame(
            id = [1, 2, 3, 4, 5],
            score = [0.9, 0.7, 0.8, 0.6, 0.95],
            qval = [0.01, 0.03, 0.02, 0.06, 0.005],
            target = [true, true, false, true, true]
        )
        input_file = joinpath(temp_dir, "input.arrow")
        Arrow.write(input_file, test_df)
        
        input_ref = PSMFileReference(input_file)
        output_file = joinpath(temp_dir, "output.arrow")
        
        # Test basic filter operation
        filter_func = row -> row.score > 0.75
        output_ref = stream_filter_rows(input_ref, output_file, filter_func)
        
        @test exists(output_ref)
        result_df = DataFrame(Arrow.Table(file_path(output_ref)))
        @test nrow(result_df) == 3  # Scores: 0.9, 0.8, 0.95
        @test all(result_df.score .> 0.75)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Memory-Efficient Batch Processing" begin
        temp_dir = mktempdir()
        
        # Create larger test dataset
        n_rows = 1000
        large_df = DataFrame(
            id = 1:n_rows,
            value = rand(n_rows),
            group = rand(["A", "B", "C"], n_rows)
        )
        large_file = joinpath(temp_dir, "large_input.arrow")
        Arrow.write(large_file, large_df)
        
        input_ref = PSMFileReference(large_file)
        output_file = joinpath(temp_dir, "large_output.arrow")
        
        # Test with small batch size for memory efficiency
        filter_func = row -> row.value > 0.5
        output_ref = stream_filter_rows(input_ref, output_file, filter_func, batch_size=100)
        
        @test exists(output_ref)
        result_df = DataFrame(Arrow.Table(file_path(output_ref)))
        @test all(result_df.value .> 0.5)
        @test nrow(result_df) < n_rows  # Should have filtered some rows
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Stream Column Operations" begin
        temp_dir = mktempdir()
        
        # Create test data
        test_df = DataFrame(
            a = [1, 2, 3, 4],
            b = [10, 20, 30, 40],
            category = ["X", "Y", "X", "Y"]
        )
        input_file = joinpath(temp_dir, "column_input.arrow")
        Arrow.write(input_file, test_df)
        
        input_ref = PSMFileReference(input_file)
        output_file = joinpath(temp_dir, "column_output.arrow")
        
        # Test adding a computed column
        add_func = row -> row.a + row.b
        output_ref = stream_add_column(input_ref, output_file, :sum, add_func)
        
        @test exists(output_ref)
        result_df = DataFrame(Arrow.Table(file_path(output_ref)))
        @test hasproperty(result_df, :sum)
        @test result_df.sum == [11, 22, 33, 44]
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Stream Select Columns" begin
        temp_dir = mktempdir()
        
        # Create test data with many columns
        test_df = DataFrame(
            id = [1, 2, 3],
            col_a = [10, 20, 30],
            col_b = [1.5, 2.5, 3.5],
            col_c = ["X", "Y", "Z"],
            col_d = [true, false, true],
            unwanted = [99, 99, 99]
        )
        input_file = joinpath(temp_dir, "select_input.arrow")
        Arrow.write(input_file, test_df)
        
        input_ref = PSMFileReference(input_file)
        output_file = joinpath(temp_dir, "select_output.arrow")
        
        # Test selecting specific columns
        selected_cols = [:id, :col_a, :col_c]
        output_ref = stream_select_columns(input_ref, output_file, selected_cols)
        
        @test exists(output_ref)
        result_df = DataFrame(Arrow.Table(file_path(output_ref)))
        @test ncol(result_df) == 3
        @test hasproperty(result_df, :id)
        @test hasproperty(result_df, :col_a)
        @test hasproperty(result_df, :col_c)
        @test !hasproperty(result_df, :col_b)
        @test !hasproperty(result_df, :unwanted)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Error Handling in Stream Operations" begin
        temp_dir = mktempdir()
        
        # Test with non-existent input file
        fake_ref = PSMFileReference(joinpath(temp_dir, "nonexistent.arrow"))
        output_file = joinpath(temp_dir, "output.arrow")
        
        @test_throws SystemError stream_filter_rows(fake_ref, output_file, row -> true)
        
        # Test with invalid column selection
        test_df = DataFrame(a = [1, 2], b = [3, 4])
        input_file = joinpath(temp_dir, "test.arrow")
        Arrow.write(input_file, test_df)
        input_ref = PSMFileReference(input_file)
        
        # This should work fine - missing columns are ignored in stream operations
        output_ref = stream_select_columns(input_ref, output_file, [:a, :nonexistent])
        @test exists(output_ref)
        result_df = DataFrame(Arrow.Table(file_path(output_ref)))
        @test hasproperty(result_df, :a)
        @test !hasproperty(result_df, :nonexistent)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end