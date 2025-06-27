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
Tests for actual streaming operations that are implemented and used.

Only tests functions that exist in the codebase:
- stream_transform (used in ProteinInference.jl)
- process_with_memory_limit (used in scoring_interface.jl)  
- estimate_batch_size (used in scoring_interface.jl)
"""

@testset "Stream Operations Tests" begin
    
    @testset "Stream Transform Function" begin
        temp_dir = mktempdir()
        
        # Create test data
        test_df = DataFrame(
            id = [1, 2, 3, 4, 5],
            score = [0.9, 0.7, 0.8, 0.6, 0.95],
            target = [true, true, false, true, true]
        )
        input_file = joinpath(temp_dir, "input.arrow")
        Arrow.write(input_file, test_df)
        
        input_ref = PSMFileReference(input_file)
        output_file = joinpath(temp_dir, "output.arrow")
        
        # Test stream_transform with a simple transformation
        transform_fn = function(df)
            # Add a computed column and filter
            df.doubled_score = df.score .* 2
            filter!(row -> row.score > 0.7, df)
            return df
        end
        output_ref = stream_transform(input_ref, output_file, transform_fn)
        
        @test exists(output_ref)
        result_df = DataFrame(Arrow.Table(file_path(output_ref)))
        @test hasproperty(result_df, :doubled_score)
        @test all(result_df.score .> 0.7)
        @test nrow(result_df) == 3  # Should have 3 rows with score > 0.7
        @test result_df.doubled_score ≈ result_df.score .* 2
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Memory Limit Processing" begin
        temp_dir = mktempdir()
        
        # Create test data
        test_df = DataFrame(
            id = 1:100,
            value = rand(100)
        )
        test_file = joinpath(temp_dir, "test.arrow")
        Arrow.write(test_file, test_df)
        
        test_ref = PSMFileReference(test_file)
        
        # Test process_with_memory_limit
        result_sum = Ref(0.0)  # Use Ref for mutable reference
        process_fn = function(df)
            result_sum[] = sum(df.value)
        end
        process_with_memory_limit(test_ref, process_fn; max_memory_mb=10)
        
        @test result_sum[] > 0.0  # Should have processed the data
        @test result_sum[] ≈ sum(test_df.value)  # Should match original sum
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Batch Size Estimation" begin
        # Create a mock schema for testing
        test_schema = FileSchema([:id, :score, :name, :target])
        
        # Test batch size estimation
        batch_size_small = estimate_batch_size(test_schema, 10)  # 10 MB limit
        batch_size_large = estimate_batch_size(test_schema, 100)  # 100 MB limit
        
        @test batch_size_small > 1000  # Should be reasonable size
        @test batch_size_large > batch_size_small  # Larger memory should allow larger batches
        @test batch_size_small isa Int
        @test batch_size_large isa Int
    end
    
    @testset "Error Handling" begin
        temp_dir = mktempdir()
        
        # Test with non-existent input file
        fake_ref = PSMFileReference(joinpath(temp_dir, "nonexistent.arrow"))
        output_file = joinpath(temp_dir, "output.arrow")
        
        identity_fn = df -> df
        @test_throws Exception stream_transform(fake_ref, output_file, identity_fn)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end