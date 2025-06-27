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


@testset "Merge Operations Tests" begin
    
    @testset "Basic Two-File Merge" begin
        temp_dir = mktempdir()
        
        # Create two sorted files
        df1 = DataFrame(
            id = UInt32[1, 3, 5],
            score = Float32[0.9, 0.7, 0.8],
            group = ["A", "B", "A"]
        )
        df2 = DataFrame(
            id = UInt32[2, 4, 6],
            score = Float32[0.85, 0.6, 0.75],
            group = ["B", "A", "B"]
        )
        
        file1 = joinpath(temp_dir, "file1.arrow")
        file2 = joinpath(temp_dir, "file2.arrow")
        Arrow.write(file1, df1)
        Arrow.write(file2, df2)
        
        ref1 = PSMFileReference(file1)
        ref2 = PSMFileReference(file2)
        mark_sorted!(ref1, :id)
        mark_sorted!(ref2, :id)
        
        # Test merge
        output_file = joinpath(temp_dir, "merged.arrow")
        merged_ref = stream_sorted_merge([ref1, ref2], output_file, :id)
        
        @test exists(merged_ref)
        @test is_sorted_by(merged_ref, :id)
        
        result_df = DataFrame(Arrow.Table(output_file))
        @test nrow(result_df) == 6
        @test issorted(result_df.id)
        @test result_df.id == UInt32[1, 2, 3, 4, 5, 6]
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Multi-File Merge Operations" begin
        temp_dir = mktempdir()
        
        # Create multiple files with overlapping ranges
        files_data = [
            DataFrame(
                id = UInt32[1, 2],
                score = Float32[0.9, 0.8],
                group = ["A", "B"]
            ),
            DataFrame(
                id = UInt32[3, 4],
                score = Float32[0.7, 0.6],
                group = ["C", "D"]
            ),
            DataFrame(
                id = UInt32[5, 6],
                score = Float32[0.5, 0.4],
                group = ["E", "F"]
            ),
            DataFrame(
                id = UInt32[7, 8],
                score = Float32[0.3, 0.2],
                group = ["G", "H"]
            )
        ]
        
        # Write files and create references
        file_refs = PSMFileReference[]
        for (i, df) in enumerate(files_data)
            filepath = joinpath(temp_dir, "merge_file_$i.arrow")
            Arrow.write(filepath, df)
            ref = PSMFileReference(filepath)
            mark_sorted!(ref, :id)
            push!(file_refs, ref)
        end
        
        # Test merge operation
        output_file = joinpath(temp_dir, "merged_output.arrow")
        merged_ref = stream_sorted_merge(file_refs, output_file, :id; batch_size=2)
        
        # Verify merge results
        @test exists(merged_ref)
        @test is_sorted_by(merged_ref, :id)
        
        merged_df = DataFrame(Arrow.Table(output_file))
        @test nrow(merged_df) == 8
        @test issorted(merged_df.id)
        @test merged_df.id == UInt32[1, 2, 3, 4, 5, 6, 7, 8]
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Multi-Key Merge Operations" begin
        temp_dir = mktempdir()
        
        # Create files sorted by multiple keys
        df1 = DataFrame(
            protein = ["P1", "P1", "P2"],
            target = Bool[true, false, true],
            score = Float32[0.9, 0.8, 0.7]
        )
        df2 = DataFrame(
            protein = ["P1", "P2", "P3"],
            target = Bool[true, true, false],
            score = Float32[0.85, 0.75, 0.65]
        )
        
        file1 = joinpath(temp_dir, "multi1.arrow")
        file2 = joinpath(temp_dir, "multi2.arrow")
        
        # Sort by protein, then target (descending), then score (descending)
        sort!(df1, [:protein, :target, :score], rev=[false, true, true])
        sort!(df2, [:protein, :target, :score], rev=[false, true, true])
        
        Arrow.write(file1, df1)
        Arrow.write(file2, df2)
        
        ref1 = PSMFileReference(file1)
        ref2 = PSMFileReference(file2)
        mark_sorted!(ref1, :protein, :target, :score)
        mark_sorted!(ref2, :protein, :target, :score)
        
        # Test multi-key merge
        output_file = joinpath(temp_dir, "multi_merged.arrow")
        merged_ref = stream_sorted_merge([ref1, ref2], output_file, :protein, :target, :score; 
                                       reverse=[false, true, true])
        
        @test exists(merged_ref)
        result_df = DataFrame(Arrow.Table(output_file))
        @test issorted(result_df, [:protein, :target, :score], rev=[false, true, true])
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Reverse Sort Merge" begin
        temp_dir = mktempdir()
        
        # Create files sorted in descending order
        df1 = DataFrame(
            score = Float32[0.9, 0.7, 0.5],
            id = [1, 2, 3]
        )
        df2 = DataFrame(
            score = Float32[0.8, 0.6, 0.4],
            id = [4, 5, 6]
        )
        
        file1 = joinpath(temp_dir, "rev1.arrow")
        file2 = joinpath(temp_dir, "rev2.arrow")
        Arrow.write(file1, df1)
        Arrow.write(file2, df2)
        
        ref1 = PSMFileReference(file1)
        ref2 = PSMFileReference(file2)
        mark_sorted!(ref1, :score)
        mark_sorted!(ref2, :score)
        
        # Test reverse merge
        output_file = joinpath(temp_dir, "rev_merged.arrow")
        merged_ref = stream_sorted_merge([ref1, ref2], output_file, :score; reverse=true)
        
        @test exists(merged_ref)
        result_df = DataFrame(Arrow.Table(output_file))
        @test issorted(result_df.score, rev=true)
        expected_scores = Float32[0.9, 0.8, 0.7, 0.6, 0.5, 0.4]
        @test result_df.score == expected_scores
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Large File Merge Performance" begin
        temp_dir = mktempdir()
        
        # Create moderately large files for performance testing
        n_per_file = 1000
        num_files = 5
        
        file_refs = PSMFileReference[]
        for i in 1:num_files
            # Create non-overlapping ranges for predictable merge
            start_id = (i-1) * n_per_file + 1
            end_id = i * n_per_file
            
            df = DataFrame(
                id = UInt32.(start_id:end_id),
                value = rand(Float32, n_per_file),
                category = rand(["A", "B", "C"], n_per_file)
            )
            
            filepath = joinpath(temp_dir, "large_file_$i.arrow")
            Arrow.write(filepath, df)
            ref = PSMFileReference(filepath)
            mark_sorted!(ref, :id)
            push!(file_refs, ref)
        end
        
        # Test merge with batch processing
        output_file = joinpath(temp_dir, "large_merged.arrow")
        @time merged_ref = stream_sorted_merge(file_refs, output_file, :id; batch_size=500)
        
        @test exists(merged_ref)
        result_df = DataFrame(Arrow.Table(output_file))
        @test nrow(result_df) == num_files * n_per_file
        @test issorted(result_df.id)
        @test result_df.id[1] == 1
        @test result_df.id[end] == num_files * n_per_file
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Merge Error Handling" begin
        temp_dir = mktempdir()
        
        # Test with unsorted files
        df1 = DataFrame(id = [3, 1, 2], value = [1, 2, 3])
        df2 = DataFrame(id = [6, 4, 5], value = [4, 5, 6])
        
        file1 = joinpath(temp_dir, "unsorted1.arrow")
        file2 = joinpath(temp_dir, "unsorted2.arrow")
        Arrow.write(file1, df1)
        Arrow.write(file2, df2)
        
        ref1 = PSMFileReference(file1)
        ref2 = PSMFileReference(file2)
        # Don't mark as sorted - should cause error
        
        output_file = joinpath(temp_dir, "error_output.arrow")
        @test_throws ErrorException stream_sorted_merge([ref1, ref2], output_file, :id)
        
        # Test with missing columns
        df3 = DataFrame(different_col = [1, 2, 3])
        file3 = joinpath(temp_dir, "missing_col.arrow")
        Arrow.write(file3, df3)
        ref3 = PSMFileReference(file3)
        mark_sorted!(ref3, :id)  # Mark as sorted by :id even though it doesn't have that column
        
        # Create a properly sorted file to test with ref3
        df4 = DataFrame(id = [10, 11, 12], value = [7, 8, 9])
        file4 = joinpath(temp_dir, "sorted.arrow")
        Arrow.write(file4, df4)
        ref4 = PSMFileReference(file4)
        mark_sorted!(ref4, :id)
        
        @test_throws BoundsError stream_sorted_merge([ref4, ref3], output_file, :id)
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
    
    @testset "Empty File Merge Handling" begin
        temp_dir = mktempdir()
        
        # Create one empty file and one normal file
        empty_df = DataFrame(id = UInt32[], value = String[])
        normal_df = DataFrame(id = UInt32[1, 2, 3], value = ["a", "b", "c"])
        
        empty_file = joinpath(temp_dir, "empty.arrow")
        normal_file = joinpath(temp_dir, "normal.arrow")
        Arrow.write(empty_file, empty_df)
        Arrow.write(normal_file, normal_df)
        
        empty_ref = PSMFileReference(empty_file)
        normal_ref = PSMFileReference(normal_file)
        mark_sorted!(empty_ref, :id)
        mark_sorted!(normal_ref, :id)
        
        # Test merge with empty file
        output_file = joinpath(temp_dir, "mixed_output.arrow")
        merged_ref = stream_sorted_merge([empty_ref, normal_ref], output_file, :id)
        
        @test exists(merged_ref)
        result_df = DataFrame(Arrow.Table(output_file))
        @test nrow(result_df) == 3  # Only rows from normal file
        @test result_df.id == UInt32[1, 2, 3]
        
        # Clean up
        rm(temp_dir, recursive=true)
    end
end