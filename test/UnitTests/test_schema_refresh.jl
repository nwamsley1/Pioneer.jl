"""
Test that FileReference schema updates correctly when columns are added.
"""

using Test
using DataFrames
using Arrow

# Include necessary modules
include("../../src/Routines/SearchDIA/SearchMethods/FileReferences.jl")
include("../../src/Routines/SearchDIA/SearchMethods/FileOperations.jl")
include("../../src/utils/writeArrow.jl")

@testset "FileReference Schema Refresh" begin
    # Create a temporary test file
    test_path = tempname() * ".arrow"
    
    # Create initial DataFrame
    df = DataFrame(
        protein_name = ["Protein1", "Protein2", "Protein3"],
        target = [true, true, false],
        entrap_id = UInt8[0, 0, 1],
        pg_score = Float32[10.0, 8.0, 5.0]
    )
    
    # Write initial file
    writeArrow(test_path, df)
    
    # Create reference
    ref = ProteinGroupFileReference(test_path)
    
    # Test initial state
    @test exists(ref)
    @test row_count(ref) == 3
    @test has_column(schema(ref), :protein_name)
    @test has_column(schema(ref), :pg_score)
    @test !has_column(schema(ref), :global_pg_score)
    
    # Add new column using FileOperations
    add_column_to_file!(ref, :global_pg_score, 
        batch -> batch.pg_score .* 2.0f0
    )
    
    # Test that schema was updated
    @test has_column(schema(ref), :global_pg_score)
    @test row_count(ref) == 3  # Row count should be unchanged
    
    # Verify we can sort by the new column
    sort_file_by_keys!(ref, :global_pg_score; reverse=true)
    @test is_sorted_by(ref, :global_pg_score)
    
    # Read back and verify data
    df_result = DataFrame(Arrow.Table(test_path))
    @test size(df_result) == (3, 5)  # 5 columns now
    @test df_result.global_pg_score == Float32[20.0, 16.0, 10.0]
    @test df_result.protein_name[1] == "Protein1"  # Should be first after sort
    
    # Test add_column_and_sort!
    test_path2 = tempname() * ".arrow"
    writeArrow(test_path2, df)
    ref2 = ProteinGroupFileReference(test_path2)
    
    add_column_and_sort!(ref2, :new_score, 
        batch -> batch.pg_score .+ 1.0f0,
        :new_score, :target;
        reverse=true
    )
    
    @test has_column(schema(ref2), :new_score)
    @test is_sorted_by(ref2, :new_score, :target)
    
    # Cleanup
    rm(test_path, force=true)
    rm(test_path2, force=true)
end

println("Schema refresh tests passed!")