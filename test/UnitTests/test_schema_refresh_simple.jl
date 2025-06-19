"""
Simple test that FileReference schema updates correctly when columns are added.
"""

using Test
using DataFrames
using Arrow

# Include necessary modules
include("../../src/Routines/SearchDIA/SearchMethods/FileReferences.jl")
include("../../src/Routines/SearchDIA/SearchMethods/SearchResultReferences.jl")
include("../../src/Routines/SearchDIA/SearchMethods/FileOperations.jl")
include("../../src/Routines/SearchDIA/SearchMethods/ScoringSearch/scoring_interface.jl")
include("../../src/utils/writeArrow.jl")

# Import ProteinKey from protein_inference_types
include("../../src/structs/protein_inference_types.jl")

@testset "FileReference Schema Refresh Simple" begin
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
    
    # Test the calculate_and_add_global_scores! function
    refs = [ref]
    acc_to_max_pg_score = calculate_and_add_global_scores!(refs)
    
    # Test that schema was updated
    @test has_column(schema(ref), :global_pg_score)
    @test row_count(ref) == 3  # Row count should be unchanged
    
    # Verify the function returned the correct scores
    @test length(acc_to_max_pg_score) == 3
    @test acc_to_max_pg_score[ProteinKey("Protein1", true, UInt8(0))] == 10.0f0
    
    # Read back and verify data
    df_result = DataFrame(Arrow.Table(test_path))
    @test size(df_result) == (3, 5)  # 5 columns now
    @test "global_pg_score" in names(df_result)
    @test df_result.global_pg_score == [10.0f0, 8.0f0, 5.0f0]
    
    # Test that the file is sorted by global_pg_score
    @test issorted(df_result.global_pg_score, rev=true)
    
    # Test write_arrow_file function
    test_path2 = tempname() * ".arrow"
    df2 = DataFrame(a=[1,2,3], b=[4,5,6])
    writeArrow(test_path2, df2)  # Create file first
    ref2 = ProteinGroupFileReference(test_path2)
    write_arrow_file(ref2, df2)  # Now update it
    
    @test exists(ref2)
    @test row_count(ref2) == 3
    @test has_column(schema(ref2), :a)
    @test has_column(schema(ref2), :b)
    
    # Cleanup
    rm(test_path, force=true)
    rm(test_path2, force=true)
end

println("Simple schema refresh tests passed!")