"""
FileOperations Test Suite Runner

This file runs all tests for the modular FileOperations system.
It serves as the entry point for testing the entire FileOperations
functionality in the new modular structure.

Usage:
```julia
include("test/utils/FileOperations/test_file_operations_suite.jl")
```

Test Structure:
- core/: Tests for FileSchema, FileReferences, SortStateManagement
- streaming/: Tests for StreamOperations, ColumnOperations, MergeOperations  
- pipeline/: Tests for PipelineAPI and transformation operations
- algorithms/: Tests for ProteinInference and PSM update algorithms
- io/: Tests for ValidationFunctions and file I/O operations
"""

@testset "FileOperations Comprehensive Test Suite" begin
    
    # Core module tests
    @testset "Core FileOperations Tests" begin
        include("core/test_file_schema.jl")
        include("core/test_file_references.jl") 
        include("core/test_sort_state_management.jl")
    end
    
    # Streaming operations tests
    @testset "Streaming Operations Tests" begin
        include("streaming/test_stream_operations.jl")
        include("streaming/test_merge_operations.jl")
    end
    
    # Pipeline API tests
    @testset "Pipeline API Tests" begin
        include("pipeline/test_pipeline_api.jl")
    end
    
    # Algorithm tests
    @testset "Algorithm Tests" begin
        include("algorithms/test_protein_inference.jl")
        include("algorithms/test_psm_updates.jl")
    end
    
    # I/O and validation tests
    @testset "I/O and Validation Tests" begin
        include("io/test_validation_functions.jl")
    end
    
    println("✅ All FileOperations tests completed successfully!")
end