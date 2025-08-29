# Test suite for Koina API integration tests
# Main entry point that includes all Koina-related tests
# Uses real Koina API calls for reliable integration testing

using Test, DataFrames, HTTP, JSON
using Pioneer

# KOINA_URLS is already defined in runtests.jl

@testset "Koina API Integration Tests" begin
    include("test_koina_api.jl") 
    include("test_koina_batch_prep.jl")
    include("test_koina_batch_parse.jl")
    include("test_fragment_predict_integration.jl")
end