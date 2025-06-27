# Comprehensive test suite for BuildSpecLib/fragments functionalit

@testset "BuildSpecLib Fragments Test Suite" begin
    # Run all fragment-related tests
    include("test_fragment_annotation.jl")
    include("test_fragment_parse.jl") 
    include("test_fragment_predict.jl")
    include("test_get_frag_bounds.jl")
end

# To run this test suite:
# julia> include("test/Routines/BuildSpecLib/fragments/test_fragments_suite.jl")