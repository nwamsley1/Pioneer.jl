# Standalone test for EFDR trait system
using Pkg
Pkg.activate(".")
using Test, DataFrames

# Load necessary functions (without the @test lines from entrapment_helper_funcs.jl)
include("entrapment_helper_funcs.jl")
include("efdr_funcs.jl")

# Run the tests
include("test_efdr_traits.jl")