"""
Response structure for Koina batch processing results
"""
struct KoinaBatchResult{T}
    fragments::DataFrame
    frags_per_precursor::Int64
    extra_data::T
end

"""
Custom error type for Koina API issues
"""
struct KoinaRequestError <: Exception
    message::String
    attempt::Int
    response::Union{Nothing, Dict{String, Any}}
end
