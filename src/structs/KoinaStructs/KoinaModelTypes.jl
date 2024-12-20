abstract type KoinaModelType end
struct InstrumentSpecificModel <: KoinaModelType
    name::String 
end

struct InstrumentAgnosticModel <: KoinaModelType
    name::String
end

struct SplineCoefficientModel <: KoinaModelType 
    name::String
end

struct RetentionTimeModel <: KoinaModelType
    name::String 
end

abstract type FragAnnotation end
struct UniSpecFragAnnotation <: FragAnnotation
    annotation::String
end

struct GenericFragAnnotation <: FragAnnotation
    annotation::String
end
getAnnotation(fa::FragAnnotation) = fa.annotation


abstract type FragRegexIterator end
struct UniSpecFragRegexIterator <: FragRegexIterator
    annotation_pieces::Base.RegexMatchIterator
end
struct GenericFragRegexIterator <: FragRegexIterator
    annotation_pieces::Base.RegexMatchIterator
end
getAnnotationPieces(fri::FragRegexIterator) = fri.annotation_pieces

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
