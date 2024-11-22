abstract type KoinaModelType end
struct InstrumentSpecificModel <: KoinaModelType
    model_name::String
end
struct InstrumentAgnosticModel <: KoinaModelType
    model_name::String
end
struct RetentionTimeModel <: KoinaModelType
    model_name::String
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
