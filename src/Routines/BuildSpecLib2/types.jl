# Project Structure:
#
# src/
#   Pioneer.jl               # Main module file
#   types.jl                # Type definitions
#   constants.jl            # Constants and configuration
#   fasta/                  # FASTA parsing and processing
#     fasta_types.jl
#     fasta_parser.jl
#     fasta_digest.jl
#   fragments/              # Fragment handling
#     fragment_types.jl  
#     fragment_parser.jl
#     fragment_index.jl
#   koina/                  # Koina API integration
#     koina_types.jl
#     koina_api.jl
#     koina_models.jl
#   utils/                  # Utility functions
#     io.jl
#     math.jl
#     modifications.jl
#   build/                  # Library building
#     build_types.jl
#     build_lib.jl
#     build_index.jl

# src/types.jl
# Core type definitions
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
Type alias for m/z to eV interpolation functions.
Uses GriddedInterpolation with linear interpolation and line extrapolation.
"""
const InterpolationTypeAlias = Interpolations.Extrapolation{
    Float32,  # Value type
    1,        # Dimension
    Interpolations.GriddedInterpolation{
        Float32,                            # Value type
        1,                                  # Dimension
        Vector{Float32},                    # Values
        Gridded{Linear{Throw{OnGrid}}},     # Method
        Tuple{Vector{Float32}}              # Grid type
    },
    Gridded{Linear{Throw{OnGrid}}},         # Method
    Line{Nothing}                           # Extrapolation
}
# src/constants.jl
const PREDICTION_MODELS = Set([
    "unispec",
    "prosit_2020_hcd", 
    "AlphaPeptDeep"
])

const MODEL_CONFIGS = Dict(
    "unispec" => (
        annotation_type = UniSpecFragAnnotation("y1^1"),
        model_type = InstrumentSpecificModel("unispec"),
        instruments = Set(["QE","QEHFX","LUMOS","ELITE","VELOS","NONE"])
    ),
    "prosit_2020_hcd" => (
        annotation_type = GenericFragAnnotation("y1+1"), 
        model_type = InstrumentAgnosticModel("prosit_2020_hcd"),
        instruments = Set([])
    ),
    "AlphaPeptDeep" => (
        annotation_type = GenericFragAnnotation("y1+1"),
        model_type = InstrumentSpecificModel("AlphaPeptDeep"),
        instruments = Set(["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"])
    )
)

const KOINA_URLS = Dict(
    "unispec" => "https://koina.wilhelmlab.org:443/v2/models/UniSpec/infer",
    "prosit_2020_hcd" => "https://koina.wilhelmlab.org:443/v2/models/Prosit_2020_intensity_HCD/infer",
    "AlphaPeptDeep" => "https://koina.wilhelmlab.org:443/v2/models/AlphaPeptDeep_ms2_generic/infer",
    "chronologer" => "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer"
)

