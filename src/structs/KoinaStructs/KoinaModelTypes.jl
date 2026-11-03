# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

abstract type KoinaModelType end

# Instrument-agnostic scalar-intensity model (Prosit). Returns a fixed vector of
# raw fragment intensities at a single collision energy — no NCE spline (contrast
# SplineCoefficientModel) and no instrument_type input. Re-added for Prosit support.
struct InstrumentAgnosticModel <: KoinaModelType
    name::String
end

struct SplineCoefficientModel <: KoinaModelType
    name::String
end

struct RetentionTimeModel <: KoinaModelType
    name::String
end

# Ion-mobility (CCS) model. Inputs are `peptide_sequences` + `precursor_charges`,
# the single output tensor is `ccs` (Å²). CCS is converted to 1/K0 downstream
# (chronologer_predict.jl) from the precursor m/z and charge.
struct IonMobilityModel <: KoinaModelType
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
