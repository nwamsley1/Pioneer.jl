# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

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
