module CompactLibraryStorage

using Arrow
import Pioneer

const P = Pioneer

export encode_index, decode_index, write_index, read_index
export encode_details, decode_details, write_details, read_details
export write_precursors, read_precursors, koina_sequence

require(condition, message) = condition || throw(ArgumentError(message))

"""Smallest unsigned integer type that preserves every value in a column."""
function unsigned_type(maximum_value::Unsigned)
    maximum_value <= typemax(UInt8) && return UInt8
    maximum_value <= typemax(UInt16) && return UInt16
    maximum_value <= typemax(UInt32) && return UInt32
    return UInt64
end

include("index.jl")
include("details.jl")
include("precursors.jl")

end
