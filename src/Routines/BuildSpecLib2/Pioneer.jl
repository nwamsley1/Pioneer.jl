module Pioneer

using Arrow
using CSV
using CodecZlib
using Combinatorics
using DataFrames
using Dictionaries
using FASTX
using Interpolations
using JSON
using ProgressBars
using Polynomials
using Plots
# Type exports
export KoinaModelType, InstrumentSpecificModel, InstrumentAgnosticModel
export SplineCoefficientModel, RetentionTimeModel
export FragAnnotation, UniSpecFragAnnotation, GenericFragAnnotation
export FastaEntry, PioneerFrag, PioneerSplineFrag

# Function exports
export parse_fasta, digest_fasta, build_library
export prepare_koina_batch, make_koina_request
export build_fragment_index, parse_fragment_annotation

# Include all module files
include("types.jl")
include("constants.jl")

# FASTA processing
include("fasta/fasta_types.jl")
include("fasta/fasta_parser.jl") 
include("fasta/fasta_digest.jl")

# Fragment handling
include("fragments/fragment_types.jl")
include("fragments/get_frag_bounds.jl")
include("fragments/fragment_parser.jl")
include("fragments/fragment_index.jl")

# Koina integration
include("koina/koina_types.jl")
include("koina/koina_api.jl")
include("koina/koina_models.jl")

# Utilities
include("utils/io.jl")
include("utils/estimate_collision_ev.jl")
include("utils/math.jl")
include("utils/modifications.jl")

# Library building
include("build/build_types.jl")
include("build/build_lib.jl")
include("build/build_index.jl")


include("chronologer/chronologer_types.jl")
include("chronologer/chronologer_prep.jl")
include("chronologer/chronologer_predict.jl")
include("chronologer/chronologer_parse.jl")

end # module