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
using JLD2

# Include all module files
include("types.jl")
include("constants.jl")
include("../../structs/Ion.jl")
include("../../structs/LibraryIon.jl")
include("../../structs/LibraryFragmentIndex.jl")
# FASTA processing
include("fasta/fasta_types.jl")
include("fasta/fasta_parser.jl") 
include("fasta/fasta_digest.jl")
include("fasta/fasta_utils.jl")
# Fragment handling
include("fragments/fragment_types.jl")
include("fragments/get_frag_bounds.jl")
include("fragments/fragment_parse.jl")
include("fragments/fragment_index.jl")
include("fragments/fragment_annotation.jl")
include("fragments/fragment_predict.jl")

# Koina integration
include("koina/koina_types.jl")
include("koina/koina_api.jl")
include("koina/koina_batch_prep.jl")
include("koina/koina_batch_parse.jl")

# Utilities
include("utils/io.jl")
include("utils/estimate_collision_ev.jl")
include("utils/math.jl")
include("utils/modifications.jl")
include("utils/get_mz.jl")
include("utils/parse_isotope_mods.jl")
include("utils/check_params.jl")

# Library building
include("build/build_types.jl")
include("build/build_lib.jl")
include("build/build_index.jl")
include("build/build_poin_lib.jl")


include("chronologer/chronologer_types.jl")
include("chronologer/chronologer_prep.jl")
include("chronologer/chronologer_predict.jl")
include("chronologer/chronologer_parse.jl")

end # module