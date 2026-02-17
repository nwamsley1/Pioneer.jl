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

function importSpecLibScripts()
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
    include("build/partition_library.jl")


    include("chronologer/chronologer_types.jl")
    include("chronologer/pair_decoys.jl")
    include("chronologer/chronologer_prep.jl")
    include("chronologer/chronologer_predict.jl")
    include("chronologer/chronologer_parse.jl")
end