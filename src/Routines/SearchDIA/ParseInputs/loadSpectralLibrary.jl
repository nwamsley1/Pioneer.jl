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


#=
function load_detailed_frags(filename::String)
    jldopen(filename, "r") do file
        data = read(file, "data")
        spline_type_name = eltype(Vector{SplineDetailedFrag{4, Float32}}(undef, 0)).name
        if eltype(data).name != spline_type_name
            return map(x -> DetailedFrag{Float32}(
                x.prec_id,
                x.mz,
                x.intensity,
                x.ion_type,
                x.is_y,
                x.is_b,
                x.is_p,
                x.is_isotope,
                x.frag_charge,
                x.ion_position,
                x.prec_charge,
                x.rank,
                x.sulfur_count
            ), data)
        else
            return map(x -> eltype(data)(
                x.prec_id,
                x.mz,
                x.intensity,
                x.ion_type,
                x.is_y,
                x.is_b,
                x.is_p,
                x.is_isotope,
                x.frag_charge,
                x.ion_position,
                x.prec_charge,
                x.rank,
                x.sulfur_count
            ), data)
        end
    end
end
=#


function loadSpectralLibrary(SPEC_LIB_DIR::String,
                             params::PioneerParameters)
    # Note: Can't use @user_info here as LoggingSystem is loaded after ParseInputs
    # This message will be captured by the logging system when SearchDIA runs
    spec_lib = Dict{String, Any}()

    # Load detailed fragments (load_detailed_frags handles backwards compatibility)
    detailed_frags = load_detailed_frags(joinpath(SPEC_LIB_DIR, "detailed_fragments.jls"))

    # Load precursor-to-fragment indices with backwards compatibility
    prec_frag_ranges = if isfile(joinpath(SPEC_LIB_DIR, "precursor_to_fragment_indices.jls"))
        deserialize_from_jls(joinpath(SPEC_LIB_DIR, "precursor_to_fragment_indices.jls"))
    elseif isfile(joinpath(SPEC_LIB_DIR, "precursor_to_fragment_indices.jld2"))
        @warn "Loading legacy JLD2 format for precursor_to_fragment_indices. Consider rebuilding library."
        load(joinpath(SPEC_LIB_DIR, "precursor_to_fragment_indices.jld2"))["pid_to_fid"]
    else
        error("precursor_to_fragment_indices file not found in $SPEC_LIB_DIR")
    end

    library_fragment_lookup_table = nothing
    if (eltype(detailed_frags).name == (eltype(Vector{SplineDetailedFrag{4, Float32}}(undef, 0)).name))
        try
            # Load spline knots with backwards compatibility
            spl_knots = if isfile(joinpath(SPEC_LIB_DIR, "spline_knots.jls"))
                deserialize_from_jls(joinpath(SPEC_LIB_DIR, "spline_knots.jls"))
            elseif isfile(joinpath(SPEC_LIB_DIR, "spline_knots.jld2"))
                @warn "Loading legacy JLD2 format for spline_knots. Consider rebuilding library."
                load(joinpath(SPEC_LIB_DIR, "spline_knots.jld2"))["spl_knots"]
            else
                error("spline_knots file not found in $SPEC_LIB_DIR")
            end
            library_fragment_lookup_table = SplineFragmentLookup(
                detailed_frags,
                prec_frag_ranges,
                Tuple(spl_knots),
                3
            )
        catch e
            @user_warn "Could not load spline_knots"
            throw(e)
        end

    else
        library_fragment_lookup_table = StandardFragmentLookup(detailed_frags, prec_frag_ranges)
    end
    #Is this still necessary?
    #last_range = library_fragment_lookup_table.prec_frag_ranges[end] #0x29004baf:(0x29004be8 - 1)
    #last_range = range(first(last_range), last(last_range) - 1)
    #library_fragment_lookup_table.prec_frag_ranges[end] = last_range
    spec_lib["f_det"] = library_fragment_lookup_table

    precursors = Arrow.Table(joinpath(SPEC_LIB_DIR, "precursors_table.arrow"))
    proteins = Arrow.Table(joinpath(SPEC_LIB_DIR, "proteins_table.arrow"))

    # Load partitioned fragment indexes
    partitioned_index = deserialize_from_jls(joinpath(SPEC_LIB_DIR, "partitioned_fragment_index.jls"))
    presearch_partitioned_index = deserialize_from_jls(joinpath(SPEC_LIB_DIR, "presearch_partitioned_fragment_index.jls"))

    if typeof(library_fragment_lookup_table) == Pioneer.StandardFragmentLookup{Float32}
        return FragmentIndexLibrary(
            presearch_partitioned_index,
            partitioned_index,
            SetPrecursors(precursors),
            SetProteins(proteins),
            spec_lib["f_det"]
        )
    else
        return SplineFragmentIndexLibrary(
            presearch_partitioned_index,
            partitioned_index,
            SetPrecursors(precursors),
            SetProteins(proteins),
            spec_lib["f_det"]
        )
    end
end