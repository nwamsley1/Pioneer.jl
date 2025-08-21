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
    f_index_fragments = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_fragments.arrow"))
    f_index_rt_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_rt_bins.arrow"))
    f_index_frag_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_fragment_bins.arrow"))


    presearch_f_index_fragments = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_fragments.arrow"))
    presearch_f_index_rt_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_rt_bins.arrow"))
    presearch_f_index_frag_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_fragment_bins.arrow"))


    # Note: Can't use @user_info here as LoggingSystem is loaded after ParseInputs
    # This message will be captured by the logging system when SearchDIA runs
    # println("Loading spectral library from $SPEC_LIB_DIR into main memory...")
    spec_lib = Dict{String, Any}()
    detailed_frags = load_detailed_frags(joinpath(SPEC_LIB_DIR,"detailed_fragments.jld2"))
    prec_frag_ranges = load(joinpath(SPEC_LIB_DIR,"precursor_to_fragment_indices.jld2"))["pid_to_fid"]
    library_fragment_lookup_table = nothing
    if (eltype(detailed_frags).name == (eltype(Vector{SplineDetailedFrag{4, Float32}}(undef, 0)).name))
        try
            #Model that encodes initial nce guess. 
            nmc = PiecewiseNceModel(Float32(params.acquisition[:nce]))

            spl_knots = load(joinpath(SPEC_LIB_DIR,"spline_knots.jld2"))["spl_knots"]
            library_fragment_lookup_table = SplineFragmentLookup(
                detailed_frags, 
                prec_frag_ranges, 
                Tuple(spl_knots), 
                Ref(nmc),
                3
            )
        catch e
            @warn "Could not load `spline_knots.jld2`"
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

    f_index = FragmentIndex(
        f_index_frag_bins[:FragIndexBin],
        f_index_rt_bins[:FragIndexBin],
        f_index_fragments[:IndexFragment],
    );
    presearch_f_index = FragmentIndex(
        presearch_f_index_frag_bins[:FragIndexBin],
        presearch_f_index_rt_bins[:FragIndexBin],
        presearch_f_index_fragments[:IndexFragment],
    );
    spec_lib["f_index"] = f_index;
    spec_lib["presearch_f_index"] = presearch_f_index;
    spec_lib["precursors"] = precursors;
    spec_lib["proteins"] = proteins;

    if typeof(library_fragment_lookup_table) == Pioneer.StandardFragmentLookup{Float32}
        return FragmentIndexLibrary(
            spec_lib["presearch_f_index"], 
            spec_lib["f_index"], 
            SetPrecursors(spec_lib["precursors"]), 
            SetProteins(spec_lib["proteins"]),
            spec_lib["f_det"]
        )
    else
        return SplineFragmentIndexLibrary(
            spec_lib["presearch_f_index"], 
            spec_lib["f_index"], 
            SetPrecursors(spec_lib["precursors"]), 
            SetProteins(spec_lib["proteins"]),
            spec_lib["f_det"]
        )
    end
end