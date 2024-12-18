
#=
function load_detailed_frags(filename::String)
    jldopen(filename, "r") do file
        data = read(file, "data")
        spline_type_name = eltype(Vector{SplineDetailedFrag{4, Float32}}(undef, 0)).name
        println("a ", eltype(data).name)
        println("b ", spline_type_name)
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
                             params::Any)
    f_index_fragments = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_fragments.arrow"))
    f_index_rt_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_rt_bins.arrow"))
    f_index_frag_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_fragment_bins.arrow"))


    presearch_f_index_fragments = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_fragments.arrow"))
    presearch_f_index_rt_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_rt_bins.arrow"))
    presearch_f_index_frag_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_fragment_bins.arrow"))


    println("Loading spectral libraries into main memory...")
    spec_lib = Dict{String, Any}()
    detailed_frags = load_detailed_frags(joinpath(SPEC_LIB_DIR,"detailed_fragments.jld2"))
    prec_frag_ranges = load(joinpath(SPEC_LIB_DIR,"precursor_to_fragment_indices.jld2"))["pid_to_fid"]
    library_fragment_lookup_table = nothing
    if (eltype(detailed_frags).name == (eltype(Vector{SplineDetailedFrag{4, Float32}}(undef, 0)).name))
        try
            #Model that encodes initial nce guess. 
            nmc = PiecewiseNceModel(Float32(params[:presearch_params]["nce_guess"]))

            spl_knots = load(joinpath(SPEC_LIB_DIR,"spline_knots.jld2"))["spl_knots"]
            library_fragment_lookup_table = SplineFragmentLookup(
                detailed_frags, 
                prec_frag_ranges, 
                spl_knots, 
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

    precursors = Arrow.Table(joinpath(SPEC_LIB_DIR, "precursors_table.arrow"))#DataFrame(precursors)
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
    return spec_lib
end