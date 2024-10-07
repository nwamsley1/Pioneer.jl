function loadSpectralLibrary(SPEC_LIB_DIR::String)
    f_index_fragments = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_fragments.arrow"))
    f_index_rt_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_rt_bins.arrow"))
    f_index_frag_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_fragment_bins.arrow"))


    presearch_f_index_fragments = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_fragments.arrow"))
    presearch_f_index_rt_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_rt_bins.arrow"))
    presearch_f_index_frag_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_fragment_bins.arrow"))


    println("Loading spectral libraries into main memory...")
    spec_lib = Dict{String, Any}()
    detailed_frags = load(joinpath(SPEC_LIB_DIR,"detailed_fragments.jld2"))["data"]
    prec_frag_ranges = load(joinpath(SPEC_LIB_DIR,"precursor_to_fragment_indices.jld2"))["pid_to_fid"]
    library_fragment_lookup_table = LibraryFragmentLookup(detailed_frags, prec_frag_ranges)
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