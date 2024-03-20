#Change frag index precursor bins from vector of vectors to single vector. 


function fillTestFrags!(test_frags::Vector{DetailedFrag{Float32}},
                        arrow_frags::Arrow.Struct{DetailedFrag, Tuple{Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{Float16, Vector{Float16}}, Arrow.BoolVector{Bool}, Arrow.BoolVector{Bool}, Vararg{Arrow.Primitive{UInt8, Vector{UInt8}}, 5}}, (:prec_id, :mz, :intensity, :is_y_ion, :is_isotope, :frag_charge, :ion_position, :prec_charge, :rank, :sulfur_count)},
                        N::Int64
    )
    n = 1
    i = 1000000
    for n in range(1, 10000)
        test_frags[n] = getindex(arrow_frags, i)
        i += 1
    end
end
N = 200000
test_frags = Vector{DetailedFrag{Float32}}(undef, N)
@time fillTestFrags!(test_frags, library_fragment_lookup_table.frags, N)


function fillTestFragsRef!(test_frags_ref::Vector{Base.RefValue{DetailedFrag{Float32}}},
                        arrow_frags::Arrow.Struct{DetailedFrag, Tuple{Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{Float16, Vector{Float16}}, Arrow.BoolVector{Bool}, Arrow.BoolVector{Bool}, Vararg{Arrow.Primitive{UInt8, Vector{UInt8}}, 5}}, (:prec_id, :mz, :intensity, :is_y_ion, :is_isotope, :frag_charge, :ion_position, :prec_charge, :rank, :sulfur_count)},
                        N::Int64
    )
    n = 1
    i = 1000000
    for n in range(1, N)
        test_frags_ref[n] = Ref(arrow_frags[i])
        i += 1
    end
end
N = 100000
test_frags_ref = Vector{Base.RefValue{DetailedFrag{Float32}}}(undef, N)
@time fillTestFragsRef!(test_frags_ref, library_fragment_lookup_table.frags, N)


function fillTestFrags!(test_frags::Vector{DetailedFrag{Float32}},
                        ref_frags::Vector{DetailedFrag{Float32}},
                        N::Int64
    )
    n = 1
    i = 1
    for n in range(1, 10000)
        test_frags[n] = ref_frags[i]
        i += 1
    end
end
N = 100000
test_frags = Vector{DetailedFrag{Float32}}(undef, N)
test_frags_ref = Vector{DetailedFrag{Float32}}(undef, N)
fillTestFrags!(test_frags_ref, library_fragment_lookup_table.frags, N)
@time fillTestFrags!(test_frags, test_frags_ref, N)

using Profile, PProf
Profile.clear()
MS_TABLE = Arrow.Table(MS_TABLE_PATH)  
@profile PSMs = vcat(mainLibrarySearch(
    MS_TABLE,
    prosit_lib["f_index"],
    prosit_lib["precursors"],
    library_fragment_lookup_table,
    RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
    UInt32(ms_file_idx), #MS_FILE_IDX
    frag_err_dist_dict[ms_file_idx],
    irt_errs[ms_file_idx],
    params_,
    ionMatches,
    ionMisses,
    all_fmatches,
    IDtoCOL,
    ionTemplates,
    iso_splines,
    scored_PSMs,
    unscored_PSMs,
    spectral_scores,
    precursor_weights,
    precs
#scan_range = (100000, 100010)
)...);
pprof(;webport=58600)
#pprof(;webport=58599)