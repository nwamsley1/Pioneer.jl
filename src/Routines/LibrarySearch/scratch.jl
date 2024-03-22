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

eight_bit_ints = rand(UInt8, 10000)
function isZeroSimple()
    @inbounds @fastmath begin
        for 
    end
end

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
pprof(;webport=58701)
#pprof(;webport=58599)

FragmentIndex

describe(precursors[:mz])

minimum([x for x in MS_TABLE[:precursorMZ] if ismissing(x)==false])
maximum([x for x in MS_TABLE[:precursorMZ] if ismissing(x)==false])


function getPrecursorBins(
    min_prec_mz::T,
    max_prec_mz::T,
    bin_width::T,
    fragment_bins::Arrow.Struct{FragIndexBin, 
        Tuple{Arrow.Primitive{T, Vector{T}}, 
        Arrow.Primitive{T, Vector{T}}, 
        Arrow.Primitive{UInt32, Vector{UInt32}}, 
        Arrow.Primitive{UInt32, Vector{UInt32}}}, (:lb, :ub, :first_bin, :last_bin)},
    fragments:: Arrow.Struct{IndexFragment, 
        Tuple{Arrow.Primitive{UInt32, Vector{UInt32}}, 
        Arrow.Primitive{T, Vector{T}}, 
        Arrow.Primitive{UInt8, Vector{UInt8}}, 
        Arrow.Primitive{UInt8, Vector{UInt8}}}, (:prec_id, :prec_mz, :score, :charge)}
    ) where {T<:AbstractFloat}
    frag_bin_mzs = collect(range(
                            min_prec_mz, 
                            max_prec_mz, 
                            step = bin_width)
                            )
    M = length(frag_bin_mzs)
    N = length(fragment_bins)
    precursor_bins = zeros(UInt32, M*N)
    i = 1
    for (n, frag_bin) in ProgressBar(enumerate(fragment_bins))
        frag_range = getSubBinRange(frag_bin)
        idx = first(frag_range)
        for (m, prec_mz) in enumerate(frag_bin_mzs)
            while idx <= last(frag_range)
                if getPrecMZ(fragments[idx]) > prec_mz
                    precursor_bins[i] = idx
                    break
                end
                idx += 1
            end
            if idx > last(frag_range)
                precursor_bins[(n - 1)*M + m] = last(frag_range)
            end
            i += 1
        end
    end
    return precursor_bins, frag_bin_mzs
end

precursor_bins, frag_bin_mzs = getPrecursorBins(396.0f0, 
                1001.0f0,
                0.5f0, 
f_index_frag_bins[:FragIndexBin],
f_index_fragments[:IndexFragment]
)

precursor_bins  = (PrecursorBin = precursor_bins ,)
frag_bin_mzs = (FragBinMZ = frag_bin_mzs ,)
Arrow.write("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/f_index_top5_7ppm_1hi_rtmajor_precursor_bins_031924.arrow", precursor_bins);
Arrow.write("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/f_index_top5_7ppm_1hi_rtmajor_frag_bin_mzs_031924.arrow", frag_bin_mzs);


#Come up with some test cases. Keep trying till they work. 
N = 1005000
M = length(frag_bin_mzs)
test_f_bin = f_index_frag_bins[:FragIndexBin][N]
f_index_fragments[:IndexFragment][getSubBinRange(test_f_bin)]



for i in range(1, length(frag_bin_mzs))
    idx = precursor_bins[(N-1)*M + i]
    if idx!=0
        println(frag_bin_mzs[i], "-", idx, "-", getPrecMZ(f_index_fragments[:IndexFragment][idx]))
    else
        println(frag_bin_mzs[i], "-", idx, "-", "null")
    end
end

hcat(frag_bin_mzs, precursor_bins[(N-1)*M:(N-1)*M + M - 1])
@time PSMs = vcat(mainLibrarySearch(
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

function branchless_binary(t::Vector{Float32},
                           x::Float32,
                           y::Float32,
                           lo::UInt32,
                           hi::UInt32)
    #hi_f = hi
    base = lo
    @inbounds @fastmath begin
        len = hi - lo + UInt32(1)

        while len > 1
            mid = len>>>0x01
            base += (t[base + mid - UInt32(1)] < x)*mid
            len -= mid# - UInt32(1)
        end
        window_start = base
        len = hi - base + UInt32(1)
        println("len $len")
        base = hi
        while len > 1
            mid = len>>>0x01
            base -= (t[base - mid + UInt32(1)] > y)*mid
            #base = (t[base - mid + UInt32(1)] > y) ? base - mid : base
            len -= mid# - UInt32(1)
        end
        window_stop = base

        if window_start === window_stop
            if t[window_start]>y
                return one(UInt32), zero(UInt32)
            end
            if t[window_stop]<x
                return one(UInt32), zero(UInt32)
            end
        end
    end
    return window_start, window_stop
end

x = 4.1f0
y = 4.5f0
test_t = 100.0f0.*sort(rand(Float32, 100000))   
test_t = test_t[(test_t .< 4.0) .| (test_t .> 6.0)]
#@btime answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(100000))
answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(100000))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]

test_t = 100.0f0.*sort(rand(Float32, 100000))   
#@btime answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(100000))
answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(100000))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]
test_t[first(answer)-1:last(answer)+1]

x = 3.5f0
y = 5.0f0
test_t = 100.0f0.*sort(rand(Float32, 100))   
#@btime answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(100))
answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(100))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]
test_t[first(answer)-1:last(answer)+1]

test_t = [20.0f0, 100.0f0, 200.0f0]
answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(3))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]

test_t = [1.0f0, 20.0f0, 100.0f0, 200.0f0, 200.0f0, 300.0f0]
answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(5))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]

test_t = [1.0f0, 4.2f0, 20.0f0, 100.0f0, 200.0f0, 200.0f0, 300.0f0]
answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(2))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]

test_t = [1.0f0, 4.2f0, 20.0f0, 100.0f0, 200.0f0, 200.0f0, 300.0f0]
answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(5))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]


test_t = [1.0f0, 4.2f0, 20.0f0, 100.0f0, 200.0f0, 200.0f0, 300.0f0]
answer = branchless_binary(test_t, x, y, UInt32(3), UInt32(5))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]


test_t = [1.0f0, 20.0f0, 100.0f0, 200.0f0, 200.0f0, 300.0f0]
answer = branchless_binary(test_t, x, y, UInt32(6), UInt32(6))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]


function branchless_binary_a(t::Vector{Float32},
                           x::Float32,
                           lo::UInt32,
                           hi::UInt32)
    #hi_f = hi
    base = lo
    @fastmath len = hi - lo
    while len > 1
        mid = len>>>0x01
        base += (t[base + mid - UInt32(1)] < x)*mid
        len -= mid# - UInt32(1)
    end
    #lo_f = lo
    return base
end
x = 3.5f0
y = 5.0f0
test_t = 100.0f0.*sort(rand(Float32, 100))   
#@btime answer = branchless_binary(test_t, x, y, UInt32(1), UInt32(100))
answer = branchless_binary_a(test_t, x, UInt32(1), UInt32(100))
println(Int64(answer), " ", Int64(answer))
test_t[answer-1:answer+1]
#test_t[first(answer)-1:last(answer)+1]
test_t = [3.6f0, 4.2f0, 20.0f0, 100.0f0, 200.0f0, 200.0f0, 300.0f0]
answer = branchless_binary_a(test_t, x, UInt32(1), UInt32(2))
println(Int64(answer), " ", Int64(answer))
test_t[answer-1:answer+1]


test_t = [1.0f0, 4.2f0, 20.0f0, 100.0f0, 200.0f0, 200.0f0, 300.0f0]
answer = branchless_binary_a(test_t, x, UInt32(1), UInt32(2))
println(Int64(answer), " ", Int64(answer))
test_t[answer-1:answer+1]

function branchless_binary(t::Vector{Float32},
                           x::Float32,
                           y::Float32,
                           lo::UInt32,
                           hi::UInt32)
    #hi_f = hi
    base = lo
    @fastmath len = hi - lo
    while len > 1
        mid = len>>>0x01
        base += (t[base + mid - UInt32(1)] < x)*mid
        len -= mid# - UInt32(1)
    end
    a = base
    @fastmath len = hi - lo
    base = hi
    while len > 1
        mid = len>>>0x01
        base -= (t[base - mid + UInt32(1)] > y)*mid
        len -= mid# - UInt32(1)
    end

    #lo_f = lo
    return a, base
end

x = 4.1f0
y = 4.5f0
test_t = 10.0f0.*sort(rand(Float32, 1000000))   
#test_t = test_t[(test_t .< 4.0) .| (test_t .> 6.0)]
answer = branchless_binary(test_t, x, y, UInt32(1000), UInt32(900000))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]


x = 7.1f0
y = 8.1f0
answer = branchless_binary(test_t, x, y, UInt32(1000), UInt32(900000))
println(Int64(answer[1]), " ", Int64(answer[2]))
test_t[first(answer):last(answer)]


for i in 10:1
    println(i)
end



function FindPrecBinRange(t::Arrow.Primitive{Float32, Vector{Float32}},
                           x::Float32,
                           y::Float32)
    #hi_f = hi
    base, hi = one(UInt32), UInt32(length(t))
    @fastmath len = hi - base
    while len > 1
        mid = len>>>0x01
        base += (t[base + mid - UInt32(1)] < x)*mid
        len -= mid# - UInt32(1)
    end

    #=
    a = base
    @fastmath len = hi - lo
    base = hi
    while len > 1
        mid = len>>>0x01
        base -= (t[base - mid + UInt32(1)] > y)*mid
        len -= mid# - UInt32(1)
    end
    =#
    a = base
    @fastmath len = hi - base
    while len > 1
        mid = len>>>0x01
        base += (t[base + mid - UInt32(1)] < y)*mid
        len -= mid# - UInt32(1)
    end

    #lo_f = lo
    return a, base
end

x = 800.9f0
y = 809.4f0

@btime answer = branchless_binary(f_index.frag_bin_mzs, x, y)
println(Int64(answer[1]), " ", Int64(answer[2]))
f_index.frag_bin_mzs[first(answer):last(answer)]
