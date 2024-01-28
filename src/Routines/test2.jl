frag_count = 0
for prec_frags in f_det
    for frag in prec_frags
        frag_count += 1
    end
end

detailed_frags = Vector{DetailedFrag{Float32}}(undef, frag_count)
prec_frag_idxs = Vector{UnitRange{UInt32}}(undef, length(precursors))

frag_idx = 0
prec_idx = 1
for prec_frags in ProgressBar(f_det)
    start = frag_idx + 1
    for frag in prec_frags
        frag_idx += 1
        detailed_frags[frag_idx] = frag
    end
    prec_frag_idxs[prec_idx] = range(UInt32(start), UInt32(frag_idx))
    prec_idx += 1
end

struct LibraryFragmentLookup{T<:AbstractFloat}
    frags::Vector{DetailedFrag{T}}
    prec_frag_ranges::Vector{UnitRange{UInt32}}
end

library_fragment_lookup_table = LibraryFragmentLookup(
    detailed_frags,
    prec_frag_idxs
)


#function Base.getindex(x::LibraryFragmentLookup{T}, k::S) where {T<:AbstractFloat,S<:Integer}
#    return x.summands[k]
#end

@save "C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\nOf3_y4b3_012724_sulfur\\HumanYeastEcoli_NCE33COR_012724_nOf3_indy4b3_ally3b2_library_fragment_lookup_table.jld2" library_fragment_lookup_table 
