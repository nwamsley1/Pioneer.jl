const NativeFragment = Union{
    P.CompactFrag{Float32}, P.SplineCompactFrag{N,Float32} where N
}

"""Disk-only columns. Intensities retain their existing scalar or tuple type."""
struct StoredDetails{C<:Unsigned,I}
    counts::Vector{C}
    mz::Vector{Float32}
    intensity::Vector{I}
    tag_low::Vector{UInt16}
    tag_high::Vector{UInt8}
end

"""
    encode_details(fragments, offsets, charges) -> StoredDetails

Replace UInt64 CSR offsets with fragment counts. Precursor IDs and charges are
implicit in each CSR range; all remaining annotation bits and prediction bits
are preserved. No fragment dictionaries or runtime representations are added.
"""
function encode_details(fragments::AbstractVector{F}, offsets, charges) where {F<:NativeFragment}
    n = length(fragments)
    require(length(offsets) == length(charges) + 1, "One CSR boundary is required per precursor plus the final boundary")
    require(first(offsets) == 1 && last(offsets) == n + 1 && issorted(offsets),
        "CSR boundaries must cover the complete fragment vector in order")
    require(all(charge -> 0 <= charge <= 7, charges), "Precursor charge exceeds its three-bit field")
    counts = UInt64.(diff(offsets))
    C = unsigned_type(maximum(counts; init = UInt64(0)))
    mz = Vector{Float32}(undef, n)
    intensity = Vector{fieldtype(F, :intensity)}(undef, n)
    tag_low = Vector{UInt16}(undef, n)
    tag_high = Vector{UInt8}(undef, n)
    for pid in eachindex(charges), i in offsets[pid]:offsets[pid + 1] - 1
        fragment = fragments[i]
        require(fragment.prec_id == pid && P.getPrecCharge(fragment) == charges[pid],
            "Fragment owner and charge must agree with the precursor table and CSR")
        require(fragment.packed_a & 0xf000 == 0 && fragment.packed_b & 0xf800 == 0,
            "Unrecognized fragment annotation bits")
        # Bits 0:5 = ion flags and fragment charge; 6:8 = sulfur;
        # 9:14 = ion position; 15:19 = intensity rank.
        tag = UInt32(fragment.packed_a & 0x3f) |
            (UInt32((fragment.packed_a >> 9) & 7) << 6) |
            (UInt32(fragment.packed_b) << 9)
        mz[i] = fragment.mz
        intensity[i] = fragment.intensity
        tag_low[i] = UInt16(tag & 0xffff)
        tag_high[i] = UInt8(tag >> 16)
    end
    return StoredDetails(C.(counts), mz, intensity, tag_low, tag_high)
end

fragment_type(::Type{Float16}) = P.CompactFrag{Float32}
fragment_type(::Type{NTuple{N,Float32}}) where {N} = P.SplineCompactFrag{N,Float32}

"""Decode both scalar and spline columns directly into the original native vector."""
function decode_details(stored::StoredDetails{C,I}, charges) where {C,I}
    n = length(stored.mz)
    require(length(stored.counts) == length(charges), "Fragment counts must match precursor charges")
    require(length(stored.intensity) == length(stored.tag_low) == length(stored.tag_high) == n,
        "Detailed-fragment columns must have equal lengths")
    require(all(tag -> tag & 0xf0 == 0, stored.tag_high), "Unrecognized stored annotation bits")
    offsets = Vector{UInt64}(undef, length(charges) + 1)
    offsets[1] = 1
    for pid in eachindex(charges)
        require(0 <= charges[pid] <= 7, "Precursor charge exceeds its three-bit field")
        offsets[pid + 1] = Base.Checked.checked_add(offsets[pid], UInt64(stored.counts[pid]))
    end
    require(last(offsets) == n + 1, "Stored fragment counts must cover the complete vector")
    F = fragment_type(I)
    fragments = Vector{F}(undef, n)
    for pid in eachindex(charges), i in offsets[pid]:offsets[pid + 1] - 1
        tag = UInt32(stored.tag_low[i]) | (UInt32(stored.tag_high[i]) << 16)
        a = UInt16(tag & 0x3f) | (UInt16(charges[pid]) << 6) |
            (UInt16((tag >> 6) & 7) << 9)
        b = UInt16(tag >> 9)
        fragments[i] = F(UInt32(pid), stored.mz[i], stored.intensity[i], a, b)
    end
    return fragments, offsets
end

write_details(path, fragments, offsets, charges) =
    P.serialize_to_jls(path, encode_details(fragments, offsets, charges))
read_details(path, charges) = decode_details(P.deserialize_from_jls(path), charges)
