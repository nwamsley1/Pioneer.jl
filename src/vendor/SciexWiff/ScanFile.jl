# Decoding one `.wiff.scan` block. Layout and evidence: notes/format.md §3–4.
#
#   [varint len][protobuf: 1 = {1: f64 a, 2: f64 b}, 2* = {1: start bin, 2: seek offset}]
#   [ff ff ff ff][u32 start bin][00][tokens ...][0xff padding to a multiple of 4]
#
# Token: delta prefix (none = 1, 80–fb = b − 0x7f, fc u8, fd u16, fe u32; the escapes store Δ − 1)
# then intensity (00–7b literal, 7c u8, 7d u16, 7e u32). bin = start + 8·Σ Δ (first Δ included).
# m/z = (a · (bin / 40 − b))².

struct ScanFormatError <: Exception
    msg::String
end
Base.showerror(io::IO, e::ScanFormatError) = print(io, "ScanFormatError: ", e.msg)

"Per-thread scratch for one decoded block."
mutable struct ScanBuffer
    n::Int
    bin::Vector{UInt32}
    intensity::Vector{UInt32}
    cal_a::Float64
    cal_b::Float64
end
ScanBuffer() = ScanBuffer(0, UInt32[], UInt32[], NaN, NaN)

@inline bin_to_mz(bin, a::Float64, b::Float64) = (v = a * (bin / 40 - b); v * v)
@inline mz_to_bin(mz, a::Float64, b::Float64) = 40 * (sqrt(mz) / a + b)
@inline mz(buf::ScanBuffer, k::Integer) = bin_to_mz(buf.bin[k], buf.cal_a, buf.cal_b)

@inline function _varint(s::AbstractVector{UInt8}, i::Int)
    v = UInt64(0); sh = 0
    @inbounds while true
        x = s[i+1]; v |= UInt64(x & 0x7f) << sh; i += 1; sh += 7
        x < 0x80 && return v, i
        sh > 63 && throw(ScanFormatError("varint too long at byte $i"))
    end
end

"Read calibration (a, b) from the block metadata occupying [meta_start, ff)."
function read_calibration(s::AbstractVector{UInt8}, meta_start::Int64, ff::Int64)
    len, i = _varint(s, Int(meta_start))
    i + len == ff || throw(ScanFormatError("metadata at $meta_start does not end at block $ff"))
    @inbounds (s[i+1] == 0x0a && s[i+2] == 0x12 && s[i+3] == 0x09 && s[i+12] == 0x11) ||
        throw(ScanFormatError("unexpected calibration message at $i"))
    _rd(Float64, s, i + 3), _rd(Float64, s, i + 12)
end

"""
    decode_block!(buf, s, meta_start, ff, size) -> buf

Decode the block whose `ffffffff` marker is at 0-based offset `ff` of the `.wiff.scan` bytes `s`.
"""
function decode_block!(buf::ScanBuffer, s::AbstractVector{UInt8}, meta_start::Int64, ff::Int64, size::Integer)
    stop = Int(ff + size)
    stop <= length(s) || throw(ScanFormatError("block at $ff runs past end of file"))
    buf.cal_a, buf.cal_b = read_calibration(s, meta_start, ff)
    @inbounds begin
        (s[ff+1] == 0xff && s[ff+2] == 0xff && s[ff+3] == 0xff && s[ff+4] == 0xff && s[ff+9] == 0x00) ||
            throw(ScanFormatError("no block marker at $ff"))
        start = _rd(UInt32, s, ff + 4)
        while stop > ff + 9 && s[stop] == 0xff
            stop -= 1
        end
    end
    # every token is at least one byte, so this bounds the peak count
    cap = stop - Int(ff) - 9
    length(buf.bin) < cap && (resize!(buf.bin, cap); resize!(buf.intensity, cap))
    bins = buf.bin; ints = buf.intensity
    i = Int(ff) + 9; n = 0; pos = UInt64(start)
    @inbounds while i < stop
        c = s[i+1]
        if c < 0x80
            d = UInt64(1)
        elseif c <= 0xfb
            d = UInt64(c - 0x7f); i += 1
        elseif c == 0xfc
            d = UInt64(s[i+2]) + 1; i += 2
        elseif c == 0xfd
            d = UInt64(s[i+2]) | UInt64(s[i+3]) << 8 + 1; i += 3
        elseif c == 0xfe
            d = (UInt64(s[i+2]) | UInt64(s[i+3]) << 8 | UInt64(s[i+4]) << 16 | UInt64(s[i+5]) << 24) + 1; i += 5
        else
            throw(ScanFormatError("0xff delta prefix inside block at $ff (byte $i)"))
        end
        i < stop || throw(ScanFormatError("block at $ff ends inside a token"))
        v = s[i+1]
        if v <= 0x7b
            iv = UInt32(v); i += 1
        elseif v == 0x7c
            iv = UInt32(s[i+2]); i += 2
        elseif v == 0x7d
            iv = UInt32(s[i+2]) | UInt32(s[i+3]) << 8; i += 3
        elseif v == 0x7e
            iv = UInt32(s[i+2]) | UInt32(s[i+3]) << 8 | UInt32(s[i+4]) << 16 | UInt32(s[i+5]) << 24; i += 5
        else
            throw(ScanFormatError("bad intensity prefix 0x$(string(v, base = 16)) in block at $ff"))
        end
        i <= stop || throw(ScanFormatError("block at $ff ends inside a token"))
        pos += 8d
        n += 1
        bins[n] = pos % UInt32; ints[n] = iv
    end
    buf.n = n
    buf
end
