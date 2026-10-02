# Pure-Julia reader for Thermo .raw files, reverse-engineered from data files only (format version 66, Orbitrap
# Exploris; verified against RawFileReader-converted .arrow). No Thermo code is used.
#
# Layout (little-endian):
#   FileHeader   UInt16 0xA101, "Finnigan" (UTF-16), ..., UInt32 version at 0x24 (66). After variable-length file info:
#                UInt64 data-stream address, then per controller (Int32 type, 0 = MS; Int32 index; UInt64 RunHeader
#                address). Located here through the RunHeader's self-pointer instead (see find_ms_run_header).
#   RunHeader    +8 Int32 first scan, +12 Int32 last scan; UInt64 addresses at +7408 scan index, +7416 data stream,
#                +7424 status log, +7432 (log), +7448 scan events, +7456 trailer, +7472 RunHeader itself.
#   Scan index   88 bytes per scan: +4 Int32 index (0-based), +8 UInt16 scan event / UInt16 segment, +12 Int32 next,
#                +16 Int32 packet type, +20 UInt32 packet size, +24 Float64 x6 (RT min, TIC, base-peak intensity,
#                base-peak m/z, low m/z, high m/z), +72 UInt64 packet offset into the data stream.
#   Packet       40-byte header: UInt32 x8 (?, profile words, peak-list words, layout, descriptor words, two more
#                stream sizes, one more size) + Float32 low/high m/z; then profile, then the peak list
#                (UInt32 count + count x (Float32 m/z, Float32 intensity)), then one UInt32 descriptor per peak
#                (flags in bits 16-23; 0x10 = reference / lock-mass peak, which RawFileReader-based converters drop).
#   Scan events  UInt32, then one record per scan: 136-byte preamble (byte 6 = MS order), UInt32 n reactions,
#                56 bytes each (Float64 precursor/isolation centre, Float64 isolation width, Float64 collision energy,
#                UInt32 x2, Float64 isolation low/high, 8 bytes), UInt32 n mass ranges x (Float64 low, high),
#                UInt32 n coefficients x Float64 (5 on Exploris, 7 on Eclipse Orbitrap scans, 0 for ion trap), 12 bytes.
module ThermoRaw
using Mmap

struct ScanInfo
    rt::Float64; tic::Float64; base_peak_intensity::Float64; base_peak_mz::Float64
    low_mz::Float64; high_mz::Float64; packet_type::Int32; packet_offset::Int64
    ms_order::UInt8; center_mz::Float64; isolation_width::Float64; collision_energy::Float64
end

struct RawFile
    bytes::Vector{UInt8}
    version::UInt32
    first_scan::Int32
    last_scan::Int32
    data_addr::Int64
    scans::Vector{ScanInfo}
end

@inline rd(b, ::Type{T}, o) where {T} = ltoh(unsafe_load(Ptr{T}(pointer(b, o + 1))))

"""
Find the mass-spectrometer RunHeader. Each acquisition device (MS, pump, UV, analog) has one; a RunHeader stores its own
address at +7472. The MS one is the one with scans and a scan-event stream. (The controller list that points at them
sits after variable-length file info, so it is located by this self-pointer rather than at a fixed offset.)
"""
function find_ms_run_header(b::Vector{UInt8})
    best = -1; best_n = -1
    GC.@preserve b for o in 7472:length(b)-8
        rd(b, UInt64, o) == o - 7472 || continue
        rh = o - 7472
        n = rd(b, Int32, rh + 12) - rd(b, Int32, rh + 8) + 1
        (rd(b, UInt64, rh + 7448) != 0 && n > best_n) && (best = rh; best_n = n)
    end
    best >= 0 || error("no mass-spectrometer RunHeader found")
    best
end

function open_raw(path::AbstractString)
    b = Mmap.mmap(path)
    rd(b, UInt16, 0) == 0xa101 || error("$path is not a Thermo .raw file")
    version = rd(b, UInt32, 0x24)
    version == 66 || @warn "untested .raw format version $version"
    rh = find_ms_run_header(b)
    first_scan, last_scan = rd(b, Int32, rh + 8), rd(b, Int32, rh + 12)
    index_addr, data_addr, events_addr = (Int(rd(b, UInt64, rh + o)) for o in (7408, 7416, 7448))
    n = last_scan - first_scan + 1
    scans = Vector{ScanInfo}(undef, n)
    r = events_addr + 4
    for k in 0:n-1
        e = index_addr + 88k
        rd(b, Int32, e + 4) == k || error("scan index entry $k is out of order")
        order = b[r+7]
        nreact = Int(rd(b, UInt32, r + 136))
        center, width, ce = nreact > 0 ? (rd(b, Float64, r + 140), rd(b, Float64, r + 148), rd(b, Float64, r + 156)) :
                                         (NaN, NaN, NaN)
        q = r + 140 + 56nreact
        nrange = Int(rd(b, UInt32, q)); q += 4 + 16nrange
        ncoef = Int(rd(b, UInt32, q)); q += 4 + 8ncoef + 12
        scans[k+1] = ScanInfo(rd(b, Float64, e + 24), rd(b, Float64, e + 32), rd(b, Float64, e + 40),
                              rd(b, Float64, e + 48), rd(b, Float64, e + 56), rd(b, Float64, e + 64),
                              rd(b, Int32, e + 16), data_addr + Int(rd(b, UInt64, e + 72)),
                              order, center, width, ce)
        r = q
    end
    RawFile(b, version, first_scan, last_scan, data_addr, scans)
end

"Centroids of scan `i` (1-based) as (m/z, intensity) Float32 vectors; reference (lock-mass) peaks are dropped unless `keep_reference`."
function centroids(f::RawFile, i::Integer; keep_reference::Bool = false)
    # Orbitrap (type 21: profile + centroids) and Astral (type 20: centroids only, empty profile) packets carry this
    # centroid list. Ion-trap packets (type 18) are laid out differently; RawFileReader-based converters write no
    # peaks for them either.
    f.scans[i].packet_type in (20, 21) || return Float32[], Float32[]
    b = f.bytes; p = f.scans[i].packet_offset
    prof, pkw = Int(rd(b, UInt32, p + 4)), Int(rd(b, UInt32, p + 8))
    q = p + 40 + 4prof
    n = pkw == 0 ? 0 : Int(rd(b, UInt32, q))
    d = q + 4pkw
    mz = Float32[]; it = Float32[]
    for j in 0:n-1
        (keep_reference || (rd(b, UInt32, d + 4j) >> 16) & 0x10 == 0) || continue
        push!(mz, rd(b, Float32, q + 4 + 8j)); push!(it, rd(b, Float32, q + 8 + 8j))
    end
    mz, it
end

end
