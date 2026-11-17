# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Little-endian reads at 0-based offset k. Bounds are always checked: a wrong layout guess must fail loudly,
# never read outside the mapped file.
@inline function _rd(::Type{T}, b, k) where T
    checkbounds(b, k + 1:k + sizeof(T))
    GC.@preserve b ltoh(unsafe_load(Ptr{T}(pointer(b, k + 1))))
end
@inline ru8(b, k)  = b[k+1]
@inline ru16(b, k) = _rd(UInt16, b, k)
@inline ru32(b, k) = _rd(UInt32, b, k)
@inline ru64(b, k) = _rd(UInt64, b, k)
@inline rf32(b, k) = reinterpret(Float32, ru32(b, k))
@inline rf64(b, k) = reinterpret(Float64, ru64(b, k))

"""Per-scan description, from the scan table and region C."""
struct ScanInfo
    scan_number::Int32
    ms_order::UInt8
    analyser::UInt8              # region C fixed part +40: 0 ion trap, 4 orbitrap, 7 astral
    retention_time::Float64      # minutes
    low_mz::Float64
    high_mz::Float64
    tic::Float64
    base_peak_mz::Float64
    base_peak_intensity::Float64
    center_mz::Float64           # isolation centre (MS2); NaN for MS1
    isolation_width::Float64     # NaN for MS1
    collision_energy::Float64    # NaN for MS1
    faims_cv::Float64            # region C n_extra list (FAIMS compensation voltage); NaN when absent
    instrument_cycle::Int32
    packet_type::UInt32          # scan table +16
    packet_offset::Int64         # absolute
    packet_size::Int32
end

struct Acquisition
    path::String
    bytes::Vector{UInt8}         # memory-mapped file
    version::Int                 # u32 at 0x24
    layout::VersionLayout
    device_block::Int64
    scan_table::Int64
    scan_record::Int                 # bytes per scan-table record (88 or 80)
    packet_start::Int64
    region_c::Int64
    region_c_fixed::Int              # bytes in the fixed part of a region C record (136 or 128)
    region_c_prec::Int               # bytes per precursor entry (56 or 32)
    first_scan::Int
    last_scan::Int
    scans::Vector{ScanInfo}
    region_c_starts::Vector{Int64}   # absolute start of each scan's region C record
end

Base.show(io::IO, f::Acquisition) = print(io, "Acquisition(", basename(f.path), ", ", length(f.scans), " scans)")
Base.length(f::Acquisition) = length(f.scans)
scans(f::Acquisition) = f.scans

"""
    device_directory(b) -> (position, [(kind, index, offset), ...])

Find the device directory: the unique position k in the head where the u32 pair at k - DEVICE_COUNT_BACK is (n, n)
and the n entries at k have index == e, kind < 16, increasing in-file offsets, each landing on the device-block
signature.
"""
function device_directory(b::AbstractVector{UInt8}, count_back::Int = layout_for(b).count_back)
    fs = length(b); found = Tuple{Int,Vector{Tuple{UInt32,UInt32,Int64}}}[]
    sig(o) = (ru32(b, o), ru32(b, o + 4), ru32(b, o + 8)) == DEVICE_BLOCK_SIGNATURE
    for k in count_back:2:min(HEAD_SEARCH_BYTES, fs) - DEVICE_DIR_ENTRY_SIZE * DEVICE_DIR_MAX
        n = Int(ru32(b, k - count_back))
        (1 <= n <= DEVICE_DIR_MAX && ru32(b, k - count_back + 4) == n) || continue
        ents = Tuple{UInt32,UInt32,Int64}[]; prev = 0
        for e in 0:n-1
            x = k + e * DEVICE_DIR_ENTRY_SIZE; off = ru64(b, x + 8) % Int64
            (ru32(b, x + 4) == e && ru32(b, x) < 16 && prev < off < fs - 12 && sig(off)) || break
            push!(ents, (ru32(b, x), ru32(b, x + 4), off)); prev = off
        end
        length(ents) == n && push!(found, (k, ents))
    end
    if isempty(found)
        # a directory whose entries point past the end of the file means the file is truncated (the device blocks
        # sit at the end of the file)
        for k in count_back:2:min(HEAD_SEARCH_BYTES, fs) - DEVICE_DIR_ENTRY_SIZE * DEVICE_DIR_MAX
            n = Int(ru32(b, k - count_back))
            (1 <= n <= DEVICE_DIR_MAX && ru32(b, k - count_back + 4) == n) || continue
            all(e -> ru32(b, k + e * DEVICE_DIR_ENTRY_SIZE + 4) == e && ru32(b, k + e * DEVICE_DIR_ENTRY_SIZE) < 16, 0:n-1) || continue
            last_off = maximum(e -> ru64(b, k + e * DEVICE_DIR_ENTRY_SIZE + 8) % Int64, 0:n-1)
            last_off >= fs && error("file is truncated: its device directory points to offset $last_off, beyond the file size $fs")
        end
    end
    length(found) == 1 || error("device directory: expected one candidate, found $(length(found))")
    found[1]
end

"""Layout for the file's format version (u32 at 0x24); unknown versions are refused."""
function layout_for(b)
    v = Int(ru32(b, VERSION_OFFSET))
    haskey(LAYOUTS, v) || error("format version $v has not been decoded; refusing to guess")
    LAYOUTS[v]
end

"""
    open_acquisition(path) -> Acquisition

Memory-map the file, locate the mass-spec device block and decode the per-scan tables.
"""
function open_acquisition(path::AbstractString)
    b = Mmap.mmap(String(path))
    ru16(b, 0) == PREAMBLE_MARKER || error("$path: unexpected leading marker")
    L = layout_for(b); version = Int(ru32(b, VERSION_OFFSET))
    _, dirs = device_directory(b, L.count_back)
    i = findfirst(d -> d[1] == KIND_MASS_SPEC, dirs)
    i === nothing && error("$path: no mass-spec device in the device directory")
    dev = dirs[i][3]
    doff(k) = L.dev_offset_width == 8 ? Int64(ru64(b, dev + k)) : Int64(ru32(b, dev + k))
    doff(L.dev_self) == dev || error("$path: device block self-offset mismatch")
    first_scan = Int(ru32(b, dev + DEV_FIRST_SCAN)); last_scan = Int(ru32(b, dev + DEV_LAST_SCAN))
    st = doff(L.dev_scan_table); p0 = doff(L.dev_packet_start); rc = doff(L.dev_region_c); rd = doff(L.dev_region_d)
    n = last_scan - first_scan + 1
    rs = L.scan_record
    rc - st == n * rs || error("$path: scan table spans $(rc - st) bytes, expected $n × $rs")
    rcf, rct, rcp = L.rc_fixed, L.rc_tail, L.rc_prec
    recs = region_c_records(b, rc, rd, n, rcf, rct, rcp)
    info = Vector{ScanInfo}(undef, n)
    for j in 0:n-1
        t = st + j * rs; c = recs[j+1]
        ms = ru8(b, c + RC_MS_ORDER)
        np = ru32(b, c + rcf)
        has_prec = ms >= 0x02 && np >= 1
        pe = c + rcf + 4
        x = pe + rcp * Int64(np); x += 4 + RC_RANGE_SIZE * Int64(ru32(b, x)); x += 4 + 8 * Int64(ru32(b, x))
        cv = ru32(b, x) == 0 ? NaN : rf64(b, x + 4)
        info[j+1] = ScanInfo(
            Int32(ru32(b, t + ST_SCAN)), ms, ru8(b, c + RC_ANALYSER),
            rf64(b, t + ST_RT), rf64(b, t + ST_LOW_MZ), rf64(b, t + ST_HIGH_MZ),
            rf64(b, t + ST_TIC), rf64(b, t + ST_BP_MZ), rf64(b, t + ST_BP_INT),
            has_prec ? rf64(b, pe + PE_CENTER) : NaN, has_prec ? rf64(b, pe + PE_WIDTH) : NaN,
            has_prec ? rf64(b, pe + PE_ENERGY) : NaN, cv,
            Int32(L.st_cycle < 0 ? 0 : ru32(b, t + L.st_cycle)), ru32(b, t + ST_PKT_TYPE),
            p0 + (L.st_pkt_off_width == 8 ? Int64(ru64(b, t + L.st_pkt_off)) : Int64(ru32(b, t + L.st_pkt_off))),
            Int32(ru32(b, t + ST_PKT_SIZE)))
    end
    f = Acquisition(String(path), b, version, L, dev, st, rs, p0, rc, rcf, rcp, first_scan, last_scan, info, recs)
    # Before version 66 the precursor entry's second f64 is 1.0, not the isolation width: the per-scan extras'
    # "MS2 Isolation Width" is used instead; NaN if unavailable.
    if version < 66
        ex = scan_extras(f)
        w = ex === nothing ? nothing : get(ex.columns, "MS2 Isolation Width", nothing)
        # fallback: a single "Isolation Width: x" statement in the readable method text of the head
        mw = method_text_widths(f); fallback = length(mw) == 1 ? mw[1] : NaN
        for j in eachindex(info)
            info[j].ms_order >= 0x02 || continue
            s = info[j]; nw = w === nothing ? fallback : Float64(w[j])
            info[j] = ScanInfo(s.scan_number, s.ms_order, s.analyser, s.retention_time, s.low_mz, s.high_mz, s.tic,
                               s.base_peak_mz, s.base_peak_intensity, s.center_mz, nw, s.collision_energy, s.faims_cv,
                               s.instrument_cycle, s.packet_type, s.packet_offset, s.packet_size)
        end
    end
    f
end

"""Distinct isolation widths stated in the readable method text of the file head ("Isolation Width: 2.00")."""
function method_text_widths(f::Acquisition)
    b = f.bytes; hi = min(f.packet_start, length(b))
    txt8 = String(map(c -> 0x20 <= c < 0x7f ? c : UInt8('\n'), b[1:hi]))
    txt16 = String(map(i -> (b[i+1] == 0 && 0x20 <= b[i] < 0x7f) ? Char(b[i]) : '\n', 1:2:hi-1))   # method text may be UTF-16
    unique([parse(Float64, m.captures[1]) for t in (txt8, txt16) for m in eachmatch(r"Isolation Width:\s*([0-9]+\.?[0-9]*)", t)])
end

"""
    region_c_records(b, rc, rd, n) -> Vector{Int64}

Walk region C (variable-length records, see layout.jl) and return each record's absolute start. The walk must
produce exactly `n` records and end exactly at region D.
"""
function region_c_records(b, rc::Int64, rd::Int64, n::Int, fixed::Int, tail::Int, prec::Int)
    out = Vector{Int64}(undef, n); k = rc + RC_LEAD
    for j in 1:n
        out[j] = k
        x = k + fixed
        x += 4 + prec * Int64(ru32(b, x))
        x += 4 + RC_RANGE_SIZE * Int64(ru32(b, x))
        x += 4 + 8 * Int64(ru32(b, x))
        x += 4 + 8 * Int64(ru32(b, x)) + tail
        x <= rd || error("region C walk overruns region D at record $j")
        k = x
    end
    k == rd || error("region C walk ends at $k, region D starts at $rd")
    out
end

"""Byte offset of the peak list within a packet at absolute offset `p`, and its peak count."""
@inline function peak_list_location(b, p::Int64)
    pk = p + PKT_HEADER + 4 * Int64(ru32(b, p + PKT_W_PROFILE))
    # an empty scan can have no peak-list section at all (W2 == 0, not even the count word)
    ru32(b, p + PKT_W_PEAKS) == 0 && return pk, 0
    pk, Int(ru32(b, pk))
end

"""
    peaks(f, i) -> (mz::Vector{Float32}, intensity::Vector{Float32})

Centroided peak list of the `i`-th scan (1-based position in `scans(f)`).
"""
function peaks(f::Acquisition, i::Integer)
    mz = Float32[]; it = Float32[]
    peaks!(mz, it, f, i)
end

function peaks!(mz::Vector{Float32}, it::Vector{Float32}, f::Acquisition, i::Integer)
    b = f.bytes
    t = f.scans[i].packet_type
    t in KNOWN_PACKET_TYPES || error("scan $(f.scans[i].scan_number): packet type $t has not been decoded; refusing to guess")
    p = f.scans[i].packet_offset
    pk, n = peak_list_location(b, p)
    pw = f.layout.peak_words
    n == 0 || ru32(b, p + PKT_W_PEAKS) == 1 + pw * n || error("scan $(f.scans[i].scan_number): peak-list size disagrees with layout")
    t == 0x00000013 && n != 0 && error("scan $(f.scans[i].scan_number): profile-only packet (type 19) with a peak list; not decoded")
    resize!(mz, n); resize!(it, n)
    if pw == 2
        @inbounds for j in 0:n-1
            mz[j+1] = rf32(b, pk + 4 + 8j); it[j+1] = rf32(b, pk + 8 + 8j)
        end
    else   # version 62: f64 m/z + f32 intensity (converted to f32 for Pioneer's schema)
        @inbounds for j in 0:n-1
            mz[j+1] = Float32(rf64(b, pk + 4 + 12j)); it[j+1] = rf32(b, pk + 12 + 12j)
        end
    end
    mz, it
end

"""
    profile(f, i) -> Union{Nothing, (mz::Vector{Float64}, intensity::Vector{Float32})}

Profile points of a scan that stores them (packet type 21; `nothing` otherwise). Runs of consecutive points:
u32 start index, u32 count, [f32 per-run value, versions >= 63], count × f32 intensity. m/z from the per-scan
calibration (region C coefficients 4 and 5 = extras "Conversion Parameter B, C"): f = A + B'·index with A, B' the
profile section's two f64, m/z = B/f² + C/f⁴. Packet type 19 (profile only) uses a linear m/z scale
(`_profile_linear`). `peaks` returns the stored centroids, not centroids of these points.
"""
function profile(f::Acquisition, i::Integer)
    s = f.scans[i]
    s.packet_type == 0x00000013 && return _profile_linear(f, i)
    s.packet_type == 0x00000015 || return nothing
    b = f.bytes; o = s.packet_offset; ru32(b, o + PKT_W_PROFILE) == 0 && return nothing
    c = f.region_c_starts[i] + f.region_c_fixed; c += 4 + f.region_c_prec * Int64(ru32(b, c)); c += 4 + RC_RANGE_SIZE * Int64(ru32(b, c))
    nc = Int(ru32(b, c)); nc >= 5 || return nothing
    B = rf64(b, c + 4 + 8 * 3); C = rf64(b, c + 4 + 8 * 4)
    ps = o + PKT_HEADER; A0 = rf64(b, ps); A1 = rf64(b, ps + 8); nr = Int(ru32(b, ps + 16)); k = ps + 24
    hdr = f.version >= 63 ? 12 : 8
    mz = Float64[]; it = Float32[]
    for _ in 1:nr
        st = Int(ru32(b, k)); n = Int(ru32(b, k + 4))
        for j in 0:n-1
            fr = A0 + A1 * (st + j); push!(mz, B / fr^2 + C / fr^4); push!(it, rf32(b, k + hdr + 4j))
        end
        k += hdr + 4n
    end
    p = sortperm(mz); mz[p], it[p]
end

"""
Profile points of an ion-trap profile-only scan (packet type 19): one run of n points after an 8-byte run header
(u32 start, u32 n); point j has m/z = A0 + A1·(start + j + 1) with A0, A1 the profile section's two f64 (A0 = the scan's
low m/z, A1 = the point spacing). Packets with more than one run are refused.
"""
_profile_linear(f::Acquisition, i::Integer) = _profile_linear!(Float64[], Float32[], f, i)

function _profile_linear!(mz::Vector{Float64}, it::Vector{Float32}, f::Acquisition, i::Integer)
    s = f.scans[i]; b = f.bytes; ps = s.packet_offset + PKT_HEADER
    A0 = rf64(b, ps); A1 = rf64(b, ps + 8); nr = Int(ru32(b, ps + 16)); k = ps + 24
    nr == 1 || error("scan $(s.scan_number): profile-only packet with $nr runs; only single runs have been decoded")
    st = Int(ru32(b, k)); n = Int(ru32(b, k + 4))
    checkbounds(b, k + 8 + 1 : k + 8 + 4n)
    resize!(mz, n); resize!(it, n)
    @inbounds for j in 0:n-1
        mz[j+1] = A0 + A1 * (st + j + 1); it[j+1] = rf32(b, k + 8 + 4j)
    end
    mz, it
end

"""
    features(f, i) -> Union{Nothing, NamedTuple}

Packet section 6 (`nothing` when absent or not of this form): the instrument's charge deconvolution. Returns
  `peak_feature` :: Vector{Int}   feature number per centroid (0-based; -1 = none), aligned with `peaks(f, i)`,
  `peak_flags`   :: Vector{UInt8} top 3 bits of each per-peak u16 (6 = 0xC000 = monoisotopic peak of a charge state),
  `mass`, `score`, `peak` :: per feature: neutral monoisotopic mass (f64), score (f32), 0-based peak index.
Layout: [4 words: tag, 0x00080000, ceil(n/2)+1, 0][n × u16][pad][3 words][K × (f64, f32, u32)];
the section size must equal 4 + ceil(n/2) + 3 + 4K words. Not written to Pioneer's `.arrow`.
"""
function features(f::Acquisition, i::Integer)
    b = f.bytes; o = f.scans[i].packet_offset
    W = [Int(ru32(b, o + q)) for q in (PKT_W_PROFILE, PKT_W_PEAKS, PKT_W_REST...)]
    W[6] == 0 && return nothing
    pk, n = peak_list_location(b, o); n == 0 && return nothing
    s6 = o + PKT_HEADER + 4 * sum(W[1:5])
    Int(ru32(b, s6 + 8)) == cld(n, 2) + 1 || return nothing          # another section-6 form (e.g. first word 0x1)
    u = [ru16(b, s6 + 16 + 2k) for k in 0:n-1]
    ids = Int.(u .& 0x1fff); K = isempty(filter(!=(0x1fff), ids)) ? 0 : maximum(filter(!=(0x1fff), ids)) + 1
    W[6] == 4 + cld(n, 2) + 3 + 4K || return nothing
    rb = s6 + 16 + 4 * cld(n, 2) + 12
    (peak_feature = [x == 0x1fff ? -1 : x for x in ids], peak_flags = UInt8.(u .>> 13),
     mass = [rf64(b, rb + 16c) for c in 0:K-1], score = [Float64(rf32(b, rb + 16c + 8)) for c in 0:K-1],
     peak = [Int(ru32(b, rb + 16c + 12)) for c in 0:K-1])
end
