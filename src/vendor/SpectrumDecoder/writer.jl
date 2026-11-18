# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Pioneer's spectrum `.arrow` input (the schema read by src/structs/MassSpecData/{types,getters}.jl).

const ARROW_BATCH_PEAKS = 4_000_000   # peaks per record batch: ~32 MB per buffer set, 3 sets

"""
    pioneer_cycle_idxs(ms_orders, centers, widths)

Cycle index with the same rule as Pioneer's own fallback (`compute_cycle_idxs`): a new cycle starts each time the
first MS2 isolation window recurs (within 1e-3 m/z), i.e. one cycle = one full pass over all windows.
"""
function pioneer_cycle_idxs(ms_orders, centers, widths)
    out = Vector{Int32}(undef, length(ms_orders))
    cyc = Int32(1); have = false; c0 = 0f0; w0 = 0f0
    for i in eachindex(ms_orders)
        if ms_orders[i] == 0x02
            c = Float32(centers[i]); w = Float32(widths[i])
            if !have
                have = true; c0 = c; w0 = w
            elseif abs(c - c0) <= 1f-3 && abs(w - w0) <= 1f-3
                cyc += Int32(1)
            end
        end
        out[i] = cyc
    end
    out
end

"""
    _peak_list_column(values, offsets)

A list column of `Union{Missing,Float32}` (Pioneer's peak-array element type) built directly as an `Arrow.List` over a
plain `Vector{Float32}` and 0-based `Int32` offsets, with no nulls. Arrow.jl passes an `Arrow.List` through unchanged
and writes its value buffer in one call, rather than element by element through a boxed `Union`.
"""
function _peak_list_column(values::Vector{Float32}, offsets::Vector{Int32})
    n = length(offsets) - 1; m = length(values)
    child = Arrow.Primitive(Union{Missing,Float32}, UInt8[], Arrow.ValidityBitmap(UInt8[], 1, m, 0), values, m, nothing)
    T = SubArray{Union{Missing,Float32}, 1, typeof(child), Tuple{UnitRange{Int64}}, true}
    Arrow.List{T, Int32, typeof(child)}(UInt8[], Arrow.ValidityBitmap(UInt8[], 1, n, 0),
                                        Arrow.Offsets{Int32}(UInt8[], offsets), child, n, nothing)
end

"""
Tell the OS it may drop the mapped pages of the input file between byte offsets `lo` and `hi` (0-based). They are clean
file pages: reading them again faults them back in. On Linux this removes them from the process's resident set, so
peak memory does not grow with the file size. Failures are ignored.
"""
function _release_pages(f::Acquisition, lo::Integer, hi::Integer)
    Sys.iswindows() && return   # no madvise on Windows
    ps = Int(ccall(:getpagesize, Cint, ())); a = (Int(lo) ÷ ps) * ps; z = min((cld(Int(hi), ps)) * ps, length(f.bytes))
    z > a || return
    GC.@preserve f ccall(:madvise, Cint, (Ptr{Cvoid}, Csize_t, Cint), pointer(f.bytes) + a, z - a, Mmap.MADV_DONTNEED)
    nothing
end
_packet_span(f::Acquisition, rows) = (f.scans[first(rows)].packet_offset,
                                      f.scans[last(rows)].packet_offset + f.scans[last(rows)].packet_size)

"""Contiguous ranges covering `rows`, at most `k`, with roughly equal total `weight` (one extra unit per scan)."""
function _split_ranges(weight, rows::UnitRange{Int}, k::Int)
    k = max(1, min(k, length(rows))); k == 1 && return [rows]
    tot = 0; for i in rows; tot += weight(i) + 1; end
    out = UnitRange{Int}[]; lo = first(rows); acc = 0
    for i in rows
        acc += weight(i) + 1
        if length(out) < k - 1 && acc * k >= tot * (length(out) + 1)
            push!(out, lo:i); lo = i + 1
        end
    end
    lo <= last(rows) && push!(out, lo:last(rows))
    out
end

"""Peak count of every scan, from the packet headers alone (no peak data read), split over `ntasks` tasks. With
`release`, each task releases the mapped pages behind every 64 MiB it has walked (a fault maps neighbouring pages
too, so a header-only pass would otherwise make most of the file resident)."""
function peak_counts(f::Acquisition; release::Bool = false, ntasks::Int = 1)
    b = f.bytes; S = scans(f); out = Vector{Int}(undef, length(S))
    isempty(S) && return out
    @sync for r in _split_ranges(_ -> 0, 1:length(S), ntasks)
        Threads.@spawn begin
            done = S[first(r)].packet_offset
            for i in r
                p = S[i].packet_offset
                out[i] = peak_list_location(b, p)[2]
                if release && p - done >= RELEASE_STEP
                    _release_pages(f, done, p); done = p
                end
            end
        end
    end
    out
end

const RELEASE_STEP = 64 * 2^20   # bytes of input walked between page releases

_ms2f32(s::ScanInfo, x::Float64) = s.ms_order >= 0x02 && !isnan(x) ? Float32(x) : missing

"""Per-task scratch: one scan's peaks, and buffers for centroiding profile-only scans."""
struct _ScanBufs
    mz::Vector{Float32}; it::Vector{Float32}
    pmz::Vector{Float64}; pit::Vector{Float32}; tmp::Vector{Float32}; cent::Vector{Centroid}
end
_ScanBufs() = _ScanBufs(Float32[], Float32[], Float64[], Float32[], Float32[], Centroid[])

"""
Peaks written for scan `i` (into `B.mz`, `B.it`, sorted by m/z): the stored centroids, or for a profile-only scan
(packet type 19) centroids of its profile points when `centroid_profile` (fitted height as intensity; the intensity
gate is half the median of the positive points, estimated from every 8th point). Stored peak lists can hold a peak
out of m/z order; every scan is written sorted by m/z.
"""
function _scan_peaks!(B::_ScanBufs, f::Acquisition, i::Int, centroid_profile::Bool)
    peaks!(B.mz, B.it, f, i)
    if centroid_profile && f.scans[i].packet_type == 0x00000013
        _profile_linear!(B.pmz, B.pit, f, i)
        centroid!(B.cent, B.pmz, B.pit; threshold = 0.5f0 * _median_positive!(B.tmp, B.pit; stride = 8))
        resize!(B.mz, length(B.cent)); resize!(B.it, length(B.cent))
        @inbounds for (j, c) in enumerate(B.cent); B.mz[j] = Float32(c.mz); B.it[j] = c.height; end
    end
    if !issorted(B.mz)
        o = sortperm(B.mz); permute!(B.mz, o); permute!(B.it, o)
    end
    nothing
end

"""
Fill one record batch's peak arrays for scans `rows`, split over `ntasks` tasks. Stored peak counts are known from
the header pass, so every scan's place in the batch is known in advance: tasks fill disjoint slices of `mz`/`it`
directly and the order is the scan order whichever task ran. Batches holding centroided profile-only scans (counts
not known in advance) are assembled per task first and copied in after.
"""
function _fill_batch!(mz::Vector{Float32}, it::Vector{Float32}, offs::Vector{Int32}, f::Acquisition, rows::UnitRange{Int},
                      counts::Vector{Int}, centroid_profile::Bool, ntasks::Int, bufs::Vector{_ScanBufs})
    S = scans(f); n = length(rows)
    grow = centroid_profile && any(i -> S[i].packet_type == 0x00000013, rows)
    ranges = _split_ranges(i -> counts[i], rows, ntasks)
    if !grow
        offs[1] = 0
        for (r, i) in enumerate(rows); offs[r+1] = offs[r] + counts[i]; end
        npk = Int(offs[n+1]); resize!(mz, npk); resize!(it, npk)
        @sync for (t, r) in enumerate(ranges)
            Threads.@spawn begin
                B = bufs[t]
                for i in r
                    _scan_peaks!(B, f, i, false)
                    o = Int(offs[i - first(rows) + 1]); m = length(B.mz)
                    m == counts[i] || error("scan $(S[i].scan_number): peak count changed between header pass and read")
                    copyto!(mz, o + 1, B.mz, 1, m); copyto!(it, o + 1, B.it, 1, m)
                end
            end
        end
    else
        parts = Vector{Tuple{Vector{Float32},Vector{Float32},Vector{Int}}}(undef, length(ranges))
        @sync for (t, r) in enumerate(ranges)
            Threads.@spawn begin
                B = bufs[t]; pm = Float32[]; pi_ = Float32[]; pc = Int[]
                for i in r
                    _scan_peaks!(B, f, i, true); append!(pm, B.mz); append!(pi_, B.it); push!(pc, length(B.mz))
                end
                parts[t] = (pm, pi_, pc)
            end
        end
        npk = sum(p -> length(p[1]), parts); npk <= typemax(Int32) || error("record batch exceeds Int32 list offsets")
        resize!(mz, npk); resize!(it, npk); k = 0; r = 1; offs[1] = 0
        for (pm, pi_, pc) in parts
            copyto!(mz, k + 1, pm, 1, length(pm)); copyto!(it, k + 1, pi_, 1, length(pi_)); k += length(pm)
            for c in pc; offs[r+1] = offs[r] + c; r += 1; end
        end
    end
    nothing
end

function _batch_table(f::Acquisition, rows::UnitRange{Int}, mz, it, offs, cycles)
    R = view(scans(f), rows); n = length(rows)
    (mz_array = _peak_list_column(mz, offs),
     intensity_array = _peak_list_column(it, offs),
     scanHeader = fill("", n),
     scanNumber = Int32[s.scan_number for s in R],
     packetType = zeros(Int32, n),
     retentionTime = Float32[s.retention_time for s in R],
     lowMz = Float32[s.low_mz for s in R], highMz = Float32[s.high_mz for s in R],
     TIC = Float32[s.tic for s in R],
     basePeakMz = Float32[s.base_peak_mz for s in R],
     basePeakIntensity = Float32[s.base_peak_intensity for s in R],
     centerMz = Union{Missing,Float32}[_ms2f32(s, s.center_mz) for s in R],
     isolationWidthMz = Union{Missing,Float32}[_ms2f32(s, s.isolation_width) for s in R],
     collisionEnergyField = Union{Missing,Float32}[_ms2f32(s, s.collision_energy) for s in R],
     collisionEnergyEvField = zeros(Float32, n),
     msOrder = UInt8[s.ms_order for s in R],
     cycle_idx = cycles[rows])
end

"""Threads used for conversion by default: at most 4, fewer when Julia was started with fewer."""
default_threads() = min(4, Threads.nthreads())

"""
    write_pioneer_arrow(f::Acquisition, out_path; threads = default_threads(), batch_peaks = ARROW_BATCH_PEAKS,
                        centroid_profile = true, release_pages = true)

Write every scan's centroided peaks and header fields in Pioneer's input schema, as record batches of at most
`batch_peaks` peaks (Int32 list offsets cap it at 2^31-1).

Pipeline (bounded memory): the header pass and the filling of each batch are split over `threads` tasks;
Arrow.jl's own writer task writes batch b while batch b+1 is filled (on its own thread when Julia has more than one).
The writer's message channel is bounded to one batch (`ntasks = 1`; Arrow.jl's default is unbounded), so once batch b
has been handed over, batch b-2 is fully written: peak buffers come from a ring of three sets and are reused, never
reallocated. Peak memory ≈ 3 × batch buffers + the mapped pages not yet released (`release_pages`, every 64 MiB).
`centroid_profile`: profile-only scans (ion-trap MS1, packet type 19) are written with centroids of their profile.
`acquisition_label`: the acquisition mode stated by the embedded method (`acquisition_metadata`) is written as schema
metadata; files whose method states none get no metadata.
"""
function write_pioneer_arrow(f::Acquisition, out_path::AbstractString; auxiliary::Bool = false,
                             batch_peaks::Int = ARROW_BATCH_PEAKS, centroid_profile::Bool = true,
                             release_pages::Bool = true, threads::Int = default_threads(),
                             acquisition_label::Bool = true)
    auxiliary && write_auxiliary(f, out_path)
    S = scans(f); nt = max(1, threads)
    cycles = pioneer_cycle_idxs([s.ms_order for s in S], [s.center_mz for s in S], [s.isolation_width for s in S])
    counts = peak_counts(f; release = release_pages, ntasks = nt)
    # profile-only scans have no stored peak count: batch sizing counts 1,024 centroids for each
    est(i) = centroid_profile && S[i].packet_type == 0x00000013 ? 1024 : counts[i]
    ring = [(Float32[], Float32[], Int32[]) for _ in 1:3]
    bufs = [_ScanBufs() for _ in 1:nt]
    meta = acquisition_label ? acquisition_metadata(f) : nothing
    writer = open(Arrow.Writer, String(out_path); ntasks = 1, metadata = meta)
    try
        lo = 1; b = 0
        while lo <= length(S) || (b == 0 && isempty(S))
            hi = lo; npk = isempty(S) ? 0 : est(lo)
            while hi < length(S) && npk + est(hi + 1) <= batch_peaks; hi += 1; npk += est(hi); end
            rows = isempty(S) ? (1:0) : (lo:hi)
            mz, it, offs = ring[b % 3 + 1]; b += 1
            resize!(offs, length(rows) + 1); offs[1] = 0
            isempty(rows) || _fill_batch!(mz, it, offs, f, rows, counts, centroid_profile, nt, bufs)
            Arrow.write(writer, _batch_table(f, rows, mz, it, offs, cycles))
            release_pages && !isempty(rows) && _release_pages(f, _packet_span(f, rows)...)
            isempty(S) && break
            lo = hi + 1
        end
    finally
        close(writer)
    end
    out_path
end

"""
    write_auxiliary(f, out_path) -> Vector{String} (warnings)

Optional sidecar files next to `out_path` (Pioneer's file is unchanged): `<base>.lc.arrow` (long table: device,
name, time_min, value for every trace and status-log field) and `<base>.scan_extras.arrow` (scanNumber plus one
column per per-scan field). Never throws; returns the warnings.
"""
function write_auxiliary(f::Acquisition, out_path::AbstractString)
    a = read_auxiliary(f); w = copy(a.warnings); base = replace(String(out_path), r"\.arrow$" => "")
    try
        dev = Int32[]; name = String[]; t = Float64[]; v = Float64[]
        for tr in a.traces
            append!(dev, fill(Int32(tr.device), length(tr.times))); append!(name, fill(tr.name, length(tr.times)))
            append!(t, tr.times); append!(v, tr.values)
        end
        for lg in a.logs, (c, lab) in enumerate(lg.labels)
            col = lg.values[:, c]; all(x -> x isa Real, col) || continue   # text columns (event log) are not numeric
            append!(dev, fill(Int32(lg.device), length(lg.times))); append!(name, fill(lab, length(lg.times)))
            append!(t, lg.times); append!(v, Float64.(col))
        end
        isempty(t) || Arrow.write(base * ".lc.arrow", (device = dev, name = name, time_min = t, value = v))
    catch e
        push!(w, "lc sidecar: " * sprint(showerror, e))
    end
    try
        if a.extras !== nothing
            cols = Pair{Symbol,Any}[:scanNumber => Int32[s.scan_number for s in scans(f)]]
            for l in a.extras.labels
                c = a.extras.columns[l]
                push!(cols, Symbol(l) => (eltype(c) <: Number || all(x -> x isa Number, c)) ? Float64.(c) : String.(string.(c)))
            end
            Arrow.write(base * ".scan_extras.arrow", NamedTuple(cols))
        end
    catch e
        push!(w, "extras sidecar: " * sprint(showerror, e))
    end
    w
end
