# The `.scxs` container: centroided SCIEX scans, one zstd block per scan (the `.tdfs` layout of TimsSlices.jl,
# with a per-scan m/z calibration). See notes/scxs_format.md.
#
#   <name>.scxs/
#     meta.json    format_version, source, converter parameters, bin_scale, int_scale, counts
#     scans.arrow  one row per scan (= one Pioneer scan), in acquisition order; empty records are not written
#     blocks.bin   one zstd block per scan, in scans.arrow order, no headers
#
# m/z of a stored bin q:        (cal_a · (q / bin_scale / 40 − cal_b))²
# intensity of a stored value v: v / int_scale

using Arrow, JSON3

const SCXS_FORMAT_VERSION = 1

"Per-scan columns of scans.arrow (struct of vectors)."
Base.@kwdef mutable struct ScanRows
    record::Vector{Int32} = Int32[]            # 1-based Idx record in the .wiff
    cycle::Vector{Int32} = Int32[]             # acquisition cycle (1-based)
    experiment::Vector{Int16} = Int16[]        # experiment within the cycle (1 = MS1)
    ms_order::Vector{UInt8} = UInt8[]
    retention_time::Vector{Float32} = Float32[]   # minutes
    low_mz::Vector{Float32} = Float32[]
    high_mz::Vector{Float32} = Float32[]
    tic::Vector{Float32} = Float32[]           # the instrument's TIC (Idx), as msConvert reports it
    base_peak_mz::Vector{Float32} = Float32[]          # vendor base peak (Idx), else the profile maximum
    base_peak_intensity::Vector{Float32} = Float32[]   # profile units, not the stored centroids' (see `base_peak`)
    center_mz::Vector{Float32} = Float32[]     # NaN for MS1
    isolation_width::Vector{Float32} = Float32[]
    cal_a::Vector{Float64} = Float64[]
    cal_b::Vector{Float64} = Float64[]
    n_peaks::Vector{Int32} = Int32[]
    block_offset::Vector{Int64} = Int64[]
    block_size::Vector{Int32} = Int32[]
end
Base.length(r::ScanRows) = length(r.record)

mutable struct ScxsWriter
    dir::String
    blocks::IOStream
    offset::Int64
    rows::ScanRows
    meta::Dict{String, Any}
end

function ScxsWriter(dir::AbstractString, meta::Dict{String, Any})
    mkpath(dir)
    ScxsWriter(String(dir), open(joinpath(dir, "blocks.bin"), "w"), 0, ScanRows(), meta)
end

"Append one scan's row (with `block_offset`/`block_size` filled in here) and its compressed block."
function write_scan!(w::ScxsWriter, row::NamedTuple, zbytes::AbstractVector{UInt8})
    r = w.rows
    push!(r.record, row.record); push!(r.cycle, row.cycle); push!(r.experiment, row.experiment)
    push!(r.ms_order, row.ms_order); push!(r.retention_time, row.retention_time)
    push!(r.low_mz, row.low_mz); push!(r.high_mz, row.high_mz); push!(r.tic, row.tic)
    push!(r.base_peak_mz, row.base_peak_mz); push!(r.base_peak_intensity, row.base_peak_intensity)
    push!(r.center_mz, row.center_mz); push!(r.isolation_width, row.isolation_width)
    push!(r.cal_a, row.cal_a); push!(r.cal_b, row.cal_b); push!(r.n_peaks, row.n_peaks)
    push!(r.block_offset, w.offset); push!(r.block_size, Int32(length(zbytes)))
    write(w.blocks, zbytes)
    w.offset += length(zbytes)
    w
end

function Base.close(w::ScxsWriter)
    close(w.blocks)
    r = w.rows
    Arrow.write(joinpath(w.dir, "scans.arrow"), NamedTuple(f => getfield(r, f) for f in fieldnames(ScanRows)))
    m = copy(w.meta)
    m["format_version"] = SCXS_FORMAT_VERSION
    m["n_scans"] = length(r); m["n_peaks"] = sum(Int, r.n_peaks; init = 0); m["blocks_bytes"] = w.offset
    open(joinpath(w.dir, "meta.json"), "w") do io
        JSON3.pretty(io, m)
    end
    nothing
end

# --- reader -------------------------------------------------------------------------------------------------

struct ScxsFile
    dir::String
    meta::Dict{String, Any}
    scans::NamedTuple
    blocks::Vector{UInt8}
    bin_scale::Int
    int_scale::Float64
end

function open_scxs(dir::AbstractString)
    meta = Dict{String, Any}(JSON3.read(read(joinpath(dir, "meta.json"), String), Dict{String, Any}))
    Int(meta["format_version"]) == SCXS_FORMAT_VERSION ||
        error("scxs format version $(meta["format_version"]) (reader is $SCXS_FORMAT_VERSION)")
    t = Arrow.Table(joinpath(dir, "scans.arrow"))
    scans = NamedTuple(k => collect(getproperty(t, k)) for k in propertynames(t))
    blocks = open(joinpath(dir, "blocks.bin"), "r") do io
        Mmap.mmap(io, Vector{UInt8}, filesize(io))
    end
    ScxsFile(String(dir), meta, scans, blocks, Int(meta["bin_scale"]), Float64(meta["int_scale"]))
end
n_scans(f::ScxsFile) = length(f.scans.record)

"Decode scan `s` (1-based row of scans.arrow) into `bins`/`ints`; returns the peak count."
function read_block!(bins::Vector{UInt32}, ints::Vector{UInt32}, c::BlockCodec, f::ScxsFile, s::Integer)
    off = f.scans.block_offset[s]; bs = f.scans.block_size[s]; np = Int(f.scans.n_peaks[s])
    off + bs <= length(f.blocks) || error("scxs scan $s: block runs beyond blocks.bin")
    decode_block!(bins, ints, c, view(f.blocks, off+1:off+bs), np)
    np
end

@inline stored_mz(f::ScxsFile, s::Integer, q::UInt32) = bin_to_mz(Float64(q) / f.bin_scale, f.scans.cal_a[s], f.scans.cal_b[s])
@inline stored_intensity(f::ScxsFile, v::UInt32) = Float64(v) / f.int_scale
