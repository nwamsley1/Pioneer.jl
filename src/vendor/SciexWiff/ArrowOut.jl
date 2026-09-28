# Pioneer's scan Arrow schema (the columns `convertMzML` writes, PioneerScanElement), in record batches.

const ARROW_BATCH_PEAKS = 50_000_000   # list offsets are Int32 per record batch; this also bounds pending RAM

mutable struct PioneerArrowWriter
    writer::Arrow.Writer
    scan_number::Int32
    mz::Vector{Union{Missing, Float32}}
    intensity::Vector{Union{Missing, Float32}}
    starts::Vector{Int}
    retention_time::Vector{Float32}; low_mz::Vector{Float32}; high_mz::Vector{Float32}; tic::Vector{Float32}
    center_mz::Vector{Union{Missing, Float32}}; isolation_width::Vector{Union{Missing, Float32}}
    ms_order::Vector{UInt8}; cycle_idx::Vector{Int32}
end

PioneerArrowWriter(path::AbstractString; metadata = Dict{String, String}()) =
    PioneerArrowWriter(open(Arrow.Writer, String(path); metadata = metadata), 0,
        Union{Missing, Float32}[], Union{Missing, Float32}[], Int[1], Float32[], Float32[], Float32[], Float32[],
        Union{Missing, Float32}[], Union{Missing, Float32}[], UInt8[], Int32[])

"Append one scan: peaks `mz[1:n]`, `intensity[1:n]` and its row (same fields as `ScanRows`)."
function write_scan!(w::PioneerArrowWriter, row::NamedTuple, mz::AbstractVector{Float32}, intensity::AbstractVector{Float32}, n::Int)
    @inbounds for k in 1:n
        push!(w.mz, mz[k]); push!(w.intensity, intensity[k])
    end
    push!(w.starts, length(w.mz) + 1)
    ms1 = row.ms_order == 0x01
    push!(w.retention_time, row.retention_time); push!(w.low_mz, row.low_mz); push!(w.high_mz, row.high_mz)
    push!(w.tic, row.tic)
    push!(w.center_mz, ms1 ? missing : row.center_mz); push!(w.isolation_width, ms1 ? missing : row.isolation_width)
    push!(w.ms_order, row.ms_order); push!(w.cycle_idx, row.cycle)
    length(w.mz) >= ARROW_BATCH_PEAKS && flush_batch!(w)
    w
end

function flush_batch!(w::PioneerArrowWriter)
    n = length(w.retention_time)
    n == 0 && return w
    tbl = (mz_array = [view(w.mz, w.starts[r]:w.starts[r+1]-1) for r in 1:n],
           intensity_array = [view(w.intensity, w.starts[r]:w.starts[r+1]-1) for r in 1:n],
           scanHeader = fill("", n), scanNumber = Int32.(w.scan_number+1:w.scan_number+n), packetType = zeros(Int32, n),
           retentionTime = w.retention_time, lowMz = w.low_mz, highMz = w.high_mz, TIC = w.tic,
           centerMz = w.center_mz, isolationWidthMz = w.isolation_width,
           collisionEnergyField = Vector{Union{Missing, Float32}}(missing, n), collisionEnergyEvField = zeros(Float32, n),
           msOrder = w.ms_order, cycle_idx = w.cycle_idx)
    Arrow.write(w.writer, tbl)
    w.scan_number += Int32(n)
    # Arrow.Writer may serialise from these columns after returning: hand them over and start new ones
    w.mz = Union{Missing, Float32}[]; w.intensity = Union{Missing, Float32}[]; w.starts = Int[1]
    for f in (:retention_time, :low_mz, :high_mz, :tic, :center_mz, :isolation_width, :ms_order, :cycle_idx)
        setfield!(w, f, similar(getfield(w, f), 0))
    end
    w
end

function Base.close(w::PioneerArrowWriter)
    flush_batch!(w)
    close(w.writer)
    nothing
end
