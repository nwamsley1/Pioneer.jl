# Copyright (C) 2026 Nathan Wamsley
#
# This file is part of TimsSlices.jl
#
# TimsSlices.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# Pioneer's slice Arrow schema (the prototype's output), written in record batches from quantised frames.
# Used for search-level validation and for parameter sweeps against today's numbers.

const ARROW_BATCH_PEAKS = 100_000_000    # Arrow list offsets are Int32 per record batch; smaller batches bound the RAM of the pending columns

mutable struct SliceArrowWriter
    path::String
    writer::Arrow.Writer
    lo_mz::Float32
    hi_mz::Float32
    cal::LinearMzCal
    bin_scale::Int
    int_scale::Float64
    scan_number::Int32
    cycle::Int32
    # pending batch
    mz::Vector{Union{Missing, Float32}}
    intensity::Vector{Union{Missing, Float32}}
    starts::Vector{Int}
    retention_time::Vector{Float32}; tic::Vector{Float32}
    center_mz::Vector{Union{Missing, Float32}}; isolation_width::Vector{Union{Missing, Float32}}
    ce_field::Vector{Union{Missing, Float32}}; ce_ev::Vector{Float32}
    ms_order::Vector{UInt8}; cycle_idx::Vector{Int32}; frame_id::Vector{Int32}; im_scan::Vector{UInt16}; window_group::Vector{UInt8}
    n_batches::Int
end

function SliceArrowWriter(path::AbstractString, meta::Dict{String, String}, cal::LinearMzCal, bin_scale::Integer, int_scale::Real, lo_mz::Real, hi_mz::Real)
    w = open(Arrow.Writer, String(path); metadata = meta)
    SliceArrowWriter(String(path), w, Float32(lo_mz), Float32(hi_mz), cal, Int(bin_scale), Float64(int_scale), 0, 0,
                     Union{Missing, Float32}[], Union{Missing, Float32}[], Int[1], Float32[], Float32[],
                     Union{Missing, Float32}[], Union{Missing, Float32}[], Union{Missing, Float32}[], Float32[],
                     UInt8[], Int32[], Int32[], UInt16[], UInt8[], 0)
end

"Append the slices of one frame (frames in order; the cycle index is derived from MS1 frames as in the prototype)."
function write_frame!(w::SliceArrowWriter, fm::FrameMeta, rows::SliceRows, blk::SliceBlock)
    fm.ms_order == 0x01 && (w.cycle += Int32(1))
    ms1 = fm.ms_order == 0x01
    @inbounds for j in 1:blk.n_slices
        r = slice_range(blk, j)
        for k in r
            push!(w.mz, Float32(bin_to_mz(w.cal, Float64(blk.bin[k]) / w.bin_scale)))
            push!(w.intensity, Float32(blk.intensity[k] / w.int_scale))   # stored integers are int_scale x intensity
        end
        push!(w.starts, length(w.mz) + 1)
        push!(w.retention_time, rows.retention_time[j]); push!(w.tic, rows.tic[j])
        push!(w.center_mz, ms1 ? missing : rows.center_mz[j]); push!(w.isolation_width, ms1 ? missing : rows.isolation_width[j])
        push!(w.ce_field, ms1 ? missing : rows.window_ce[j]); push!(w.ce_ev, rows.collision_energy_ev[j])
        push!(w.ms_order, ms1 ? 0x01 : 0x02); push!(w.cycle_idx, w.cycle); push!(w.frame_id, fm.frame_id)
        push!(w.im_scan, rows.im_scan[j]); push!(w.window_group, ms1 ? 0x00 : fm.window_group)
    end
    length(w.mz) >= ARROW_BATCH_PEAKS && flush_batch!(w)
    w
end

function flush_batch!(w::SliceArrowWriter)
    n = length(w.retention_time)
    n == 0 && return w
    mz_views = [view(w.mz, w.starts[r]:w.starts[r+1]-1) for r in 1:n]
    int_views = [view(w.intensity, w.starts[r]:w.starts[r+1]-1) for r in 1:n]
    tbl = (mz_array = mz_views, intensity_array = int_views,
           scanHeader = fill("", n), scanNumber = Int32.(w.scan_number+1:w.scan_number+n), packetType = zeros(Int32, n),
           retentionTime = w.retention_time, lowMz = fill(w.lo_mz, n), highMz = fill(w.hi_mz, n), TIC = w.tic,
           centerMz = w.center_mz, isolationWidthMz = w.isolation_width,
           collisionEnergyField = w.ce_field, collisionEnergyEvField = w.ce_ev,
           msOrder = w.ms_order, cycle_idx = w.cycle_idx, frameId = w.frame_id, imScan = w.im_scan, windowGroup = w.window_group)
    Arrow.write(w.writer, tbl)
    w.n_batches += 1
    w.scan_number += Int32(n)
    # Arrow.Writer serialises asynchronously from the columns it was handed: give it ownership and start new ones
    w.mz = Union{Missing, Float32}[]; w.intensity = Union{Missing, Float32}[]; w.starts = Int[1]
    for f in (:retention_time, :tic, :center_mz, :isolation_width, :ce_field, :ce_ev, :ms_order, :cycle_idx, :frame_id, :im_scan, :window_group)
        setfield!(w, f, similar(getfield(w, f), 0))
    end
    w
end

function Base.close(w::SliceArrowWriter)
    flush_batch!(w)
    close(w.writer)
    nothing
end
