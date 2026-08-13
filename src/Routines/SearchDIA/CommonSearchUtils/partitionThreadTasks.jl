# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
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

struct ThreadTask
    range::UnitRange{Int64}
    thread_id::Int64
end

getRange(tt::ThreadTask) = tt.range
getThreadID(tt::ThreadTask) = tt.thread_id

function partitionThreadTasks(n_tasks::Int, tasks_per_thread::Int, n_threads::Int)
    chunk_size = max(1, n_tasks ÷ (tasks_per_thread*n_threads))
    return partition(1:n_tasks, chunk_size)
end

#############################################################################
# Internal helpers
#############################################################################

@inline _round_float32(x::Float32, decimals::Int)::Float32 = Float32(round(x; digits=decimals))
@inline _round_float32(::Missing, decimals::Int)::Float32 = zero(Float32)

"""
Sort scan indices in-place by precursor m/z within 1-minute RT bins.
Used for MS2 scans to improve cache locality during thread processing.
"""
function _sort_scans_by_mz_in_rt_bins!(scan_indices::Vector{Int64},
                                        rt::AbstractVector{Float32},
                                        prec_mz::AbstractVector{Union{Missing, Float32}})
    length(scan_indices) <= 1 && return
    bin_start = 1
    for i in 2:length(scan_indices)
        if rt[scan_indices[i]] - rt[scan_indices[bin_start]] > 1.0f0
            sort!(@view(scan_indices[bin_start:i-1]),
                  by = x -> _round_float32(prec_mz[x], 6))
            bin_start = i
        end
    end
    sort!(@view(scan_indices[bin_start:end]),
          by = x -> _round_float32(prec_mz[x], 6))
end

"""
Distribute scan indices to threads using round-robin batches of 10.
Returns Vector of [thread_id, scan_ids_vector] pairs.
"""
function _distribute_scans_to_threads(scan_indices::Vector{Int64}, n_threads::Int)
    batch_sz = 10
    spectra_count = length(scan_indices)
    if spectra_count == 0
        return [[i, Int64[]] for i in 1:n_threads]
    end

    scans_per_thread = max(1, spectra_count ÷ n_threads) + max(n_threads, batch_sz) + 1
    thread_tasks = [[0, zeros(Int64, scans_per_thread)] for _ in 1:n_threads]

    for (thread_id, task) in enumerate(thread_tasks)
        task[1] = thread_id
        n = 1
        for i in 1:length(task[2])
            id = n_threads * (i - 1) * batch_sz + thread_id * batch_sz - batch_sz + 1
            for j in id:(id + batch_sz - 1)
                if j <= spectra_count && n <= length(task[2])
                    task[2][n] = scan_indices[j]
                    n += 1
                else
                    break
                end
            end
            if id > spectra_count
                break
            end
        end
        task[2] = task[2][1:(n-1)]
    end

    return thread_tasks
end

#############################################################################
# Public API
#############################################################################

function partitionScansToThreads(spectra::AbstractArray,
                                rt::AbstractVector{Float32},
                                prec_mz::AbstractVector{Union{Missing, Float32}},
                                ms_order::AbstractVector{UInt8},
                                n_threads::Int,
                                tasks_per_thread::Int)
    total_peaks = sum(length.(spectra))
    scan_indices = Int64[x for x in 1:length(spectra) if ms_order[x] == 2]
    _sort_scans_by_mz_in_rt_bins!(scan_indices, rt, prec_mz)
    return _distribute_scans_to_threads(scan_indices, n_threads), total_peaks
end

function partitionScansToThreadsMS1(spectra::AbstractArray,
                                rt::AbstractVector{Float32},
                                prec_mz::AbstractVector{Union{Missing, Float32}},
                                ms_order::AbstractVector{UInt8},
                                n_threads::Int,
                                tasks_per_thread::Int)
    total_peaks = sum(length.(spectra))
    scan_indices = Int64[x for x in 1:length(spectra) if ms_order[x] == 1]
    return _distribute_scans_to_threads(scan_indices, n_threads), total_peaks
end

"""
Specialized partitioning for IndexedMassSpecData MS2 scans.
Returns virtual indices (1, 2, 3...) properly distributed to threads based on
the underlying scan properties (RT, m/z) but respecting the IndexedMassSpecData view.
"""
function partitionScansToThreadsIndexed(
    spectra::AbstractArray,
    rt::AbstractVector{Float32},
    prec_mz::AbstractVector{Union{Missing, Float32}},
    ms_order::AbstractVector{UInt8},
    actual_scan_indices::AbstractVector{Int32},
    n_threads::Int,
    tasks_per_thread::Int
)
    total_peaks = sum(length.(spectra))
    scan_indices = Int64[i for i in eachindex(ms_order) if ms_order[i] == 2]
    @debug_l2 "partitionScansToThreadsIndexed: Found $(length(scan_indices)) MS2 scans out of $(length(ms_order)) total"
    _sort_scans_by_mz_in_rt_bins!(scan_indices, rt, prec_mz)
    thread_tasks = _distribute_scans_to_threads(scan_indices, n_threads)
    @debug_l2 "partitionScansToThreadsIndexed: Distributed $(length(scan_indices)) scans across $(n_threads) threads"
    return thread_tasks, total_peaks
end

"""
Specialized partitioning for IndexedMassSpecData MS1 scans.
"""
function partitionScansToThreadsMS1Indexed(
    spectra::AbstractArray,
    rt::AbstractVector{Float32},
    prec_mz::AbstractVector{Union{Missing, Float32}},
    ms_order::AbstractVector{UInt8},
    actual_scan_indices::AbstractVector{Int32},
    n_threads::Int,
    tasks_per_thread::Int
)
    total_peaks = sum(length.(spectra))
    scan_indices = Int64[i for i in eachindex(ms_order) if ms_order[i] == 1]
    return _distribute_scans_to_threads(scan_indices, n_threads), total_peaks
end
