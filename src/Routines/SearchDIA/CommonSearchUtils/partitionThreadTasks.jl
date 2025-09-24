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


function partitionScansToThreads(spectra::AbstractArray,
                                rt::AbstractVector{Float32},
                                prec_mz::AbstractVector{Union{Missing, Float32}},
                                ms_order::AbstractVector{UInt8},
                                n_threads::Int,
                                tasks_per_thread::Int)
    total_peaks = sum(length.(spectra))
    n_tasks = n_threads*tasks_per_thread
    peaks_per_task = total_peaks÷(n_tasks)
    function round_float32_alt(x::Float32, decimals::Int)::Float32
        Float32(round(x; digits=decimals))
    end
    function round_float32_alt(x::Missing, decimals::Int)::Float32
        zero(Float32)
    end

    spectra_ids = collect([x for x in range(1, length(spectra)) if ms_order[x]==2])
    bin_start, bin_stop = 1, 1
    for i in range(1, length(spectra_ids))
        if rt[i] - rt[bin_start] > 1.0f0
            bin_stop = i - 1
            sort!(@view(spectra_ids[bin_start:bin_stop]), by = x->round_float32_alt(prec_mz[x],6))
            bin_start, bin_stop = i, i
        end
    end
    #sort final bin
    bin_stop = length(spectra_ids)
    sort!(@view(spectra_ids[bin_start:bin_stop]), by = x->round_float32_alt(prec_mz[x],6))

    spectra_count = length(spectra_ids)
    scans_per_thread = spectra_count÷n_threads + max(n_threads, 10) + 1
    thread_tasks = [[0, zeros(Int64, scans_per_thread)] for _ in range(1, n_threads)]
    for (thread_id, task) in enumerate(thread_tasks)
        task[1] = thread_id
        n = 1
        for i in range(1, length(task[2]))
            id = n_threads*(i - 1)*10 + thread_id*10 - 10 + 1
            for j in range(id, id + 10 - 1)
                if j <= length(spectra_ids)
                    task[2][n] = spectra_ids[j]
                    n += 1
                else
                    break
                end
            end
        end
    end
    return thread_tasks, total_peaks
end

function partitionScansToThreadsMS1(spectra::AbstractArray,
                                rt::AbstractVector{Float32},
                                prec_mz::AbstractVector{Union{Missing, Float32}},
                                ms_order::AbstractVector{UInt8},
                                n_threads::Int,
                                tasks_per_thread::Int)
    total_peaks = sum(length.(spectra))
    n_tasks = n_threads*tasks_per_thread
    peaks_per_task = total_peaks÷(n_tasks)
    spectra_ids = collect([x for x in range(1, length(spectra)) if ms_order[x]==1])
    bin_start, bin_stop = 1, 1
    #sort!(spectra_ids, by = x->)

    spectra_count = length(spectra_ids)
    scans_per_thread = spectra_count÷n_threads + max(n_threads, 10) + 1
    thread_tasks = [[0, zeros(Int64, scans_per_thread)] for _ in range(1, n_threads)]
    for (thread_id, task) in enumerate(thread_tasks)
        task[1] = thread_id
        n = 1
        for i in range(1, length(task[2]))
            id = n_threads*(i - 1)*10 + thread_id*10 - 10 + 1
            for j in range(id, id + 10 - 1)
                if j <= length(spectra_ids)
                    task[2][n] = spectra_ids[j]
                    n += 1
                else
                    break
                end
            end
        end
    end
    return thread_tasks, total_peaks
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
    n_tasks = n_threads * tasks_per_thread
    peaks_per_task = total_peaks ÷ n_tasks

    function round_float32_alt(x::Float32, decimals::Int)::Float32
        Float32(round(x; digits=decimals))
    end
    function round_float32_alt(x::Missing, decimals::Int)::Float32
        zero(Float32)
    end

    # Create virtual indices for MS2 scans only
    virtual_ms2_indices = Int64[]
    for i in eachindex(ms_order)
        if ms_order[i] == 2
            push!(virtual_ms2_indices, i)  # Virtual index, not actual scan index
        end
    end

    @debug_l2 "partitionScansToThreadsIndexed: Found $(length(virtual_ms2_indices)) MS2 scans out of $(length(ms_order)) total"

    if isempty(virtual_ms2_indices)
        # Return empty thread tasks if no MS2 scans
        empty_tasks = [[i, Int64[]] for i in 1:n_threads]
        return empty_tasks, total_peaks
    end

    # Sort virtual indices by RT, then by m/z within RT bins
    bin_start, bin_stop = 1, 1
    for i in 2:length(virtual_ms2_indices)
        virtual_idx = virtual_ms2_indices[i]
        start_virtual_idx = virtual_ms2_indices[bin_start]
        if rt[virtual_idx] - rt[start_virtual_idx] > 1.0f0
            bin_stop = i - 1
            # Sort by m/z within this RT bin
            sort!(@view(virtual_ms2_indices[bin_start:bin_stop]),
                  by = x -> round_float32_alt(prec_mz[x], 6))
            bin_start = i
        end
    end
    # Sort final bin
    bin_stop = length(virtual_ms2_indices)
    sort!(@view(virtual_ms2_indices[bin_start:bin_stop]),
          by = x -> round_float32_alt(prec_mz[x], 6))

    # Distribute virtual indices to threads
    spectra_count = length(virtual_ms2_indices)
    if spectra_count == 0
        empty_tasks = [[i, Int64[]] for i in 1:n_threads]
        return empty_tasks, total_peaks
    end

    scans_per_thread = max(1, spectra_count ÷ n_threads) + max(n_threads, 10) + 1
    thread_tasks = [[0, zeros(Int64, scans_per_thread)] for _ in 1:n_threads]

    for (thread_id, task) in enumerate(thread_tasks)
        task[1] = thread_id
        n = 1
        for i in 1:length(task[2])
            id = n_threads * (i - 1) * 10 + thread_id * 10 - 10 + 1
            for j in id:(id + 10 - 1)
                if j <= length(virtual_ms2_indices) && n <= length(task[2])
                    task[2][n] = virtual_ms2_indices[j]  # These are virtual indices
                    n += 1
                else
                    break
                end
            end
            if id > length(virtual_ms2_indices)
                break
            end
        end
        # Resize to actual used length
        task[2] = task[2][1:(n-1)]
    end

    @debug_l2 "partitionScansToThreadsIndexed: Distributed $(spectra_count) scans across $(n_threads) threads"
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

    # Create virtual indices for MS1 scans only
    virtual_ms1_indices = Int64[]
    for i in eachindex(ms_order)
        if ms_order[i] == 1
            push!(virtual_ms1_indices, i)  # Virtual index
        end
    end

    if isempty(virtual_ms1_indices)
        empty_tasks = [[i, Int64[]] for i in 1:n_threads]
        return empty_tasks, total_peaks
    end

    # Simple distribution for MS1 scans
    spectra_count = length(virtual_ms1_indices)
    scans_per_thread = max(1, spectra_count ÷ n_threads) + max(n_threads, 10) + 1
    thread_tasks = [[0, zeros(Int64, scans_per_thread)] for _ in 1:n_threads]

    for (thread_id, task) in enumerate(thread_tasks)
        task[1] = thread_id
        n = 1
        for i in 1:length(task[2])
            id = n_threads * (i - 1) * 10 + thread_id * 10 - 10 + 1
            for j in id:(id + 10 - 1)
                if j <= length(virtual_ms1_indices) && n <= length(task[2])
                    task[2][n] = virtual_ms1_indices[j]  # Virtual indices
                    n += 1
                else
                    break
                end
            end
            if id > length(virtual_ms1_indices)
                break
            end
        end
        # Resize to actual used length
        task[2] = task[2][1:(n-1)]
    end

    return thread_tasks, total_peaks
end
