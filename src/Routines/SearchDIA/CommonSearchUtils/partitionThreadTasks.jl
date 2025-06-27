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
    scans_per_thread = spectra_count÷n_threads + n_threads + 1
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
    scans_per_thread = spectra_count÷n_threads + n_threads + 1
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