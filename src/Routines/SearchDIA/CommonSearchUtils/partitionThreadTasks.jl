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
    scans_per_thread = spectra_count÷n_threads + n_threads
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