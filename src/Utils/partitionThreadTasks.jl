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

function partitionScansToThreads(spectra::Arrow.List{Union{Missing, SubArray{Union{Missing, Float32}, 1, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}, Tuple{UnitRange{Int64}}, true}}, Int64, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}},
                                n_threads::Int,
                                tasks_per_thread::Int)
    total_peaks = sum(length.(spectra))
    n_tasks = n_threads*tasks_per_thread
    peaks_per_task = total_peaks÷(n_tasks)
    thread_tasks = Vector{ThreadTask}(undef, n_tasks)
    
    n = 0
    start = 1
    thread_id = 1
    task = 1
    for (i, spectrum) in enumerate(spectra)
        n += length(spectrum)
        if (n > peaks_per_task) | ((i + 1) == length(spectra))

            thread_tasks[task] = ThreadTask(
                                                UnitRange(start, i), 
                                                thread_id
                                            )
            thread_id += 1
            if thread_id > n_threads
                thread_id = 1
            end
            task += 1
            start = i + 1
            n = 0
        end
    end
    return thread_tasks, total_peaks
end

function partitionScansToThreads(spectra::Arrow.List{Union{Missing, SubArray{Union{Missing, Float32}, 1, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}, Tuple{UnitRange{Int64}}, true}}, Int64, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}},
                                rt::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}},
                                prec_mz::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}},
                                ms_order::Arrow.Primitive{Union{Missing, Int32}, Vector{Int32}},
                                n_threads::Int,
                                tasks_per_thread::Int)
    total_peaks = sum(length.(spectra))
    n_tasks = n_threads*tasks_per_thread
    peaks_per_task = total_peaks÷(n_tasks)

    spectra_ids = collect([x for x in range(1, length(spectra)) if ms_order[x]==2])
    bin_start, bin_stop = 1, 1
    for i in range(1, length(spectra_ids))
        if rt[i] - rt[bin_start] > 1.0f0
            bin_stop = i - 1
            sort!(@view(spectra_ids[bin_start:bin_stop]), by = x->Int64(round(prec_mz[x])))
            bin_start, bin_stop = i, i
        end
    end
    bin_stop = length(spectra_ids)
    sort!(@view(spectra_ids[bin_start:bin_stop]), by = x->Int64(round(prec_mz[x])))
    spectra_count = length(spectra_ids)
    scans_per_thread = spectra_count÷n_threads + n_threads
    thread_tasks = [[0, zeros(Int64, scans_per_thread)] for _ in range(1, n_threads)]
    for (thread_id, task) in enumerate(thread_tasks)
        task[1] = thread_id
        for i in range(1, length(task[2]))
            id = n_threads*(i - 1) - thread_id + n_threads + 1
            if id <= length(spectra_ids)
            task[2][i] = spectra_ids[id]
            end
        end
    end
    return thread_tasks, total_peaks
end