struct ThreadTask
    range::UnitRange{Int64}
    thread_id::Int64
end

getRange(tt::ThreadTask) = tt.range
getThreadID(tt::ThreadTask) = tt.thread_id

function partitionThreadTasks(n_tasks::Int, tasks_per_thread::Int, n_threads::Int)
    chunk_size = max(1, n_tasks ÷ (tasks_per_thread*n_threads))
    return partition(1:n_chroms, chunk_size)
end

function partitionScansToThreads(spectra::AbstractArray{Vector{Union{Missing, Float32}}},
                                n_threads::Int,
                                tasks_per_thread::Int)
    total_peaks = sum(length.(spectra))
    n_tasks = n_threads*tasks_per_thread
    peaks_per_task = total_peaks÷(n_tasks)
    thread_tasks = Vector{ThreadTask}(undef, n_tasks)
    
    n = 0
    start = 1
    thread_id = 1
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
            start = i + 1
            n = 0
        end
    end
    return thread_tasks, total_peaks
end