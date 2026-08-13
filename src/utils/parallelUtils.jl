"""
    parallel_foreach!(f::Function, n::Int;
                      tasks_per_thread::Int=5,
                      n_threads::Int=Threads.nthreads())

Partition `1:n` into chunks and execute `f(chunk)` in parallel using
`Threads.@spawn`. Blocks until all tasks complete.

# Example
```julia
results = Vector{Float64}(undef, n)
parallel_foreach!(n) do chunk
    for i in chunk
        results[i] = expensive_work(data[i])
    end
end
```
"""
function parallel_foreach!(f::F, n::Int;
                           tasks_per_thread::Int=5,
                           n_threads::Int=Threads.nthreads()) where {F<:Function}
    chunk_size = max(1, n ÷ (tasks_per_thread * n_threads))
    tasks = map(partition(1:n, chunk_size)) do chunk
        Threads.@spawn f(chunk)
    end
    fetch.(tasks)
    return nothing
end
