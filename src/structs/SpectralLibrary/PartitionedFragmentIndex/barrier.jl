"""
    PartitionBarrier(n_workers::Int)

Reusable rendezvous for `n_workers` logical index-search workers. Waiting releases
the condition lock and yields the task, so workers may outnumber runtime threads.
Every worker must participate at each partition, including partitions with no work.
"""
mutable struct PartitionBarrier
    condition::Threads.Condition
    parties::Int
    arrived::Int
    generation::Int
    broken::Bool
    cause::Any

    function PartitionBarrier(n_workers::Int)
        n_workers > 0 || throw(ArgumentError("Partition barrier needs at least one worker"))
        new(Threads.Condition(), n_workers, 0, 0, false, nothing)
    end
end

"""Wait for every worker to reach the next partition; throw if a peer failed."""
function wait_partition!(barrier::PartitionBarrier)
    lock(barrier.condition) do
        barrier.broken && throw(barrier.cause)
        generation = barrier.generation
        barrier.arrived += 1
        if barrier.arrived == barrier.parties
            barrier.arrived = 0
            barrier.generation += 1
            notify(barrier.condition; all=true)
        else
            while barrier.generation == generation && !barrier.broken
                wait(barrier.condition)
            end
            barrier.broken && throw(barrier.cause)
        end
    end
    return nothing
end

"""Release waiting peers after a worker fails; the original worker rethrows its error."""
function abort_partition!(barrier::PartitionBarrier,
        cause=ErrorException("Partition barrier aborted after a worker failed"))
    lock(barrier.condition) do
        if !barrier.broken
            barrier.cause = cause
            barrier.broken = true
        end
        notify(barrier.condition; all=true)
    end
    return nothing
end

# Calibration searches and single-worker calls do not need synchronization.
wait_partition!(::Nothing) = nothing
abort_partition!(::Nothing, cause) = nothing
