using Test
using Pioneer: PartitionBarrier, wait_partition!, abort_partition!

"""On timeout, release waiters and avoid fetching any unfinished task."""
function join_partition_test_tasks(tasks, barrier)
    status = timedwait(() -> all(istaskdone, tasks), 10.0; pollint=0.01)
    if status != :ok
        cleanup = Threads.@spawn abort_partition!(barrier)
        timedwait(() -> istaskdone(cleanup) && all(istaskdone, tasks),
                  2.0; pollint=0.01)
    end
    @test status == :ok
    all(istaskdone, tasks) || return nothing
    return fetch.(tasks)
end

@testset "Partition barrier generations" begin
    @test_throws ArgumentError PartitionBarrier(0)
    @test_throws ArgumentError PartitionBarrier(-1)

    single = PartitionBarrier(1)
    for _ in 1:6
        @test wait_partition!(single) === nothing
    end
    abort_partition!(single)
    @test_throws ErrorException wait_partition!(single)

    # Oversubscription must work even with a single Julia runtime thread.
    nworkers = Threads.nthreads() + 3
    barrier = PartitionBarrier(nworkers)
    work_counts = [2nworkers + 3, 0, 1, nworkers - 1, 3nworkers + 1, 0, 2]
    events = Tuple{Symbol,Int,Int}[]
    event_lock = ReentrantLock()
    output = [Tuple{Int,Int}[] for _ in 1:nworkers]
    tasks = map(1:nworkers) do worker
        Threads.@spawn try
            for partition in eachindex(work_counts)
                wait_partition!(barrier)
                lock(event_lock) do
                    push!(events, (:start, worker, partition))
                end
                # Some workers have no scans, and two whole partitions are empty.
                for scan in worker:nworkers:work_counts[partition]
                    for _ in 1:mod(worker + scan, 4)
                        yield()
                    end
                    push!(output[worker], (partition, scan))
                end
                lock(event_lock) do
                    push!(events, (:finish, worker, partition))
                end
            end
            return nothing
        catch
            abort_partition!(barrier)
            rethrow()
        end
    end
    results = join_partition_test_tasks(tasks, barrier)
    if results !== nothing
        @test all(isnothing, results)
        finished = zeros(Int, nworkers)
        for (kind, worker, partition) in events
            if kind == :start
                # No next-partition work may start before all previous work ends.
                @test all(>=(partition - 1), finished)
            else
                finished[worker] = partition
            end
        end
        @test all(==(length(work_counts)), finished)
        for worker in 1:nworkers
            expected = [(partition, scan) for partition in eachindex(work_counts)
                        for scan in worker:nworkers:work_counts[partition]]
            @test output[worker] == expected
        end
    end
end

@testset "Partition worker failure releases peers" begin
    nworkers = Threads.nthreads() + 3
    barrier = PartitionBarrier(nworkers)
    entered = Threads.Atomic{Int}(0)
    failure = ErrorException("intentional partition worker failure")
    tasks = map(1:nworkers) do worker
        Threads.@spawn try
            if worker == nworkers
                # This participant never arrives, so peers cannot pass normally.
                ready = timedwait(() -> entered[] == nworkers - 1,
                                  5.0; pollint=0.01)
                ready == :ok || error("Partition test peers did not enter")
                throw(failure)
            end
            Threads.atomic_add!(entered, 1)
            wait_partition!(barrier)
            return :unexpected_success
        catch error
            abort_partition!(barrier, error)
            return error
        end
    end
    results = join_partition_test_tasks(tasks, barrier)
    if results !== nothing
        @test results[end] === failure
        @test all(result -> result === failure, results)

        # Later aborts must preserve the first cause for future arrivals as well.
        abort_partition!(barrier, ErrorException("later partition worker failure"))
        late_task = Threads.@spawn try
            wait_partition!(barrier)
            :unexpected_success
        catch error
            error
        end
        late_results = join_partition_test_tasks([late_task], barrier)
        if late_results !== nothing
            @test only(late_results) === failure
        end
    end
end
