using DataFrames
using Test

# Test: Does concurrent DataFrame column access cause issues?

@testset "DataFrame concurrent column access patterns" begin
    # Create test DataFrame
    test_df = DataFrame(precursor_idx = rand(UInt32(1):UInt32(100000), 50000))

    println("Testing Pattern 1: Direct access inside @spawn (BUGGY PATTERN)")
    @testset "Pattern 1: Direct DataFrame access in @spawn" begin
        n_tasks = 20
        iterations = 100000

        tasks = map(1:n_tasks) do task_id
            Threads.@spawn begin
                # Access DataFrame and create Set inside @spawn
                precursors_passing = Set(test_df[!, :precursor_idx])

                hits = 0
                for i in 1:iterations
                    test_idx = UInt32(i % 150000)
                    if test_idx ∈ precursors_passing
                        hits += 1
                    end
                end
                (task_id, hits, objectid(precursors_passing))
            end
        end

        results = fetch.(tasks)

        println("Pattern 1 Results:")
        println("  Task hits: ", [r[2] for r in results])
        println("  Unique Set object IDs: ", length(unique([r[3] for r in results])))

        @test all(r -> r[2] > 0, results)
        # Check if each thread got a different Set object
        @test length(unique([r[3] for r in results])) == n_tasks
    end

    println("\nTesting Pattern 2: Pre-extracted vector (PROPOSED FIX)")
    @testset "Pattern 2: Pre-extracted vector before @spawn" begin
        n_tasks = 20
        iterations = 100000

        # Extract column BEFORE spawning
        precursor_indices = test_df[!, :precursor_idx]
        println("  Pre-extracted vector object ID: ", objectid(precursor_indices))

        tasks = map(1:n_tasks) do task_id
            Threads.@spawn begin
                # Create Set from pre-extracted vector
                thread_local_precursors = Set(precursor_indices)

                hits = 0
                for i in 1:iterations
                    test_idx = UInt32(i % 150000)
                    if test_idx ∈ thread_local_precursors
                        hits += 1
                    end
                end
                (task_id, hits, objectid(thread_local_precursors), objectid(precursor_indices))
            end
        end

        results = fetch.(tasks)

        println("Pattern 2 Results:")
        println("  Task hits: ", [r[2] for r in results])
        println("  Unique Set object IDs: ", length(unique([r[3] for r in results])))
        println("  All threads see same Vector ID: ", allequal([r[4] for r in results]))

        @test all(r -> r[2] > 0, results)
        @test length(unique([r[3] for r in results])) == n_tasks
        @test allequal([r[4] for r in results])  # All see same Vector
    end

    println("\nTesting Pattern 3: Check if DataFrame access returns same object")
    @testset "Pattern 3: DataFrame column access consistency" begin
        # Does test_df[!, :precursor_idx] return the same object every time?
        vec1 = test_df[!, :precursor_idx]
        vec2 = test_df[!, :precursor_idx]

        println("  Same object returned: ", vec1 === vec2)
        println("  Object IDs: ", objectid(vec1), " vs ", objectid(vec2))

        @test vec1 === vec2  # Should be same object (no copy)
    end
end
