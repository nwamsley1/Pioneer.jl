# Run explicitly after functional tests pass; never included in the test suite.
using Pioneer
using Arrow, DataFrames, Tables, Statistics, Printf

const MERGE_BENCH_BATCH = DataFrame([Symbol("x$i") => Int32.(1:1000) for i in 1:12])

function legacy_merge_batches(path, n_batches)
    open(path, "w") do io
        Arrow.write(io, MERGE_BENCH_BATCH; file=false)
    end
    for _ in 2:n_batches
        Arrow.append(path, MERGE_BENCH_BATCH)
    end
end

function persistent_merge_batches(path, n_batches)
    Pioneer._with_merge_output(path) do output
        for _ in 1:n_batches
            Pioneer._write_batch_typed(output, MERGE_BENCH_BATCH, nrow(MERGE_BENCH_BATCH))
        end
    end
end

function measure_merge_writer(f, n_batches)
    mktempdir() do dir
        try
            path = joinpath(dir, "merged.arrow")
            sample = @timed f(path, n_batches)
            # Byte-backed readers avoid retaining output mappings during cleanup.
            rows = sum(length(Tables.getcolumn(table, 1)) for table in Arrow.Stream(read(path)))
            rows == n_batches * nrow(MERGE_BENCH_BATCH) || error("Incorrect output row count")
            return (seconds=sample.time, allocated_bytes=sample.bytes, file_bytes=filesize(path))
        finally
            # The legacy Arrow.append baseline can leave unreachable mappings.
            # Release those outside the measured interval before Windows cleanup.
            if Sys.iswindows()
                lock(Pioneer._WINDOWS_DELETE_GC_LOCK) do
                    GC.gc(true)
                end
            end
        end
    end
end

function main()
    for f in (legacy_merge_batches, persistent_merge_batches)
        measure_merge_writer(f, 4)
    end
    println("# Julia $(VERSION); Arrow $(pkgversion(Arrow)); median of three runs")
    println("method\tbatches\tseconds\tallocated_bytes\tfile_bytes")
    for n_batches in (16, 32, 64, 128, 256, 512)
        for (label, f) in (("legacy_append", legacy_merge_batches), ("persistent_output", persistent_merge_batches))
            samples = [measure_merge_writer(f, n_batches) for _ in 1:3]
            @printf("%s\t%d\t%.6f\t%d\t%d\n", label, n_batches,
                    median(getproperty.(samples, :seconds)),
                    Int(median(getproperty.(samples, :allocated_bytes))),
                    Int(median(getproperty.(samples, :file_bytes))))
        end
    end
end

main()
