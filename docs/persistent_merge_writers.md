# Persistent merge writers: implementation and pending validation

This branch implements the first stage of the persistent merge writer plan.
Each ordinary merge, hierarchical staging merge, and output chunk uses one
`Arrow.Writer` instead of reopening and inspecting prior batches through
`Arrow.append`.

Tests and benchmarks have **not been run** for this implementation, at the
user's request. The package was loaded before that request; that does not
validate these edits. No branch has been pushed and no CI run has been started.
Static Julia parsing and `git diff --check` passed.

## Lifecycle and Windows behavior

`_MergeOutput` owns a unique temporary sibling file, its `IOStream`, and the
Arrow writer. Output uses the existing uncompressed Arrow stream format
(`file=false`) with an unbuffered message channel (`ntasks=0`). Each submitted
batch copies its filled column vectors because Arrow's consumer can still be
writing when the submission returns. Nested values from input Arrow tables
remain read-only; the merge never mutates their contents.

Arrow 2.x does not bind its consumer task to its message channel. The helper
uses `Base.bind(writer.msgs, writer.task)` so a consumer IO error closes the
channel and wakes blocked producers. This is a deliberate, small dependency on
Arrow's writer fields and needs regression coverage when upgrading Arrow.
The implementation does not alter Arrow's block arrays or encoding internals.

Finalization closes the writer, waits for its consumer through Arrow's close
implementation, and closes the underlying IO in a `finally` block. Only then
does the helper replace or rename the destination. Windows replacements use
the existing `safeRm` policy and verify that the old destination path is gone.
There are no new unsynchronized production `GC.gc()` calls.

Merge failures close the IO and discard the partial temporary file. A failure
to publish a completed stream retains the closed temporary file and reports
its absolute path for recovery. Cleanup warnings do not replace the original
failure. Replacements keep the old file until the new stream has finalized;
this requires temporary disk space for both versions. Replacement on Windows
is not atomic. `safeRm` can leave a renamed backup under its existing policy.
If a later chunk fails, earlier published chunks remain on disk; the error
reports their count and directory. They do not represent a complete merge.

Writer ownership and input memory mappings are different concerns. A reader
held by a caller can still block replacement/deletion; closing the writer
cannot release that reader. The existing whole-table input readers and
hierarchical staging-directory cleanup remain separate work.

## Behavior changes and memory limits

All-empty inputs now produce a readable zero-row stream with the input schema.
Batch sizes must be positive, fan-in must be at least two, and chunk byte
targets must be positive. Aliased input/output paths are rejected before
writing; standard chunk-name aliases are checked before hierarchical staging.
Each active chunk also checks its destination against its current inputs.

Chunk rotation uses the source file bytes-per-row estimate throughout a chunk.
It no longer resets that estimate by probing an output file while the Arrow
consumer could still be writing. Consequently, **chunk boundaries and file
sizes can differ from develop**, although row order and complete-group
boundaries should remain unchanged. `max_chunk_bytes` is a soft estimated
target. Complete groups larger than the target remain intact. The group key
must be a leading sort key for the complete-group guarantee.

The live batch payload and queued writes are bounded at a fixed batch size.
**This is not yet a constant-memory merge.** Arrow 2.8 still retains metadata
for every written record batch, including stream writers. Whole-file input
tables and growing chunk-reference lists also remain. The next stage requires
an Arrow stream-writer metadata fix and separate bounded reader/reference
work; this branch makes no claim that those problems are solved.

## Prepared tests

`test/utils/FileOperations/streaming/test_persistent_merge_writer.jl` covers
immediate buffer reuse, values and nullable schemas, list columns, exact and
partial batches, empty sources, mixed sort directions and ties, replacement,
whole-group chunk rotation, oversized/missing groups, multiple staging levels,
alias rejection, and invalid bounds.

Lifecycle tests keep the output session reachable while checking immediate
rename/delete, inject merge and IO failures, and check recovery after a failed
publication. A native Windows test opens a destination with `CreateFileW`
without file sharing and verifies that an external lock preserves the old
destination and leaves a closed recoverable temporary output. Test readers
use `Arrow.Table(read(path))` so they do not create output memory mappings.

The new file is registered in `test/runtests_part3_units.jl`. A focused Windows
job in `.github/workflows/tests.yml` is prepared to run it and the existing
basic merge tests with one and four Julia threads once the branch reaches CI.
Native Windows results are still required before claiming compatibility.

## Commands to run later

Run from this worktree:

```sh
cd /Users/nathanwamsley/Projects/Pioneer-persistent-merge-writers

julia --startup-file=no --threads=1 --project=. -e 'using Pioneer; include("test/utils/FileOperations/streaming/test_stream_sorted_merge_basic.jl"); include("test/utils/FileOperations/streaming/test_persistent_merge_writer.jl")'

julia --startup-file=no --threads=4 --project=. -e 'using Pioneer; include("test/utils/FileOperations/streaming/test_stream_sorted_merge_basic.jl"); include("test/utils/FileOperations/streaming/test_persistent_merge_writer.jl")'

julia --startup-file=no --threads=4 --project=. -e 'using Pioneer; include("test/utils/FileOperations/io/test_safe_file_ops.jl"); include("test/utils/FileOperations/io/test_arrow_operations_basic.jl"); include("test/utils/FileOperations/core/test_core_references_basic.jl")'

julia --startup-file=no --threads=4 --project=. scripts/benchmark_merge_writer.jl
```

The benchmark measures the actual new helper (including batch snapshots,
finalization, publication, and row-count verification) against the old append
strategy. Its allocations are cumulative, not peak memory. Run broader
FileOperations/quantification tests after the focused checks pass. Do not use
the earlier audit's writer-only timings as performance results for this branch.
