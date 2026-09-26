# Fragment index: UInt32 partition-local precursor IDs

Branch `feat/frag-index-uint32` (stacked on the timsTOF branch). Status: opt-in via
`library_params.frag_index_local_id_type = "UInt32"`; default stays `"UInt16"`.

## The constraint

The partitioned fragment index stores partition-local precursor IDs as UInt16, so a partition can hold at most
65,535 precursors. `build_partitioned_index_from_lib` first groups precursors into `prec_partition_width` bins of
precursor m/z, then halves any bin over the limit (by count, recursively). On a dense library the ID width, not
`prec_partition_width`, sets the layout:

| HYE timsTOF library, 10.66 M precursors (m/z 359-1034) | partitions | median width |
|---|---|---|
| 5 Da, UInt16 (94 of 135 bins split) | 233 | 2.5 Da |
| 10 Da, UInt16 (would split 67 of 68) | ~231 | ~2.5 Da |
| 10 Da, UInt32 | 68 | 10 Da |

A 25 Da diaPASEF window therefore overlaps ~10 UInt16 partitions. The search is partition-major, so each scan's
peaks are fetched (for `.tdfs`, decoded) and its mass windows rebuilt once per overlapping partition.

## Design

- New sibling types `LocalFragment32` (8 bytes), `LocalPartition32{T}`, `LocalPartitionedFragmentIndex32{T}` under
  abstract supertypes `AbstractLocalFragment`, `AbstractLocalPartition{T}`, `AbstractLocalPartitionedFragmentIndex{T}`.
- The original UInt16 types are unchanged. Julia `Serialization` records concrete type names and parameters, so
  making the existing types parametric would break every serialized `.poin` (including the prebuilt catalog).
- The search is written against the abstract types and `Counter{I,UInt8}`; `I = local_id_type(pfi)` is a
  compile-time constant for a concrete index, so each variant compiles to its own specialized code.
  `FragIndexScratch` carries a second set of UInt32 counters (`scratch_counters(s, I)`).
- Builder: `build_partitioned_index_from_lib(...; id_type)`; `max_local_precs(UInt32)` removes the split.

## No regression for UInt16

- Canonicalized LLVM IR of 20 hot-path method instances (`_run_thread`, `_score_partition_hinted!` for `.tdfs`
  and Arrow peaks x 3 mass error models, `queryFragmentHinted!`, `emit_candidates!` for both filters with/without
  the IM gate, the accumulator, `get_partition_range`) is identical before and after. Raw IR differs run to run
  in SSA numbering, GC-frame slots and JIT symbol counters, so those are canonicalized; two baseline runs match.
- The top-level `searchFragmentIndexPartitionMajorHinted` is type-stable for both index types.
- keap1 synthetic build: identical UInt16 index content hash from the old and new code. An existing library
  (E. coli timsTOF, built before the change) loads unchanged.
- Tests: existing 966 pass; `partitionedFragmentIndex32.jl` (52) and the UInt32 case in
  `test_build_determinism.jl` (15).

## Result: HYE timsTOF, 5 Da UInt16 vs 10 Da UInt32 (2026-09-25)

Two HYE libraries built from Koina with the paper HYE config + `im_model = alphapept_ccs`, seed 1844, differing
only in width / ID type. Searched two HYE 50 ng diaPASEF runs (Ultra 2), warm start (a warm-up search in the
same process first), two alternating reps per library; reps agree within ~1%.

| | 5 Da / UInt16 | 10 Da / UInt32 | change |
|---|---|---|---|
| fragment-index time, 2 files | 205.9 / 204.4 s | 125.3 / 125.9 s | -39% |
| Main Search | 344.6 / 338.2 s | 265.7 / 265.8 s | -22% |
| total search | 498.9 / 494.2 s | 422.5 / 422.9 s | -15% |
| Parameter Tuning | 21.0 / 20.6 s | 22.4 / 22.4 s | +1.5 s |
| allocated / peak RSS | 144 GiB / ~43 GiB | 141 GiB / ~45 GiB | ~same |
| precursors (A, B) | 85,573 / 87,821 | 86,822 / 87,220 | +0.4% total |
| protein groups (A, B) | 9,856 / 9,940 | 9,909 / 9,881 | ~same |

Fewer partitions per scan (3-4 instead of ~10) outweigh the larger counter (~230k slots) and 8-byte fragments.

Caveat: the two libraries are not byte-identical even with the same config and seed, because Koina predictions are
not bit-reproducible. Sequences, m/z, iRT, decoys and order match; 1/K0 differs by <= 6e-7; fragment intensities by
median 5e-7 relative, with 0.03% of fragments > 1e-3 and 5,106 rank swaps. The +/-1% per-file ID differences may be
this library noise; not yet isolated.

### Size

| HYE timsTOF | 5 Da / UInt16 | 10 Da / UInt32 |
|---|---|---|
| `partitioned_fragment_index.jls` on disk | 658.8 MB | 485.2 MB (-26%) |
| `presearch_partitioned_fragment_index.jls` on disk | 417.3 MB | 359.4 MB (-14%) |
| whole library on disk | 4.7 GB | 4.5 GB |
| main index in memory (fragment array) | 1.17 GB (0.34) | 1.20 GB (0.68) |
| presearch index in memory | 0.70 GB | 0.86 GB |

Fragment entries double (4 -> 8 bytes), but 3.4x fewer partitions means fewer fragment m/z bins (bins are built per
partition and RT bin), which nearly offsets it in the main index; the files compress better (zero upper ID bytes).

## Result: E. coli timsTOF, 10 Da UInt16 vs 10 Da UInt32 (2026-09-26)

Same config as the HYE pair but the E. coli canonical FASTA: 1,016,212 precursors, at most 20,754 per 10 Da bin,
so UInt16 does not split and both indexes have the same 68 partitions (74.8 vs 73.9 MB). E. coli 50 ng 5 min
diaPASEF (Ultra 2, TimsSlices 0.3.0), warm start, two alternating reps:

| | 10 Da / UInt16 | 10 Da / UInt32 |
|---|---|---|
| fragment-index time | 11.0 / 10.5 s | 10.8 / 10.6 s |
| Main Search | 25.3 / 24.7 s | 25.3 / 24.9 s |
| total search | 45.3 / 44.7 s | 45.2 / 45.0 s |
| allocated / peak RSS | 18.1 GiB / 11.1 GiB | 18.2 GiB / 11.1 GiB |
| precursors / protein groups | 16,416 / 1,781 | 16,357 / 1,788 |

- With the same partition layout the ID width has no measurable cost at ~21k precursors per partition (the counter
  stays cache-resident either way). The HYE gain comes from wider partitions, which UInt32 makes possible.
- The two libraries can only differ in search results through Koina's prediction noise: -0.4% precursors / +0.4%
  protein groups. That is the scale of the HYE differences, so those are most likely library noise, not the index.

## Open: how to expose it

Proposed: `frag_index_local_id_type = "auto"`. After precursor m/z are known, use UInt32 only when some
`prec_partition_width` bin would exceed 65,535 (the ID width is binding), else UInt16. Not yet tested on narrow-window
(Thermo / Astral) data, where a wider effective partition also admits more out-of-window candidates.
