# Compact library storage — handoff for Nathan

This branch is based on develop `226abc0b7b15aa5e6ab2c563837228376ad09f9a`.
It contains an executable reference for the accepted storage changes, focused
tests, and measured results. Nathan can reimplement these codecs in his branches
without carrying over the research harness or its experimental loader overrides.

The reference readers accept only the new payload types and precursor schema.
There is no legacy decoding, schema negotiation, or runtime method replacement.
The normal build/load entry points are intentionally left for Nathan to wire
within his branches. The only production-source adjustment here removes five
concrete `Arrow.List` return annotations so native precursor accessors can retain
mapped Arrow dictionaries. Numeric getters and search kernels are unchanged.

## Files

| File | Responsibility |
| --- | --- |
| `CompactLibraryStorage.jl` | Reference module and shared validation helpers |
| `index.jl` | Disk-only index types, exact encoding, native reconstruction |
| `details.jl` | Scalar/spline columns and count-based CSR reconstruction |
| `precursors.jl` | String dictionaries and on-demand prediction-input strings |
| `verify_reference.jl` | Offline verification against saved develop libraries |
| `index-byte-layout.svg` | Byte and bit allocation diagram |
| `measurements.csv` | Actual sizes from the cleaned reference implementation |
| `validation.json` | Test and earlier search-validation evidence |

## Format decisions

### Indexes

Each `StoredPartition` contains native RT bins, skip hints, and the local/global
precursor mapping. Fragment entries become separate UInt16/UInt32 ID and UInt8
rank-mask columns. All eight mask bits are retained.

Fragment-m/z bins store Float32 lows, UInt32 range ends, and unsigned high-bound
bit deltas. The delta is `bits(high) - bits(low)`, not a rounded numerical mass
difference. Its type is selected independently for each partition: UInt8 through
255, UInt16 through 65,535, otherwise UInt32. The first entry of each bin is the
previous bin's end plus one. SIMD infinity padding is recreated from its UInt8
count. Precursor partition bounds remain exact Float32 pairs.

`decode_index` returns the same native index and partition types the current
search uses. There are no packed runtime entries or extra search indirections.

### Detailed fragments and CSR

`StoredDetails` stores fragment counts, Float32 masses, the original intensity
type, and two annotation columns (UInt16 + UInt8). Scalar intensities remain
Float16; spline intensities remain `NTuple{N,Float32}`. No prediction is rounded,
recomputed, or deduplicated.

The annotation's low 20 bits are:

| Bits | Meaning |
| --- | --- |
| 0–3 | y, b, precursor-ion, isotope flags |
| 4–5 | Fragment charge |
| 6–8 | Sulfur count |
| 9–14 | Ion position |
| 15–19 | Intensity rank |

Precursor ID follows from the CSR range and precursor charge comes from the
precursor table. Fragment counts use the smallest sufficient unsigned type;
they reconstruct the original UInt64 offsets, including empty ranges. Counts
reside in `detailed_fragments.jls`, replacing the separate
`precursor_to_fragment_indices.jls` file. Fragment order is preserved exactly.

### Precursors

String columns use Arrow dictionary encoding. Numeric and list columns retain
their existing values and types. `koina_sequence` is not stored. If supplied to
the writer, it must be exactly derivable from sequence and structural mods;
otherwise the writer fails instead of silently discarding information.
`koina_sequence(sequence, mods)` is a cold helper for prediction tools. It is
not used during search; model-specific IM formatting remains in prediction code.

`read_precursors` returns mapped dictionary columns, not an expanded copy of all
strings. Standard precursor CV folds and text IDs are checked for equality.
The new precursor schema has one fixed metadata identifier; it has no old-schema
fallback or optional stored Koina column.

## Integration points for the port

1. Include the disk-only types after the native fragment/index definitions.
2. In `write_fragment_indexes`, serialize `encode_index(index)` for every main
   and presearch index width. Decode once when loading each selected index.
3. After the builder's existing detailed-fragment m/z sort, write one encoded
   details payload. Read it with precursor charges to reconstruct native
   fragments and UInt64 CSR before constructing the normal lookup object.
4. Finalize the dictionary precursor table **after** intensity, RT, and IM
   prediction steps. Prediction intermediates can still contain Koina strings.
5. Update all final-library readers/writers together: normal loading, index
   rebuilding, DIA-NN decoy postprocessing, fixture/precompile helpers, and any
   tool that currently reads the separate CSR file. Remove the old storage paths
   rather than maintaining a dual reader. Prediction tools that reopen a final
   library must explicitly regenerate their prediction-input strings.
6. Keep proteins, spline knots, output schemas, IDs, CV assignment, ranking,
   and runtime fragment/index types unchanged. Do not infer rank from row order:
   detailed fragments remain sorted by m/z during search.

The normal entry points in this branch do not consume the reference payloads.
This is a porting handoff, not a deployable library-format migration. A complete
integration must update those entry points and rebuild its test libraries.

## Validation and reproduction

From the repository root:

```sh
julia --compiled-modules=existing --project test/UnitTests/test_compact_library_storage.jl
julia --compiled-modules=existing --project reference/library_storage/verify_reference.jl BASELINE_LIBRARY EMPTY_OUTPUT_DIRECTORY
```

The offline verifier reads a develop benchmark library, writes the reference
payloads to a new directory, and checks every precursor column other than the
verified-derived Koina string, all fragment values/annotations, all CSR offsets,
all index entries/bin metadata, CV folds, and precursor text IDs. It does not
contact prediction services or modify its source library. Its baseline read is
test code, not a legacy path in any new reader.

Tests exercise both local-ID widths, one/two/four-byte delta boundaries, full
rank masks, empty indexes and precursor ranges, adaptive count widths, annotation
extremes, signed zero and NaN payload preservation, scalar and spline file
round trips, dictionaries, derived strings, and invalid inputs. The test is also
included in the repository's third test partition.

The cleaned reference is rechecked on the 18 human narrow-window cases and the
54,410-precursor E. coli spline reference. `measurements.csv` measures all four
indexes, details **plus CSR**, and the precursor table. Fixed protein/config
files are excluded. Unlike the earlier adaptive benchmark, this implementation
always writes the dictionary precursor schema, even for tiny controls where
the Arrow dictionary overhead makes it larger.

Earlier prototype evidence: real-DIA candidate masks matched for all matrix
cases, and full develop/current/experimental E. coli searches matched every
field in eight intermediate/final tables. Single cold-JIT runs took about
297/296/298 seconds respectively; these are correctness checks and rough timing,
not a statistically powered speedup claim. That full-search test used the
earlier prototype loader, not this cleaned, unwired reference. The port should
repeat full-search comparisons after integration.

The human sample uses 64 full proteins and 500–501, 700–701, and 900–901 m/z
precursor windows, with tryptic/semi/nonenzymatic digestion and 0–3 STY phospho
modifications. Predictions are actual Prosit outputs. Narrow unmodified tryptic
controls contain only 16–34 precursors; their fixed overhead is not representative
of large libraries. Neither RAM reduction nor billion-precursor scaling is
claimed: fragments and indexes expand to native arrays for computation.

The cleaned reference saves 42.9–43.8% across the three components in the
nonenzymatic/phospho cases and 29.1% in the wider spline reference. The smallest
16-precursor tryptic control grows 0.3% with the fixed new schema; the benchmark
CSV retains that result rather than selecting a legacy passthrough format.

## Deliberately excluded

Shared fragment templates, mass dictionaries, packed runtime entries, protein-
coordinate sequence reconstruction, generic compression changes, prediction
cache work, compatibility shims, and runtime loader overrides are not part of
this handoff. The research sources, predictions, full logs, and generated
libraries are preserved locally outside the repository.
