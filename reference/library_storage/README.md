# Compact library storage

Reference codecs for reducing library size while reconstructing the existing
native fragment and index arrays for search. Based on develop `226abc0b7`.

## Storage changes

| Component | Encoding |
| --- | --- |
| Index entries | Separate UInt16/UInt32 precursor-ID and UInt8 rank-mask arrays, removing record padding. |
| Fragment-m/z bins | Store lows, range ends, and exact high-bound bit deltas. Reconstruct range starts and SIMD padding during loading. |
| Detailed fragments | Separate mass, intensity, and annotation columns. Derive precursor ID from its range and charge from the precursor table. |
| Fragment offsets | Store counts using the smallest sufficient unsigned type; reconstruct UInt64 CSR offsets. |
| Precursor table | Dictionary-encode strings and derive `koina_sequence` when prediction tools need it. |

High-bound deltas are `bits(high) - bits(low)`. Each partition uses UInt8,
UInt16, or UInt32 according to its largest delta. Reconstruction is exact.
All eight rank-mask bits, fragment order, Float16 intensities, and Float32
spline coefficients are preserved.

Detailed-fragment annotations occupy 20 bits in UInt16 + UInt8 columns:

| Bits | Field |
| --- | --- |
| 0–3 | y, b, precursor-ion, isotope flags |
| 4–5 | Fragment charge |
| 6–8 | Sulfur count |
| 9–14 | Ion position |
| 15–19 | Intensity rank |

See [the byte-layout diagram](index-byte-layout.svg) and the implementations in
[index.jl](index.jl), [details.jl](details.jl), and [precursors.jl](precursors.jl).

## Integration

The codecs are in `CompactLibraryStorage.jl`; production build/load wiring
remains to be added.

1. Encode every main/presearch index in `write_fragment_indexes`; decode once
   when loading the selected indexes.
2. Write detailed columns and counts after the existing per-precursor m/z sort.
   Reconstruct native fragments and CSR before constructing the lookup object.
3. Finalize precursor dictionaries after intensity, RT, and IM prediction.
   Keep them mapped during search. The five string getters in
   `src/structs/SpectralLibrary/precursors.jl` already accept dictionary columns.
4. Update index rebuilding, DIA-NN decoy processing, fixture/precompile helpers,
   and tools that read the separate CSR file. Counts now live inside
   `detailed_fragments.jls`. Prediction tools can regenerate Koina strings.
