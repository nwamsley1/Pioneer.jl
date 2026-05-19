# Duplicate Precursor Attribution and Fragment Index Check

Updated: 2026-05-08 15:35:43 CDT

## Score Attribution Concern

The fused scan path gets a per-scan precursor range in `process_scans_fused!`:

- `src/Routines/SearchDIA/process_scans_fused.jl:89`
  - `prec_range = get_prec_range(prec_index, scan_idx)`
- `src/Routines/SearchDIA/process_scans_fused.jl:90`
  - `precs_vec = get_precursors(prec_index)`
- `src/Routines/SearchDIA/process_scans_fused.jl:98`
  - `precs_vec, prec_range` are passed into `run_fused!` as `precursors_passed, prec_range`.

Inside `run_fused!`, each candidate precursor ID is read here:

- `src/Routines/SearchDIA/CommonSearchUtils/fusedMatch.jl:645-646`
  - `for i in prec_range`
  - `prec_idx = precursors_passed[i]`

For normal `FusedStandard`, the column map key is just the precursor ID:

- `src/Routines/SearchDIA/CommonSearchUtils/fusedMatch.jl:95`
  - `match_column_id(::FusedSearchKind, prec_idx, iso_pass) = prec_idx`

When a precursor gets its first matched fragment, a new design-matrix column is allocated and the precursor-to-column map is written:

- `src/Routines/SearchDIA/CommonSearchUtils/fusedMatch.jl:725-729`
  - `col += UInt16(1)`
  - `this_col = col`
  - `id_to_col[match_column_id(kind, prec_idx, iso_pass)] = this_col`

Later, scoring loops over column-positioned `unscored_PSMs`, but it re-fetches the column index from `IDtoCOL`:

- `src/Routines/SearchDIA/PSMs/ScoredPSMs.jl:98-101`
  - `for i in range(1, n_vals)`
  - `precursor_idx = UInt32(unscored_PSMs[i].precursor_idx)`
  - `scores_idx = IDtoCOL[precursor_idx]`

If a scan candidate list contained the same precursor twice, for example:

```julia
precursors_passed[prec_range] == UInt32[42, 42]
```

and both entries matched, then `run_fused!` could create two columns but overwrite the map:

```julia
# first matching occurrence
this_col = 1
id_to_col[42] = 1
unscored_psms[1].precursor_idx = 42

# second matching occurrence
this_col = 2
id_to_col[42] = 2
unscored_psms[2].precursor_idx = 42
```

Then both scored rows would look up `IDtoCOL[42] == 2`. The first scored row would combine fragment-count/error features from `unscored_psms[1]` with deconvolution weight and spectral metrics from column 2.

Current assessment: this is a real attribution risk only if a per-scan candidate list contains duplicate precursor IDs. The main partitioned fragment-index path appears to satisfy uniqueness for the checked library below.

## Serialized Fragment Index Check

Checked file:

```text
/Users/nathanwamsley/Data/SPEC_LIBS/altimeter_3P_len7o40_ch2o3_mc1_OlsenExploris_mzsorted.poin/partitioned_fragment_index.jls
```

Inspection method:

- Deserialized the `LocalPartitionedFragmentIndex`.
- For each partition, mapped each `LocalFragment.local_id` through `partition.local_to_global`.
- Traversed RT bins correctly as `RT bin -> fragment m/z bin range -> fragment entry range`.
- Recorded each global precursor's unique `(partition_idx, rt_bin_idx)` locations.

Results:

```text
n_partitions: 157
total_rt_bins: 1571
total_fragment_entries: 54195771
unique_precursors_seen: 6799457
empty_partitions: 0
precursors_in_multiple_partitions: 0
precursors_in_multiple_partition_rt_bin_locations: 0
```

Conclusion: for this serialized index, every precursor appears in exactly one precursor partition and exactly one RT bin location. This means the specific duplicate-precursor attribution failure should not arise from this library's main partitioned fragment index.

Note: every RT bin contains repeated precursor IDs at the fragment-entry level because each precursor contributes multiple fragment ions. That is expected and is not the duplicate-candidate condition described above.

Repro script:

```text
investigation_notes/check_fragment_index_uniqueness.jl
```

