# Parameter Configuration

Pioneer uses JSON configuration files to control analysis. This guide documents the parameters for both `SearchDIA` and `BuildSpecLib`.

The schema below mirrors [`assets/example_config/defaultSearchParams.json`](https://github.com/nwamsley1/Pioneer.jl/blob/main/assets/example_config/defaultSearchParams.json) and [`defaultBuildLibParams.json`](https://github.com/nwamsley1/Pioneer.jl/blob/main/assets/example_config/defaultBuildLibParams.json). Any field not listed here has been hardcoded in source.

A simplified config (`defaultSearchParamsSimplified.json`) is also available for users who only need to specify paths plus a few high-level toggles; everything else falls back to the defaults documented below.

## SearchDIA Configuration

### Frequently Modified Parameters

Most parameters can stay at their defaults. The few worth tuning per-experiment:

* **`global.q_value_threshold`** — final FDR cutoff for output (default `0.01`). Loosen to `0.05` for exploratory work; tighten only with a strong rationale.

* **`search.n_isotopes`** — number of fragment isotopes to consider (default `2`, includes M and M+1). If using a non-Altimeter library (Prosit, UniSpec) where M+1 intensities are not modeled, set to `1`.

* **`search.global_prescore_qvalue_threshold`** — intermediate FDR cutoff between MainSearch's first LightGBM (15 features) and ScoringSearch's second LightGBM (20 features). Default `0.02`. Looser values (e.g. `0.05`) carry more borderline candidates into the second model — better recall at the cost of more compute. Tighter values speed up ScoringSearch.

* **`acquisition.nce`** — initial collision-energy guess for the pre-search before NCE tuning runs (default `26`, suitable for Thermo Orbitrap/Astral instruments). Inspect the QC plots after a run; if the auto-fitted NCE is far from this guess, re-running with a closer initial guess can improve results.

* **`optimization.machine_learning.max_psm_memory_mb`** — memory budget for in-memory LightGBM training (default `2000` = 2 GB). Workstations with more RAM benefit from raising this; lower it on memory-constrained machines to force the out-of-memory path earlier.

* **`maxLFQ.run_to_run_normalization`** — apply the run-to-run median-spline normalization to peak areas (default `false`). Turn on for multi-run experiments where systemic between-run intensity bias is expected.

### Global

| Parameter | Type | Description |
|-----------|------|-------------|
| `global.q_value_threshold` | Float | Final FDR threshold applied to ScoringSearch output (default: `0.01`). |

### Search

| Parameter | Type | Description |
|-----------|------|-------------|
| `search.n_isotopes` | Int | Number of fragment isotopes to consider in matching (default: `2`). Set `1` for non-Altimeter libraries. |
| `search.global_prescore_qvalue_threshold` | Float | FDR cutoff between the first and second LightGBM models (default: `0.02`). |

### Acquisition

| Parameter | Type | Description |
|-----------|------|-------------|
| `acquisition.nce` | Int | Initial NCE guess used during the pre-search before NCE tuning (default: `26`). |

### Optimization (Machine Learning)

| Parameter | Type | Description |
|-----------|------|-------------|
| `optimization.machine_learning.max_psm_memory_mb` | Real | Memory budget (MB) for in-memory PSM scoring (default: `2000`). When PSM rows × bytes-per-row exceeds this, ScoringSearch switches to the out-of-memory path. |
| `optimization.machine_learning.pep_bin_size` | Int | Number of PSMs per bin in the empirical q-value/PEP histogram (default: `10`). Smaller = finer-grained but noisier per-bin FDR estimate; larger = smoother but coarser. |

### Protein Inference

| Parameter | Type | Description |
|-----------|------|-------------|
| `proteinInference.min_peptides` | Int | Minimum unique peptides required for a protein group to be reported (default: `1`). |

### MaxLFQ

| Parameter | Type | Description |
|-----------|------|-------------|
| `maxLFQ.run_to_run_normalization` | Bool | Apply between-run median-spline normalization to peak areas (default: `false`). |
| `maxLFQ.max_chunk_size_mb` | Int | Maximum chunk size (MB) for the chunked merge during MaxLFQ (default: `1024`). |

### Output

| Parameter | Type | Description |
|-----------|------|-------------|
| `output.write_csv` | Bool | Write CSV copies of the precursor/protein output tables alongside the Arrow files (default: `true`). |
| `output.write_decoys` | Bool | Include decoy precursors in the output tables (default: `false`). |
| `output.delete_temp` | Bool | Delete the `temp_data/` scratch directory at the end of the run (default: `true`). |

### Logging

| Parameter | Type | Description |
|-----------|------|-------------|
| `logging.debug_console_level` | Int | Console verbosity. `0` (default) = user-facing only; higher values progressively expose internal `@debug_lN` messages. |
| `logging.max_message_bytes` | Int | Maximum bytes per log line before truncation (default: `4096`). Truncation preserves valid UTF-8 and appends a `… [truncated N bytes]` suffix. Override at runtime with the `PIONEER_MAX_LOG_MSG_BYTES` env var (clamped to `[1024, 1048576]`). |

### Paths

| Parameter | Type | Description |
|-----------|------|-------------|
| `paths.library` | String | Path to the `.poin` spectral library directory (built by `BuildSpecLib`). |
| `paths.ms_data` | String | Directory of converted MS Arrow files. |
| `paths.results` | String | Output directory. |

### Constants (formerly user-tunable, now hardcoded)

A number of knobs were previously surfaced in the JSON schema but are now constants in source. They are listed here for reference; users do not need to touch them.

| Constant | Value | Where |
|----------|-------|-------|
| `min_fraction_transmitted` | `0.25` | `MainSearch/types.jl`, `IntegrateChromatogramsSearch.jl`, `ParameterTuningSearch/types.jl` |
| `prec_estimation` | `PartialPrecCapture()` | same |
| `isotope_trace_type` | `SeperateTraces()` | same |
| `isotope_err_bounds` | `(1, 0)` | same |
| `max_frag_rank` (search-time) | `255` (= `typemax(UInt8)`) | `MainSearch/types.jl`, `IntegrateChromatogramsSearch.jl` |
| `MAXLFQ_NORM_N_RT_BINS` | `100` | `MaxLFQSearch.jl` |
| `MAXLFQ_NORM_SPLINE_N_KNOTS` | `7` | `MaxLFQSearch.jl` |
| `n_files_per_plot` | `12` | `WriteOutputs/qcPlots.jl` |
| `WIDE_SCOUT_TOL_PPM` | `100.0` | `ParameterTuningSearch/constants.jl` |
| `BITVEC_MIN_EXCESS_RATE_FILTER` | `0.05` | `BitVecCalibration/BitVecCalibration.jl` |
| `DEFAULT_INDEX_SEARCH_MIN_SCORE` | `15` | `MainSearch/types.jl` (CountFilter fallback when BitVec hasn't trained a LUT) |
| `FORCE_OOM` | `false` | `ScoringSearch/ScoringSearch.jl` (developer toggle) |

To change any of these, edit the constant in source. The historical JSON paths (`global.isotope_settings.*`, `chromatogram.*`, `parameter_tuning.*`, `acquisition.quad_transmission.*`, `rt_alignment.*`, `output.plots_per_page`, `fragment_index_search.*`, etc.) are no longer accepted — old configs containing them still parse (extra keys are ignored), but the values are not consumed.

## BuildSpecLib Configuration

### FASTA Input and Regex Mapping

`BuildSpecLib` accepts FASTA inputs in three forms via `GetBuildLibParams`:

1. **Single directory** — scans for all `.fasta` and `.fasta.gz` files
2. **Single file** — uses the specified file directly
3. **Mixed array** — any combination of directories and files

The regex patterns for parsing FASTA headers can be configured in three ways:

1. **Default regex set for all files**:
   ```julia
   GetBuildLibParams(out_dir, lib_name, [dir1, dir2, file1])
   ```

2. **Custom single regex set** applied to every file:
   ```julia
   GetBuildLibParams(out_dir, lib_name, [dir1, file1],
       regex_codes = Dict(
           "accessions" => "^>(\\S+)",
           "genes"      => "GN=(\\S+)",
           "proteins"   => "\\s+(.+?)\\s+OS=",
           "organisms"  => "OS=(.+?)\\s+GN="
       ))
   ```

3. **Positional mapping** (one regex set per FASTA input):
   ```julia
   GetBuildLibParams(out_dir, lib_name, [uniprot_dir, custom_file],
       regex_codes = [
           Dict("accessions" => "^\\w+\\|(\\w+)\\|", ...),  # uniprot_dir
           Dict("accessions" => "^>(\\S+)",          ...),  # custom_file
       ])
   ```

### FASTA Digest

| Parameter | Type | Description |
|-----------|------|-------------|
| `fasta_digest_params.min_length` | Int | Minimum peptide length (default: `7`). |
| `fasta_digest_params.max_length` | Int | Maximum peptide length (default: `30`). |
| `fasta_digest_params.min_charge` | Int | Minimum charge state (default: `2`). |
| `fasta_digest_params.max_charge` | Int | Maximum charge state (default: `4`). |
| `fasta_digest_params.cleavage_regex` | String | Cleavage rule (default: `[KR][^_\|$]`). To exclude cleavage before proline: `[KR][^P\|$]`. |
| `fasta_digest_params.missed_cleavages` | Int | Maximum missed cleavages (default: `1`). |
| `fasta_digest_params.max_var_mods` | Int | Maximum variable modifications per peptide (default: `1`). |
| `fasta_digest_params.add_decoys` | Bool | Generate decoy sequences (default: `true`). |
| `fasta_digest_params.entrapment_r` | Float | Entrapment-sequence ratio (default: `0`). |
| `fasta_digest_params.decoy_method` | String | `"shuffle"`, `"reverse"`, or `"diann_mutation"` (default: `"shuffle"`). |
| `fasta_digest_params.entrapment_method` | String | `"shuffle"` or `"reverse"` (default: `"shuffle"`). |

### Modifications

| Parameter | Type | Description |
|-----------|------|-------------|
| `variable_mods.{pattern, mass, name}` | [String], [Float], [String] | Variable modifications. Default: oxidation on M (`Unimod:35`, +15.99491 Da). |
| `fixed_mods.{pattern, mass, name}` | [String], [Float], [String] | Fixed modifications. Default: carbamidomethyl on C (`Unimod:4`, +57.021464 Da). |
| `isotope_mod_groups` | [Object] | Multiplexed labelling channels (default: `[]`). |

### NCE

| Parameter | Type | Description |
|-----------|------|-------------|
| `nce_params.nce` | Float | Base NCE for fragment prediction (default: `26.0`). |
| `nce_params.default_charge` | Int | Charge state for NCE calculations (default: `2`). |
| `nce_params.dynamic_nce` | Bool | Charge-dependent NCE adjustment (default: `true`). |

### Library

| Parameter | Type | Description |
|-----------|------|-------------|
| `library_params.auto_detect_frag_bounds` | Bool | Auto-detect fragment m/z bounds from the calibration RAW (default: `true`). When `false`, manual `frag_mz_min/max` are used. |
| `library_params.frag_mz_min` / `frag_mz_max` | Float | Manual fragment m/z bounds (default: `150.0` / `2020.0`). |
| `library_params.prec_mz_min` / `prec_mz_max` | Float | Precursor m/z bounds (default: `390.0` / `1010.0`). |
| `library_params.instrument_type` | String | Instrument hint (default: `"NONE"`). |
| `library_params.prediction_model` | String | Fragment-prediction model (default: `"altimeter"`). |

### Top-level

| Parameter | Type | Description |
|-----------|------|-------------|
| `library_path` | String | Output directory for the `.poin` library. |
| `fasta_paths` | [String] | FASTA files or directories. |
| `fasta_names` | [String] | Per-FASTA proteome label (e.g. `"HUMAN"`). |
| `fasta_header_regex_accessions` | [String] | Per-FASTA accession regex (capture group). |
| `fasta_header_regex_genes` | [String] | Per-FASTA gene regex. |
| `fasta_header_regex_proteins` | [String] | Per-FASTA protein-name regex. |
| `fasta_header_regex_organisms` | [String] | Per-FASTA organism regex. |
| `calibration_raw_file` | String | Path to a representative Arrow MS file for `auto_detect_frag_bounds`. Optional. |
| `include_contaminants` | Bool | Append the bundled contaminants FASTA (default: `true`). |
| `predict_fragments` | Bool | Run Koina prediction (default: `true`). Set `false` to use library-fed intensities only. |
| `match_lib_build_batch` | Int | Batch size for Koina prediction calls (default: `100000`). |
| `rank_to_score` | [Int] | Intensity-rank weights for fragment-index scoring (default: `[8,4,4,2,2,1,1]`). |

!!! note "Koina API retry behavior"
    Koina retry warnings log at debug level 2 (set `logging.debug_console_level: 2` in the search config to see them). The build only fails if all retry attempts are exhausted.

### Output structure

A successful `BuildSpecLib` run produces a `.poin` directory containing:

| File | Purpose |
|------|---------|
| `precursors.arrow` | Precursor table (sequence, charge, m/z, iRT, decoy flag, …) |
| `proteins_table.arrow` | Protein metadata |
| `detailed_fragments.jls` | Per-precursor fragment ions, m/z-sorted within each precursor |
| `precursor_to_fragment_indices.jls` | Per-precursor fragment range pointers |
| `partitioned_fragment_index.jls` | MainSearch partitioned fragment index (bitmask-encoded) |
| `presearch_partitioned_fragment_index.jls` | Pre-search partitioned fragment index |
| `spline_knots.jls` | Spline knots for `SplineCompactFrag` libraries (Altimeter) |
| `config.json` | Snapshot of the validated build parameters |
| `build_log.txt` | Build log |
