# Parameter Configuration

Pioneer reads JSON configuration files for both `SearchDIA` and `BuildSpecLib`. The schemas below mirror [`defaultSearchParams.json`](https://github.com/nwamsley1/Pioneer.jl/blob/main/assets/example_config/defaultSearchParams.json) and [`defaultBuildLibParams.json`](https://github.com/nwamsley1/Pioneer.jl/blob/main/assets/example_config/defaultBuildLibParams.json). Only fields documented here are read; any other key is ignored.

A simplified search config (`defaultSearchParamsSimplified.json`) covers paths plus a handful of high-level toggles. Anything not specified inherits the defaults below.

## SearchDIA Configuration

Most parameters work at their defaults. The few worth tuning per experiment:

* **`global.q_value_threshold`** — final FDR cutoff for output (default `0.01`). Loosen to `0.05` for exploratory work.
* **`search.n_isotopes`** — number of fragment isotopes used in matching (default `2`, M and M+1). Set `1` for non-Altimeter libraries that do not model M+1 intensities (Prosit, UniSpec).
* **`acquisition.nce`** — initial NCE guess for the pre-search before NCE tuning (default `26`, suitable for Thermo Orbitrap/Astral). If the auto-fitted NCE in the QC plot is far from this value, re-run with a closer guess.
* **`optimization.machine_learning.max_psm_memory_mb`** — memory budget for in-memory LightGBM training (default `2000` MB). Raise on workstations with more RAM; lower to force the out-of-memory path earlier.
* **`maxLFQ.run_to_run_normalization`** — apply between-run median-spline normalization to peak areas (default `false`). Turn on when systemic between-run intensity bias is expected.

### Global

| Parameter | Type | Default | Description |
|---|---|---|---|
| `global.q_value_threshold` | Float | `0.01` | Final FDR threshold applied to ScoringSearch output. |

### Search

| Parameter | Type | Default | Description |
|---|---|---|---|
| `search.n_isotopes` | Int | `2` | Fragment isotopes considered in matching. Use `1` for non-Altimeter libraries. |

### Acquisition

| Parameter | Type | Default | Description |
|---|---|---|---|
| `acquisition.nce` | Int | `26` | Initial NCE guess used during the pre-search before NCE tuning. |

### Optimization (Machine Learning)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `optimization.machine_learning.max_psm_memory_mb` | Real | `2000` | Memory budget (MB) for in-memory PSM scoring. Above this, ScoringSearch switches to the out-of-memory path. |
| `optimization.machine_learning.pep_bin_size` | Int | `10` | PSMs per bin in the empirical q-value/PEP histogram. Smaller is finer-grained but noisier; larger is smoother but coarser. |

### Optimization (Chromatogram Integration)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `optimization.chromatogram_integration.trace_mode` | String | `"combined"` | `"combined"` integrates all isotope traces of a precursor as a single chromatogram; `"separate"` integrates each trace independently. |
| `optimization.chromatogram_integration.deconvolution_solver` | String | `"huber"` | Solver used for chromatogram deconvolution. `"huber"` is the robust default; `"pmm"` selects the Poisson MM solver. |

### Protein Scoring

| Parameter | Type | Default | Description |
|---|---|---|---|
| `proteinScoring.min_peptides` | Int | `1` | Minimum unique peptides required for a protein group to be reported. |
| `proteinScoring.global_protein_inference` | Bool | `true` | Run protein inference once across the union of passing PSMs from every file. Set `false` for the legacy per-file path. |
| `proteinScoring.write_qc_plots` | Bool | `false` | Emit protein-scoring QC plots. |

### MaxLFQ

| Parameter | Type | Default | Description |
|---|---|---|---|
| `maxLFQ.run_to_run_normalization` | Bool | `false` | Apply between-run median-spline normalization to peak areas. |
| `maxLFQ.max_chunk_size_mb` | Int | `1024` | Maximum chunk size (MB) for the chunked merge during MaxLFQ. |

### Output

| Parameter | Type | Default | Description |
|---|---|---|---|
| `output.write_csv` | Bool | `true` | Write CSV copies of the precursor and protein output tables alongside the Arrow files. |
| `output.write_decoys` | Bool | `false` | Include decoy precursors in the output tables. |
| `output.delete_temp` | Bool | `true` | Delete the `temp_data/` scratch directory at the end of the run. |

### Logging

| Parameter | Type | Default | Description |
|---|---|---|---|
| `logging.debug_console_level` | Int | `0` | Console verbosity. `0` shows user-facing messages only; higher values progressively expose internal `@debug_lN` messages. Level 1 includes protein probit coefficients and standardized feature importances. |
| `logging.max_message_bytes` | Int | `4096` | Maximum bytes per log line before truncation. Truncation preserves valid UTF-8 and appends `… [truncated N bytes]`. The `PIONEER_MAX_LOG_MSG_BYTES` env var overrides at runtime, clamped to `[1024, 1048576]`. |

### Paths

| Parameter | Type | Description |
|---|---|---|
| `paths.library` | String | Path to the `.poin` spectral library directory built by `BuildSpecLib`. |
| `paths.ms_data` | String | Directory of converted MS Arrow files. |
| `paths.results` | String | Output directory for results. |

## BuildSpecLib Configuration

### FASTA Inputs and Regex Mapping

`BuildSpecLib` accepts FASTA inputs in three forms via `GetBuildLibParams`:

1. **Single directory** — scans for all `.fasta` and `.fasta.gz` files.
2. **Single file** — uses the specified file directly.
3. **Mixed array** — any combination of directories and files.

Header-parsing regex patterns can be configured three ways:

1. Default regex set, applied to all files:
   ```julia
   GetBuildLibParams(out_dir, lib_name, [dir1, dir2, file1])
   ```
2. Custom single regex set, applied to every file:
   ```julia
   GetBuildLibParams(out_dir, lib_name, [dir1, file1],
       regex_codes = Dict(
           "accessions" => "^>(\\S+)",
           "genes"      => "GN=(\\S+)",
           "proteins"   => "\\s+(.+?)\\s+OS=",
           "organisms"  => "OS=(.+?)\\s+GN="
       ))
   ```
3. Positional mapping — one regex set per FASTA input:
   ```julia
   GetBuildLibParams(out_dir, lib_name, [uniprot_dir, custom_file],
       regex_codes = [
           Dict("accessions" => "^\\w+\\|(\\w+)\\|", ...),  # uniprot_dir
           Dict("accessions" => "^>(\\S+)",          ...),  # custom_file
       ])
   ```

### FASTA Digest

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fasta_digest_params.min_length` | Int | `7` | Minimum peptide length. |
| `fasta_digest_params.max_length` | Int | `30` | Maximum peptide length. |
| `fasta_digest_params.min_charge` | Int | `2` | Minimum charge state. |
| `fasta_digest_params.max_charge` | Int | `4` | Maximum charge state. |
| `fasta_digest_params.cleavage_regex` | String | `[KR][^_\|$]` | Cleavage rule. To exclude cleavage before proline use `[KR][^P\|$]`. |
| `fasta_digest_params.missed_cleavages` | Int | `1` | Maximum missed cleavages. |
| `fasta_digest_params.max_var_mods` | Int | `1` | Maximum variable modifications per peptide. |
| `fasta_digest_params.add_decoys` | Bool | `true` | Generate decoy sequences. |
| `fasta_digest_params.entrapment_r` | Float | `0` | Entrapment-sequence ratio. |
| `fasta_digest_params.decoy_method` | String | `"shuffle"` | One of `"shuffle"`, `"reverse"`, `"diann_mutation"`. |
| `fasta_digest_params.entrapment_method` | String | `"shuffle"` | One of `"shuffle"` or `"reverse"`. |

### Modifications

| Parameter | Type | Default | Description |
|---|---|---|---|
| `variable_mods.{pattern, mass, name}` | [String], [Float], [String] | Met oxidation (`Unimod:35`, +15.99491 Da) | Variable modifications. |
| `fixed_mods.{pattern, mass, name}` | [String], [Float], [String] | Cys carbamidomethyl (`Unimod:4`, +57.021464 Da) | Fixed modifications. |
| `isotope_mod_groups` | [Object] | `[]` | Multiplexed labelling channels. |

### Collision Energy

| Parameter | Type | Default | Description |
|---|---|---|---|
| `nce_params.nce` | Float | `26.0` | Base NCE used by the fragment-prediction model. |

### Library m/z Bounds

| Parameter | Type | Default | Description |
|---|---|---|---|
| `library_params.auto_detect_frag_bounds` | Bool | `true` | Detect fragment m/z bounds from the calibration RAW. When `false`, manual `frag_mz_min/max` are used. |
| `library_params.frag_mz_min` | Float | `150.0` | Manual lower fragment m/z bound. |
| `library_params.frag_mz_max` | Float | `2020.0` | Manual upper fragment m/z bound. |
| `library_params.prec_mz_min` | Float | `390.0` | Lower precursor m/z bound. |
| `library_params.prec_mz_max` | Float | `1010.0` | Upper precursor m/z bound. |

### Top-level

| Parameter | Type | Default | Description |
|---|---|---|---|
| `library_path` | String | — | Output directory for the `.poin` library. |
| `fasta_paths` | [String] | — | FASTA files or directories. |
| `fasta_names` | [String] | — | Per-FASTA proteome label (e.g. `"HUMAN"`). |
| `fasta_header_regex_accessions` | [String] | UniProt | Per-FASTA accession capture regex. |
| `fasta_header_regex_genes` | [String] | UniProt | Per-FASTA gene regex. |
| `fasta_header_regex_proteins` | [String] | UniProt | Per-FASTA protein-name regex. |
| `fasta_header_regex_organisms` | [String] | UniProt | Per-FASTA organism regex. |
| `calibration_raw_file` | String | — | Optional. Path to a representative MS Arrow file used by `auto_detect_frag_bounds`. |
| `include_contaminants` | Bool | `true` | Append the bundled contaminants FASTA. |
| `predict_fragments` | Bool | `true` | Run Koina fragment-intensity prediction. Set `false` to use library-fed intensities only. |
| `match_lib_build_batch` | Int | `100000` | Batch size for Koina prediction calls. |

!!! note "Koina API retry behavior"
    Koina retry warnings log at debug level 2. Set `logging.debug_console_level: 2` in the search config to see them. The build only fails if all retry attempts are exhausted.

### Output Structure

A successful `BuildSpecLib` run writes a `.poin` directory containing:

| File | Purpose |
|---|---|
| `precursors.arrow` | Precursor table (sequence, charge, m/z, iRT, decoy flag, …). |
| `proteins_table.arrow` | Protein metadata. |
| `detailed_fragments.jls` | Per-precursor fragment ions, m/z-sorted within each precursor. |
| `precursor_to_fragment_indices.jls` | Per-precursor fragment range pointers. |
| `partitioned_fragment_index.jls` | MainSearch partitioned fragment index. |
| `presearch_partitioned_fragment_index.jls` | Pre-search partitioned fragment index. |
| `spline_knots.jls` | Spline knots for `SplineCompactFrag` libraries (Altimeter). |
| `config.json` | Snapshot of the validated build parameters. |
| `build_log.txt` | Build log. |
