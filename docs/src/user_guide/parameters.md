# Parameter Configuration

Pioneer reads JSON configuration files for both `SearchDIA` and `BuildSpecLib`. The schemas below mirror [`defaultSearchParams.json`](https://github.com/nwamsley1/Pioneer.jl/blob/main/assets/example_config/defaultSearchParams.json) and [`defaultBuildLibParams.json`](https://github.com/nwamsley1/Pioneer.jl/blob/main/assets/example_config/defaultBuildLibParams.json). Only fields documented here are read; any other key is ignored.

A simplified search config (`defaultSearchParamsSimplified.json`) covers paths plus a handful of high-level toggles. Anything not specified inherits the defaults below.

## SearchDIA Configuration

Most parameters work at their defaults. The few worth tuning per experiment:

* **`global.q_value_threshold`** — final FDR cutoff for output (default `0.01`). Loosen to `0.05` for exploratory work.
* **`search.n_isotopes`** — number of fragment isotopes used in matching (default `2`, M and M+1). Set `1` for non-Altimeter libraries that do not model M+1 intensities (Prosit, UniSpec).
* **`acquisition.nce`** — initial NCE guess for the pre-search before NCE tuning (default `26`, suitable for Thermo Orbitrap/Astral). If the auto-fitted NCE in the QC plot is far from this value, re-run with a closer guess.
* **`optimization.machine_learning.max_psm_memory_mb`** — table-size budget for run-level protein training (default `2000` MB). Larger tables use a representative training pool. This is not a total search RAM limit; precursor training uses its own fixed pool cap.
* **`maxLFQ.run_to_run_normalization`** — apply between-run median-spline normalization to peak areas (default `true`). Turn off when between-run intensity differences are biological rather than systematic.

### Global

| Parameter | Type | Default | Description |
|---|---|---|---|
| `global.q_value_threshold` | Float | `0.01` | Final FDR threshold applied to ScoringSearch output. |
| `global.im_refinement` | String | `"auto"` | Ion-mobility prediction residual correction: `"none"`, `"charge"`, `"composition"`, `"charge_composition"`, or `"auto"`. Requires library mobility predictions and observed mobility data. |

Mobility refinement fits charge-specific offset/slopes and optionally residue,
modification and terminal-residue counts. `composition` shares composition
coefficients between charges; `charge_composition` fits them separately.
`auto` selects among these and no correction on an inner peptide-sequence holdout,
requiring a reduction in median-plus-tail error without materially worsening the
95th-percentile error. All corrections are applied on five outer sequence folds,
keeping charge states and modifications of a peptide together. Training uses
distinct target precursors with first-pass probability above 0.9 and q-value at
most 0.01; charges need at least 100 training anchors.

Refinement updates `im_error` before experiment-wide precursor and MBR scoring.
It preserves observed mobility (`im_obs`), scan calibration and extraction windows.
Output also records `im_pred`, `im_pred_refined` and `im_error_uncorrected` when
refinement runs. Residue tokens currently distinguish carbamidomethyl C and
oxidized M; other modifications do not receive separate coefficients.

For a reproducible prediction benchmark, run a baseline with
`global.im_refinement = "none"`, then use
`scripts/evaluate_im_refinement.jl anchors library.poin output_dir [tdfs_dir]`.
Prefer the baseline's retained `temp_data/main_search_psms` directory as `anchors`
(`output.delete_temp = false`), selecting confident targets before IM-dependent
experiment-wide scoring. A final `precursors_long.arrow` table is also supported.
Run the benchmark after the baseline search finishes; its fold tables are merged
and replaced during scoring.
The optional `.tdfs` directory uses the instrument's recorded calibration as the
observation, avoiding any peptide-fitted alignment in the evaluation targets.
The benchmark excludes MBR recoveries and measures accuracy conditional on
identification; search identification counts require a separate full search comparison.

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
| `optimization.machine_learning.max_psm_memory_mb` | Real | `2000` | Table-size budget (MB) for run-level protein training; larger tables use a representative training pool. Does not limit total search RAM. |

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

### Protein Quantification

| Parameter | Type | Default | Description |
|---|---|---|---|
| `maxLFQ.quantification_method` | String | `"sparsemaxlfq"` | Protein quantification method: `"sparsemaxlfq"` or `"maxlfq"`. |
| `maxLFQ.run_to_run_normalization` | Bool | `true` | Apply between-run median-spline normalization to peak areas. |
| `maxLFQ.max_chunk_size_mb` | Int | `1024` | Maximum chunk size (MB) for the chunked merge during protein quantification. |

Protein quantification defaults to sparse MaxLFQ with 16 partner proposals per
run and a fixed seed, plus connections that preserve the full overlap graph’s
connected components. It uses shared-precursor ratios on those selected run pairs;
results can differ from full MaxLFQ. For comparisons, select `"maxlfq"`.
The selected method and its settings are saved in `protein_quantification.json`.

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
| `fasta_digest_params.specificity` | String | `"full"` | Digestion specificity: `"full"`, `"semi"` (either terminus), `"semi-n"` (C terminus required), or `"semi-c"` (N terminus required). Protein termini count as enzymatic. |
| `fasta_digest_params.nterm_met_excision` | Bool | `true` | N-terminal Met excision. Each protein N-terminal peptide is emitted both with and without its initiator Met (`MPEPTIDEK` and `PEPTIDEK`); the excised form is length-filtered on its own and costs no missed cleavage. |
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

Modification names must be UNIMOD accessions (`Unimod:<id>`); they are sent to
Koina verbatim. Before any prediction is requested, every fixed and variable
modification (accession and residue) is checked against what the selected
`prediction_model` **and** `rt_model` were trained on, and a build whose
modifications either model cannot predict is refused — the error names the
offending modifications and the models that would accept the whole selection.
Leaving cysteine without a fixed modification means *unmodified cysteine*, which
only `prosit_2025_40ptm` (fragments) and both retention-time models can
predict; every other fragment model assumes carbamidomethyl-C.

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
| `library_params.im_model` | String | `""` | Koina ion-mobility model for timsTOF libraries (`"alphapept_ccs"` or `"im2deep"`); adds `ccs` / `inv_ion_mobility` precursor columns. Empty skips it. |
| `library_params.frag_index_local_id_type` | String | `"auto"` | Width of the fragment index's partition-local precursor IDs: `"auto"`, `"UInt16"` or `"UInt32"`. UInt16 partitions hold at most 65,535 precursors and a denser partition is split, so on large libraries the effective partition width drops below its nominal width (about 2.5 Da at 5 Da for a 10 M-precursor library). `"auto"` picks UInt32 only in that case. The choice is logged and recorded in the library's `config.json`. |

### Prediction Models

| Parameter | Type | Default | Description |
|---|---|---|---|
| `library_params.prediction_model` | String | `"altimeter"` | Koina fragment-intensity model: `altimeter`, `prosit_2020_hcd`, `prosit_2024_ptm`, or `prosit_2025_40ptm`. |
| `library_params.rt_model` | String | `"chronologer"` | Koina retention-time model: `chronologer` (hydrophobic index, %ACN) or `prosit_2024_irt_ptm` (Prosit iRT, the sibling of the Prosit PTM fragment models). Either scale works for the search, which calibrates RT↔iRT per file. The choice is recorded in the library's `config.json`. |

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
| `partitioned_fragment_index.jls` | MainSearch partitioned fragment index (5 Da precursor partitions). |
| `presearch_partitioned_fragment_index.jls` | Pre-search partitioned fragment index (5 Da). |
| `partitioned_fragment_index_w10.jls`, `presearch_partitioned_fragment_index_w10.jls` | The same two indexes with 10 Da partitions. SearchDIA loads one pair per search, chosen from the data's MS2 isolation windows (10 Da for windows of 7.5 m/z or wider, such as timsTOF diaPASEF). |
| `fragment_indices.json` | Lists the fragment indexes and their partition widths. Libraries built before it existed have only the 5 Da pair. |
| `spline_knots.jls` | Spline knots for `SplineCompactFrag` libraries (Altimeter). |
| `config.json` | Snapshot of the validated build parameters. |
| `build_log.txt` | Build log. |
