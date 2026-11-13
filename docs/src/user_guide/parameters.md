# Parameter Configuration

Pioneer reads JSON configuration files for both `SearchDIA` and `BuildSpecLib`. The schemas below mirror [`defaultSearchParams.json`](https://github.com/nwamsley1/Pioneer.jl/blob/main/assets/example_config/defaultSearchParams.json) and [`defaultBuildLibParams.json`](https://github.com/nwamsley1/Pioneer.jl/blob/main/assets/example_config/defaultBuildLibParams.json). Parameters left out of a file keep their defaults.

A simplified search config (`defaultSearchParamsSimplified.json`) covers paths plus a handful of high-level toggles. Anything not specified inherits the defaults below.

## Frequently Modified Parameters

Most parameters should remain at their defaults for most experiments, but a few are more frequently modified.

| Parameter | Default | Description |
|---|---|---|
| `global.q_value_threshold` | `0.01` | FDR threshold for the reported precursors and protein groups. |
| `global.match_between_runs` | `true` | Transfer identifications to runs where a precursor was not identified on its own. |
| `maxLFQ.run_to_run_normalization` | `true` | Normalize peak areas between runs. Turn off when between-run intensity differences are biological rather than systematic. |
| `output.write_csv` | `true` | Write `.tsv` copies of the result tables. Turning it off is recommended for large searches. |
| `proteinScoring.write_qc_plots` | `false` | Write additional protein-scoring QC plots. |

## SearchDIA Configuration

### Global

| Parameter | Type | Default | Description |
|---|---|---|---|
| `global.q_value_threshold` | Float | `0.01` | FDR threshold for the precursors and protein groups reported in the output. |
| `global.match_between_runs` | Bool | `true` | Match-between-runs: transfer identifications to runs where a precursor was not identified on its own. Not in the default file; omitted means `true`. |

### Search

| Parameter | Type | Default | Description |
|---|---|---|---|
| `search.n_isotopes` | Int | `2` | Fragment isotopes considered in matching. Use `1` for non-Altimeter libraries. |

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
run and a fixed seed, plus connections that preserve the full overlap graph's
connected components. It uses shared-precursor ratios on those selected run pairs;
results can differ from full MaxLFQ. For comparisons, select `"maxlfq"`.

### Output

| Parameter | Type | Default | Description |
|---|---|---|---|
| `output.write_csv` | Bool | `true` | Write tab-separated (`.tsv`) copies of the precursor and protein group tables alongside the Arrow files. |
| `output.write_decoys` | Bool | `false` | Include decoy precursors in the output tables. |
| `output.delete_temp` | Bool | `true` | Delete the `temp_data/` scratch directory at the end of the run. |

### Logging

| Parameter | Type | Default | Description |
|---|---|---|---|
| `logging.debug_console_level` | Int | `0` | Console verbosity. `0` shows user-facing messages only; higher values progressively expose internal `@debug_lN` messages. Level 1 includes the protein-scoring model's feature gains. |
| `logging.debug_file_level` | Int | `1` | Verbosity of `pioneer_search_debug.log`, on the same scale as `debug_console_level`. |
| `logging.max_message_bytes` | Int | `4096` | Maximum bytes per log line before truncation. Truncation preserves valid UTF-8 and appends `… [truncated N bytes]`. The `PIONEER_MAX_LOG_MSG_BYTES` environment variable overrides it at runtime. Either value is clamped to `[1024, 1048576]`. |

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

Header-parsing regex patterns can be configured three ways. The accession pattern is matched against the FASTA identifier, the first word of the header without the leading `>`; a pattern that does not match falls back to the whole identifier.

1. Default regex set, applied to all files:
   ```julia
   GetBuildLibParams(out_dir, lib_name, [dir1, dir2, file1])
   ```
2. Custom single regex set, applied to every file:
   ```julia
   GetBuildLibParams(out_dir, lib_name, [dir1, file1],
       regex_codes = Dict(
           "accessions" => "^(\\S+)",
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
           Dict("accessions" => "^(\\S+)",           ...),  # custom_file
       ])
   ```

### FASTA Digest

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fasta_digest_params.min_length` | Int | `7` | Minimum peptide length. |
| `fasta_digest_params.max_length` | Int | `30` | Maximum peptide length (`40` in the simplified template). |
| `fasta_digest_params.min_charge` | Int | `2` | Minimum charge state. |
| `fasta_digest_params.max_charge` | Int | `4` | Maximum charge state (`3` in the simplified template). |
| `fasta_digest_params.cleavage_regex` | String | `[KR][^_\|$]` | Cleavage rule. To exclude cleavage before proline use `[KR][^P\|$]`. |
| `fasta_digest_params.missed_cleavages` | Int | `1` | Maximum missed cleavages. |
| `fasta_digest_params.specificity` | String | `"full"` | Digestion specificity: `"full"`, `"semi"` (either terminus), `"semi-n"` (C terminus required), or `"semi-c"` (N terminus required). Protein termini count as enzymatic. |
| `fasta_digest_params.nterm_met_excision` | Bool | `true` | N-terminal Met excision. Each protein N-terminal peptide is emitted both with and without its initiator Met (`MPEPTIDEK` and `PEPTIDEK`); the excised form is length-filtered on its own and costs no missed cleavage. |
| `fasta_digest_params.max_var_mods` | Int | `1` | Maximum variable modifications per peptide. |
| `fasta_digest_params.add_decoys` | Bool | `true` | Generate decoy sequences. |
| `fasta_digest_params.entrapment_r` | Int | `0` | Number of entrapment sequences generated per target sequence; `0` disables entrapment. |
| `fasta_digest_params.decoy_method` | String | `"shuffle"` | One of `"shuffle"`, `"reverse"`, `"diann_mutation"`. |
| `fasta_digest_params.entrapment_method` | String | `"shuffle"` | `"shuffle"` or `"reverse"`. `"reverse"` cannot be combined with `decoy_method = "reverse"`. |

### Modifications

| Parameter | Type | Default | Description |
|---|---|---|---|
| `variable_mods.{pattern, mass, name}` | [String], [Float], [String] | Met oxidation (`Unimod:35`, +15.99491 Da) | Variable modifications. |
| `fixed_mods.{pattern, mass, name}` | [String], [Float], [String] | Cys carbamidomethyl (`Unimod:4`, +57.021464 Da) | Fixed modifications. |
| `isotope_mod_groups` | [Object] | `[]` | Multiplexed labeling channels. Each object requires `pattern`, `mass`, `mod_names` (arrays) and `label_name` (string). |

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
| `nce_params.nce` | Float | `26.0` | Collision energy sent to Koina for the Prosit fragment models. Altimeter, the default, does not use it. |

### Library m/z Bounds

| Parameter | Type | Default | Description |
|---|---|---|---|
| `library_params.auto_detect_frag_bounds` | Bool | `true` | Detect the fragment and precursor m/z bounds from `calibration_raw_file` (precursor bounds padded by 1 m/z). When `false`, or without a calibration file, the manual bounds below are used. |
| `library_params.frag_mz_min` | Float | `150.0` | Manual lower fragment m/z bound. |
| `library_params.frag_mz_max` | Float | `2020.0` | Manual upper fragment m/z bound. |
| `library_params.prec_mz_min` | Float | `390.0` | Manual lower precursor m/z bound. |
| `library_params.prec_mz_max` | Float | `1010.0` | Manual upper precursor m/z bound. |
| `library_params.im_model` | String | `"alphapept_ccs"` | Koina ion-mobility model for timsTOF libraries (`"alphapept_ccs"` or `"im2deep"`); adds `ccs` / `inv_ion_mobility` precursor columns. Empty skips it. |
| `library_params.frag_index_local_id_type` | String | `"auto"` | Width of the fragment index's partition-local precursor IDs: `"auto"`, `"UInt16"` or `"UInt32"`. UInt16 partitions hold at most 65,535 precursors and a denser partition is split, so on large libraries the effective partition width drops below its nominal width (about 2.5 Da at 5 Da for a 10-million-precursor library). `"auto"` picks UInt32 only in that case. The choice is logged and recorded in the library's `config.json`. |

### Prediction Models

| Parameter | Type | Default | Description |
|---|---|---|---|
| `library_params.prediction_model` | String | `"altimeter"` | Koina fragment-intensity model: `altimeter`, `prosit_2020_hcd`, `prosit_2024_ptm`, or `prosit_2025_40ptm`. Not in the default file; omitted means `altimeter`. |
| `library_params.rt_model` | String | `"chronologer"` | Koina retention-time model: `chronologer` (hydrophobic index, %ACN) or `prosit_2024_irt_ptm` (Prosit iRT, the sibling of the Prosit PTM fragment models). Either scale works for the search, which calibrates RT↔iRT per file. The choice is recorded in the library's `config.json`. Not in the default file; omitted means `chronologer`. |

### Top-level

| Parameter | Type | Default | Description |
|---|---|---|---|
| `library_path` | String | — | Output directory for the `.poin` library. |
| `fasta_paths` | [String] | — | FASTA files or directories. |
| `fasta_names` | [String] | — | Per-FASTA proteome label (e.g. `"HUMAN"`), reported in the `species` column of the search results. |
| `fasta_header_regex_accessions` | [String] | UniProt | Per-FASTA accession capture regex. |
| `fasta_header_regex_genes` | [String] | UniProt | Per-FASTA gene regex. |
| `fasta_header_regex_proteins` | [String] | UniProt | Per-FASTA protein-name regex. |
| `fasta_header_regex_organisms` | [String] | UniProt | Per-FASTA organism regex. |
| `calibration_raw_file` | String | — | Optional. A representative converted run (`.arrow`, timsTOF `.tdfs`, or SCIEX `.scxs`) used by `auto_detect_frag_bounds`. |
| `include_contaminants` | Bool | `true` | Append the bundled contaminants FASTA. |
| `predict_fragments` | Bool | `true` | Run Koina fragment-intensity prediction. `false` skips prediction and resumes from the prediction files already in the library folder; the build stops if they are missing. |
| `seed` | Int | `1844` | Random seed for decoy and entrapment shuffling. |

!!! note "Koina API retry behavior"
    Koina retry warnings log at debug level 2. Set `logging.debug_console_level` to `2` in the BuildSpecLib parameter file to see them. The build only fails if all retry attempts are exhausted.

### Output Structure

A successful `BuildSpecLib` run writes a `.poin` directory containing:

| File | Purpose |
|---|---|
| `precursors_table.arrow` | Precursor table (sequence, charge, m/z, iRT, decoy flag, …). |
| `proteins_table.arrow` | Protein metadata. |
| `detailed_fragments.jls` | Per-precursor fragment ions, m/z-sorted within each precursor. |
| `precursor_to_fragment_indices.jls` | Per-precursor fragment range pointers. |
| `partitioned_fragment_index.jls` | MainSearch partitioned fragment index (5 Da precursor partitions). |
| `presearch_partitioned_fragment_index.jls` | Pre-search partitioned fragment index (5 Da). |
| `partitioned_fragment_index_w10.jls`, `presearch_partitioned_fragment_index_w10.jls` | The same two indexes with 10 Da partitions. SearchDIA loads one pair per search, chosen from the data's MS2 isolation windows (10 Da for windows of 7.5 m/z or wider, such as timsTOF diaPASEF). |
| `fragment_indices.json` | Lists the fragment indexes and their partition widths. Libraries built before it existed have only the 5 Da pair. |
| `spline_knots.jls` | Spline knots for `SplineCompactFrag` libraries (Altimeter). |
| `frag_name_to_idx.jls`, `ion_annotations.jls` | Fragment ion annotation tables. |
| `config.json` | Snapshot of the validated build parameters, plus `pioneer_version`, `build_date` and `fasta_provenance` (name, size, SHA-256 and bundled path of each input FASTA). |
| `fasta/` | Copies of every input FASTA, including contaminants. A `<fasta>.provenance.json` placed next to an input FASTA (for example with its UniProt release, source URL and download time) is copied too, and its fields are merged into that FASTA's `fasta_provenance` record. |
| `pioneer_search_report.txt`, `pioneer_search_log.log`, `pioneer_search_debug.log`, `pioneer_warnings.log` | Build summary and logs. |

Downloaded prebuilt libraries also include `SHA256SUMS`, which `DownloadSpecLib` uses to verify the download.
