# Search Results

SearchDIA writes its results to the folder chosen under **Results folder** (`paths.results` in the parameter file).

## Output Files

| File | Contents |
|---|---|
| `precursors_long.tsv` / `.arrow` | Identified precursors, one row per precursor per run |
| `precursors_wide.tsv` / `.arrow` | Identified precursors, one row per precursor with a quantity column for each run |
| `protein_groups_long.tsv` / `.arrow` | Protein groups, one row per protein group per run |
| `protein_groups_wide.tsv` / `.arrow` | Protein groups, one row per protein group with an abundance column for each run |
| `run_summary.tsv` | One row of quality-control metrics per run |
| `qc_plots/` | Retention time alignment, mass error, collision energy, and quadrupole transmission plots |
| `pioneer_search_report.txt` | Summary of the search |
| `pioneer_search_log.log`, `pioneer_search_debug.log`, `pioneer_warnings.log` | Full log, debug log, and warnings |
| `config.json` | All parameters used, including defaults |

Only target identifications are written unless `output.write_decoys` is on. Every table row passes the q-value threshold (`global.q_value_threshold`, 1% FDR by default).

## Per-Run Summary (`run_summary.tsv`)

Every search writes `run_summary.tsv` next to the result tables: one row per
raw file with QC metrics in the spirit of DIA-NN's `report.stats.tsv`.
*Identified* counts target IDs regardless of abundance, including missing values.
*Quantified* counts the subset with positive abundance. Missed cleavages are
reported as 100 × their mean over identified target precursors.

| Column | Description |
|---|---|
| `file_name` | Run name (raw file name without extension) |
| `precursors_identified` / `precursors_quantified` | Target precursors passing FDR / with a positive `peak_area` |
| `precursors_mbr` | Identified precursors recovered by match-between-runs |
| `peptides_identified` | Distinct stripped sequences among identified precursors |
| `protein_groups_identified` / `protein_groups_quantified` | Target protein groups for the run / with a positive abundance |
| `total_peak_area`, `median_peak_area` | Sum and median of quantified precursor areas |
| `median_normalization_factor` | Median `peak_area_normalized / peak_area` (empty when normalization is off) |
| `median_irt_error` | Median absolute iRT prediction error |
| `median_rt_fwhm` | Median RT span (min) of the scans at or above half the apex intensity; 0 when only the apex scan qualifies (fast gradients / few points per peak) |
| `median_points_integrated` | Median number of scans integrated per peak |
| `median_peptide_length`, `median_charge` | Peptide-property medians |
| `missed_cleavage_percentage` | 100 × mean missed cleavages among identified target precursors; can exceed 100 |
| `ms2_mass_tol_low/high/unit`, `ms1_mass_tol_low/high/unit` | Fragment / precursor mass tolerance the search settled on (`ppm` or `Da` depending on the fitted model; MS1 empty when no MS1 model was fit) |
| `gradient_length_min`, `n_ms1_scans`, `n_ms2_scans` | Run length (last retention time, minutes) and scan counts |
| `<stage>_calibration_qc` | Calibration check for each stage: `rt`, `ms2_mass`, `ms1_mass`, `quadrupole`, `nce`, and `ion_mobility` for ion mobility data. One of `normal`, `suspicious`, `warning`, `failed`, or `not_assessed` (stage skipped or not reached) |
| `<stage>_calibration_qc_reason` | Why the check was not `normal`, for example too few identifications to fit the model or a fallback calibration; empty when `normal` |

## Long and Wide Formats

The precursor and protein group results are each written in two layouts.

* **Long:** one row per identification per run. A precursor identified in three runs has three rows, distinguished by `file_name`. Long tables carry every per-run value: scores, q-values, peak areas, and retention times.
* **Wide:** one row per precursor or protein group across the whole experiment, followed by one column per run, named after the run. Each run column holds the quantity in that run and is empty where the precursor or protein group was not identified or not quantified. Wide tables carry only experiment-level descriptors (sequence, protein, global score and q-value), since per-run scores do not fit in a single row.

For example, a precursor identified in runs `A` and `B`:

| Long: `file_name` | `sequence` | `peak_area` |
|---|---|---|
| A | PEPTIDEK | 1.2e7 |
| B | PEPTIDEK | 9.8e6 |

| Wide: `sequence` | `A` | `B` | `C` |
|---|---|---|---|
| PEPTIDEK | 1.2e7 | 9.8e6 | |

With run-to-run normalization on, the wide run columns hold the normalized areas (`peak_area_normalized`).

## TSV and Arrow Files

* **`.tsv`** files are tab-separated text that opens in Excel, R, Python, or any text editor. They are written when `output.write_csv` is on (the default); turning it off is recommended for large searches.
* **`.arrow`** files are binary [Apache Arrow](https://arrow.apache.org/) tables, always written. They load much faster than text for large experiments (`Arrow.Table` in Julia, `pyarrow.feather.read_table` or `polars.read_ipc` in Python, `arrow::read_feather` in R).

* `precursors_long.arrow` keeps internal column names (`prec_prob` for `score`, `global_prob` for `global_score`, `new_best_scan` for `apex_scan`), omits `gene_names`, and includes additional internal scoring columns.
* `protein_groups_long.arrow` omits `gene_names` and `protein_names`, and stores `peptides` as library precursor indices instead of sequences.

## Precursor Tables

A precursor is a peptide sequence with a specific set of modifications and a specific charge state.

### `precursors_long`

| Column | Description |
|---|---|
| `file_name` | Run name (raw file name without extension) |
| `species` | Species of the matching proteins, from the FASTA files used to build the library; `;`-separated when several |
| `gene_names` | Gene names of the matching proteins; `;`-separated |
| `inferred_protein_group` | Protein group assigned by protein inference: `;`-separated accessions |
| `accession_numbers` | All protein accessions containing the peptide; `;`-separated |
| `sequence` | Unmodified peptide sequence |
| `charge` | Precursor charge |
| `structural_mods` | Modifications as `(position,residue,Unimod:id)` entries, for example `(4,C,Unimod:4)` for carbamidomethyl cysteine at position 4; empty when unmodified |
| `isotopic_mods` | Isotope label modifications; present only for libraries built with isotope labels |
| `prec_mz` | Precursor m/z |
| `peptide_start_positions` | Start position(s) of the peptide within the matching protein sequence(s) |
| `missed_cleavage` | Number of missed cleavages |
| `num_enzymatic_termini` | Number of termini (0–2) consistent with the digestion enzyme |
| `global_score` | Experiment-wide precursor score; the same in every run |
| `score` | Precursor score in this run |
| `global_qval` | Experiment-wide precursor q-value |
| `qval` | Precursor q-value in this run |
| `pep` | Posterior error probability in this run |
| `peak_area` | Integrated chromatographic peak area; empty when the precursor was identified but not quantified |
| `peak_area_normalized` | Peak area after run-to-run normalization; present only when normalization is on (`maxLFQ.run_to_run_normalization`) |
| `points_integrated` | Number of scans integrated across the peak |
| `quant_withheld` | `true` when the peak area was withheld because less than a fifth of the signal survived baseline subtraction; explains an empty `peak_area` |
| `precursor_fraction_transmitted` | Estimated fraction of the precursor isotope envelope transmitted by the quadrupole isolation window |
| `rt` | Retention time (minutes) of the best-scoring scan |
| `apex_scan` | Scan index of the chromatographic peak apex |
| `global_pg_score`, `pg_score` | Score of the inferred protein group: experiment-wide / in this run |
| `global_pg_qval`, `pg_qval` | q-value of the inferred protein group: experiment-wide / in this run |
| `use_for_protein_quant` | `true` when the peptide is unique to its protein group and contributes to the group's abundance; shared peptides count toward protein scoring only |
| `precursor_idx` | Precursor index in the spectral library |
| `target` | `true` for targets, `false` for decoys |
| `entrapment_group_id` | Entrapment group (0 = original sequences); present only for libraries built with entrapment sequences |

### `precursors_wide`

The precursor descriptors from `species` through `global_qval` (above, without the per-run columns), then `use_for_protein_quant`, `precursor_idx`, `target`, and `entrapment_group_id` when present, followed by one column per run. Each run column holds `peak_area_normalized`, or `peak_area` when normalization is off.

## Protein Group Tables

A protein group is the set of proteins that the identified peptides cannot distinguish; it is reported as `;`-separated accessions.

### `protein_groups_long`

| Column | Description |
|---|---|
| `file_name` | Run name |
| `species` | Species of the proteins in the group; `;`-separated when several |
| `gene_names` | Gene names of the proteins in the group; `;`-separated |
| `protein_names` | Protein names of the proteins in the group; `;`-separated |
| `protein` | Protein group accessions; `;`-separated |
| `target` | `true` for targets, `false` for decoys |
| `entrap_id` | Entrapment group (0 = original sequences); present only for libraries built with entrapment sequences |
| `peptides` | Precursors identified for the group in this run, as `_SEQUENCE_.charge` with modifications; `;`-separated |
| `n_precursors` | Number of precursors identified for the group in this run |
| `n_precursors_quantified` | Number of those precursors with a peak area used for quantification |
| `n_modified_peptides` | Number of distinct modified peptide sequences |
| `n_peptides` | Number of distinct unmodified peptide sequences |
| `global_pg_score`, `pg_score` | Protein group score: experiment-wide / in this run |
| `global_qval`, `qval` | Protein group q-value: experiment-wide / in this run |
| `pg_pep` | Protein group posterior error probability in this run |
| `total_peak_area` | Sum of the quantified precursor peak areas (normalized when normalization is on) in this run |
| `abundance` | Protein group abundance in this run, estimated by MaxLFQ across runs (`maxLFQ.quantification_method`). Equals `total_peak_area` when the group is quantified in only one run |

### `protein_groups_wide`

`species`, `gene_names`, `protein_names`, `protein`, `target`, `entrap_id` when present, `global_pg_score`, and `global_qval`, followed by one column per run holding `abundance`.

## Quality-Control Plots (`qc_plots/`)

| Folder | Contents |
|---|---|
| `rt_alignment_plots/` | Library retention time versus observed retention time, per run |
| `mass_error_plots/` | Fragment mass errors, per run |
| `ms1_mass_error_plots/` | Precursor mass errors, per run |
| `collision_energy_alignment/` | Fitted collision energy, per run |
| `quad_transmission_model/` | Fitted quadrupole transmission, per run |
| `ion_mobility_model/` | Ion mobility calibration (timsTOF data only) |
