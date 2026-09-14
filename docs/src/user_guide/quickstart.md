# Quick Start Tutorial

## Basic Workflow
Pioneer performs three major steps:
1. Convert vendor MS files into the Arrow format using [PioneerConverter](https://github.com/nwamsley1/PioneerConverter).
2. Build in silico spectral libraries using FASTA files and the [Koina](https://koina.wilhelmlab.org/) server.
3. Search DIA experiments using a spectral library and the MS data files.

## Pioneer Converter
Pioneer operates on MS/MS data stored in the [Apache Arrow IPC format](https://arrow.apache.org/docs/python/ipc.html).
Use the bundled [PioneerConverter](https://github.com/nwamsley1/PioneerConverter) via the CLI to convert Thermo RAW files:

```bash
pioneer convert-raw /path/to/raw/or/folder --output-dir /path/to/arrow --skip-existing
```

This subcommand accepts either a single `.raw` file or a directory of files. Common options are `--output-dir`, `--skip-existing`, `--concurrent-files`, and `--threads-per-file`. For all options, run `pioneer convert-raw --help`.

## MzML to Arrow IPC (Sciex)
For mzML-formatted data, use:

```bash
pioneer convert-mzml /path/to/mzml/or/folder
```

## Starting Pioneer
After installation, Pioneer is accessed from the command line. Running `pioneer --help` displays available subcommands:

```bash
pioneer [options] <subcommand> [subcommand-args...]
```

Subcommands include `search`, `predict`, `params-search`, `params-predict`, `convert-raw`, and `convert-mzml`. On the first launch macOS performs a one-time Gatekeeper check.

A minimal end-to-end workflow is:

```bash
# Generate library build parameters (library name derived from output path)
pioneer params-predict lib_dir fasta_dir --params-path=predict_params.json

# Edit predict_params.json to customize:
# - Digestion parameters (missed cleavages, modifications, etc.)
# - For multiple FASTA sources, edit fasta_paths array to include directories and/or files
# - Set calibration file if available for automatic m/z range detection

pioneer predict predict_params.json
pioneer convert-raw raw_dir --output-dir arrow_dir --skip-existing
pioneer params-search library.poin ms_data_dir results_dir --params-path=search_params.json
pioneer search search_params.json
```

This sequence builds a predicted spectral library, converts vendor files to Arrow, generates search parameters, and searches the experiment.

!!! tip "Advanced FASTA Input"
    The CLI takes a single path (file or directory). For multiple FASTA sources (directories and/or files),
    edit the `fasta_paths` array in the generated JSON parameter file.

    When using the Julia API directly, you can specify mixed sources at creation:
    ```julia
    params = GetBuildLibParams(out_dir, lib_name,
        ["/path/to/uniprot/", "/custom/proteins.fasta"])
    ```
    Note: CLI users don't specify `lib_name` - it's automatically derived from `out_dir`.

`params-predict` and `params-search` create template JSON files. Edit these
configurations to suit your experiment before running `predict` or `search`.
See [Parameter Configuration](parameters.md) for a description of each option.

## Per-run summary (`run_summary.tsv`)

Every search writes `run_summary.tsv` next to the result tables: one row per
raw file with QC metrics in the spirit of DIA-NN's `report.stats.tsv`.
*Identified* counts everything passing the q-value threshold; *quantified* is
the subset with a positive peak area (the same rule that blanks
`peak_area` in the tables). All distribution statistics are medians over the
identified precursors of that run.

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
| `median_peptide_length`, `median_charge`, `median_missed_cleavages` | Peptide-property medians |
| `ms2_mass_tol_low/high/unit`, `ms1_mass_tol_low/high/unit` | Fragment / precursor mass tolerance the search settled on (`ppm` or `Da` depending on the fitted model; MS1 empty when no MS1 model was fit) |
| `gradient_length_min`, `n_ms1_scans`, `n_ms2_scans` | Run length (last retention time, minutes) and scan counts |
