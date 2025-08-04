# Quick Start Tutorial

## Basic Workflow
Pioneer performs two major steps:
1. Build in silico spectral libraries using FASTA files and the [Koina](https://koina.wilhelmlab.org/) server.
2. Search DIA experiments using a spectral library and the MS data files.

## Pioneer Converter
Pioneer operates on MS/MS data stored in the [Apache Arrow IPC format](https://arrow.apache.org/docs/python/ipc.html).
Use the bundled PioneerConverter via the CLI to convert Thermo RAW files:

```bash
pioneer convert-raw /path/to/raw/or/folder
```

This subcommand accepts either a single `.raw` file or a directory of files. See the PioneerConverter repository for additional options such as thread count and output paths.

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

Subcommands include `search`, `predict`, `params-search`, `params-predict`, `convert-raw`, and `convert-mzml`. The first launch may download dependencies on Windows or Linux, while macOS performs a one-time Gatekeeper check.

A minimal end-to-end workflow is:

```bash
pioneer params-predict out_dir lib_name fasta_dir
pioneer predict buildspeclib_params.json
pioneer convert-raw raw_dir
pioneer params-search library.poin ms_data results
pioneer search search_parameters.json
```

This sequence builds a predicted spectral library, converts vendor files to Arrow, generates search parameters, and searches the experiment.

`params-predict` and `params-search` create template JSON files. Edit these
configurations to suit your experiment before running `predict` or `search`.
See [Parameter Configuration](parameters.md) for a description of each option.
