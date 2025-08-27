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
pioneer convert-raw /path/to/raw/or/folder
```

This subcommand accepts either a single `.raw` file or a directory of files. See the [PioneerConverter repository](https://github.com/nwamsley1/PioneerConverter) for additional options such as thread count and output paths.

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

Subcommands include `search`, `predict`, `params-search`, `params-predict`, `convert-raw`, and `convert-mzml`. The first launch macOS performs a one-time Gatekeeper check.

A minimal end-to-end workflow is:

```bash
# Traditional approach with single FASTA directory
pioneer params-predict lib_dir lib_name fasta_dir --params-path=predict_params.json

# Or with new flexible FASTA input (mixing directories and files)
# Note: CLI support for mixed inputs requires editing the JSON parameter file
pioneer params-predict lib_dir lib_name fasta_dir --params-path=predict_params.json
# Then edit predict_params.json to set fasta_paths to include specific files

pioneer predict predict_params.json
pioneer convert-raw raw_dir
pioneer params-search library.poin ms_data_dir results_dir --params-path=search_params.json
pioneer search search_params.json
```

This sequence builds a predicted spectral library, converts vendor files to Arrow, generates search parameters, and searches the experiment.

!!! tip "Advanced FASTA Input"
    When using the Julia API directly, you can specify mixed FASTA sources:
    ```julia
    params = GetBuildLibParams(out_dir, lib_name, 
        ["/path/to/uniprot/", "/custom/proteins.fasta"])
    ```

`params-predict` and `params-search` create template JSON files. Edit these
configurations to suit your experiment before running `predict` or `search`.
See [Parameter Configuration](parameters.md) for a description of each option.
