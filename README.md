<img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/PIONEER_LOGO.jpg" align="right" width="150px"/>
<h1>Pioneer: Fast and Open-Source Analysis of Data-Independent Acquisition Proteomics Experiments

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://nwamsley1.github.io/Pioneer.jl/dev)
[![Build Status](https://github.com/nwamsley1/Pioneer.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nwamsley1/Pioneer.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/nwamsley1/Pioneer.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nwamsley1/Pioneer.jl)
</h1>

## Documentation
See [documentation for full installation and usage instructions.](https://nwamsley1.github.io/Pioneer.jl/dev)

## Installation and Usage
Download the installer for your operating system from the [releases page](https://github.com/nwamsley1/Pioneer.jl/releases). The installer adds a `pioneer` command to your `PATH`. On the first run, Windows and Linux download IntelOpenMP and MKL; macOS performs a Gatekeeper security check.

```bash
pioneer --help
```
Lists subcommands such as `predict`, `params-predict`, `search`, `params-search`, `convert-raw`, and `convert-mzml`.
A typical workflow is:

```bash
pioneer params-predict out_dir lib_name fasta_dir
pioneer predict buildspeclib_params.json
pioneer convert-raw raw_dir
pioneer params-search library.poin ms_data results
pioneer search search_parameters.json
```

`params-predict` and `params-search` write template JSON files. Review and edit these
configurations before running `predict` or `search`. See the
[Parameter Configuration](https://nwamsley1.github.io/Pioneer.jl/dev/user_guide/parameters/)
guide for available options.

## Development
To work on Pioneer itself:

```bash
git clone https://github.com/nwamsley1/Pioneer.jl.git
cd Pioneer.jl
julia --project=dev
```

```julia
pkg> develop ./
julia> using Revise, Pioneer
```

Install [PioneerConverter](https://github.com/nwamsley1/PioneerConverter) to convert Thermo RAW files to Arrow format before calling functions like `GetBuildLibParams`, `BuildSpecLib`, `GetSearchParams`, and `SearchDIA` directly.

| Subcommand       | Julia function   |
|------------------|------------------|
| `params-predict` | `GetBuildLibParams` |
| `predict`        | `BuildSpecLib`     |
| `params-search`  | `GetSearchParams`  |
| `search`         | `SearchDIA`        |
| `convert-raw`    | `PioneerConverter` |
| `convert-mzml`   | `convertMzML`      |

Developers can invoke these functions directly:

```julia
params = GetBuildLibParams(out_dir, lib_name, fasta_dir)
BuildSpecLib(params)
params = GetSearchParams("library.poin", "ms_data", "results")
SearchDIA(params)
```

## Introduction

Pioneer and its companion tool, Altimeter, are together an open-source and performant solution for analysis of protein MS data acquired by data-independent acquisition (DIA). Poineer includes routines for searching DIA experments from Thermo and Sciex instruments and for building spectral libraries using the [Koina](https://koina.wilhelmlab.org/) interface. Given a spectral library of precursor fragment ion intensities and retention time estimates, Pioneer identifies and quantifies peptides and protein groups from the library in the data. 

## Design Goals

- **Open-Source:** Pioneer is completely open source. 
- **Cross-Platform:** Pioneer and the vendor-specific file conversion tool run on Linux, MacOS, and Windows
- **High-Performance:** Pioneer achieves high sensitivity, FDR control, and both quantitative precision and accuracy on benhcmark datat-sets 
- **Scalability:** Memory consumption and speed should remain constant as the number of raw files in an anslysis grows. Pioneer should scale to very large experiments with hundreds to thousands of raw files (experimental)
- **Fast:** Pioneer searches data several times faster than it can be aquired and faster than state-of-the-art search tools.
<h1>Goldfarb Lab </h1>
 Pioneer is developed in the Goldfarb Lab: https://goldfarblab.wustl.edu   <img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/goldfarb.png" align="left" width="125px"/> 
<br><br><br><br><br>

## ASMS 2025 
<img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/Pioneer.jpg"/>

## ASMS 2024
<img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/asms_2024_image.jpg"/>
