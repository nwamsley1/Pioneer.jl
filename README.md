<img src="figures/PIONEER_LOGO.svg" align="right" width="200px"/>
<h1>Pioneer: Fast and Open-Source Analysis of Data-Independent Acquisition Proteomics Experiments

[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue)](LICENSE)
[![Docs: dev](https://img.shields.io/badge/docs-dev-blue)](https://nwamsley1.github.io/Pioneer.jl/dev)
[![Main branch tests](https://img.shields.io/github/actions/workflow/status/nwamsley1/Pioneer.jl/tests.yml?branch=main&label=Main%20tests)](https://github.com/nwamsley1/Pioneer.jl/actions/workflows/tests.yml?query=branch%3Amain)
[![Develop branch tests](https://img.shields.io/github/actions/workflow/status/nwamsley1/Pioneer.jl/tests.yml?branch=develop&label=Develop%20tests)](https://github.com/nwamsley1/Pioneer.jl/actions/workflows/tests.yml?query=branch%3Adevelop)
[![Coverage](https://codecov.io/gh/nwamsley1/Pioneer.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nwamsley1/Pioneer.jl)
</h1>

## Introduction

Pioneer and its companion tool, Altimeter, are together an open-source and performant solution for analysis of protein MS data acquired by data-independent acquisition (DIA). Poineer includes routines for searching DIA experments from Thermo and Sciex instruments and for building spectral libraries using the [Koina](https://koina.wilhelmlab.org/) interface. Given a spectral library of precursor fragment ion intensities and retention time estimates, Pioneer identifies and quantifies peptides and protein groups from the library in the data. 

## Design Goals

- **Open-Source:** Pioneer is completely open source. 
- **Cross-Platform:** Pioneer and the .raw file conversion tool run on Linux, MacOS, and Windows
- **High-Performance:** Pioneer achieves high sensitivity, FDR control, and both quantitative precision and accuracy on benhcmark datat-sets 
- **Scalability:** Memory consumption and speed should remain constant as the number of raw files in an analysis grows. Pioneer should scale to very large experiments with hundreds to thousands of raw files (experimental)
- **Fast:** Pioneer searches data several times faster than it can be aquired and faster than state-of-the-art search tools.

## Documentation
See [documentation for full installation and usage instructions.](https://nwamsley1.github.io/Pioneer.jl/dev)

## Installation
Download the installer for your operating system from the [releases page](https://github.com/nwamsley1/Pioneer.jl/releases). The installer adds a `pioneer` command to your `PATH`.

```bash
pioneer --help
```
Lists subcommands such as `predict`, `params-predict`, `search`, `params-search`, `convert-raw`, and `convert-mzml`.

On the first run macOS performs a Gatekeeper security check.

## Quick Start
A minimal end-to-end workflow is:

```bash
pioneer params-predict lib_dir fasta_dir --params-path=predict_params.json
pioneer predict predict_params.json
pioneer convert-raw raw_dir
pioneer params-search library.poin ms_data_dir results_dir --params-path=search_params.json
pioneer search search_params.json
```

`params-predict` and `params-search` write template JSON files. Review and edit these configurations before running `predict` or `search`. See the
[Parameter Configuration](https://nwamsley1.github.io/Pioneer.jl/dev/user_guide/parameters/)
guide for available options.


### Docker

Pioneer can also be run from a Docker container:

```bash
docker pull dennisgoldfarb/pioneer:latest
```

Run Pioneer inside the container, mounting a host directory to access your data:

```bash
docker run --rm -it -v /path/on/host:/work dennisgoldfarb/pioneer:latest pioneer --help
```

Replace `/path/on/host` with the directory containing your data and `pioneer --help` with any Pioneer subcommand. The repository includes a `Dockerfile` for building the image locally:

```bash
docker build -t pioneer .
```


## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details on our Git Flow workflow and development process.

<h1>Goldfarb Lab </h1>
 Pioneer is developed in the Goldfarb Lab: https://goldfarblab.wustl.edu   <img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/goldfarb.png" align="left" width="125px"/> 
<br><br><br><br><br>

## ASMS 2025 
<img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/Pioneer.jpg"/>

## ASMS 2024
<img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/asms_2024_image.jpg"/>
