<img src="figures/PIONEER_LOGO.svg" align="right" width="200px"/>

# Pioneer: Fast and Open-Source Analysis of Data-Independent Acquisition Proteomics Experiments

[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-stable-blue)](https://nwamsley1.github.io/Pioneer.jl/docs/stable/)
[![Regression Reports](https://img.shields.io/badge/regression-reports-orange)](https://nwamsley1.github.io/Pioneer.jl/reports/)
[![Main branch tests](https://img.shields.io/github/actions/workflow/status/nwamsley1/Pioneer.jl/tests.yml?branch=main&label=Main%20tests)](https://github.com/nwamsley1/Pioneer.jl/actions/workflows/tests.yml?query=branch%3Amain)
[![Develop branch tests](https://img.shields.io/github/actions/workflow/status/nwamsley1/Pioneer.jl/tests.yml?branch=develop&label=Develop%20tests)](https://github.com/nwamsley1/Pioneer.jl/actions/workflows/tests.yml?query=branch%3Adevelop)
[![Coverage](https://codecov.io/gh/nwamsley1/Pioneer.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nwamsley1/Pioneer.jl)

Pioneer and its companion tool, Altimeter, are an open-source and performant solution for analysis of protein MS data acquired by data-independent acquisition (DIA). Pioneer includes routines for searching DIA experiments from Thermo and Sciex instruments and for building spectral libraries using the [Koina](https://koina.wilhelmlab.org/) interface. Given a spectral library of precursor fragment ion intensities and retention time estimates, Pioneer identifies and quantifies peptides and protein groups from the library in the data.

**[Documentation](https://nwamsley1.github.io/Pioneer.jl/docs/stable/)** &bull; **[Regression Reports](https://nwamsley1.github.io/Pioneer.jl/reports/)** &bull; **[Landing Page](https://nwamsley1.github.io/Pioneer.jl/)**

## Design Goals

- **Open-Source:** Pioneer is completely open source.
- **Cross-Platform:** Pioneer and the .raw file conversion tool run on Linux, macOS, and Windows
- **High-Performance:** Pioneer achieves high sensitivity, FDR control, and both quantitative precision and accuracy on benchmark datasets
- **Scalability:** Memory consumption and speed remain constant as the number of raw files in an analysis grows. Pioneer scales to very large experiments with hundreds to thousands of raw files (experimental)
- **Fast:** Pioneer searches data several times faster than it can be acquired and faster than state-of-the-art search tools.

## Installation

Download the installer for your operating system from the [releases page](https://github.com/nwamsley1/Pioneer.jl/releases):

| Platform | Installer |
|----------|-----------|
| Windows | `.msi` |
| macOS | `.pkg` (Intel and Apple Silicon) |
| Linux | `.deb` |

The installer adds a `pioneer` command to your `PATH`. On the first run macOS performs a Gatekeeper security check.

```bash
pioneer --help
```

See the [Installation Guide](https://nwamsley1.github.io/Pioneer.jl/docs/stable/user_guide/installation/) for Docker, development setup, and more details.

## Quick Start

```bash
pioneer params-predict lib_dir fasta_dir --params-path=predict_params.json
pioneer predict predict_params.json
pioneer convert-raw raw_dir
pioneer params-search library.poin ms_data_dir results_dir --params-path=search_params.json
pioneer search search_params.json
```

`params-predict` and `params-search` write template JSON files. Review and edit these configurations before running `predict` or `search`. See the [Parameter Configuration](https://nwamsley1.github.io/Pioneer.jl/docs/stable/user_guide/parameters/) guide for available options.

## Current Limitations

- **Variable modifications:** Only oxidation of methionine (Unimod:35) is currently supported as a variable PTM
- **Digestion:** Fully enzymatic digestion only (no semi-enzymatic or non-specific searches)
- **Interface:** Command-line only; no graphical user interface yet

## Regression Tests

Pioneer maintains a public suite of regression tests across diverse instruments and acquisition schemes. Each tagged release triggers end-to-end searches on real datasets and the results are published as browsable reports.

**[Browse regression reports](https://nwamsley1.github.io/Pioneer.jl/reports/)**

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details on our Git Flow workflow and development process.

## Labs
Pioneer is developed in the [Major Lab](https://majorlab.wustl.edu/) and [Goldfarb Lab](https://goldfarblab.wustl.edu/) at Washington University in St. Louis.

<a href="https://majorlab.wustl.edu/"><img src="figures/majorlab.png" width="125px"/></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="https://goldfarblab.wustl.edu/"><img src="figures/goldfarb.png" width="125px"/></a>
<br>

## ASMS 2025
<img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/Pioneer.jpg"/>

## ASMS 2024
<img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/asms_2024_image.jpg"/>
