<img src="figures/PIONEER_LOGO.svg" align="right" width="200px"/>

# Pioneer: Fast and Scalable Analysis of Protein Mass Spectrometry Data

[![License: WashU non-commercial](https://img.shields.io/badge/License-WashU%20non--commercial-blue)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-stable-blue)](https://nwamsley1.github.io/Pioneer.jl/docs/stable/)
[![Dev Docs](https://img.shields.io/badge/docs-dev-cyan)](https://nwamsley1.github.io/Pioneer.jl/docs/dev/)
[![Regression Reports](https://img.shields.io/badge/regression-reports-orange)](https://nwamsley1.github.io/Pioneer.jl/reports/)
[![Main branch tests](https://img.shields.io/github/actions/workflow/status/nwamsley1/Pioneer.jl/tests.yml?branch=main&label=Main%20tests)](https://github.com/nwamsley1/Pioneer.jl/actions/workflows/tests.yml?query=branch%3Amain)
[![Develop branch tests](https://img.shields.io/github/actions/workflow/status/nwamsley1/Pioneer.jl/tests.yml?branch=develop&label=Develop%20tests)](https://github.com/nwamsley1/Pioneer.jl/actions/workflows/tests.yml?query=branch%3Adevelop)
[![Coverage](https://codecov.io/gh/nwamsley1/Pioneer.jl/branch/develop/graph/badge.svg)](https://codecov.io/gh/nwamsley1/Pioneer.jl/branch/develop)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2026.02.16.706201-B31B1B)](https://www.biorxiv.org/content/10.64898/2026.02.16.706201)

Pioneer is a fast and scalable solution for the analysis of protein mass spectrometry data. Pioneer includes routines for searching data-independent acquisition (DIA) experiments from Thermo Fisher, Sciex, and Bruker instruments and for building in silico spectral libraries using the [Koina](https://koina.wilhelmlab.org/) interface. Pre-built libraries for several proteomes are downloadable through Pioneer. Given a spectral library of precursor fragment ion intensities and retention time estimates, Pioneer identifies and quantifies peptides and protein groups from the library in the data. 

**[Documentation](https://nwamsley1.github.io/Pioneer.jl/docs/stable/)** &bull; **[Regression Reports](https://nwamsley1.github.io/Pioneer.jl/reports/)** &bull; **[Landing Page](https://nwamsley1.github.io/Pioneer.jl/)**

## Paper

The latest Pioneer manuscript is available on bioRxiv:

*Pioneer and Altimeter: Fast Analysis of DIA Proteomics Data Optimized for Narrow Isolation Windows*. bioRxiv (2026).

[https://www.biorxiv.org/content/10.64898/2026.02.16.706201](https://www.biorxiv.org/content/10.64898/2026.02.16.706201)

## Features

- **Cross-Platform:** Pioneer and the file conversion tool run on Linux, macOS, and Windows.
- **Multi-Vendor:** Supports Thermo Fisher, Sciex, and Bruker instruments.
- **GUI:** Releases include a graphical user interface.
- **Fast:** Pioneer searches files many times faster than data acquisition with whole-proteome, in silico libraries even with the latest, fast-scanning TOF instruments. 
- **Scalable:** Pioneer may process experiments with hundreds to thousands of runs using standard desktop hardware (32-64GB RAM). Search time grows linearly with the number of MS data files. 
- **Deep Proteome Coverage:** Pioneer achieves high sensitivity, FDR control, and both quantitative precision and accuracy on benchmark datasets covering a variety of instruments and parameters.
- **Regression Testing:** Pioneer maintains a public suite of regression tests across diverse instruments and acquisition schemes. Each tagged release triggers end-to-end searches on the benchmark datasets. **[Browse regression reports](https://nwamsley1.github.io/Pioneer.jl/reports/)**.

## License 

Licenses for commercial, human subjects, or clinical use are available upon request. Please contact Washington University's Office of Technology Management to request a license agreement: otm@wustl.edu. 

## Installation

Download the installer for your operating system from the [releases page](https://github.com/nwamsley1/Pioneer.jl/releases):

| Platform | Installer |
|----------|-----------|
| Windows | `.exe` |
| macOS | `.pkg` (Intel and Apple Silicon) |
| Linux | `.deb` |

The installer adds a `pioneer` command to your `PATH`. On the first run macOS performs a Gatekeeper security check.

```bash
pioneer --help
```

See the [Installation Guide](https://nwamsley1.github.io/Pioneer.jl/docs/stable/user_guide/installation/) for Docker, development setup, and more details.

## Labs
Pioneer is developed in the [Major Lab](https://majorlab.wustl.edu/) and [Goldfarb Lab](https://goldfarblab.wustl.edu/) at Washington University in St. Louis.

<a href="https://majorlab.wustl.edu/"><img src="figures/majorlab.png" width="125px"/></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="https://goldfarblab.wustl.edu/"><img src="figures/goldfarb.png" width="125px"/></a>
<br>

