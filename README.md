<div align="center">
  <img src="figures/PIONEER_LOGO.svg" alt="Pioneer logo" width="200">

  <h3>Pioneer: Fast and Scalable Analysis of Protein Mass Spectrometry Data</h3>

  <p>
    <a href="https://www.biorxiv.org/content/10.64898/2026.02.16.706201">Paper</a>
    &middot;
    <a href="https://github.com/nwamsley1/Pioneer.jl/releases/latest">Download</a>
    &middot;
    <a href="https://nwamsley1.github.io/Pioneer.jl/docs/stable/">Documentation</a>
    &middot;
    <a href="https://nwamsley1.github.io/Pioneer.jl/reports/">Regression Reports</a>
    &middot;
    <a href="https://github.com/nwamsley1/Pioneer.jl/discussions">Discussions</a>
    &middot;
    <a href="https://nwamsley1.github.io/Pioneer.jl/">Landing Page</a>
  </p>

  <p>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-WashU%20non--commercial-blue" alt="License: WashU non-commercial"></a>
  <a href="https://nwamsley1.github.io/Pioneer.jl/docs/stable/"><img src="https://img.shields.io/badge/docs-stable-blue" alt="Docs"></a>
  <a href="https://nwamsley1.github.io/Pioneer.jl/docs/dev/"><img src="https://img.shields.io/badge/docs-dev-cyan" alt="Dev Docs"></a>
  <a href="https://nwamsley1.github.io/Pioneer.jl/reports/"><img src="https://img.shields.io/badge/regression-reports-orange" alt="Regression Reports"></a>
  <a href="https://github.com/nwamsley1/Pioneer.jl/actions/workflows/tests.yml?query=branch%3Amain"><img src="https://img.shields.io/github/actions/workflow/status/nwamsley1/Pioneer.jl/tests.yml?branch=main&label=Main%20tests" alt="Main branch tests"></a>
  <a href="https://github.com/nwamsley1/Pioneer.jl/actions/workflows/tests.yml?query=branch%3Adevelop"><img src="https://img.shields.io/github/actions/workflow/status/nwamsley1/Pioneer.jl/tests.yml?branch=develop&label=Develop%20tests" alt="Develop branch tests"></a>
  <a href="https://codecov.io/gh/nwamsley1/Pioneer.jl/branch/develop"><img src="https://codecov.io/gh/nwamsley1/Pioneer.jl/branch/develop/graph/badge.svg" alt="Coverage"></a>
  <a href="https://www.biorxiv.org/content/10.64898/2026.02.16.706201"><img src="https://img.shields.io/badge/bioRxiv-2026.02.16.706201-B31B1B" alt="bioRxiv"></a>
  </p>
</div>

Pioneer is a fast and scalable solution for the analysis of protein mass spectrometry data. Pioneer runs through both a graphical and a command-line interface, and it scales to thousands of runs on standard desktop hardware. Pioneer includes routines for searching data-independent acquisition (DIA) experiments from Thermo Fisher, SCIEX, and Bruker instruments and for building *in silico* spectral libraries using the [Koina](https://koina.wilhelmlab.org/) interface. Prebuilt libraries for several proteomes are downloadable through Pioneer. Given a spectral library of precursor fragment ion intensities and retention time estimates, Pioneer identifies and quantifies peptides and protein groups from the library in the data. 

## Features

- **Cross-Platform:** Pioneer runs on Linux, macOS, and Windows.
- **Multi-Vendor:** Supports Thermo Fisher, SCIEX, and Bruker instruments. See [Tested Instruments](https://nwamsley1.github.io/Pioneer.jl/docs/stable/user_guide/instruments/) for the instruments and mass analyzers tested with Pioneer.
- **Fast:** Pioneer searches files many times faster than data acquisition with whole-proteome, *in silico* libraries even with the latest, fast-scanning TOF instruments. 
- **Scalable:** Pioneer may process experiments with hundreds to thousands of runs using standard desktop hardware (32-64GB RAM). Search time grows linearly with the number of MS data files. 
- **Deep Proteome Coverage:** Pioneer achieves high sensitivity and both quantitative precision and accuracy on benchmark datasets covering a variety of instruments and parameters.
- **Regression Testing and FDR Control:** Pioneer maintains a public suite of regression tests across diverse instruments and acquisition schemes. Each tagged release triggers searches on the benchmark datasets, which include entrapment false discovery rate and false transfer rate analyses. **[Browse regression reports](https://nwamsley1.github.io/Pioneer.jl/reports/)**.

## Graphical User Interface

<p align="center">
  <img src="figures/gui_searchdia.png" alt="The Pioneer graphical interface, SearchDIA screen" width="800">
</p>

## Paper

The latest Pioneer manuscript is available on bioRxiv:

*Pioneer and Altimeter: Fast Analysis of DIA Proteomics Data Optimized for Narrow Isolation Windows*. bioRxiv (2026).

[https://www.biorxiv.org/content/10.64898/2026.02.16.706201](https://www.biorxiv.org/content/10.64898/2026.02.16.706201)

## License

Non-commercial research use that is non-clinical and does not involve human subjects is permitted under the [Washington University license](LICENSE). Licenses for commercial, clinical, or human subjects use are available upon request. Please contact Washington University's Office of Technology Management to request a license agreement: otm@wustl.edu.

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

See the [Installation Guide](https://nwamsley1.github.io/Pioneer.jl/docs/stable/user_guide/installation/) for details and development setup, and the [Quick Start Tutorial](https://nwamsley1.github.io/Pioneer.jl/docs/stable/user_guide/quickstart/) for the command line and Docker.

## Labs
Pioneer is developed in the [Major Lab](https://majorlab.wustl.edu/) and [Goldfarb Lab](https://goldfarblab.wustl.edu/) at Washington University in St. Louis.

<a href="https://majorlab.wustl.edu/"><img src="figures/majorlab.png" width="125px"/></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="https://goldfarblab.wustl.edu/"><img src="figures/goldfarb.png" width="125px"/></a>
<br>

