<div align="center">
  <img src="figures/PIONEER_LOGO.svg" alt="Pioneer logo" width="200">

  <h3>Pioneer: Fast and Scalable Analysis of Protein Mass Spectrometry Data</h3>

  <p>
    <a href="https://nwamsley1.github.io/Pioneer.jl/"><img src="figures/button-website-download.svg" alt="Website &amp; Download" height="28"></a>
    <a href="https://nwamsley1.github.io/Pioneer.jl/docs/stable/"><img src="figures/button-documentation.svg" alt="Documentation" height="28"></a>
  </p>
</div>

Pioneer is a fast solution for the quantitative analysis of data-independent acquisition (DIA) proteomics data from Thermo Fisher, SCIEX, and Bruker instruments. The Pioneer software can process large experiments of thousands of MS data files on desktop hardware. Each release is tested on a suite of benchmark datasets and tests, which include protein group and precursor identification counts, quantitative precision and accuracy analysis, and entrapment false discovery rate and false transfer rate experiments.

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

