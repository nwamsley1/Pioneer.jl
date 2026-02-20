```@meta
CurrentModule = Pioneer
```

## Introduction

Pioneer and its companion tool Altimeter are an open-source and performant solution for analysis of protein MS data acquired by data-independent acquisition (DIA). Poineer includes routines for searching DIA experments from Thermo and Sciex instruments and for building spectral libraries using the [Koina](https://koina.wilhelmlab.org/) interface. Given a spectral library of precursor fragment ion intensities and retention time estimates, Pioneer identifies and quantifies peptides from the library in the data.

## Key Features

* **Isotope-Aware DIA Analysis**: Narrow isolation windows distort fragment ion isotope distributions because the quadrupole partially transmits precursor isotopic envelopes. Pioneer addresses this by estimating a quadrupole transmission efficiency function for each scan and re-isotoping library spectra accordingly, using methods from [Goldfarb et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6166224/). This correction is critical for accurate matching and quantification in narrow-window DIA.

* **Altimeter: Collision Energy-Independent Spectral Libraries**: Altimeter predicts coefficients for B-splines that model total rather than monoisotopic fragment ion intensities as a function of normalized collision energy (NCE). Evaluating the splines at a given NCE produces a complete spectrum, so a single library works across different instruments and acquisition settings. Pioneer calibrates the optimal NCE per data file automatically.

* **Intensity-Aware Fragment Index**: Pioneer implements a fast fragment index search inspired by [MSFragger](https://pubmed.ncbi.nlm.nih.gov/28394336/) and [Sage](https://pubmed.ncbi.nlm.nih.gov/37819886/). Pioneer's implementation uniquely leverages accurate fragment intensity predictions from in silico libraries—indexing only the highest-ranked fragments—to improve both speed and specificity of candidate identification.

* **Spectral Deconvolution with Robust Regression**: Pioneer explains each observed mass spectrum as a linear combination of template spectra from the library. To reduce quantitative bias from interfering signals in chimeric spectra, Pioneer minimizes the [pseudo-Huber](https://en.wikipedia.org/wiki/Huber_loss) loss rather than squared error. For other examples of linear regression applied to DIA analyses, see [Specter](https://pubmed.ncbi.nlm.nih.gov/29608554/) and [Chimerys](https://www.biorxiv.org/content/10.1101/2024.05.27.596040v2).

* **Dual-Window Quantification**: In narrow-window DIA, a precursor's isotopic envelope is split across adjacent windows. Pioneer normalizes quantification by the isolated precursor fraction and combines signal from adjacent windows for denser chromatographic sampling and improved quantitative accuracy.

* **Match Between Runs**: Pioneer transfers peptide identifications across runs with false transfer rate (FTR) control, increasing coverage in large-scale experiments.

* **Spectral Library Prediction via Koina**: Using [Koina](https://koina.wilhelmlab.org/), Pioneer constructs fully predicted spectral libraries from a FASTA file and an internet connection. Pioneer uses Chronologer for retention time prediction and Altimeter for fragment ion intensity prediction.

## Performance

- **Speed**: 2–6x faster than DIA-NN and AlphaDIA on benchmark datasets
- **FDR Control**: Conservative false discovery rate control validated by entrapment analysis
- **Scalability**: Memory consumption remains constant as the number of raw files grows, scaling to experiments with hundreds of runs

## Current Limitations

- **Variable modifications:** Only oxidation of methionine (Unimod:35) is currently supported as a variable PTM
- **Digestion:** Fully enzymatic digestion only (no semi-enzymatic or non-specific searches)
- **Interface:** Command-line only; no graphical user interface yet

## Quick Links

- [Installation Guide](@ref)
- [Quick Start Tutorial](@ref)

## Authors and Development
Pioneer is developed and maintained by:
- Nathan Wamsley ([Major Lab](https://majorlab.wustl.edu/)/[Goldfarb Lab](https://goldfarblab.wustl.edu/), Washington University)
- Dennis Goldfarb ([Goldfarb Lab](https://goldfarblab.wustl.edu/), Washington University)

## Citation
If you use Pioneer or Altimeter in your research, please cite:

> Wamsley, N. T., Wilkerson, E. M., Major, M., & Goldfarb, D. "Pioneer and Altimeter: Fast Analysis of DIA Proteomics Data Optimized for Narrow Isolation Windows." *bioRxiv* (2025). DOI: [forthcoming]

## Contact
For questions about Pioneer or to collaborate, please contact:
- Nathan Wamsley (wamsleynathan@gmail.com)
- Dennis Goldfarb (dennis.goldfarb@wustl.edu)

For troubleshooting use the [Issues page](https://github.com/nwamsley1/Pioneer.jl/issues) on GitHub. To critique methods or propose features use the [Discussions page](https://github.com/nwamsley1/Pioneer.jl/discussions).

## Exported Methods
```@index
```

```@docs
SearchDIA
```
```@docs
GetSearchParams
```
```@docs
BuildSpecLib
```
```@docs
GetBuildLibParams
```
```@docs
Pioneer.convertMzML
```
