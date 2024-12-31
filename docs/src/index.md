```@meta
CurrentModule = Pioneer
```

## Introduction

Pioneer.jl is a search engine for processing data-independent aquisition (DIA) proteomics data. Poineer includes routines for searching DIA experments from Thermo and Sciex instruments and for building spectral libraries using the [Koina](https://koina.wilhelmlab.org/) interface. Given a spectral library of precursor fragment ion intensities and retention time estimates, Pioneer identifies and quantifies peptides from the library in the data. 

## Design Goals
Pioneer was desinged to be an open source and performant solution for DIA analyses that third-parties can easily adapt to special use cases and experimental techniques. 

- **Open-Source:** Pioneer is completely open source. 
- **Cross-Platform:** Pioneer and the vendor-specific file conversion tool run on Linux, MacOS, and Windows
- **High-Performance:** Pioneer achieves high sensitivity, FDR control, quantitative precision and accuracy on benhcmark datat-sets 
- **Scalability:** Memory consumption and speed should remain constant as the number of raw files in an anslysis grows. Pioneer should scale to very large experiments with hundreds to thousands of raw files (experimental)
- **Fast:** Pioneer searches data several times faster than it can be aquired and faster than state-of-the-art search tools in our benchmarks. 

## Features
Pioneer builds on previous search engines and introduces several new concepts:

* **Intensity-Aware Fragment Index Search** Pioneer implements a fragment index search (see [MSFragger](https://pubmed.ncbi.nlm.nih.gov/28394336/) and [Sage](https://pubmed.ncbi.nlm.nih.gov/37819886/)) to quickly limit the precursor search space for more intensive computation. The Pioneer implementation leverages fragment intensity information from in silico spectral libraries for faster and more specific indexing. 
* **Linear Regrssion onto Library Templates** Pioneer regresses observed mass spectra onto template fragment spectra from in an silico spectral library. To reduce quantitative bias from interfering signals, Pioneer minimizes the [pseudo-Huber](https://en.wikipedia.org/wiki/Huber_loss) loss function rather than squared error. The pseudo-Huber loss behaves like absolute error for large errors and like squared error for small errors. For other examples of linear regression as applied to DIA analyses, see [Specter](https://pubmed.ncbi.nlm.nih.gov/29608554/) and [Chimerys](https://www.biorxiv.org/content/10.1101/2024.05.27.596040v2).
* **Fragment Isotope Correction** Fragment isotope distributions are conditional on the precursor isotope distributions as distorted by a quadrupole mass filter. To account for this, Pioneer uses a spectral library tool, Altimeter, that predicts total fragment ion abundances rather than mono-isotopic abundances. Pioneer uses methods from [Goldfarb et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6166224/) to re-isotope these library spectra. For narrow-window data, Pioneer can optionally estimate a quadrupole-transmission efficiency function for more accurate re-isotoping. 


## Quick Links

- [Installation Guide](@ref)
- [Quick Start Tutorial](@ref)

## Authors and Development
Pioneer is developed and maintained by:
- Nathan Wamsley ([Major Lab](https://majorlab.wustl.edu/)/[Goldfarb Lab](https://goldfarblab.wustl.edu/), Washington University)
- Dennis Goldfarb ([Goldfarb Lab](https://goldfarblab.wustl.edu/), Washington University)

## Contact
For questions about Pioneer or to collaborate, please contact:
- Nathan Wamsley (wamsleynathan@gmail.com)
- Dennis Goldfarb (dennis.goldfarb@wustl.edu)

For toubleshooting use the [Issues page](https://github.com/nwamsley1/Pioneer.jl/issues) on github. To critique methods or propose features use the [Discussions page](https://github.com/nwamsley1/Pioneer.jl/discussions).

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
