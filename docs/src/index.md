```@meta
CurrentModule = Pioneer
```

## Introduction

Pioneer and its companion tool Altimeter are an open-source and performant solution for analysis of protein MS data acquired by data-independent acquisition (DIA). Poineer includes routines for searching DIA experments from Thermo and Sciex instruments and for building spectral libraries using the [Koina](https://koina.wilhelmlab.org/) interface. Given a spectral library of precursor fragment ion intensities and retention time estimates, Pioneer identifies and quantifies peptides from the library in the data. 

## Design Goals

- **Open-Source:** Pioneer is completely open source. 
- **Cross-Platform:** Pioneer and the vendor-specific file conversion tool run on Linux, MacOS, and Windows
- **High-Performance:** Pioneer achieves high sensitivity, FDR control, quantitative precision and accuracy on benhcmark datat-sets 
- **Scalability:** Memory consumption and speed should remain constant as the number of raw files in an anslysis grows. Pioneer should scale to very large experiments with hundreds to thousands of raw files (experimental)
- **Fast:** Pioneer searches data several times faster than it can be aquired and faster than state-of-the-art search tools.

## Features
Pioneer and Altimeter build on previous search engines and introduce several new concepts:

* **Spectral Library Prediction with Koina**: Using [Koina](https://koina.wilhelmlab.org/) Pioneer can construct fully predicted spectral libraries given an internet connection and a FASTA file with protein sequences. Pioneer uses Chronologer to predict peptide retention times and is optimized to use Altimeter for fragment ion intensity predictions.

* **Collision Energy Independent Spectral Libraries**: Rather than predicting a single intensity value for each fragment ion, Altimeter predicts 4 B-spline coefficients. Evaluating the fragment splines at a given collision energy gives a fragment ion intensity. Pioneer calibrates the library to find the optimal collision energy value to use for each MS data file in an experiment. In this way, it is possible to use a single spectral library for different instruments and scan settings. 

* **Fragment Isotope Correction**: Fragment isotope distributions depend on precursor isotope distributions as distorted by quadrupole mass filtering. Altimeter addresses this by predicting total fragment ion intensities rather than monoisotopic ones. Pioneer then accurately re-isotopes these library spectra using methods from [Goldfarb et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6166224/). This is particularly important for narrow-window DIA methods where precursor isotopic envelopes frequently straddle multiple windows.

* **Qaudrupole Transmission Modeling**: For narrow-window data, Pioneer can optionally estimate a quadrupole-transmission efficiency function for more accurate re-isotoping. 

* **Intensity-Aware Fragment Index Search**: Pioneer implements a fast fragment index search inspired by [MSFragger](https://pubmed.ncbi.nlm.nih.gov/28394336/) and [Sage](https://pubmed.ncbi.nlm.nih.gov/37819886/). Pioneer's implementation uniquely leverages accurate fragment intensity predictions from in silico libraries to improve both speed and specificity of the search.

* **Linear Regression onto Library Templates**: Pioneer explains each observed mass spectrum as a linear combination of template spectra from the library. To reduce quantitative bias from interfering signals, Pioneer minimizes the [pseudo-Huber](https://en.wikipedia.org/wiki/Huber_loss) loss rather than squared error. This provides robust quantification even in complex spectra. For other examples of linear regression applied to DIA analyses, see [Specter](https://pubmed.ncbi.nlm.nih.gov/29608554/) and [Chimerys](https://www.biorxiv.org/content/10.1101/2024.05.27.596040v2).

* **Scalability**: Pioneer was designed to scale to large experiments with many MS data files. Memory consumption remains constant as the number of data files in an experiment grows large. 

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
```@docs
ParseSpecLib
```
```@docs
GetParseSpecLibParams
```
```@docs
Pioneer.convertMzML
```
