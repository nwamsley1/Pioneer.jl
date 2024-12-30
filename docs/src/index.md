```@meta
CurrentModule = Pioneer
```

# Pioneer.jl

Pioneer.jl is a Julia package for analyzing Data-Independent Acquisition (DIA) mass spectrometry data. Poineer includes routines for searching DIA experments from Thermo and Sciex instruments and for building spectral libraries using the [Koina](https://koina.wilhelmlab.org/) interface. Pioneer supports the Altimeter spectral library generation tool. 

## Features and Aims

- **Open-Source:** Pioneer is completely open source
- **Cross-Platform:** Pioneer and the vendor-specific file conversion tool run on Linux, MacOS, and Windows
- **High-Performance:** Pioneer achieves high sensitivity, FDR control, quantitative precision and accuracy on benhcmark datat-sets 
- **Scalability:** Should scale to very large experiments with hundreds to thousands of raw files
- **Fast:** Pioneer searches data several times faster than it can be aquired. 

## Quick Links

- [Installation Guide](@ref)
- [Quick Start Tutorial](@ref)
- [API Reference](@ref api-reference)
- [Algorithm Documentation](@ref algorithm-documentation)

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
