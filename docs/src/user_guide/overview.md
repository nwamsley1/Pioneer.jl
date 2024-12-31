# [User Guide](@id user-guide)

## Introduction 
Pioneer.jl is a search engine for processing data-independent aquisition proteomics data. given a spectral library of precursor fragment ion intensities and retention time estimates, Pioneer identifies and quantifies peptides from the library in the data. 

## Design Goals
Pioneer was designed with the following goals: 

- **Open-Source:** Methods should be understood and open to scrutiny by users
- **Cross-Platform:** All steps of analysis, including vendor file conversion, should run on all major operating systems
- **High-Performance:** The sensitivity, FDR control, and quantitative precision and accuracy should be competitive with state-of-the-art commercial software packages
- **Scalability:** Should scale to very large experiments with hundreds to thousands of raw files
- **Fast:** Use of simple heuristics and carefully implemented, efficient algorithms should ensure that data can be analyzed many times faster than it is aquired for typical experiments

## Features 
Pioneer combines elements from previous sear. A fragment-index search (MSFragger and Sage) that 