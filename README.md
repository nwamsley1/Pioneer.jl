
<img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/PIONEER_LOGO.jpg" align="right" width="150px"/>
<h1>Pioneer: Fast and Open-Source Analysis of Data-Indepdendent Aquisition Proteomics Experiments<br><br><h1>

  
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://nwamsley1.github.io/Titus.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nwamsley1.github.io/Titus.jl/dev/)
[![Coverage](https://codecov.io/gh/nwamsley1/Titus.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nwamsley1/Titus.jl)
<!---
[![Build Status](https://github.com/nwamsley1/Titus.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nwamsley1/Titus.jl/actions/workflows/CI.yml?query=branch%3Amain)
--->
##  Development Aims
  Pioneer is a cross-platform and open-source tool fully implemented 
in Juilia that identifies and quantifies proteins and peptides from data independent acquisition (DIA) experiments. Given a 
spectral library of fragment ion intensities and retention time estimates on an arbitrary scale, Pioneer employs a spectrum-centric 
algorithm and heuristics to statistically infer the identification status and abundance of each library precursor in the data.  

Pioneer :    
- **Open-Source** Operating under an MIT license all source code should be freely available and sufficiently documented. Methods should be understood and open to scrutiny by the proteomics community
- **Cross-Platform** All steps of analysis, including vendor file conversion, should run all major operating systems.  
- **** The number of identifications, proper FDR, control as state-of-the-art commercial software packages
- **Fast** 

## Status
- We are exited to present preiminary results at US HUPO . 
- Accepts raw ms data in the Apache Arrow format. See the following for cross-platform conversion of Thermo .raw files to the arrow format https://github.com/nwamsley1/ThermoRawFileToParquetConverter.
- Analysis of survey methods. Given a list table of protein-peptide pairs, identifies the best charge state for each precursor (by XTandem hyperscore), the best transitions, and the MS1 peak height. If the survey analyses are split accross multiple experiments, these can be analyzed at once and combined. In addition, can run survey analyses at multiple collision energies/FAIMS CV's to identify the optimum for each analyte. Output is given in a format freindly to XCalibur method editor for Thermo Tribrid instruments.
- Supports variable and fixed modifications defined by regular expressions and includes examples. 
- Estimates peak area ratios using an MM-Estimator (https://github.com/getzze/RobustModels.jl/blob/main/docs/make.jl). Enables accurate par estimation in the pressence of noisy or interfered transitions. High uncertainty in estimation can be used as grounds for exclusion. 
- Summarizaiton of peptide-level quantitation to protein-level quantitation using the MaxLFQ Algorithm without normalization [2].
- Generates a multi-page pdf for each experiment file including chromatogram plots.

## References
<a id="1">[1]</a> 
Gallien S, Kim SY, Domon B. Large-Scale Targeted Proteomics Using Internal Standard Triggered-Parallel Reaction Monitoring (IS-PRM). Mol Cell Proteomics. 2015 Jun;14(6):1630-44. doi: 10.1074/mcp.O114.043968. Epub 2015 Mar 9. PMID: 25755295; PMCID: PMC4458725.
<br>
<a id="1">[2]</a> 
Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep;13(9):2513-26. doi: 10.1074/mcp.M113.031591. Epub 2014 Jun 17. PMID: 24942700; PMCID: PMC4159666
<br>
<a id="1">[3]</a> 
Stopfer LE, Flower CT, Gajadhar AS, Patel B, Gallien S, Lopez-Ferrer D, White FM. High-Density, Targeted Monitoring of Tyrosine Phosphorylation Reveals Activated Signaling Networks in Human Tumors. Cancer Res. 2021 May 1;81(9):2495-2509. doi: 10.1158/0008-5472.CAN-20-3804. Epub 2021 Jan 28. PMID: 33509940; PMCID: PMC8137532.
<br>
<a id="1">[4]</a> 
Wamsley et al. Targeted proteomic quantitation of NRF2 signaling and predictive biomarkers in HNSCC. https://doi.org/10.1101/2023.03.13.532474 
