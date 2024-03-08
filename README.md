
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
- **Open-Source:** Operating under an MIT license all source code should be freely available and sufficiently documented. Methods should be understood and open to scrutiny by the proteomics community
- **Cross-Platform:** All steps of analysis, including vendor file conversion, should run all major operating systems.  
- **High-Performance:** The by sensitivity, FDR control, and quantitative precision and accuracy, should be competitive with state-of-the-art commercial software packages
- **Fast:** Use of simple heuristics and carefully implemented, efficient algorithms should ensure that data can be analyzed several times faster than it is aquired for the great majority of experiments. 

## Status
- We are exited to present preiminary results at US HUPO 2024 in Portland, Oregon! See a copy of the poster below.
- Pioneer is in an early stage of development and not yet ready for If curious, please contact us at n.t.wamsley@wustl.edu.
- Updates will be continuously added to github as the project progresses.
- Cross-platgorm onversion of Thermo RAW files to Pioneer compatible Apache Arrow tables. https://github.com/nwamsley1/ThermoRawFileToParquetConverter

## US HUPO 2024
<img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/HUPO_POSTER_2024_FORFEDEX.jpg"/>
