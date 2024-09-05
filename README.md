<img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/PIONEER_LOGO.jpg" align="right" width="150px"/>
<h1>Pioneer: Fast and Open-Source Analysis of Data-Indepdendent Aquisition Proteomics Experiments

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://nwamsley1.github.io/Pioneer.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nwamsley1.github.io/Pioneer.jl/dev/)
[![Build Status](https://github.com/nwamsley1/Pioneer.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nwamsley1/Pioneer.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/nwamsley1/Pioneer.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nwamsley1/Pioneer.jl)
</h1>


##  Development Aims
  Pioneer is a cross-platform and open-source tool fully implemented 
in Juilia that identifies and quantifies proteins and peptides from data independent acquisition (DIA) experiments. Given a 
spectral library of fragment ion intensities and retention time estimates on an arbitrary scale, Pioneer employs a spectrum-centric 
algorithm and heuristics to statistically infer the identification status and abundance of each library precursor in the data. We develop Pioneer with the following goals:

- **Open-Source:** Methods should be understood and open to scrutiny by users. 
- **Cross-Platform:** All steps of analysis, including vendor file conversion, should run on all major operating systems.  
- **High-Performance:** The sensitivity, FDR control, and quantitative precision and accuracy should be competitive with state-of-the-art commercial software packages.
- **Fast:** Use of simple heuristics and carefully implemented, efficient algorithms should ensure that data can be analyzed many times faster than it is aquired for typical experiments. 

## Instalation
1) Pioneer requires Julia 1.10. Download [julia](https://pages.github.com/) and add it to the PATH. 
2) Open an instance of the julia REPL
3) Type ";" to activate the shell from within the REPL. Then, navigate to the desired directory and clone the Pioneer.jl repository.
```
shell> git clone https://github.com/nwamsley1/Pioneer.jl.git
```
and not move into the package directory
```
shell> cd Pioneer.jl
```
4) Return to julia by hitting the backspace key. Activate the julia package manager by typing "]" into the REPL and enter the following:
```
(@v1.10) pkg> activate
(@v1.10) pkg> develop ./
(@v1.10) pkg> add ./
```

## Usage
Pioneer exports three the "SearchDIA" method, which takes a single argument, that is a file path to a .json parameters files (examples included). To access these methods
import the Pioneer package from within the Julia REPL. 
```
julia> using Pioneer
```

### File Conversion
`Pioneer` requires Thermo .raw files be converted to an Apache .arrow format with a specific column specification. Use [PioneerConverter](https://github.com/nwamsley1/PioneerConverter) to convert .raw files. 

### SearchDIA

###### .pion Spectral Library
`SearchDIA` requires a properly formated spectral library. The contents include
- detailed-fragments.jld2
- f_index_fragment_bins.arrow
- f_index_fragments.arrow
- f_index_rt_bins.arrow
- precursors_table.arrow
- precursor_to_fragment_indices.jld2
- presearch_f_index_fragment_bins.arrow
- presearch_f_index_fragments.arrow
- presearch_f_index_rt_bins.arrow
- simple_fragments.arrow
```
╭─n.t.wamsley@3225-AD-00020.local ~/RIS_temp/ASMS_2024/ASTRAL_THREE_PROTEOME/unispec_chronologer_1mc_1var_by_052724/spec_lib/pioneer_lib.pion  
╰─➤  ls
config.json                           f_index_rt_bins.arrow                 presearch_f_index_fragments.arrow
detailed_fragments.jld2               precursor_table.arrow                 presearch_f_index_rt_bins.arrow
f_index_fragment_bins.arrow           precursor_to_fragment_indices.jld2    simple_fragments.arrow
f_index_fragments.arrow               presearch_f_index_fragment_bins.arrow
```

## Status
- We are excited to present preiminary results at ASMS 2024 in Anaheim, California! See a copy of the poster below.
- Pioneer is in an early stage of development and not yet ready for use in research. If curious, please contact us at n.t.wamsley@wustl.edu.
- Updates will be continuously added to github as the project progresses.
- Cross-platform conversion of Thermo RAW files to Pioneer compatible Apache Arrow tables. https://github.com/nwamsley1/ThermoRawFileToParquetConverter
  
<h1>Goldfarb Lab </h1>
 Pioneer is developed in the Goldfarb Lab: https://goldfarblab.wustl.edu   <img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/goldfarb.png" align="left" width="125px"/> 
<br><br><br><br><br>

## ASMS 2024
<img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/asms_2024_image.jpg"/>

## US HUPO 2024
<img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/HUPO_POSTER_2024_FORFEDEX.jpg"/>
