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
`SearchDIA` requires a properly formated spectral library. Spectral libraries are contained in folders with the `.pion` extension. The contents include the following. 
```
╭─n.t.wamsley@3225-AD-00020.local ~/RIS_temp/ASMS_2024/ASTRAL_THREE_PROTEOME/unispec_chronologer_1mc_1var_by_052724/spec_lib/pioneer_lib.pion  
╰─➤  ls
config.json                           f_index_rt_bins.arrow                 presearch_f_index_fragments.arrow
detailed_fragments.jld2               precursor_table.arrow                 presearch_f_index_rt_bins.arrow
f_index_fragment_bins.arrow           precursor_to_fragment_indices.jld2    simple_fragments.arrow
f_index_fragments.arrow               presearch_f_index_fragment_bins.arrow
```
- detailed_fragments.jld2
```
julia> spec_lib["f_det"].frags[1:5]
5-element Vector{DetailedFrag{Float32}}:
DetailedFrag{Float32}(0x00000001, 329.6432f0, Float16(1.0), 0x02, false, 0x02, 0x06, 0x02, 0x01, 0x00)
DetailedFrag{Float32}(0x00000001, 543.2522f0, Float16(0.976), 0x02, false, 0x01, 0x05, 0x02, 0x02, 0x00)
DetailedFrag{Float32}(0x00000001, 358.1539f0, Float16(0.809), 0x02, false, 0x02, 0x07, 0x02, 0x03, 0x00)
DetailedFrag{Float32}(0x00000001, 658.2791f0, Float16(0.4597), 0x02, false, 0x01, 0x06, 0x02, 0x04, 0x00)
DetailedFrag{Float32}(0x00000001, 456.2201f0, Float16(0.4324), 0x02, false, 0x01, 0x04, 0x02, 0x05, 0x00)

julia> spec_lib["f_det"].prec_frag_ranges[1:5]
5-element Vector{UnitRange{UInt32}}:
0x00000001:0x00000013
0x00000014:0x00000027
0x00000028:0x00000039
0x0000003a:0x0000004d
0x0000004e:0x0000005c
```
- f_index_fragment_bins.arrow
```
Julia> spec_lib["f_index"].fragment_bins
1187912-element Arrow.Struct{FragIndexBin, Tuple{Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{UInt32, Vector{UInt32}}}, (:lb, :ub, :first_bin, :last_bin)}:
FragIndexBin{Float32}(172.0717f0, 172.0717f0, 0x00000001, 0x00000017)
FragIndexBin{Float32}(173.5953f0, 173.5953f0, 0x00000018, 0x00000018)
FragIndexBin{Float32}(179.6079f0, 179.6079f0, 0x00000019, 0x00000019)
FragIndexBin{Float32}(181.0713f0, 181.0713f0, 0x0000001a, 0x0000001a)
```
- f_index_fragments.arrow
```
julia> spec_lib["f_index"].fragments
62231372-element Arrow.Struct{IndexFragment, Tuple{Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{UInt8, Vector{UInt8}}, Arrow.Primitive{UInt8, Vector{UInt8}}}, (:prec_id, :prec_mz, :score, :charge)}:
IndexFragment{Float32}(0x00007028, 379.6988f0, 0x02, 0x02)
IndexFragment{Float32}(0x000009f8, 380.6941f0, 0x01, 0x02)
IndexFragment{Float32}(0x00000a69, 381.6818f0, 0x01, 0x02)
IndexFragment{Float32}(0x00000b20, 383.6736f0, 0x01, 0x02)
```
- f_index_rt_bins.arrow
These make up the fragment index for the initial search. 
```
julia> spec_lib["f_index"].rt_bins
20-element Arrow.Struct{FragIndexBin, Tuple{Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{UInt32, Vector{UInt32}}}, (:lb, :ub, :first_bin, :last_bin)}:
FragIndexBin{Float32}(-1.3982847f0, 0.60170865f0, 0x00000001, 0x00006f72)
FragIndexBin{Float32}(0.6017635f0, 2.6017632f0, 0x00006f73, 0x000162e7)
FragIndexBin{Float32}(2.6017723f0, 4.6017694f0, 0x000162e8, 0x00027f1d)
FragIndexBin{Float32}(4.6017733f0, 6.601769f0, 0x00027f1e, 0x0003b678)
```
- precursors_table.arrow
A table with one row per precursor in the library. Each precursor has a unique id that corresponds to a row in this table. 
```
julia> precursors
Arrow.Table with 8893001 rows, 12 columns, and schema:
:irt                   Float32
:mz                    Float32
:is_decoy              Bool
:proteome_identifiers  String
:accession_numbers     String
:sequence              String
:structural_mods       String
:isotopic_mods         String
:prec_charge           UInt8
:missed_cleavages      UInt8
:length                UInt8
:sulfur_count          UInt8
```
- precursor_to_fragment_indices.jld2
- presearch_f_index_fragment_bins.arrow
- presearch_f_index_fragments.arrow
- presearch_f_index_rt_bins.arrow
- simple_fragments.arrow

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
