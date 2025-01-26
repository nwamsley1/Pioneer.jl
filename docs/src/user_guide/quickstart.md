# Quick Start Tutorial

## Basic Workflow
Pioneer.jl performs two functions
1. Builds in silico spectral libraries using FASTA files and the [Koina](https://koina.wilhelmlab.org/) server.
2. Searches RAW data-independent aquisition experiments given a spectral library  


## Pioneer Converter
Pioneer.jl cannot search Thermo RAW files directly, and instead searches MS/MS data in from the [Apache Arrow IPC format](https://arrow.apache.org/docs/python/ipc.html). [Pioneer Converter](https://github.com/nwamsley1/PioneerConverter/releases/tag/v0.1.0) converts the vendor files to the the IPC format. The conversion tool accepts either a path to a single .raw file or a directory containing .raw files. In the later case Pioneer Converter converts all .raw files in the directory. To convert .raw files, navigate to the PioneerConverter folder and use the following commands. See [Installation Guide](@ref) for instructions on downloading the Pioneer Converter tool. 
###### POSIX
```
bin/Release/net8.0/PioneerConverter /path/to/raw/file.raw -b 5000

bin/Release/net8.0/PioneerConverter /directory/containing/raw/files/ -b 5000
```

###### Windows
```
cmd bin\Release\net8.0\PioneerConverter.exe \path\to\raw\file.raw -b 5000

cmd bin\Release\net8.0\PioneerConverter.exe \directory\containing\raw\files\ -b 5000
```
!!! info "options"
     * **-n --threads** flag specifies the number of threads to use. Defaults to 2.
     * **-b --batch-size** batch size. Number of scans to convert per-batch. It is not recommended to change this parameter. Setting too high will cause significant memory allocation. Defaults to 10000
     * **-o** path to folder where the converted files will be saved. Defaults to 'arrow_out' directory within the directory of the input.
## Starting Julia
Pioneer runs from within the julia [REPL](https://docs.julialang.org/en/v1/stdlib/REPL/). For optimal performance, it is important to start an instance of julia with multiple threads. It is recommended to set the total number of threads to one fewer than the number of threads available and to then set the number of threads for garbage collection to half of that number. For example, on a desktop computer with 16 threads, the REPL should be oppened as follows: 
```
julia --threads 15 --gcthreads 7,1
```

## Starting Pioneer
If Pioneer has already been installed, then open the REPL and enter the following command.
```@julia
julia> using Pioneer
```  
The Pioneer.jl package exports four methods, `GetBuildLibParams`, `BuildSpecLib`, `GetSearchParams`, and `SearchDIA`. The first two methods build the in silico spectral libraries. The later two search data-independent aquisition (DIA) proteomics data given a spectral library and raw data. 

## Spectral Library Building
Pioneer.jl includes two methods for building spectral libraries. These are `GetBuildLibParams` and `BuildSpecLib`.
### Set Library Building Parameters
`GetBuildLibParams` generates a `.json` formated config file with default parameters for running `BuildSpecLib`. `GetBuildLibParams` requires an output directory path,
a library name, and a path to a folder containing one or more `.fasta` or `fasta.gz` formated files of protein sequences from [UniProt](https://www.uniprot.org/). It inserts these options into the template configuration file, `buildspeclib_parameters.json`.
```@julia
params_path = GetBuildLibParams(
    "/path/to/output",     # Output directory
    "keap1_library",       # Library name
    "/path/to/fasta_files" # FASTA directory
    params_path = "./config.json" #Optional parameter specifying where to write the template configuration file. Defaults to `config.json1` in the working directory
)
```
!!! info "config"
     See [Parameter Configuration](@ref) for a description of each option


### Build Spectral Library

To build a spectral library, simply pass a path to a valid `.json` parameters file to the `BuildSpecLib` function from within the REPL. 
```@julia
# Build the library
julia> BuildSpecLib(params)
```

Given a `json` configuration file, `BuildSpecLib` generates a spectral library in the `.pion` format.  `BuildSpecLib` builds spectral libraries in three steps:

1. In silico digestion the all `.fasta` or `.fasta.gz` files in the directory specified by the configuration file according to a user-defined enzymatic cleavage rule. 
2. Retention time and fragment ion prediction using `chronologer` and a user-specified library prediction tool. Both `chronologer` and the predictions tools are hosted on [Koina](https://koina.wilhelmlab.org/). `BuildSpecLib` currently supports the following spectral library prediction tools:
    * Altimeter (recommended)
    * Prosit 2020 HCD
    * AlphaPeptDeep MS2 Generic
    * Unispec
3. Organization of the predictions into a properly formatted `.poin` spectral library. The output is a folder with the `.poin` extention, which includes the following contents:
```
ExampleLibrary.poin/
├── build_log.txt
├── chronologer_temp
├── config.json
├── detailed_fragments.jld2
├── f_index_fragment_bins.arrow
├── f_index_fragments.arrow
├── f_index_rt_bins.arrow
├── frag_name_to_idx.jld2
├── fragments_table.arrow
├── ion_annotations.jld2
├── prec_to_frag.arrow
├── precursor_to_fragment_indices.jld2
├── precursors.arrow
├── precursors_table.arrow
├── presearch_f_index_fragment_bins.arrow
├── presearch_f_index_fragments.arrow
├── presearch_f_index_rt_bins.arrow
├── raw_fragments.arrow
└── spline_knots.jld2
```

!!! note "'note'"
     `BuildSpecLib` does not generate spectral libraries locally. It uses the Koina servers and will fail without an internet connection. 
!!! note "'note'"
     `BuildSpecLib` performs field and type checking on the paramters file and returns warnings in case of missing parameters or invalid values.


## Search Raw Files
Pioneer.jl includes two methods for searching data-independent aquisition DIA data from .raw files. These are `GetSearchParams` and `SearchDIA`. 

### Set Library Building Parameters
`GetSearchParams` generates a `.json` formated config file with default parameters for running `SearchDIA`. `GetSearchParams` requires three inputs:
1. A path to the `.poin` formatted spectral library. 
2. A path to a folder containing one or more `.arrow` formatted MS/MS runs to search. 
3. The output directory.
GetSearchParams inserts these options into the template configuration file, `search_parameters.json`.


```@julia
#Generate Search Parameters 
params_path = GetSearchParams(
    "/path/to/speclib.poin",
    "/path/to/ms/data/dir",
    "/path/to/output/dir",
    params_path = "./config.json" #Optional parameter specifying where to write the template configuration file. Defaults to `config.json1` in the working directory
)
```
!!! info "config"
     See [Parameter Configuration](@ref) for a description of each option

### Search DIA
To search raw files, pass a path to a valid `.json` paramters file to the SearchDIA method from within the REPL. 

```@julia
# Run the search
SearchDIA(search_params)
```
The output folder contains the following:
```
results_dir/
├── config.json
├── pioneer_search_log.txt
├── qc_plots/
│   ├── collision_energy_alignment/
│   │   └── nce_alignment_plots.pdf
│   ├── quad_transmission_model/
│   │   ├── quad_data
│   │   │   └── quad_data_plots.pdf
│   │   └── quad_models
│   │       └── quad_model_plots.pdf
│   ├── rt_alignment_plots/
│   │   └── rt_alignment_plots.pdf
│   ├── mass_error_plots/
│   │   └── mass_error_plots.pdf
│   └── QC_PLOTS.pdf
├── precursors_long.arrow
├── precursors_long.tsv
├── precursors_wide.arrow
├── precurosrs_wide.tsv
├── protein_groups_long.arrow
├── protein_groups_long.tsv
├── protein_groups_wide.arrow
└── protein_groups_wide.tsv
```
## Example

In a completed example workflow, the user provides a relative path to a folder, `./test_data/fasta`, containing a fasta file(s) for the E. coli proteome. `BuildSpecLib` builds the spectral library `ecoli_test_lib.pion` in the `./libraries` folder. The user searches the raw files in `./test_data/ms_data` using the `ecoli_test_lib.pion` library and saves the output to the `.results` folder. 

```@julia
julia> using Pioneer

julia> # Build spectral library
julia> params_path = GetBuildLibParams(
    joinpath(pwd(), "libraries"),
    "ecoli_test_lib",
    joinpath(pwd(), "test_data", "fasta")
)
julia> #Edit paramters file in a text editor if needed 
julia> BuildSpecLib(params_path)

julia> # Run DIA search
julia> params_path = GetSearchParams(
    joinpath(pwd(), "libraries", "ecoli_test_lib.pion"),
    joinpath(pwd(), "test_data", "ms_data"),
    joinpath(pwd(), "results")
)
julia> #Edit paramters file if needed 
julia> SearchDIA(params_path)
```