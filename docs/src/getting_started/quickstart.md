# Quick Start Tutorial

## Basic Workflow
Pioneer takes a spectral library, raw files, and a configuration file as inputs. The workflow has three steps. 

1. Convert raw files
2. Generate a spectral library
3. Run DIA search


## Pioneer Converter
[Pioneer Converter](https://github.com/nwamsley1/PioneerConverter/releases/tag/v0.1.0) accepts either be a path to a single .raw file or a directory containing .raw files. In the later case, Pioneer Converter converts all .raw files in the directory. The output files
are written into a new directory, 'arrow_out', created within the input directory. 

* **-n --threads** flag specifies the number of threads to use. Defaults to 2.
* **-b --batch-size** batch size. Number of scans to convert per-batch. Setting to high will cause significant memory allocation. Defaults to 10000
* **-o** path to folder where the converted files will be saved. Defaults to 'arrow_out' directory within the directory of the input.

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


## Spectral Library Building

### Set Library Building Parameters
Pioneer.jl includes two methods for building spectral libraries, `GetBuildLibParams` and `BuildSpecLib`. `GetBuildLibParams` generates
a `.json` formated configuration with default parameters for running `BuildSpecLib`. `GetBuildLibParams` requires an output directory path,
a library name, and a path to a folder containing one or more `.fasta` or `fasta.gz` formated files of protein sequences from [UniProt](https://www.uniprot.org/). GetBuildLibParams
writes the `build_parameters.json` template configuration file to the current working directory.

```@julia
params_path = GetBuildLibParams(
    "/path/to/output",     # Output directory
    "keap1_library",       # Library name
    "/path/to/fasta_files" # FASTA directory
)
```
An example `build_parameters.json` is given below:
```
{
    "fasta_digest_params":
    {
        "min_length": 7,
        "max_length": 30,
        "min_charge": 2,
        "max_charge": 4,
        "cleavage_regex": "[KR][^_|$]",
        "missed_cleavages": 1,
        "max_var_mods": 1,
        "add_decoys": true,
        "entrapment_r": 0
    },
    "nce_params":
    {
        "nce": 25.0,
        "default_charge": 2,
        "dynamic_nce": true
    },
    "library_params":
    {
    "rt_bin_tol": 1.0,
    "frag_bin_tol_ppm": 10.0,
    "rank_to_score": [8, 4, 4, 2, 2, 1, 1],
    "y_start_index": 4,
    "b_start_index": 3,
    "y_start": 3,
    "b_start": 2,
    "include_p_index": false,
    "include_p": false,
    "auto_detect_frag_bounds": true,
    "calibration_raw_file": "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/OlsenMixedSpeciesAstral200ng/20230324_OLEP08_200ng_30min_E5H50Y45_180K_2Th3p5ms_01.arrow",
    "frag_mz_min":150.0, 
    "frag_mz_max":2020.0,
    "prec_mz_min":390.0,
    "prec_mz_max":1010.0,
    "max_frag_charge": 3,
    "max_frag_rank": 255,
    "min_frag_intensity": 0.00,
    "include_isotope": false,
    "include_internal": false,
    "include_immonium": false,
    "include_neutral_diff": true,
    "instrument_type": "NONE",
    "prediction_model": "altimeter"
    },
    "variable_mods": 
    {
            "pattern": ["M"],
            "mass": [15.99491],
            "name": ["Unimod:35"]
    },
    "fixed_mods": 
    {
        "pattern": ["C"],
        "mass": [57.021464],
        "name": ["Unimod:4"]
    },
    "isotope_mod_groups":
    [
    ]
    ,
    "max_koina_requests":24,
    "max_koina_batch": 1000,
    "match_lib_build_batch": 100000,
    "fasta_paths":
        [
         "/path/to/fasta/files/keap1.fasta"
        ],
    "fasta_names": ["HUMAN"],
    "out_dir": "/path/to/output",
    "lib_name": "/path/to/output/keap1_library",
    "new_lib_name":  "/path/to/output/keap1_library",
    "out_name": "keap1_library.tsv",
    "predict_fragments": true
}   
```
See [Parameter Configuration](@ref) for details on each parameter. 


### Build Spectral Library
Once the `json` configuration file is finished, `BuildSpecLib` generates a spectral library in the `pion` format. `BuildSpecLib` does this in three steps. 
First, `BuildSpecLib` in silico digests the all `.fasta` or `.fasta.gz` files in the directory specified by the configuration file according to a user-defined 
enzymatic cleavage rule. Next, `BuildSpecLib`. 
!!! note "'note'"
     `BuildSpecLib` does not generate spectral libraries locally. It uses the Koina servers and will fail without an internet connection. 

```@julia
# Build the library
BuildSpecLib(params)
```

## Running a Search

```julia
# Generate search parameters
search_params = GetSearchParams(
    "/path/to/human_library.pion",  # Library
    "/path/to/ms_data",            # MS data directory
    "/path/to/results"             # Results directory
)

# Run the search
SearchDIA(search_params)
```

## Example Dataset

Here's a complete example using the provided test data:

```julia
using Pioneer

# Build spectral library
build_params = GetBuildLibParams(
    joinpath(pwd(), "libraries"),
    "ecoli_test_lib",
    joinpath(pwd(), "test_data", "fasta")
)
BuildSpecLib(build_params)

# Run DIA search
search_params = GetSearchParams(
    joinpath(pwd(), "libraries", "ecoli_test_lib.pion"),
    joinpath(pwd(), "test_data", "ms_data"),
    joinpath(pwd(), "results")
)
SearchDIA(search_params)
```

## Understanding Results

The search produces several output files:

1. **Identifications**: Peptide/protein identifications with scores
2. **Quantification**: Abundance measurements
3. **Quality Control**: Search performance metrics
4. **Visualizations**: Diagnostic plots

## Next Steps

- See the [Parameter Configuration](@ref "Parameter Configuration")
- Review [Search Modes](@ref "Search Modes")
- Check the [Output Files](@ref "Output Files") documentation
- Learn about [Performance Tuning](@ref "Performance Tuning")