<img src="https://github.com/nwamsley1/Pioneer.jl/blob/main/figures/PIONEER_LOGO.jpg" align="right" width="150px"/>
<h1>Pioneer: Fast and Open-Source Analysis of Data-Independent Aquisition Proteomics Experiments

[![Build Status](https://github.com/nwamsley1/Pioneer.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nwamsley1/Pioneer.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/nwamsley1/Pioneer.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nwamsley1/Pioneer.jl)
</h1>


##  Development Aims
  Pioneer is a cross-platform and open-source tool fully implemented 
in Julia that identifies and quantifies proteins and peptides from data independent acquisition (DIA) experiments. Given a 
spectral library of fragment ion intensities and retention time estimates on an arbitrary scale, Pioneer employs a spectrum-centric 
algorithm and heuristics to statistically infer the identification status and abundance of each library precursor in the data. We develop Pioneer with the following goals:

- **Open-Source:** Methods should be understood and open to scrutiny by users
- **Cross-Platform:** All steps of analysis, including vendor file conversion, should run on all major operating systems
- **High-Performance:** The sensitivity, FDR control, and quantitative precision and accuracy should be competitive with state-of-the-art commercial software packages
- **Scalability:** Should scale to very large experiments with hundreds to thousands of raw files
- **Fast:** Use of simple heuristics and carefully implemented, efficient algorithms should ensure that data can be analyzed many times faster than it is aquired for typical experiments

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
(@v1.10) pkg> activate .
(@v1.10) pkg> develop ./
(@v1.10) pkg> add ./
```

## Usage

Pioneer exports two methods `SearchDIA` and `BuildSpecLib`. 

1) `SearchDIA` performs the Pioneer search algorithm on mass specrometry raw data given a spectral library in the `.pion` format.
2) `BuildSpecLib` generates a `.pion` spectral library given a FASTA file of protein sequences using [Chronologer](https://github.com/searlelab/chronologer) and the [Koina project](https://koina.wilhelmlab.org/.).

Each method takes a single argument, that is a `.json` formatted parameter file. Example parameter files for both methods are included below in addition to a description of each paramter. 

## File Conversion
`Pioneer` requires Thermo .raw files be converted to an Apache .arrow format with a specific column specification. Use [PioneerConverter](https://github.com/nwamsley1/PioneerConverter) to convert .raw files. 

## SearchDIA
```
julia> using Pioneer
julia> SearchDIA("/Users/n.t.wamsley/Projects/Pioneer.jl/data/ecoli_test/ecoli_test_params.json")
```
# Parameters JSON
Pioneer accepts parameters in a .json file format. Paramter descriptions and a .json template are provided. In most casees, the defaults are recommended. 

|Name|Type|Description|
|---|---|---|
|`isotope_settings.err_bounds`|[Int, Int]|Precursor monoisotope may lie this many isotopes (m/z) outside the quadrupole isolation window  [before, after]|
|`isotope_settings.combine_traces`|Boolean|Whether to combine precursor isotope traces in quantification. Precursor chromatograms may be split between differing quadrupole isolation windows.|
|`scoring.q_value_threshold`|Float|Global q-value threshold for filtering results|
|`normalization.n_rt_bins`|Int|Number of retention time bins for quant normalization|
|`normalization.spline_n_knots`|Int|Number of knots in quant normalization spline|

### Parameter Tuning Parameters
|Name|Type|Description|
|---|---|---|
|`fragment_settings.min_count`|Int|Minimum number of matching fragments required|
|`fragment_settings.max_rank`|Int|Maximum rank of fragments to consider|
|`fragment_settings.tol_ppm`|Float|Fragment mass tolerance in parts per million|
|`fragment_settings.min_score`|Int|Minimum score threshold for fragment matches|
|`fragment_settings.min_spectral_contrast`|Float|Minimum spectral contrast score|
|`fragment_settings.min_log2_ratio`|Float|Minimum log2 ratio of matched intensities|
|`fragment_settings.min_top_n`|[Int, Int]|Minimum number of top N matches [requirement, denominator]|
|`fragment_settings.n_isotopes`|Int|Number of isotopes to consider in matching|
|`search_settings.sample_rate`|Float|Fraction of spectra to sample during parameter tuning|
|`search_settings.min_samples`|Int|Minimum number of samples required for tuning|
|`search_settings.max_presearch_iters`|Int|Maximum number of parameter tuning iterations|
|`search_settings.frag_err_quantile`|Float|Quantile for fragment error estimation|

### First Search Parameters
|Name|Type|Description|
|---|---|---|
|`fragment_settings.min_count`|Int|Minimum number of matching fragments|
|`fragment_settings.max_rank`|Int|Maximum fragment rank to consider|
|`fragment_settings.min_score`|Int|Minimum score for fragment matches|
|`fragment_settings.min_spectral_contrast`|Float|Minimum spectral contrast required|
|`fragment_settings.min_log2_ratio`|Float|Minimum log2 ratio of matched intensities|
|`fragment_settings.min_top_n`|[Int, Int]|Minimum top N matches [requirement, denominator]|
|`fragment_settings.n_isotopes`|Int|Number of isotopes to consider|
|`scoring_settings.n_train_rounds`|Int|Number of training rounds for scoring model|
|`scoring_settings.max_iterations`|Int|Maximum iterations for scoring optimization|
|`scoring_settings.max_q_value`|Float|Maximum q-value threshold|
|`scoring_settings.max_precursors`|Int|Maximum number of precursors to consider|

### Quantification Search Parameters
|Name|Type|Description|
|---|---|---|
|`fragment_settings.min_count`|Int|Minimum fragment count for quantification|
|`fragment_settings.min_y_count`|Int|Minimum number of y-ions required|
|`fragment_settings.max_rank`|Int|Maximum fragment rank|
|`fragment_settings.min_spectral_contrast`|Float|Minimum spectral contrast score|
|`fragment_settings.min_log2_ratio`|Float|Minimum log2 ratio of intensities|
|`fragment_settings.min_top_n`|[Int, Int]|Minimum top N matches [requirement, denominator]|
|`fragment_settings.n_isotopes`|Int|Number of isotopes for quantification|
|`chromatogram.smoothing_strength`|Float|Strength of chromatogram smoothing (Whittaker-Henderson)|
|`chromatogram.padding`|Int|Number of points to pad chromatograms|
|`chromatogram.max_apex_offset`|Int|Maximum allowed apex offset|

### Acquisition Parameters
|Name|Type|Description|
|---|---|---|
|`nce`|Int|Normalized collision energy|
|`quad_transmission.fit_from_data`|Boolean|Whether to fit quadrupole transmission from data|
|`quad_transmission.overhang`|Float|Quadrupole transmission curve overhang|
|`quad_transmission.smoothness`|Float|Smoothness parameter for transmission curve|

### RT Alignment Parameters
|Name|Type|Description|
|---|---|---|
|`n_bins`|Int|Number of retention time bins|
|`bandwidth`|Float|Kernel bandwidth for RT alignment|
|`sigma_tolerance`|Int|Number of standard deviations for tolerance|
|`min_probability`|Float|Minimum probability for alignment matches|

### Optimization Parameters
|Name|Type|Description|
|---|---|---|
|`deconvolution.lambda`|Float|Regularization parameter for deconvolution|
|`deconvolution.huber_delta`|Float|Delta parameter for Huber loss|
|`deconvolution.huber_exp`|Float|Exponent for Huber delta progression|
|`deconvolution.huber_iters`|Int|Number of Huber iterations|
|`deconvolution.newton_iters`|Int|Maximum Newton iterations|
|`deconvolution.newton_accuracy`|Float|Convergence threshold for Newton method|
|`deconvolution.max_diff`|Float|Maximum allowed difference in optimization|
|`machine_learning.max_samples`|Int|Maximum number of samples for ML training|
|`machine_learning.min_trace_prob`|Float|Minimum trace probability threshold|
|`machine_learning.spline_points`|Int|Number of points for probability spline|
|`machine_learning.interpolation_points`|Int|Number of interpolation points|

### Output Parameters
|Name|Type|Description|
|---|---|---|
|`write_csv`|Boolean|Whether to write results to CSV|
|`delete_temp`|Boolean|Whether to delete temporary files|
|`plots_per_page`|Int|Number of plots per page in reports|

### Path Parameters
|Name|Type|Description|
|---|---|---|
|`library`|String|Path to spectral library file|
|`ms_data`|String|Path to mass spectrometry data directory|
|`results`|String|Path to output results directory|

### Example SearchDIA Parameters .json
```
{
    "global": {
        "isotope_settings": {
            "err_bounds": [1, 0],
            "combine_traces": false
        },
        "scoring": {
            "q_value_threshold": 0.01
        },
        "normalization": {
            "n_rt_bins": 100,
            "spline_n_knots": 7
        }
    },
    "parameter_tuning": {
        "fragment_settings": {
            "min_count": 7,
            "max_rank": 25,
            "tol_ppm": 20.0,
            "min_score": 22,
            "min_spectral_contrast": 0.9,
            "min_log2_ratio": 1.5,
            "min_top_n": [3, 3],
            "n_isotopes": 1
        },
        "search_settings": {
            "sample_rate": 0.02,
            "min_samples": 3500,
            "max_presearch_iters": 10,
            "frag_err_quantile": 0.01
        }
    },
    "first_search": {
        "fragment_settings": {
            "min_count": 4,
            "max_rank": 25,
            "min_score": 15,
            "min_spectral_contrast": 0.5,
            "min_log2_ratio": 0.0,
            "min_top_n": [2, 3],
            "n_isotopes": 1
        },
        "scoring_settings": {
            "n_train_rounds": 2,
            "max_iterations": 20,
            "max_q_value": 0.01,
            "max_precursors": 200000
        }
    },
    "quant_search": {
        "fragment_settings": {
            "min_count": 3,
            "min_y_count": 2,
            "max_rank": 255,
            "min_spectral_contrast": 0.0,
            "min_log2_ratio": -1.7,
            "min_top_n": [2, 3],
            "n_isotopes": 2
        },
        "chromatogram": {
            "smoothing_strength": 1.0,
            "padding": 20,
            "max_apex_offset": 2
        }
    },
    "acquisition": {
        "nce": 25,
        "quad_transmission": {
            "fit_from_data": false,
            "overhang": 0.25,
            "smoothness": 5.0
        }
    },
    "rt_alignment": {
        "n_bins": 200,
        "bandwidth": 0.25,
        "sigma_tolerance": 4,
        "min_probability": 0.95
    },
    "optimization": {
        "deconvolution": {
            "lambda": 0.0,
            "huber_delta": 300,
            "huber_exp": 1.5,
            "huber_iters": 15,
            "newton_iters": 100,
            "newton_accuracy": 10,
            "max_diff": 0.01
        },
        "machine_learning": {
            "max_samples": 10000000,
            "min_trace_prob": 0.75,
            "spline_points": 500,
            "interpolation_points": 10
        }
    },
    "output": {
        "write_csv": true,
        "delete_temp": true,
        "plots_per_page": 12
    },
    "paths": {
        "library": "/abs/path/to/lib",
        "ms_data": "/abs/path/to/ms_data",
        "results": "/abs/path/to/results"
    }
}
```

## BuildSpecLib
```
julia> using Pioneer
julia> SearchDIA("/Users/n.t.wamsley/Projects/Pioneer.jl/data/ecoli_test/ecoli_test_params.json")
```

### FASTA Digest Parameters
Parameters controlling protein digestion and peptide generation:

| Parameter | Type | Description |
|-----------|------|-------------|
| `min_length` | integer | Minimum peptide length (7 amino acids) |
| `max_length` | integer | Maximum peptide length (30 amino acids) |
| `min_charge` | integer | Minimum charge state to consider (2+) |
| `max_charge` | integer | Maximum charge state to consider (4+) |
| `cleavage_regex` | string | Regular expression defining cleavage sites ("[KR][^_|$]" for trypsin) |
| `missed_cleavages` | integer | Maximum number of allowed missed cleavages (1) |
| `max_var_mods` | integer | Maximum number of variable modifications per peptide (1) |
| `add_decoys` | boolean | Whether to generate decoy sequences |
| `entrapment_r` | float | Ratio of entrapment sequences to add (0 = disabled) |

### NCE Parameters 
Parameters for collision energy settings:

| Parameter | Type | Description |
|-----------|------|-------------|
| `nce` | float | Base normalized collision energy value (25.0) |
| `default_charge` | integer | Default charge state for NCE calculations (2) |
| `dynamic_nce` | boolean | Whether to use charge-dependent NCE adjustments |

### Library Parameters
Settings for spectral library generation and processing:

| Parameter | Type | Description |
|-----------|------|-------------|
| `rt_bin_tol` | float | Retention time binning tolerance in minutes (1.0) |
| `frag_bin_tol_ppm` | float | Fragment mass tolerance in PPM (10.0) |
| `rank_to_score` | array | Intensity multipliers for ranked peaks [8,4,4,2,2,1,1] |
| `y_start_index` | integer | Starting index for y-ion series annotation (4) |
| `b_start_index` | integer | Starting index for b-ion series annotation (3) |
| `y_start` | integer | Minimum y-ion number to consider (3) |
| `b_start` | integer | Minimum b-ion number to consider (2) |
| `include_p_index` | boolean | Include proline-containing index fragments |
| `include_p` | boolean | Include proline-containing fragments |
| `auto_detect_frag_bounds` | boolean | Automatically detect fragment mass bounds |
| `calibration_raw_file` | string | Path to raw file for calibration |
| `frag_mz_min` | float | Minimum fragment m/z (150.0) |
| `frag_mz_max` | float | Maximum fragment m/z (2020.0) |
| `prec_mz_min` | float | Minimum precursor m/z (390.0) |
| `prec_mz_max` | float | Maximum precursor m/z (1010.0) |
| `max_frag_charge` | integer | Maximum fragment ion charge (3) |
| `max_frag_rank` | integer | Maximum fragment rank to consider (255) |
| `min_frag_intensity` | float | Minimum relative fragment intensity (0.00) |
| `include_isotope` | boolean | Include isotope peak annotations |
| `include_internal` | boolean | Include internal fragment annotations |
| `include_immonium` | boolean | Include immonium ion annotations |
| `include_neutral_diff` | boolean | Include neutral loss annotations |
| `instrument_type` | string | Instrument type for predictions ("NONE") |
| `prediction_model` | string | Model name for fragment predictions ("prosit_2020_hcd") |

### Variable Modifications
| Parameter | Type | Description |
|-----------|------|-------------|
| `pattern` | array | Amino acids to modify (["M"] = methionine) |
| `mass` | array | Modification masses ([15.99491] = oxidation) |
| `name` | array | Modification identifiers (["Unimod:35"] = oxidation) |

### Fixed Modifications
| Parameter | Type | Description |
|-----------|------|-------------|
| `pattern` | array | Amino acids to modify (["C"] = cysteine) |
| `mass` | array | Modification masses ([57.021464] = carbamidomethyl) |
| `name` | array | Modification identifiers (["Unimod:4"] = carbamidomethyl) |

### Processing Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `max_koina_requests` | integer | Maximum concurrent Prosit API requests (24) |
| `max_koina_batch` | integer | Maximum batch size for API requests (1000) |
| `match_lib_build_batch` | integer | Batch size for library building (100000) |

### Input/Output Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `fasta_paths` | array | List of FASTA file paths to process |
| `fasta_names` | array | Names to assign to each FASTA file |
| `out_dir` | string | Output directory path |
| `lib_name` | string | Base name for library files |
| `new_lib_name` | string | Name for updated library files |
| `out_name` | string | Output filename |
| `predict_fragments` | boolean | Whether to predict fragment intensities |

### Example SearchDIA Parameters .json
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
    "include_neutral_diff": false,
    "instrument_type": "NONE",
    "prediction_model": "prosit_2020_hcd"
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
         "/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/PROTEOMES/UP000005640_9606_human.fasta.gz",
         "/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/PROTEOMES/UP000002311_559292_Saccharomyces_cerevisiae.fasta.gz",
         "/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/PROTEOMES/UP000000625_83333_Escherichia_coli.fasta.gz"
        ],
    "fasta_names": ["HUMAN","YEAST","ECOLI"],
    "out_dir": "/Users/n.t.wamsley/RIS_temp/koina_testing",
    "lib_name": "/Users/n.t.wamsley/RIS_temp/koina_testing/prosit_2020_hcd_threeproteome_122024",
    "new_lib_name": "/Users/n.t.wamsley/RIS_temp/koina_testing/prosit_2020_hcd_threeproteome_122024",
    "out_name": "prosit_2020_hcd_threeproteome_122024.tsv",
    "predict_fragments": true
}   
```

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
- f_index_fragment_bins.arrow
- f_index_fragments.arrow
- f_index_rt_bins.arrow
These make up the fragment index for the initial search. 
- precursors_table.arrow
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
