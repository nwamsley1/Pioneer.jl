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
Pioneer exports three `SearchDIA` methods, each of which takes a single argument, that is a file path to a .json parameters files (examples included). To access these methods
import the Pioneer package from within the Julia REPL. 
```
julia> using Pioneer
julia> SearchDIA("/Users/n.t.wamsley/Projects/Pioneer.jl/data/ecoli_test/ecoli_test_params.json")
[ Info: Loading Parameters...
[ Info: Loading Spectral Library...
Loading spectral libraries into main memory...
[ Info: Initializing Search Context...
[ Info: Executing Parameter Tuning...
100.0%┣████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████┫ 3/3 [00:03<00:00, 1s/it]
[ Info: Merging QC plots...
[ Info: QC plot merging complete
[ Info: Executing NCE Tuning...
100.0%┣██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████┫ 3/3 [00:00<00:00, 372it/s]
[ Info: Executing Quadrupole Tuning...
100.0%┣██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████┫ 3/3 [00:00<00:00, 712it/s]
[ Info: Executing First Pass Search...
100.0%┣████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████┫ 3/3 [00:02<00:00, 1it/s]
[ Info: Summarizing first pass search results...
[ Info: Mapping library to empirical retention times...
[ Info: Finding best precursors across runs...
[ Info: Calculating iRT errors...
[ Info: Creating RT indices...
[ Info: Search results summarization complete
[ Info: Executing Huber Tuning...
100.0%┣████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████┫ 3/3 [00:00<00:00, 5it/s]
[ Info: Executing Second Pass Search...
100.0%┣████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████┫ 3/3 [00:01<00:00, 2it/s]
[ Info: Executing Scoring...
100.0%┣██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████┫ 3/3 [00:00<00:00, 769it/s]
[ Info: Training XGBoost models...
100.0%┣████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████┫ 6/6 [00:01<00:00, 7it/s]
100.0%┣███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████┫ 6/6 [00:00<00:00, 26it/s]
[ Info: Finding best traces...
[ Info: Processing quantification results...
[ Info: Merging PSM scores...
[ Info: Calculating error probabilities...
┌ Warning: Less than 20 bins to estimate PEP. PEP results suspect...
└ @ Pioneer ~/Projects/Pioneer.jl/src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils.jl:432
┌ Warning: Failed to estimate PEP spline
└ @ Pioneer ~/Projects/Pioneer.jl/src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils.jl:436
[ Info: Filtering passing PSMs...
[ Info: Scoring protein groups...
[ Info: Executing Chromatogram Integration...
100.0%┣████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████┫ 3/3 [00:02<00:00, 1s/it]
[ Info: Executing MaxLFQ...
100.0%┣██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████┫ 3/3 [00:00<00:00, 715it/s]
[ Info: Performing intensity normalization...
[ Info: Merging quantification tables...
[ Info: Writing precursor results...
[ Info: Performing MaxLFQ...
[ Info: Writing protein group results...
[ Info: Creating QC plots...
[ Info: Generating final QC plots

==========================================================================================
DIA Search Performance Report
==========================================================================================

Detailed Step Analysis:
------------------------------------------------------------------------------------------
Step                           Time (s)     Memory (GB)  GC Time (s)  GC %        
------------------------------------------------------------------------------------------
Parameter Loading              0.00         0.00         0.00         0.0         
NCE Tuning                     0.06         0.00         0.00         0.0         
Spectral Library Loading       0.24         0.02         0.00         0.0         
Quadrupole Tuning              0.46         0.02         0.09         20.5        
Second Pass Search             1.27         0.26         0.04         3.4         
Huber Tuning                   1.28         0.15         0.03         2.4         
Scoring                        1.30         0.35         0.22         16.6        
MaxLFQ                         1.43         0.24         0.20         14.3        
Search Context Initialization  1.88         0.70         1.73         91.8        
First Pass Search              1.93         0.42         0.13         6.9         
Chromatogram Integration       2.43         2.84         0.87         35.7        
Parameter Tuning               3.21         1.25         0.17         5.3         
------------------------------------------------------------------------------------------
TOTAL                          15.50        6.24         3.49         22.5        

Memory Usage Summary:
------------------------------------------------------------------------------------------
Total Memory Allocated: 6.24 GB
Total Available Memory: 64.0 GB

Runtime Summary:
------------------------------------------------------------------------------------------
Total Runtime: 0.26 minutes
Average Runtime per Step: 1.29 seconds
Average Runtime per Raw File: 5.17 seconds

==========================================================================================
```


### File Conversion
`Pioneer` requires Thermo .raw files be converted to an Apache .arrow format with a specific column specification. Use [PioneerConverter](https://github.com/nwamsley1/PioneerConverter) to convert .raw files. 

### Parameters JSON
Pioneer accepts parameters in a .json file format. Paramter descriptions and a .json template are provided. In most casees, the defaults are recommended. 
|Name|Type|Description|
|---|---|---|
|isotope_settings.err_bounds|[Int, Int]|Precursor monoisotope may lie this many isotopes (m/z) outside the quadrupole isolation window  [before, after]|
|isotope_settings.combine_traces|Boolean|Whether to combine precursor isotope traces in quantification. Precursor chromatograms may be split between differing quadrupole isolation windows.|
|scoring.q_value_threshold|Float|Global q-value threshold for filtering results|
|normalization.n_rt_bins|Int|Number of retention time bins for quant normalization|
|normalization.spline_n_knots|Int|Number of knots in quant normalization spline|

### Parameter Tuning Parameters
|Name|Type|Description|
|---|---|---|
|fragment_settings.min_count|Int|Minimum number of matching fragments required|
|fragment_settings.max_rank|Int|Maximum rank of fragments to consider|
|fragment_settings.tol_ppm|Float|Fragment mass tolerance in parts per million|
|fragment_settings.min_score|Int|Minimum score threshold for fragment matches|
|fragment_settings.min_spectral_contrast|Float|Minimum spectral contrast score|
|fragment_settings.min_log2_ratio|Float|Minimum log2 ratio of matched intensities|
|fragment_settings.min_top_n|[Int, Int]|Minimum number of top N matches [requirement, denominator]|
|fragment_settings.n_isotopes|Int|Number of isotopes to consider in matching|
|search_settings.sample_rate|Float|Fraction of spectra to sample during parameter tuning|
|search_settings.min_samples|Int|Minimum number of samples required for tuning|
|search_settings.max_presearch_iters|Int|Maximum number of parameter tuning iterations|
|search_settings.frag_err_quantile|Float|Quantile for fragment error estimation|

### First Search Parameters
|Name|Type|Description|
|---|---|---|
|fragment_settings.min_count|Int|Minimum number of matching fragments|
|fragment_settings.max_rank|Int|Maximum fragment rank to consider|
|fragment_settings.min_score|Int|Minimum score for fragment matches|
|fragment_settings.min_spectral_contrast|Float|Minimum spectral contrast required|
|fragment_settings.min_log2_ratio|Float|Minimum log2 ratio of matched intensities|
|fragment_settings.min_top_n|[Int, Int]|Minimum top N matches [requirement, denominator]|
|fragment_settings.n_isotopes|Int|Number of isotopes to consider|
|scoring_settings.n_train_rounds|Int|Number of training rounds for scoring model|
|scoring_settings.max_iterations|Int|Maximum iterations for scoring optimization|
|scoring_settings.max_q_value|Float|Maximum q-value threshold|
|scoring_settings.max_precursors|Int|Maximum number of precursors to consider|

### Quantification Search Parameters
|Name|Type|Description|
|---|---|---|
|fragment_settings.min_count|Int|Minimum fragment count for quantification|
|fragment_settings.min_y_count|Int|Minimum number of y-ions required|
|fragment_settings.max_rank|Int|Maximum fragment rank|
|fragment_settings.min_spectral_contrast|Float|Minimum spectral contrast score|
|fragment_settings.min_log2_ratio|Float|Minimum log2 ratio of intensities|
|fragment_settings.min_top_n|[Int, Int]|Minimum top N matches [requirement, denominator]|
|fragment_settings.n_isotopes|Int|Number of isotopes for quantification|
|chromatogram.smoothing_strength|Float|Strength of chromatogram smoothing (Whittaker-Henderson)|
|chromatogram.padding|Int|Number of points to pad chromatograms|
|chromatogram.max_apex_offset|Int|Maximum allowed apex offset|

### Acquisition Parameters
|Name|Type|Description|
|---|---|---|
|nce|Int|Normalized collision energy|
|quad_transmission.fit_from_data|Boolean|Whether to fit quadrupole transmission from data|
|quad_transmission.overhang|Float|Quadrupole transmission curve overhang|
|quad_transmission.smoothness|Float|Smoothness parameter for transmission curve|

### RT Alignment Parameters
|Name|Type|Description|
|---|---|---|
|n_bins|Int|Number of retention time bins|
|bandwidth|Float|Kernel bandwidth for RT alignment|
|sigma_tolerance|Int|Number of standard deviations for tolerance|
|min_probability|Float|Minimum probability for alignment matches|

### Optimization Parameters
|Name|Type|Description|
|---|---|---|
|deconvolution.lambda|Float|Regularization parameter for deconvolution|
|deconvolution.huber_delta|Float|Delta parameter for Huber loss|
|deconvolution.huber_exp|Float|Exponent for Huber delta progression|
|deconvolution.huber_iters|Int|Number of Huber iterations|
|deconvolution.newton_iters|Int|Maximum Newton iterations|
|deconvolution.newton_accuracy|Float|Convergence threshold for Newton method|
|deconvolution.max_diff|Float|Maximum allowed difference in optimization|
|machine_learning.max_samples|Int|Maximum number of samples for ML training|
|machine_learning.min_trace_prob|Float|Minimum trace probability threshold|
|machine_learning.spline_points|Int|Number of points for probability spline|
|machine_learning.interpolation_points|Int|Number of interpolation points|

### Output Parameters
|Name|Type|Description|
|---|---|---|
|write_csv|Boolean|Whether to write results to CSV|
|delete_temp|Boolean|Whether to delete temporary files|
|plots_per_page|Int|Number of plots per page in reports|

### Path Parameters
|Name|Type|Description|
|---|---|---|
|library|String|Path to spectral library file|
|ms_data|String|Path to mass spectrometry data directory|
|results|String|Path to output results directory|

Example json
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
