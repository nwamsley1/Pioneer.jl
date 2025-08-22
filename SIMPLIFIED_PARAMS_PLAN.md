# Plan: Implement Simplified and Full Parameter Templates for SearchDIA

## Overview
Implement a system that allows users to provide a simplified parameters.json file with only essential settings, while automatically applying sensible defaults for advanced parameters. The `GetSearchParams` function will support generating either simplified or full parameter templates.

## Current State Analysis

### Existing Implementation
1. **GetSearchParams** (`src/Routines/GenerateParams.jl`):
   - Currently reads from `defaultSearchParams.json` template
   - Updates only the paths section
   - No concept of simplified vs full parameters

2. **Parameter Parsing** (`src/Routines/SearchDIA/ParseInputs/parseParams.jl`):
   - Converts JSON to `PioneerParameters` struct with NamedTuples
   - No default value merging
   - Expects all parameters to be present

3. **Parameter Validation** (`src/Routines/SearchDIA/ParseInputs/paramsChecks.jl`):
   - Strict validation requiring all parameters
   - Throws errors for missing required fields
   - Some handling for optional parameters (e.g., logging section)

## Simplified vs Full Parameters

### Simplified Parameters (Essential Only)
```json
{
    "logging": {
        "debug_console_level": 0
    },
    "global": {
        "isotope_settings": {
            "combine_traces": true,
            "partial_capture": true,
            "min_fraction_transmitted": 0.25
        },
        "scoring": {
            "q_value_threshold": 0.01
        },
        "match_between_runs": true,
        "ms1_quant": false,
        "ms1_scoring": true
    },
    "parameter_tuning": {
        "fragment_settings": {
            "min_count": 7,
            "min_score": [22,17,15]
        },
        "search_settings": {   
            "topn_peaks": 200
        },
        "iteration_settings": {
            "init_mass_tol_ppm": 10.0,
            "mass_tolerance_scale_factor": 1.5,
            "iterations_per_phase": 3
        }
    },
    "first_search": {
        "fragment_settings": {
            "min_score": 15,
            "n_isotopes": 1
        }
    },
    "quant_search": {
        "fragment_settings": {
            "max_rank": 255,
            "n_isotopes": 2
        }
    },
    "acquisition": {
        "nce": 25,
        "quad_transmission": {
            "fit_from_data": false
        }
    },
    "proteinInference": {
        "min_peptides": 1
    },
    "maxLFQ": {
        "run_to_run_normalization": false
    },
    "output": {
        "write_csv": true,
        "write_decoys": false,
        "delete_temp": true,
        "plots_per_page": 12
    },
    "paths": {
        "ms_data": "/path/to/ms/data/folder",
        "library": "/path/to/ms/data/library.pion",
        "results": "/path/to/results"
    }
}
```

### Additional Parameters in Full Version (with Defaults)

#### Global Settings
- **global.isotope_settings**: 
  - `err_bounds_first_pass: [1,0]`
  - `err_bounds_quant_search: [3,0]`
- **global.normalization**: 
  - `n_rt_bins: 100`
  - `spline_n_knots: 7`
- **global.huber_override**: 
  - `override_huber_delta_fit: false`
  - `huber_delta: 1055`

#### Parameter Tuning
- **parameter_tuning.fragment_settings**: 
  - `max_rank: 25`
  - `min_spectral_contrast: 0.5`
  - `relative_improvement_threshold: 1.5`
  - `min_log2_ratio: 1.5`
  - `min_top_n: [3,3]`
  - `n_isotopes: 1`
  - `intensity_filter_quantile: 0.50`
- **parameter_tuning.search_settings**: 
  - `initial_scan_count: 10000`
  - `max_parameter_tuning_scans: 40000`
  - `max_q_value: 0.01`
  - `min_samples: 1000`
  - `max_frags_for_mass_err_estimation: 5`
  - `min_quad_tuning_psms: 5000`
  - `min_quad_tuning_fragments: 3`
  - `max_presearch_iters: 10`
  - `frag_err_quantile: 0.0025`
- **parameter_tuning.iteration_settings**:
  - `scan_scale_factor: 4`

#### First Search
- **first_search.fragment_settings**:
  - `min_count: 4`
  - `max_rank: 25`
  - `min_spectral_contrast: 0.5`
  - `relative_improvement_threshold: 1.25`
  - `min_log2_ratio: 0.0`
  - `min_top_n: [2,3]`
- **first_search.scoring_settings**:
  - `n_train_rounds: 2`
  - `max_iterations: 20`
  - `max_q_value_probit_rescore: 0.05`
  - `max_PEP: 0.9`
- **first_search.irt_mapping**:
  - `max_prob_to_impute_irt: 0.75`
  - `fwhm_nstd: 4`
  - `irt_nstd: 4`

#### Quant Search
- **quant_search.fragment_settings**:
  - `min_count: 3`
  - `min_y_count: 2`
  - `min_spectral_contrast: 0.0`
  - `min_log2_ratio: -1.7`
  - `min_top_n: [2,3]`
- **quant_search.chromatogram**:
  - `smoothing_strength: 1e-6`
  - `padding: 0`
  - `max_apex_offset: 2`

#### Acquisition
- **acquisition.quad_transmission**:
  - `overhang: 0.25`
  - `smoothness: 5.0`

#### RT Alignment (entire section with defaults)
- **rt_alignment**:
  - `n_bins: 200`
  - `bandwidth: 0.25`
  - `sigma_tolerance: 4`
  - `min_probability: 0.95`

#### Optimization (entire section with defaults)
- **optimization.deconvolution**:
  - `lambda: 0.0`
  - `reg_type: "none"`
  - `huber_delta: 300`
  - `huber_exp: 1.5`
  - `huber_iters: 15`
  - `newton_iters: 50`
  - `bisection_iters: 100`
  - `outer_iters: 1000`
  - `newton_accuracy: 10`
  - `max_diff: 0.01`
- **optimization.machine_learning**:
  - `max_psms_in_memory: 50000000`
  - `min_trace_prob: 0.75`
  - `max_q_value_xgboost_mbr_rescore: 0.20`
  - `min_PEP_neg_threshold_xgboost_rescore: 0.90`
  - `spline_points: 500`
  - `interpolation_points: 10`

### Complete Full Parameters JSON Specification
```json
{
    "logging": {
        "debug_console_level": 0
    },
    "global": {
        "isotope_settings": {
            "err_bounds_first_pass": [1, 0],
            "err_bounds_quant_search": [3, 0],
            "combine_traces": true,
            "partial_capture": true,
            "min_fraction_transmitted": 0.25
        },
        "scoring": {
            "q_value_threshold": 0.01
        },
        "normalization": {
            "n_rt_bins": 100,
            "spline_n_knots": 7
        },
        "huber_override": {
            "override_huber_delta_fit": false,
            "huber_delta": 1055
        },
        "match_between_runs": true,
        "ms1_quant": false,
        "ms1_scoring": true
    },
    "parameter_tuning": {
        "fragment_settings": {
            "min_count": 7,
            "max_rank": 25,
            "min_score": [22, 17, 15],
            "min_spectral_contrast": 0.5,
            "relative_improvement_threshold": 1.5,
            "min_log2_ratio": 1.5,
            "min_top_n": [3, 3],
            "n_isotopes": 1,
            "intensity_filter_quantile": 0.50
        },
        "search_settings": {
            "initial_scan_count": 10000,
            "max_parameter_tuning_scans": 40000,
            "max_q_value": 0.01,
            "min_samples": 1000,
            "topn_peaks": 200,
            "max_frags_for_mass_err_estimation": 5,
            "min_quad_tuning_psms": 5000,
            "min_quad_tuning_fragments": 3,
            "max_presearch_iters": 10,
            "frag_err_quantile": 0.0025
        },
        "iteration_settings": {
            "init_mass_tol_ppm": 10.0,
            "mass_tolerance_scale_factor": 1.5,
            "iterations_per_phase": 3,
            "scan_scale_factor": 4
        }
    },
    "first_search": {
        "fragment_settings": {
            "min_count": 4,
            "max_rank": 25,
            "min_score": 12,
            "min_spectral_contrast": 0.5,
            "relative_improvement_threshold": 1.25,
            "min_log2_ratio": 0.0,
            "min_top_n": [2, 3],
            "n_isotopes": 1
        },
        "scoring_settings": {
            "n_train_rounds": 2,
            "max_iterations": 20,
            "max_q_value_probit_rescore": 0.05,
            "max_PEP": 0.9
        },
        "irt_mapping": {
            "max_prob_to_impute_irt": 0.75,
            "fwhm_nstd": 4,
            "irt_nstd": 4
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
            "smoothing_strength": 1e-6,
            "padding": 0,
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
            "reg_type": "none",
            "huber_delta": 300,
            "huber_exp": 1.5,
            "huber_iters": 15,
            "newton_iters": 50,
            "bisection_iters": 100,
            "outer_iters": 1000,
            "newton_accuracy": 10,
            "max_diff": 0.01
        },
        "machine_learning": {
            "max_psms_in_memory": 50000000,
            "min_trace_prob": 0.75,
            "max_q_value_xgboost_mbr_rescore": 0.20,
            "min_PEP_neg_threshold_xgboost_rescore": 0.90,
            "spline_points": 500,
            "interpolation_points": 10
        }
    },
    "proteinInference": {
        "min_peptides": 1
    },
    "maxLFQ": {
        "run_to_run_normalization": false
    },
    "output": {
        "write_csv": true,
        "write_decoys": false,
        "delete_temp": true,
        "plots_per_page": 12
    },
    "paths": {
        "ms_data": "/path/to/ms/data/folder",
        "library": "/path/to/ms/data/library.pion",
        "results": "/path/to/results"
    }
}
```

## Implementation Plan

### Phase 1: Create Infrastructure for Defaults

#### 1.1 Create Default Parameters Module
**File**: `src/Routines/SearchDIA/ParseInputs/paramDefaults.jl`

```julia
module ParamDefaults

export get_default_parameters, merge_with_defaults

"""
    get_default_parameters()

Returns the complete default parameter structure for SearchDIA.
All parameters have sensible defaults that work for most experiments.
Users only need to override the values they want to change.
"""
function get_default_parameters()
    return Dict(
        "logging" => Dict(
            "debug_console_level" => 0
        ),
        "global" => Dict(
            "isotope_settings" => Dict(
                "err_bounds_first_pass" => [1, 0],
                "err_bounds_quant_search" => [3, 0],
                "combine_traces" => true,
                "partial_capture" => true,
                "min_fraction_transmitted" => 0.25
            ),
            "scoring" => Dict(
                "q_value_threshold" => 0.01
            ),
            "normalization" => Dict(
                "n_rt_bins" => 100,
                "spline_n_knots" => 7
            ),
            "huber_override" => Dict(
                "override_huber_delta_fit" => false,
                "huber_delta" => 1055
            ),
            "match_between_runs" => true,
            "ms1_quant" => false,
            "ms1_scoring" => true
        ),
        # ... continue with all sections matching the full JSON spec above
    )
end

"""
    merge_with_defaults(user_params::Dict, defaults::Dict)

Recursively merges user parameters over default parameters.
User values override defaults at any nesting level.
Missing sections or parameters are filled from defaults.
"""
function merge_with_defaults(user_params::Dict, defaults::Dict)
    result = copy(defaults)
    
    function recursive_merge!(target::Dict, source::Dict)
        for (key, value) in source
            if haskey(target, key) && isa(target[key], Dict) && isa(value, Dict)
                # Both are dicts, merge recursively
                recursive_merge!(target[key], value)
            else
                # Override with user value
                target[key] = value
            end
        end
    end
    
    recursive_merge!(result, user_params)
    return result
end

end # module
```

#### 1.2 Create Template JSON Files
- **`assets/example_config/defaultSearchParamsSimplified.json`**: Simplified template
- **Keep `assets/example_config/defaultSearchParams.json`**: Full template (update with all defaults)

### Phase 2: Update Core Functions

#### 2.1 Modify Parameter Parsing
**File**: `src/Routines/SearchDIA/ParseInputs/parseParams.jl`

```julia
function parse_pioneer_parameters(json_path::String; apply_defaults::Bool = true)
    # Read user JSON
    user_params = JSON.parsefile(json_path)
    
    if apply_defaults
        # Get defaults
        defaults = ParamDefaults.get_default_parameters()
        # Merge user params over defaults
        params = ParamDefaults.merge_with_defaults(user_params, defaults)
    else
        params = user_params
    end
    
    # Convert to NamedTuples as before
    # ... existing conversion code
end
```

#### 2.2 Update GetSearchParams
**File**: `src/Routines/GenerateParams.jl`

```julia
function GetSearchParams(lib_path::String, ms_data_path::String, results_path::String; 
                        params_path::Union{String, Missing} = missing,
                        simplified::Bool = true)
    # ... path handling code
    
    # Choose template based on simplified flag
    template_name = simplified ? "defaultSearchParamsSimplified.json" : "defaultSearchParams.json"
    config_text = read(asset_path("example_config", template_name), String)
    
    # ... rest of function
end
```

### Phase 3: Update Validation

#### 3.1 Enhanced Parameter Validation
**File**: `src/Routines/SearchDIA/ParseInputs/paramsChecks.jl`

```julia
function checkParams(json_path::String)
    params = JSON.parsefile(json_path)
    
    # Apply defaults before validation
    defaults = ParamDefaults.get_default_parameters()
    params = ParamDefaults.merge_with_defaults(params, defaults)
    
    # Enhanced validation with clear comments
    # Check type correctness for user-supplied values
    # even after defaults are applied
    
    # ====== GLOBAL SETTINGS VALIDATION ======
    # Validate isotope settings
    if haskey(params["global"]["isotope_settings"], "min_fraction_transmitted")
        val = params["global"]["isotope_settings"]["min_fraction_transmitted"]
        if !(0.0 <= val <= 1.0)
            error("min_fraction_transmitted must be between 0 and 1, got $val")
        end
    end
    
    # Validate scoring threshold
    if haskey(params["global"]["scoring"], "q_value_threshold")
        val = params["global"]["scoring"]["q_value_threshold"]
        if !(0.0 < val <= 1.0)
            error("q_value_threshold must be between 0 and 1, got $val")
        end
    end
    
    # ====== PARAMETER TUNING VALIDATION ======
    # Validate min_score can be integer or array of integers
    if haskey(params["parameter_tuning"]["fragment_settings"], "min_score")
        min_score = params["parameter_tuning"]["fragment_settings"]["min_score"]
        if min_score isa Vector
            for (i, val) in enumerate(min_score)
                if !(val isa Integer) || val < 0
                    error("min_score[$i] must be non-negative integer, got $val")
                end
            end
        elseif !(min_score isa Integer) || min_score < 0
            error("min_score must be non-negative integer or array, got $min_score")
        end
    end
    
    # Validate iteration settings
    iter_settings = params["parameter_tuning"]["iteration_settings"]
    if iter_settings["mass_tolerance_scale_factor"] <= 1.0
        error("mass_tolerance_scale_factor must be > 1.0, got $(iter_settings["mass_tolerance_scale_factor"])")
    end
    if iter_settings["init_mass_tol_ppm"] <= 0.0
        error("init_mass_tol_ppm must be positive, got $(iter_settings["init_mass_tol_ppm"])")
    end
    if iter_settings["iterations_per_phase"] <= 0
        error("iterations_per_phase must be positive integer, got $(iter_settings["iterations_per_phase"])")
    end
    
    # Validate intensity_filter_quantile if present
    if haskey(params["parameter_tuning"]["fragment_settings"], "intensity_filter_quantile")
        val = params["parameter_tuning"]["fragment_settings"]["intensity_filter_quantile"]
        if !(0.0 <= val < 1.0)
            error("intensity_filter_quantile must be in [0, 1), got $val")
        end
    end
    
    # Continue with existing validation for required fields...
    # All validation messages should be clear about what's wrong
end
```

### Phase 4: CLI and User Interface

#### 4.1 Update CLI Entry Point
**File**: `src/Routines/GenerateParams.jl`

```julia
function main_GetSearchParams(argv=ARGS)::Cint
    settings = ArgParseSettings(; autofix_names = true)
    @add_arg_table! settings begin
        # ... existing arguments
        "--full"
            help = "Generate full parameter template with all advanced options"
            action = :store_true
        "--simplified"
            help = "Generate simplified parameter template (default)"
            action = :store_true
    end
    
    # Determine template type
    simplified = !parsed_args[:full]
    
    # Call GetSearchParams with simplified flag
    GetSearchParams(...; simplified=simplified)
end
```

### Phase 5: Update Integration Tests

#### 5.1 Update E. coli Test Parameters
**File**: `data/ecoli_test/ecoli_test_params.json`
- Ensure the test parameters are compatible with new defaults system
- Verify all required fields are present
- Test should continue to pass with the new parameter handling

#### 5.2 Verify RunTests.jl Compatibility
**File**: `test/runtests.jl`
- Ensure integration tests work with new parameter system
- No changes needed if tests use complete parameter files
- The default merging should be transparent to existing tests

### Phase 6: Documentation

#### 6.1 Update User Documentation
- Update CLAUDE.md with simplified/full parameter explanation
- Add examples to docs/
- Update installation guide with new GetSearchParams usage

#### 6.2 Migration Guide
- Document changes for existing users
- Provide conversion script if needed

## Risk Mitigation

### Backwards Compatibility
- All existing parameter files must continue to work
- Add `apply_defaults` flag to optionally disable new behavior
- Extensive testing with existing parameter files

### Default Value Management
- Single source of truth in `paramDefaults.jl`
- Version control for default changes
- Clear documentation of default values

### Validation Complexity
- Maintain strict validation for required parameters
- Clear error messages when required params missing
- Warnings for deprecated parameters

## Benefits

1. **User Experience**:
   - New users can start with minimal configuration
   - Reduced cognitive load for basic usage
   - Clear separation of essential vs advanced settings

2. **Maintainability**:
   - Single source of truth for defaults
   - Easier to update default values
   - Cleaner parameter files in examples

3. **Flexibility**:
   - Power users retain full control
   - Gradual learning curve from simple to advanced
   - Backwards compatible with existing workflows

## Implementation Order

1. Create `paramDefaults.jl` with default structure
2. Implement deep merge function
3. Update `parse_pioneer_parameters` to use defaults
4. Create simplified template JSON
5. Update `GetSearchParams` with simplified flag
6. Update validation to work with merged params
7. Add CLI support for --full flag
8. Write comprehensive tests
9. Update documentation

## Important Implementation Notes

### Parameter Value Validation
- **Type checking remains strict**: Even with defaults, user-supplied values must have correct types
- **Range validation**: Parameters with physical constraints (e.g., probabilities between 0-1) are validated
- **Current validation preserved**: All existing validation in `paramsChecks.jl` continues to work
- **Clear error messages**: Validation errors should clearly state what parameter failed and why

### Complete Parameter Validation Specification

#### Logging
- `debug_console_level`: Integer, [0, 3]

#### Global Settings
- `isotope_settings.err_bounds_first_pass`: Array of 2 Integers, each >= 0
- `isotope_settings.err_bounds_quant_search`: Array of 2 Integers, each >= 0
- `isotope_settings.combine_traces`: Boolean
- `isotope_settings.partial_capture`: Boolean
- `isotope_settings.min_fraction_transmitted`: Real, [0.0, 1.0]
- `scoring.q_value_threshold`: Real, (0.0, 1.0]
- `normalization.n_rt_bins`: Integer, > 0
- `normalization.spline_n_knots`: Integer, > 0
- `huber_override.override_huber_delta_fit`: Boolean
- `huber_override.huber_delta`: Real, > 0
- `match_between_runs`: Boolean
- `ms1_quant`: Boolean
- `ms1_scoring`: Boolean

#### Parameter Tuning
- `fragment_settings.min_count`: Integer, > 0
- `fragment_settings.max_rank`: Integer, [1, 255]
- `fragment_settings.min_score`: Integer > 0 OR Array of Integers each > 0
- `fragment_settings.min_spectral_contrast`: Real, [0.0, 1.0]
- `fragment_settings.relative_improvement_threshold`: Real, >= 1.0
- `fragment_settings.min_log2_ratio`: Real (can be negative)
- `fragment_settings.min_top_n`: Array of 2 Integers, each > 0
- `fragment_settings.n_isotopes`: Integer, [1, 3]
- `fragment_settings.intensity_filter_quantile`: Real, [0.0, 1.0)
- `search_settings.initial_scan_count`: Integer, > 0
- `search_settings.max_parameter_tuning_scans`: Integer, > initial_scan_count
- `search_settings.max_q_value`: Real, (0.0, 1.0]
- `search_settings.min_samples`: Integer, > 0
- `search_settings.topn_peaks`: Integer, > 0
- `search_settings.max_frags_for_mass_err_estimation`: Integer, > 0
- `search_settings.min_quad_tuning_psms`: Integer, > 0
- `search_settings.min_quad_tuning_fragments`: Integer, > 0
- `search_settings.max_presearch_iters`: Integer, > 0
- `search_settings.frag_err_quantile`: Real, [0.0, 1.0)
- `iteration_settings.init_mass_tol_ppm`: Real, > 0
- `iteration_settings.mass_tolerance_scale_factor`: Real, > 1.0
- `iteration_settings.iterations_per_phase`: Integer, > 0
- `iteration_settings.scan_scale_factor`: Real, > 1.0

#### First Search
- `fragment_settings.min_count`: Integer, > 0
- `fragment_settings.max_rank`: Integer, [1, 255]
- `fragment_settings.min_score`: Integer, > 0
- `fragment_settings.min_spectral_contrast`: Real, [0.0, 1.0]
- `fragment_settings.relative_improvement_threshold`: Real, >= 1.0
- `fragment_settings.min_log2_ratio`: Real (can be negative)
- `fragment_settings.min_top_n`: Array of 2 Integers, each > 0
- `fragment_settings.n_isotopes`: Integer, [1, 3]
- `scoring_settings.n_train_rounds`: Integer, > 0
- `scoring_settings.max_iterations`: Integer, > 0
- `scoring_settings.max_q_value_probit_rescore`: Real, (0.0, 1.0]
- `scoring_settings.max_PEP`: Real, (0.0, 1.0]
- `irt_mapping.max_prob_to_impute_irt`: Real, (0.0, 1.0]
- `irt_mapping.fwhm_nstd`: Real, > 0
- `irt_mapping.irt_nstd`: Real, > 0

#### Quant Search
- `fragment_settings.min_count`: Integer, > 0
- `fragment_settings.min_y_count`: Integer, >= 0
- `fragment_settings.max_rank`: Integer, [1, 255]
- `fragment_settings.min_spectral_contrast`: Real, [0.0, 1.0]
- `fragment_settings.min_log2_ratio`: Real (can be negative)
- `fragment_settings.min_top_n`: Array of 2 Integers, each > 0
- `fragment_settings.n_isotopes`: Integer, [1, 3]
- `chromatogram.smoothing_strength`: Real, >= 0
- `chromatogram.padding`: Integer, >= 0
- `chromatogram.max_apex_offset`: Integer, >= 0

#### Acquisition
- `nce`: Integer, [10, 50] (typical range)
- `quad_transmission.fit_from_data`: Boolean
- `quad_transmission.overhang`: Real, [0.0, 1.0]
- `quad_transmission.smoothness`: Real, > 0

#### RT Alignment
- `n_bins`: Integer, > 0
- `bandwidth`: Real, (0.0, 1.0]
- `sigma_tolerance`: Integer, > 0
- `min_probability`: Real, (0.0, 1.0]

#### Optimization
- `deconvolution.lambda`: Real, >= 0
- `deconvolution.reg_type`: String, one of ["none", "L1", "L2"]
- `deconvolution.huber_delta`: Real, > 0
- `deconvolution.huber_exp`: Real, > 0
- `deconvolution.huber_iters`: Integer, > 0
- `deconvolution.newton_iters`: Integer, > 0
- `deconvolution.bisection_iters`: Integer, > 0
- `deconvolution.outer_iters`: Integer, > 0
- `deconvolution.newton_accuracy`: Real, > 0
- `deconvolution.max_diff`: Real, > 0
- `machine_learning.max_psms_in_memory`: Integer, > 0
- `machine_learning.min_trace_prob`: Real, (0.0, 1.0]
- `machine_learning.max_q_value_xgboost_mbr_rescore`: Real, (0.0, 1.0]
- `machine_learning.min_PEP_neg_threshold_xgboost_rescore`: Real, (0.0, 1.0]
- `machine_learning.spline_points`: Integer, > 0
- `machine_learning.interpolation_points`: Integer, > 0

#### Protein Inference
- `min_peptides`: Integer, >= 1

#### MaxLFQ
- `run_to_run_normalization`: Boolean

#### Output
- `write_csv`: Boolean
- `write_decoys`: Boolean
- `delete_temp`: Boolean
- `plots_per_page`: Integer, > 0

#### Paths
- `ms_data`: String, must be valid directory path
- `library`: String, must be valid file path ending in .pion
- `results`: String, valid directory path (will be created if doesn't exist)

## Timeline Estimate

- Phase 1 (Infrastructure): 2-3 hours
- Phase 2 (Core Updates): 2-3 hours
- Phase 3 (Validation): 1-2 hours
- Phase 4 (CLI): 1 hour
- Phase 5 (Integration Test Updates): 30 minutes
- Phase 6 (Documentation): 1-2 hours

Total: ~7-11 hours of development time

## Success Criteria

1. Users can run SearchDIA with a minimal 50-line parameter file
2. All existing parameter files continue to work
3. No performance impact from default merging
4. Clear documentation and examples for both modes
5. Comprehensive test coverage for new functionality