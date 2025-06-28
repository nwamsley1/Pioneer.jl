# EntrapmentAnalysis Module

A Julia module for performing empirical false discovery rate (EFDR) analysis on proteomics data using entrapment sequences.

## Installation

The module is designed to be used within the Pioneer.jl project. Add the module directory to your Julia load path:

```julia
push!(LOAD_PATH, "/path/to/Pioneer.jl/test/entrapment_analyses")
using EntrapmentAnalysis
```

## Quick Start

```julia
using EntrapmentAnalysis

# Run EFDR analysis
results = run_efdr_analysis(
    "path/to/precursors_long.arrow",
    "path/to/library_precursors.arrow";
    output_dir = "efdr_output"
)
```

## Main Features

- **Entrapment Pairing**: Automatically pairs target sequences with entrapment sequences
- **EFDR Calculation**: Implements both Combined and Paired EFDR methods
- **Comprehensive Analysis**: Generates comparison tables, calibration metrics, and plots
- **Markdown Reports**: Creates detailed analysis reports with embedded visualizations

## Module Structure

```
entrapment_analyses/
├── EntrapmentAnalysis.jl    # Main module file
├── src/
│   ├── core/                # Core functionality
│   │   ├── efdr_methods.jl  # EFDR calculation methods
│   │   ├── entrapment_pairing.jl  # Pairing algorithms
│   │   └── scoring.jl       # Score manipulation functions
│   ├── analysis/            # Analysis functions
│   │   ├── efdr_analysis.jl # EFDR comparison and evaluation
│   │   └── calibration.jl   # Calibration error calculation
│   ├── plotting/            # Visualization
│   │   └── efdr_plots.jl    # EFDR plotting functions
│   └── api.jl              # Main API function
├── test/                    # Test files
│   ├── runtests.jl         # Test runner
│   └── test_entrapment_pairing.jl
└── examples/               # Usage examples
    ├── basic_usage.jl
    └── yeast_analysis.jl
```

## API Reference

### Main Function

```julia
run_efdr_analysis(prec_results_path, library_precursors_path; kwargs...)
```

**Arguments:**
- `prec_results_path`: Path to Arrow file with precursor results
- `library_precursors_path`: Path to Arrow file with library precursors

**Keyword Arguments:**
- `output_dir`: Output directory (default: "efdr_out")
- `method_types`: EFDR methods to use (default: [CombinedEFDR, PairedEFDR])
- `score_qval_pairs`: Score/q-value column pairs to analyze
- `r_lib`: Ratio of library to real entrapments (default: 1.0)
- `plot_formats`: Plot file formats (default: [:png, :pdf])
- `verbose`: Print progress messages (default: true)

### Core Functions

- `assign_entrapment_pairs!(df)`: Assign pair IDs to entrapment groups
- `add_efdr_columns!(df, library; kwargs...)`: Calculate and add EFDR columns
- `plot_efdr_comparison(df, score_col, qval_col; kwargs...)`: Create EFDR plots

## Examples

See the `examples/` directory for complete working examples:
- `basic_usage.jl`: Simple examples of module usage
- `yeast_analysis.jl`: Real-world analysis with yeast data

## Testing

Run the test suite:

```julia
cd("test/entrapment_analyses/test")
include("runtests.jl")
```

## Output

The analysis generates:
- Markdown report with comprehensive analysis summary
- EFDR comparison plots in multiple formats
- Processed data file with EFDR columns added
- Calibration analysis results

All outputs are saved to the specified output directory.