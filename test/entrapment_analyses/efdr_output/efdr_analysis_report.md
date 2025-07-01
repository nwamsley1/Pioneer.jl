# Empirical FDR Analysis Report

Generated: 2025-07-01T12:06:04.455

## Data Summary

- Original precursor results: 124737
- Analyzed precursor results: 124737
- Library precursors: 3509386

### Library Composition

| Entrapment Group | Count | Percentage |
|-----------------|-------|------------|
| Target (0) | 1754698 | 50.00% |
| Entrapment (1) | 1754688 | 50.00% |

### Results Composition

| Entrapment Group | Count | Percentage |
|-----------------|-------|------------|
| Target (0) | 123959 | 99.38% |
| Entrapment (1) | 778 | 0.62% |

## EFDR Method Comparison

### global_prob / global_qval

| Threshold | Q-val IDs | Actual FDR | Combined IDs | Combined EFDR | Paired IDs | Paired EFDR |
|-----------|-----------|------------|--------------|---------------|------------|-------------|
| 0.001 | 36236 | 0.0011 | 36236 | 0.0022 | 36236 | 0.0019 |
| 0.010 | 46206 | 0.0068 | 46206 | 0.0135 | 46206 | 0.0114 |
| 0.050 | 48849 | 0.0117 | 48849 | 0.0233 | 48849 | 0.0195 |
| 0.100 | 48849 | 0.0117 | 48849 | 0.0233 | 48849 | 0.0195 |

### prec_prob / qval

| Threshold | Q-val IDs | Actual FDR | Combined IDs | Combined EFDR | Paired IDs | Paired EFDR |
|-----------|-----------|------------|--------------|---------------|------------|-------------|
| 0.001 | 92668 | 0.0009 | 92668 | 0.0017 | 92668 | 0.0016 |
| 0.010 | 124737 | 0.0062 | 124737 | 0.0125 | 124737 | 0.0110 |
| 0.050 | 124737 | 0.0062 | 124737 | 0.0125 | 124737 | 0.0110 |
| 0.100 | 124737 | 0.0062 | 124737 | 0.0125 | 124737 | 0.0110 |

## Calibration Analysis

Mean absolute calibration errors:

| Method | Mean Calibration Error |
|--------|----------------------|
| global_prob_combined_efdr | 0.0009 |
| prec_prob_combined_efdr | 0.0007 |
| global_prob_paired_efdr | 0.0005 |
| prec_prob_paired_efdr | 0.0006 |

## Plots

The following plots have been generated:

- ![EFDR Comparison for global_prob](efdr_comparison_global_prob.png)
- ![EFDR Comparison for prec_prob](efdr_comparison_prec_prob.png)
- ![All EFDR Comparisons](efdr_comparison_all.png)

## Analysis Parameters

- EFDR Methods: Combined, Paired
- Score/Q-value pairs analyzed: global_prob/global_qval, prec_prob/qval
- Output directory: `efdr_output`
