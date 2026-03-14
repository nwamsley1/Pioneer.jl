# FDR Calibration Simulation Findings

## Summary

`fdr_calibration_simulation.jl` tests whether cross-validation and multi-file aggregation methods preserve calibrated FDR in a target-decoy framework using synthetic data from two multivariate Gaussian distributions.

## Experimental Design

- **Distribution A**: MVN(+μ/2, I) — all labeled "target" (true positives)
- **Distribution B**: MVN(-μ/2, I) — half labeled "target" (false targets), half labeled "decoy"
- **Metric**: inflation = false_target_pass_rate / decoy_pass_rate. If calibrated, inflation = 1.0 (both from distribution B, should pass at equal rates). If > 1.0, FDR is anti-conservative.
- **LightGBM 2-fold CV**: Train on fold 1 predict fold 0, train on fold 0 predict fold 1.

## Key Finding 1: CV Prevents LightGBM Overfitting

100-seed test, single-file, n=150k, depth=3:

| Method | Mean inflation | >1.0 |
|--------|---------------|------|
| WITH CV | 0.961 | 7/20 |
| NO CV | 1.058 | 15/20 |

Without CV, deep trees + small samples cause catastrophic FDR inflation:

| Condition | WITH CV inflation | NO CV inflation | NO CV true FDR |
|-----------|-------------------|-----------------|----------------|
| baseline (n=150k, depth=3) | 0.83x | 0.94x | 0.94% |
| deep trees (depth=10) | 0.93x | 1.64x | 1.63% |
| small+deep (n=6k, depth=10) | 1.00x | 47.5x | 47% |
| worst case (n=6k, 50feat, depth=10, 500iter) | 0.00x | 50.0x | 50% |

**Conclusion**: CV is essential for LightGBM scoring. Pioneer's FirstPass and ScoringSearch both use 2-fold CV correctly.

## Key Finding 2: PEP Calibration Introduces Systematic Anti-Conservative Bias

The isotonic regression PEP calibration step uses each sample's own target/decoy label to estimate its own PEP — a form of label leakage analogous to training without CV.

100-seed test, multi-file (5 files), n=15k per file, depth=3, CV enabled:

| Aggregation | Mean inflation | Median | Std | >1.0 | True FDR |
|-------------|---------------|--------|-----|------|----------|
| MAX | 0.987 | 1.000 | 0.272 | 49/100 | 0.99% |
| MEAN | 1.011 | 1.020 | 0.172 | 51/100 | 1.01% |
| RAW LOGODDS | 0.999 | 1.000 | 0.196 | 47/100 | 1.00% |
| PEP+LOGODDS | 1.288 | 1.280 | 0.210 | 91/100 | 1.29% |

- **MAX, MEAN, RAW LOGODDS**: All perfectly calibrated. Inflation goes both directions randomly, mean ≈ 1.0.
- **PEP+LOGODDS**: Systematically anti-conservative. 91/100 seeds show inflation > 1.0. True FDR 1.29% when reporting 1%.

**RAW LOGODDS is optimal**: perfectly calibrated (mean 0.999) with lower variance than MAX (0.196 vs 0.272), meaning more statistical power. It combines evidence across files via log-odds of the raw LightGBM probabilities without any extra fitting step that could leak labels.

**Mechanism**: Isotonic regression PEP calibration maps scores to PEP using target/decoy labels. A false target at score level s contributes label=0 (target), lowering the PEP estimate at that score level, giving itself a better PEP. A decoy contributes label=1, raising its own PEP. This creates systematic differential treatment of false targets vs decoys — the same circular dependency that CV fixes for LightGBM, but in the calibration step.

## Recommendation for Pioneer

Replace `PEPCalibratedAggregation` (isotonic regression PEP + log-odds averaging) with raw log-odds aggregation of LightGBM probabilities in `aggregate_prescore_globally!`. This eliminates the only remaining source of label leakage in the pipeline while maintaining statistical power from multi-file evidence combination.
