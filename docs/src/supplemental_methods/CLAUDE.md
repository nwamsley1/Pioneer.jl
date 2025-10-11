# Supplemental Methods Documentation

## Overview

`all_methods.tex` is a comprehensive LaTeX draft of supplemental methods for a manuscript describing **Pioneer and Altimeter**, a DIA (Data-Independent Acquisition) proteomics analysis software suite optimized for narrow isolation windows.

## Document Structure

The document is formatted for Nature journal submission (`sn-nature` document class) with ~1079 lines and includes detailed mathematical formulations for all algorithms.

### Main Sections

1. **Altimeter Training Data & Model** (Lines 51-149)
   - Training dataset from ProteomeTools project (PRIDE archives)
   - Database searching with sage v0.14.7
   - Spectrum filtering and quality control
   - Fragment ion annotation and deisotoping
   - NCE (normalized collision energy) alignment
   - Transformer-based model architecture (11.8M parameters)
   - Cubic B-spline coefficients for fragment intensity prediction
   - Koina framework deployment with 4 inference variants

2. **Pioneer - File Conversion** (Lines 168-170)
   - Cross-platform MS file conversion to Apache Arrow IPC format
   - Supports Thermo .raw and .mzML files
   - PioneerConverter tool (GitHub: nwamsley1/PioneerConverter)

3. **Spectral Library Generation** (Lines 175-184)
   - FASTA protein sequence digestion
   - Target-decoy sequence generation (reverse/shuffle)
   - Entrapment sequence integration for FDR calibration
   - Protein group definition (mathematical formulation)

4. **Intensity-Aware Fragment Index Search** (Lines 188-283)
   - Modified fragment-index algorithm accounting for library intensities
   - Hierarchical bin structure (retention time → fragment m/z)
   - Score counter data structure
   - Binary search-based MS/MS scan queries

5. **Parameter Tuning** (Lines 293-381)
   - Pre-search for run-specific parameter estimation
   - Retention time alignment via B-splines
   - Mass error/tolerance estimation (exponential distributions)
   - NCE alignment with piecewise linear model

6. **First Pass Search** (Lines 386-526)
   - Iterative training procedure based on Percolator
   - Probit regression model for PSM scoring
   - PEP (Posterior Error Probability) estimation via wPAVA
   - PSM aggregation across runs
   - Refined RT alignment and tolerance estimation

7. **Spectral Deconvolution** (Lines 534-598)
   - Matrix representation of library/empirical spectra
   - Sparse column-major layout for efficiency
   - Pseudo-Huber loss minimization with non-negativity constraints
   - Coordinate descent optimization with Newton-Raphson solver
   - Hot-start initialization

8. **Target-Decoy Model & Match-Between-Runs** (Lines 602-860)
   - LightGBM models with cross-validation
   - Iterative training with negative mining
   - MBR (Match-Between-Runs) features:
     - RV coefficient for chromatographic similarity
     - Retention time differences
     - Intensity ratios
     - Cross-run evidence
   - False Transfer Rate (FTR) filtering
   - Two-stage FDR control (global + experiment-wide q-values)
   - Probability aggregation (trace → precursor-run → global)

9. **Chromatogram Quantification** (Lines 864-906)
   - Whittaker-Henderson smoothing
   - Apex refinement
   - Peak boundary detection via second derivative
   - Baseline subtraction
   - Trapezoidal integration

10. **Protein Inference & Quantification** (Lines 911-1067)
    - Parsimony-based inference algorithm
    - Two-phase approach: unique peptides → greedy set cover
    - Bipartite graph decomposition via DFS
    - LightGBM-based protein scoring (optional)
    - MaxLFQ algorithm for label-free quantification

11. **Fragment Isotope Correction** (Lines 1073-end)
    - Conditional fragment isotope probabilities
    - Quadrupole-filtered precursor isotope distribution
    - Re-isotoping of library spectra

## Key Features

### Mathematical Rigor
- Formal definitions using set notation
- Detailed algorithm pseudocode (e.g., protein inference algorithm)
- Loss functions, optimization objectives clearly stated
- Statistical models (probit regression, isotonic regression, exponential fits)

### Cross-References to Code
The document is designed to map to the Julia implementation in Pioneer.jl. Line numbers in the protein inference algorithm (lines 117-380) appear to reference actual source code.

### FDR Control Strategy
Multi-layered approach:
- Spectrum-level (Sage searches)
- PSM-level (Percolator-style iterative training)
- Precursor-level (two-stage global + experiment-wide)
- Protein-level (target-decoy competition)
- Transfer-level (FTR for match-between-runs)

### Novel Contributions
1. **Intensity-aware fragment indexing** - accounts for predicted intensities in scoring
2. **NCE alignment** - maps vendor collision energies to PROCAL Lumos scale
3. **Conditional fragment isotopes** - models isotopes based on quadrupole transmission
4. **FTR filtering** - controls false transfers in match-between-runs
5. **Integrated workflow** - spectral prediction (Altimeter) → search → quantification (Pioneer)

## Technical Details

### Authors
- Nathan T. Wamsley (lead)
- Emily M. Wilkerson
- Ben Major
- Dennis Goldfarb (corresponding)

Washington University School of Medicine

### Software Stack
- **Altimeter**: PyTorch transformer model (UniSpec architecture)
- **Pioneer**: Julia implementation
- **Search**: Sage v0.14.7
- **ML Models**: LightGBM, probit regression
- **Quantification**: MaxLFQ
- **File I/O**: Apache Arrow IPC format

### Data Sources
- ProteomeTools PRIDE datasets: PXD021013, PXD010595, PXD004732, PXD006832
- UniProt human reference proteome (2024-06-04)

## Implementation Notes

When working with Pioneer.jl code:
- Protein inference logic is in `src/utils/proteinInference.jl`
- Line numbers in the algorithm pseudocode (lines 117-380) map to source
- The two-phase approach (unique peptides first, then greedy set cover) is the core logic
- Connected components are discovered via DFS before inference
- Quantification flags (`use_for_quant`) distinguish unique vs. shared peptides

## Document Status

**DRAFT** - This is supplemental methods for an in-preparation manuscript. Expect:
- Missing cross-references (marked with `\ref{REFHERE}`, `\ref{ref here}`)
- Placeholder citations (e.g., "LightGBM paper...")
- TODO markers in red text (e.g., `\textcolor{red}{...}`)
- Potential formatting/equation adjustments before submission
