# Missing Sections and Content Gaps: all_methods.tex

**Document Reviewed:** `all_methods.tex`
**Review Date:** 2025-10-30
**Review Type:** Content Completeness Analysis

## Purpose

This document identifies **missing content** and **gaps** in the supplemental methods documentation. It complements `methods_review.md` (which focuses on fixing existing content) by cataloging content that should be added to make the methods section complete and reproducible.

---

## Table of Contents

1. [Critical Missing Content](#critical-missing-content) (3 items)
2. [Recommended Additional Sections](#recommended-additional-sections) (6 items)
3. [Document Organization Issues](#document-organization-issues) (1 item)
4. [Action Items Summary](#action-items-summary)

---

## Critical Missing Content

These are essential components that are referenced but not explained, or are necessary for reproducibility.

### 1. Chronologer / Retention Time Prediction Model

**Priority:** 🔴 **CRITICAL**

**Problem:**
Chronologer is mentioned twice in the document but never explained:
- **Line 66:** "retention times were aligned to Chronologer predictions using spline fitting"
- **Line 177:** "Pioneer constructs target-decoy sequence libraries as inputs to Altimeter and Chronologer for *in silico* spectrum prediction"

Throughout the document, "library retention time (iRT)" and "indexed retention time (iRT)" are used extensively (50+ references), but their initial source is never explained.

**What's Missing:**
1. **What is Chronologer?**
   - Is it a separate tool, or part of Pioneer?
   - What is its relationship to Pioneer and Altimeter?

2. **How does it work?**
   - What algorithm/model does it use for retention time prediction?
   - What are its inputs and outputs?
   - What training data was used?
   - Model architecture (if machine learning-based)

3. **How are iRT values assigned?**
   - How does Pioneer obtain initial iRT predictions for library peptides?
   - Are these predictions part of the spectral library generation?
   - What happens if Chronologer predictions are unavailable?

4. **Validation and accuracy**
   - What is the typical prediction accuracy?
   - How robust is it across different instrument types or gradients?

**Where to Add:**
- Option A: New subsection after Altimeter (between lines 149-169): `\subsection{Chronologer Retention Time Prediction}`
- Option B: Clarify in the introduction that Chronologer is an external tool and cite appropriately
- Option C: Add as part of "Spectral Library Sequence Generation" (line 176)

**Impact:** HIGH - Without this information, readers cannot understand:
- Where iRT values come from
- How to reproduce the analysis
- Why retention time alignment works

**Status:** [ ] Needs to be added

---

### 2. Fragment Scoring Function Table (Missing Table)

**Priority:** 🔴 **CRITICAL**

**Problem:**
Referenced at **line 204** but table does not exist:
```latex
Table~\ref{tab:frag_score} gives the scoring function used throughout this work.
```

**What's Missing:**
A table showing the fragment scoring function that assigns weights based on fragment rank. The text says:
> "Assign each fragment a score, $s \in \mathbb{N}$, based on its rank, $r$. The scoring function is user defined and ought to assign greater weights to the more probable fragment ions."

**Required Table Content:**
| Rank ($r$) | Score ($s$) | Description |
|------------|-------------|-------------|
| 1 | ? | Top-ranked (most intense) fragment |
| 2 | ? | Second-ranked fragment |
| 3 | ? | Third-ranked fragment |
| ... | ... | ... |
| N | ? | Nth-ranked fragment |

Or, if formulaic:
- Provide the mathematical formula: $s(r) = f(r)$
- Explain the rationale (e.g., exponential decay, linear, etc.)
- State default values used in Pioneer

**Where to Add:**
- In Section 4 "Intensity-Aware Fragment Index Search" (around line 205)
- Add after the text that references it
- Label as `\label{tab:frag_score}`

**Impact:** HIGH - This directly affects:
- Fragment index scoring
- PSM ranking
- Reproducibility of search results

**Status:** [ ] Needs to be added

---

### 3. Spectral Similarity Scores Table (Missing Table)

**Priority:** 🔴 **CRITICAL**

**Problem:**
Referenced at **line 389** but table does not exist:
```latex
These include the spectral similarity scores listed in Table~\ref{tab:simple_scores}
```

The text describes PSM features but doesn't list the spectral similarity scores.

**What's Missing:**
A table listing all spectral similarity metrics calculated by Pioneer for PSM scoring. Based on line 389, these are used as features in the first-pass search alongside:
- Charge state
- Number of missed cleavages
- Number of variably oxidized methionines
- Total ion current
- Number of matched y ions
- Mean fragment error
- Number of peaks in the spectrum

**Required Table Content:**
| Feature Name | Description | Formula (if applicable) | Range |
|--------------|-------------|------------------------|-------|
| Scribe Score | Intensity-weighted spectral similarity | (cited at line 354) | [0, 1] |
| Cosine Similarity | ? | ? | [0, 1] |
| Spectral Angle | ? | ? | [0, π/2] |
| ... | ... | ... | ... |

**Where to Add:**
- In Section 6 "First Pass Search" (around line 389)
- Add immediately after or before the text that references it
- Label as `\label{tab:simple_scores}`

**Note:** Line 171 also references `tab:simple_scores` for "Arrow table schema" which seems incorrect - this may need to be a different table (`tab:arrow_schema`).

**Impact:** HIGH - These scores are critical features for:
- PSM discrimination
- Probit regression model training
- FDR estimation

**Status:** [ ] Needs to be added

---

## Recommended Additional Sections

These sections would improve completeness and usability but are not critical for understanding the core algorithms.

### 4. Software and Computational Requirements

**Priority:** 🟡 **MEDIUM**

**What's Missing:**
- **Hardware requirements**
  - Minimum RAM (likely substantial for LightGBM models and large datasets)
  - Recommended CPU cores
  - GPU support (if any for Altimeter inference)
  - Disk space requirements

- **Software dependencies**
  - Julia version requirements
  - Required Julia packages
  - Python version (for Altimeter)
  - External tools (e.g., ThermoRawFileParser version)

- **Operating system support**
  - Linux, macOS, Windows compatibility
  - Known platform-specific issues

- **Runtime expectations**
  - Typical runtime for different dataset sizes
  - Scalability considerations
  - Parallelization capabilities

**Where to Add:**
- New subsection at the beginning or end of methods
- Could be combined with "Cross-Platform Mass Spectrometry File Conversion" section

**Impact:** MEDIUM - Helps users determine if they have adequate resources

**Status:** [ ] Consider adding

---

### 5. Output File Formats and Interpretation

**Priority:** 🟡 **MEDIUM**

**What's Missing:**
- **Output file types**
  - What files does Pioneer produce?
  - File naming conventions
  - Directory structure

- **File format specifications**
  - CSV/TSV column descriptions
  - Arrow table schemas for outputs
  - HDF5/Binary format specifications (if any)

- **Result interpretation**
  - How to interpret q-values
  - Understanding protein group IDs
  - Interpreting entrapment results
  - Quality metrics in output

**Where to Add:**
- New subsection: `\subsection{Pioneer Output Files and Results Interpretation}`
- Could be placed before or after "Protein Inference and Quantification"

**Impact:** MEDIUM - Important for users to understand and use results

**Status:** [ ] Consider adding

---

### 6. Supported Post-Translational Modifications

**Priority:** 🟡 **MEDIUM**

**What's Missing:**
- **Fixed modifications**
  - Which modifications are applied by default?
  - Carbamidomethylation of cysteine?

- **Variable modifications**
  - Which PTMs are supported?
  - Maximum number of variable modifications per peptide
  - How modifications affect:
    - Library generation
    - Fragment ion prediction (Altimeter)
    - Mass calculations
    - Retention time prediction (Chronologer)

- **Custom modifications**
  - Can users define custom modifications?
  - What parameters are needed?

- **Localization**
  - Does Pioneer perform PTM localization?
  - If so, what algorithm/scoring is used?

**Current mentions in document:**
- Line 389: "number of variably oxidized methionines" - suggests methionine oxidation is supported
- No other specific modifications mentioned

**Where to Add:**
- New subsection in or near "Spectral Library Sequence Generation" (line 176)
- `\subsubsection{Supported Modifications}`

**Impact:** MEDIUM - Critical for users working with modified peptides

**Status:** [ ] Consider adding

---

### 7. Supported Instruments and Data Acquisition Parameters

**Priority:** 🟠 **MEDIUM-LOW**

**What's Missing:**
- **Instrument compatibility**
  - Thermo instruments: which models?
  - Other vendors through mzML conversion?
  - Instrument-specific considerations

- **Fragmentation methods**
  - HCD (clearly supported - line 66)
  - CID, ETD, EThcD?
  - Does Altimeter support non-HCD fragmentation?

- **Acquisition parameters**
  - DIA window schemes supported
  - Overlap vs. non-overlap windows
  - Narrow windows: how narrow? (title mentions "optimized for narrow isolation windows")
  - Wide DIA compatibility?

- **Resolution requirements**
  - MS1 resolution
  - MS2 resolution
  - Low-resolution compatibility?

**Where to Add:**
- Could expand "Cross-Platform Mass Spectrometry File Conversion" (line 169)
- Or add new subsection: `\subsection{Supported Instruments and Acquisition Methods}`

**Impact:** MEDIUM-LOW - Helps users determine if Pioneer is appropriate for their data

**Status:** [ ] Consider adding

---

### 8. Quality Control Metrics and Diagnostics

**Priority:** 🟠 **MEDIUM-LOW**

**What's Missing:**
- **QC metrics reported**
  - What QC metrics does Pioneer calculate?
  - Where are they reported?

- **Data quality assessment**
  - How to assess if data quality is sufficient?
  - What distributions should be examined?
  - Diagnostic plots or outputs

- **Recommended thresholds**
  - FDR thresholds used in benchmarks
  - When to use stricter/looser thresholds
  - QC metric thresholds for filtering

- **Troubleshooting guidance**
  - What to do if few PSMs are identified?
  - How to diagnose poor RT alignment?
  - When to adjust search parameters?

**Where to Add:**
- New section after methods: `\section{Quality Control and Diagnostics}`
- Or integrate into existing sections where relevant

**Impact:** MEDIUM-LOW - Helps users ensure good results

**Status:** [ ] Consider adding

---

### 9. Limitations and Known Issues

**Priority:** 🟢 **LOW** (but good practice for transparency)

**What's Missing:**
- **Method limitations**
  - Dataset size limitations
  - Computational complexity scaling
  - When Pioneer might not be optimal

- **Known issues**
  - Edge cases
  - Known failure modes
  - Workarounds for common problems

- **Future improvements**
  - Planned enhancements
  - Current development limitations

- **Comparison to other methods**
  - When to use Pioneer vs. other DIA tools?
  - Specific advantages for narrow windows

**Where to Add:**
- Optional section at end of methods: `\subsection{Limitations and Considerations}`
- Or in main manuscript discussion

**Impact:** LOW - Improves transparency and manages expectations

**Status:** [ ] Optional

---

## Document Organization Issues

### 10. Incorrect Table Reference

**Priority:** 🟡 **MEDIUM**

**Problem:**
**Line 171** references `tab:simple_scores` for Arrow table schema:
```latex
Pioneer-compatible arrow tables follow the schema in Table~\ref{tab:simple_scores}.
```

However, `tab:simple_scores` is supposed to contain spectral similarity scores (line 389), not the Arrow schema.

**Solutions:**
1. **Option A:** Create two separate tables
   - `tab:arrow_schema` - Arrow IPC format schema (referenced at line 171)
   - `tab:simple_scores` - Spectral similarity scores (referenced at line 389)

2. **Option B:** Fix the reference at line 171
   - Change to reference a new table: `Table~\ref{tab:arrow_schema}`
   - Add the Arrow schema table in the "File Conversion" section

**Recommended:** Option A - create both tables

**Impact:** MEDIUM - Causes confusion about what the table should contain

**Status:** [ ] Needs correction

---

## Action Items Summary

### Phase 1: Critical Additions (Do First)
- [ ] **1.1** Add Chronologer explanation or clarify it's external and cite
- [ ] **1.2** Create Fragment Scoring Function Table (`tab:frag_score`)
- [ ] **1.3** Create Spectral Similarity Scores Table (`tab:simple_scores`)
- [ ] **1.4** Fix incorrect table reference at line 171 (create `tab:arrow_schema`)

### Phase 2: High-Value Additions (Strongly Recommended)
- [ ] **2.1** Add Software/Computational Requirements section
- [ ] **2.2** Add Output File Formats section
- [ ] **2.3** Add Supported Modifications section

### Phase 3: Completeness Improvements (Optional)
- [ ] **3.1** Add Supported Instruments section
- [ ] **3.2** Add Quality Control Metrics section
- [ ] **3.3** Add Limitations section

---

## Priority Legend

| Symbol | Priority | Description |
|--------|----------|-------------|
| 🔴 | CRITICAL | Referenced but missing, or essential for reproducibility |
| 🟡 | MEDIUM | Important for usability and completeness |
| 🟠 | MEDIUM-LOW | Helpful but not essential |
| 🟢 | LOW | Nice to have for transparency |

---

## Notes for Authors

1. **Chronologer** - This is the #1 priority. Either:
   - Add a full methods section if it's part of your software suite
   - Add a brief explanation and citation if it's external/published
   - Explain if iRT values are obtained differently (e.g., from external tools, manual annotation, etc.)

2. **Missing Tables** - These are straightforward to add and are directly referenced in the text. Creating them will fix broken references and improve clarity.

3. **Additional Sections** - These would make the methods more complete, but prioritize based on:
   - Your target audience
   - Journal requirements
   - Space constraints
   - Whether this information will be in supplementary materials vs. main text

4. **Related Documents**
   - See `methods_review.md` for issues with existing content (typos, references, etc.)
   - This document (`missing_sections.md`) focuses on content that should be added

---

**Review completed:** 2025-10-30
**Reviewer notes:** The algorithmic and mathematical content in the methods is excellent and comprehensive. The main gaps are practical/implementation details and the Chronologer component which is referenced but never explained.