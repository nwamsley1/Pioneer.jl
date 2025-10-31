# Supplemental Methods Review: all_methods.tex

**Document Reviewed:** `all_methods.tex`
**Review Date:** 2025-10-30
**Total Issues Found:** 28

## Table of Contents

1. [Missing or Incomplete Cross-References](#1-missing-or-incomplete-cross-references) (6 issues)
2. [Hard-Coded Section/Equation Numbers](#2-hard-coded-sectionequation-numbers) (5 issues)
3. [Typos and Duplicate Words](#3-typos-and-duplicate-words) (4 issues)
4. [Incomplete Citations](#4-incomplete-citations) (1 issue)
5. [Terminology Inconsistency](#5-terminology-inconsistency) (1 category)
6. [Potential Conceptual Confusion](#6-potential-conceptual-confusion) (1 major issue)
7. [Mathematical Notation Issues](#7-mathematical-notation-issues) (2 issues)
8. [Unclear Explanations](#8-unclear-explanations) (2 issues)
9. [Index Inconsistency](#9-index-inconsistency) (2 issues)

---

## 1. Missing or Incomplete Cross-References

### Issue 1.1: Missing reference to linear regression section
- **Line:** 867
- **Priority:** HIGH
- **Current text:**
  ```latex
  Pioneer does this by first repeating the linear regression of mass spectra onto the library spectra \ref{REFHERE}
  ```
- **Problem:** Placeholder reference `\ref{REFHERE}`
- **Recommended fix:** Replace with `\ref{sec:linear_regression}` (defined at line 573)
- **Status:** [x] Fixed

### Issue 1.2: Missing reference to isotope trace selection
- **Line:** 873
- **Priority:** HIGH
- **Current text:**
  ```latex
  The isotope trace pattern to integrate for each precursor is selected as described in \ref{ref here}.
  ```
- **Problem:** Placeholder reference `\ref{ref here}`
- **Recommended fix:** Replace with `\ref{subsubsec:best-isotope-trace-selection}` (defined at line 825)
- **Status:** [x] Fixed

### Issue 1.3: Empty reference to isotope correction
- **Line:** 875
- **Priority:** HIGH
- **Current text:**
  ```latex
  See section \ref{} for documentation on correcting the abundance of separate isotope traces.
  ```
- **Problem:** Empty reference `\ref{}`
- **Recommended fix:** Replace with `\ref{subsec:fragment-isotope-correction}` (defined at line 1078)
- **Status:** [x] Fixed

### Issue 1.4: Outdated label reference in protein inference input
- **Line:** 919
- **Priority:** HIGH
- **Current text:**
  ```latex
  from the two-stage FDR control described in Section~\ref{subsec2}
  ```
- **Problem:** References old label `subsec2` that was changed during label standardization
- **Recommended fix:** Replace with `\ref{subsec:target-decoy-mbr-fdr}` (defined at line 603)
- **Status:** [x] Fixed

### Issue 1.5: Outdated label reference in protein-peptide bipartite graph
- **Line:** 932
- **Priority:** HIGH
- **Current text:**
  ```latex
  Let $\mathcal{A}$ denote the set of all protein accession numbers in the spectral library (as defined in Section~\ref{subsec2})
  ```
- **Problem:** References old label `subsec2` - ambiguous which section is intended
- **Recommended fix:**
  - If referring to protein groups definition: Replace with `\ref{subsec:spectral-library-sequence-generation}` (line 176)
  - If referring to the set $\mathcal{A}$ definition: Add explicit label where $\mathcal{A}$ is first defined
- **Status:** [x] Fixed

### Issue 1.6: Missing reference to MBR iteration explanation
- **Line:** 754
- **Priority:** MEDIUM
- **Current text:**
  ```latex
  Let $\text{prob}_{p,i,r}^{\text{test},(2)}$ be the out-of-fold probabilities from iteration 2 (before full MBR feature refinement)
  ```
- **Problem:** "iteration 2" is mentioned without reference to where iterations are explained
- **Recommended fix:** Add reference: `from iteration 2 (before full MBR feature refinement, see Section~\ref{subsubsec:iterative-model-training})`
- **Status:** [x] Fixed

---

## 2. Hard-Coded Section/Equation Numbers

### Issue 2.1: Red placeholder text for precursor isotope estimation
- **Line:** 80
- **Priority:** HIGH
- **Current text:**
  ```latex
  precursor isotope contributions were estimated using the approach described \textbf{\textcolor{red}{elsewhere in the Methods}}
  ```
- **Problem:** Red placeholder text indicating missing reference
- **Recommended fix:** Replace with `described in Section~\ref{subsubsec:estimating-precursor-isotope-abundances}`
- **Status:** [x] Fixed

### Issue 2.2: Hard-coded equation number reference
- **Line:** 1100
- **Priority:** MEDIUM
- **Current text:**
  ```latex
  In practice, Pioneer truncates the sum in equation 37 at $p=5$
  ```
- **Problem:** Hard-coded equation number "37" instead of using LaTeX reference
- **Recommended fix:** Replace with `in Equation~\eqref{eq:fragment-isotope-conditional}` and add label to equation at line 1097
- **Status:** [x] Fixed

### Issue 2.3: Hard-coded subsection numbers
- **Line:** 1129
- **Priority:** MEDIUM
- **Current text:**
  ```latex
  We propose a solution in subsection 4.14.2. Subsection 4.14.1 describes how to estimate each $x_{i}^{k}$.
  ```
- **Problem:** Hard-coded subsection numbers "4.14.2" and "4.14.1"
- **Recommended fix:** Replace with:
  ```latex
  We propose a solution in Section~\ref{subsubsec:asymmetric-generalized-bell-function}. Section~\ref{subsubsec:estimating-precursor-isotope-abundances} describes how to estimate each $x_{i}^{k}$.
  ```
- **Status:** [x] Fixed

### Issue 2.4: Hard-coded equation number in quadrupole transmission
- **Line:** 1188
- **Priority:** MEDIUM
- **Current text:**
  ```latex
  Pioneer minimizes the squared error between the respective natural logarithms of $R(x_0, x_1)$ and the data as in equation (33).
  ```
- **Problem:** Hard-coded equation number "(33)"
- **Recommended fix:** Add label `\label{eq:quadrupole-ratio}` to equation at line 1126, then replace with `as in Equation~\eqref{eq:quadrupole-ratio}`
- **Status:** [x] Fixed

### Issue 2.5: External equation reference
- **Line:** 1094
- **Priority:** LOW (external reference)
- **Current text:**
  ```latex
  and given equation 3.5 from Goldfarb, D. we have the following:
  ```
- **Problem:** References external equation from citation
- **Recommended fix:** Keep as-is (external reference is acceptable) OR add citation: `equation 3.5 from \citet{Goldfarb2018-ai}`
- **Status:** [x] Fixed

---

## 3. Typos and Duplicate Words

### Issue 3.1: Duplicate word "isotopes"
- **Line:** 1100
- **Priority:** HIGH
- **Current text:**
  ```latex
  since the first 5 isotopes isotopes account for
  ```
- **Problem:** Duplicate word "isotopes isotopes"
- **Recommended fix:** Remove one instance: `since the first 5 isotopes account for`
- **Status:** [x] Fixed

### Issue 3.2: Duplicate word "the"
- **Line:** 1137
- **Priority:** HIGH
- **Current text:**
  ```latex
  conditional abundance of $l$-th isotope of the the $j$-th fragment
  ```
- **Problem:** Duplicate word "the the"
- **Recommended fix:** Remove one instance: `conditional abundance of $l$-th isotope of the $j$-th fragment`
- **Status:** [x] Fixed

### Issue 3.3: Grammar error "for and observed"
- **Line:** 1149
- **Priority:** HIGH
- **Current text:**
  ```latex
  Now for and observed MS/MS spectrum
  ```
- **Problem:** Grammar error "for and" should be "for an"
- **Recommended fix:** `Now for an observed MS/MS spectrum`
- **Status:** [x] Fixed

### Issue 3.4: Duplicate word "using"
- **Line:** 1086
- **Priority:** HIGH
- **Current text:**
  ```latex
  using using sulfur-specific splines
  ```
- **Problem:** Duplicate word "using using"
- **Recommended fix:** Remove one instance: `using sulfur-specific splines`
- **Status:** [x] Fixed

---

## 4. Incomplete Citations

### Issue 4.1: Placeholder citation for LightGBM
- **Line:** 605
- **Priority:** HIGH
- **Current text:**
  ```latex
  Pioneer trains LightGBM target-decoy discrimination models using cross-validation to score each precursor isotope trace \cite{LightGBM paper...}.
  ```
- **Problem:** Incomplete placeholder citation `\cite{LightGBM paper...}`
- **Recommended fix:** Replace with proper citation. The official LightGBM paper is:
  ```latex
  \cite{Ke2017-lightgbm}
  ```
  Citation: Ke, G., Meng, Q., Finley, T., Wang, T., Chen, W., Ma, W., Ye, Q., & Liu, T. Y. (2017). LightGBM: A highly efficient gradient boosting decision tree. In Advances in Neural Information Processing Systems (pp. 3146-3154).
- **Status:** [x] Fixed

---

## 5. Terminology Inconsistency

### Issue 5.1: "Parent ion" vs "Precursor" terminology
- **Lines:** 211-213, 537, 539, 544
- **Priority:** MEDIUM
- **Current text:** Uses "parent ion" in several locations:
  - Line 211: `$P_k \in \mathbb{N}$ uniquely identifies the parent ion`
  - Line 212: `$M_k \in \mathbb{R}^+$ is the parent ion m/z`
  - Line 213: `$Z_k \in  \mathbb{N}^+$ is the parent ion charge`
  - Line 537: `for parent ions within a specified mass-to-charge`
  - Line 539: `$j_k^{(F)}$ identifies the parent ion for each fragment`
  - Line 544: `$j^{(F_{m})}, j_k^{(F_{u})}$ identifies the parent ion`
- **Problem:** Rest of document consistently uses "precursor" terminology
- **Recommended fix:** Replace all instances of "parent ion" with "precursor" for consistency
- **Impact:** This affects Fragment Index Search (lines 189-284) and Matrix Representation sections (lines 535-599)
- **Status:** [ ] Fixed

---

## 6. Potential Conceptual Confusion

### Issue 6.1: Dual definition of "protein group"
- **Lines:** 179 (first definition) and 913+ (second context)
- **Priority:** HIGH (conceptual clarity)
- **Problem:** "Protein group" is defined in two different contexts with potentially different meanings:

  **Definition 1 (Line 179)** - During spectral library construction:
  ```latex
  Pioneer defines a protein group as the maximal set of UniProt accession numbers
  in the sequence database that uniquely correspond to a maximal set of peptide
  sequences through the sequence digestion process.
  ```
  This is a **pre-computed** grouping based on theoretical digestion of the FASTA database.

  **Definition 2 (Lines 913-1034)** - During protein inference:
  ```latex
  Pioneer implements a two-phase parsimony-based inference algorithm to identify
  a minimal set of protein groups that explains all observed peptides.
  ```
  This is an **inferred** grouping based on observed peptides after FDR control.

- **Conceptual issue:** These appear to be different concepts:
  - Library protein groups: Pre-defined during library construction (all theoretically possible peptides)
  - Inferred protein groups: Discovered from actually observed peptides

- **Questions needing clarification:**
  1. Are library protein groups used to constrain the inference process?
  2. Can inferred protein groups differ from library protein groups?
  3. What happens if observed peptides map to proteins in different library groups?

- **Recommended fix:**
  - Option A: Use distinct terminology (e.g., "library protein group" vs "inferred protein group")
  - Option B: Add a paragraph explicitly describing the relationship between these two definitions
  - Option C: Clarify that line 179 defines groups for CV fold assignment (line 611) but inference independently discovers groups

- **Status:** [ ] Needs author clarification and revision

---

## 7. Mathematical Notation Issues

### Issue 7.1: Undefined parameters in quadrupole transmission function
- **Lines:** 1088-1094
- **Priority:** MEDIUM
- **Current text:**
  ```latex
  Pioneer models the transmission efficiency of each precursor isotope with which
  the quadrupole transmits each precursor isotope. The quadrupole transmission
  function is defined as the following:
  \begin{equation}
      Q(z^{p}; c, w) = \Pr\bigl( P=p \mid Q)
  \end{equation}
  ```
- **Problem:**
  1. Parameters $c$ and $w$ are used in function signature but never defined in the text
  2. Redundant phrasing: "transmission efficiency...with which the quadrupole transmits"
  3. Later (line 1114) they are mentioned as "centered at $c$, with an isolation width, $w$" but not when first introduced

- **Recommended fix:**
  1. Define parameters when first introducing the function:
     ```latex
     The quadrupole transmission function, $Q(z^{p}; c, w)$, is defined as follows,
     where $c$ is the center of the isolation window and $w$ is the isolation width:
     ```
  2. Simplify redundant text to: "Pioneer models the transmission efficiency with which the quadrupole transmits each precursor isotope."

- **Status:** [ ] Fixed

### Issue 7.2: Missing word in sentence about bad transfers
- **Line:** 769
- **Priority:** HIGH
- **Current text:**
  ```latex
  First, target isotope traces transferring identifications their paired decoy isotope traces.
  ```
- **Problem:** Missing the word "to" between "identifications" and "their"
- **Recommended fix:** `target isotope traces transferring identifications to their paired decoy isotope traces`
- **Status:** [ ] Fixed

---

## 8. Unclear Explanations

### Issue 8.1: Undefined statistical notation
- **Line:** 334
- **Priority:** MEDIUM
- **Current text:**
  ```latex
  \delta_{\text{RT}} = \frac{4}{\Phi(3/4)^{-1}} \cdot \text{MAD}
  ```
- **Problem:** $\Phi$ is not defined in the text
- **Clarification needed:** $\Phi$ is the cumulative distribution function (CDF) of the standard normal distribution, and $\Phi^{-1}$ is its inverse (the quantile function)
- **Recommended fix:** Add explanation:
  ```latex
  where $\Phi^{-1}$ is the inverse cumulative distribution function (quantile function)
  of the standard normal distribution.
  ```
- **Status:** [ ] Fixed

### Issue 8.2: Incomplete explanation of slope property
- **Line:** 1176
- **Priority:** LOW
- **Current text:**
  ```latex
  The shape parameters $b_l$ and $b_r$ determine the smoothness of the transition,
  such that the slope $-b_l/(2a_l)$ at $x = a_l$ and $-b_r/(2a_r)$ at $x = a_r$.
  ```
- **Problem:** Sentence is incomplete - "such that the slope..." doesn't finish the thought. What about the slope? Is this the slope at the half-maximum point?
- **Recommended fix:** Complete the sentence or clarify what property of the slope is being described:
  ```latex
  The shape parameters $b_l$ and $b_r$ determine the smoothness of the transition,
  with the derivative at the half-maximum points ($x = \pm a_{l/r}$) being
  $-b_l/(2a_l)$ and $-b_r/(2a_r)$ respectively.
  ```
- **Status:** [ ] Consider revision for clarity

---

## 9. Index Inconsistency

### Issue 9.1: Inconsistent index variable in set definition
- **Line:** 539
- **Priority:** MEDIUM
- **Current text:**
  ```latex
  Pioneer combines the library fragmentation spectra, $L_j$, into a single set,
  $F = \{(I_i^{(F)}, Z_i^{(F)}, j_i^{(F)})\}_{k=1}^{T}$
  ```
- **Problem:** Set elements use index $i$ in superscripts $(I_i^{(F)}, Z_i^{(F)}, j_i^{(F)})$ but the set range uses index $k$ with $\}_{k=1}^{T}$
- **Recommended fix:** Make consistent - either:
  - Change to: `$F = \{(I_k^{(F)}, Z_k^{(F)}, j_k^{(F)})\}_{k=1}^{T}$`, OR
  - Change to: `$F = \{(I_i^{(F)}, Z_i^{(F)}, j_i^{(F)})\}_{i=1}^{T}$`
- **Note:** Check subsequent usage to determine which index is used consistently
- **Status:** [ ] Fixed

### Issue 9.2: Missing subscript in notation
- **Line:** 544
- **Priority:** MEDIUM
- **Current text:**
  ```latex
  \item $j^{(F_{m})}, j_k^{(F_{u})}$ identifies the parent ion
  ```
- **Problem:** First notation `$j^{(F_{m})}$` is missing the subscript $k$ to match the pattern
- **Recommended fix:** `$j_k^{(F_{m})}, j_k^{(F_{u})}$ identifies the parent ion`
- **Status:** [ ] Fixed

---

## Summary Statistics

| Category | Count | High Priority | Medium Priority | Low Priority |
|----------|-------|---------------|-----------------|--------------|
| Missing Cross-References | 6 | 5 | 1 | 0 |
| Hard-Coded Numbers | 5 | 1 | 3 | 1 |
| Typos/Duplicates | 4 | 4 | 0 | 0 |
| Incomplete Citations | 1 | 1 | 0 | 0 |
| Terminology Issues | 1 | 0 | 1 | 0 |
| Conceptual Confusion | 1 | 1 | 0 | 0 |
| Notation Issues | 2 | 1 | 1 | 0 |
| Unclear Explanations | 2 | 0 | 1 | 1 |
| Index Inconsistency | 2 | 0 | 2 | 0 |
| **TOTAL** | **28** | **13** | **9** | **2** |

---

## Recommended Order of Fixes

### Phase 1: Critical Fixes (Do First)
1. Fix all typos and duplicate words (Issues 3.1-3.4)
2. Fix missing cross-references (Issues 1.1-1.6)
3. Add LightGBM citation (Issue 4.1)
4. Fix missing word in line 769 (Issue 7.2)

### Phase 2: Important Clarifications
5. Clarify protein group definition confusion (Issue 6.1) - **requires author input**
6. Replace hard-coded equation/section numbers (Issues 2.1-2.4)
7. Define quadrupole transmission parameters (Issue 7.1)
8. Fix index inconsistencies (Issues 9.1-9.2)

### Phase 3: Polish
9. Standardize terminology (Issue 5.1)
10. Add statistical notation clarification (Issue 8.1)
11. Improve incomplete explanation (Issue 8.2)
12. Consider external citation addition (Issue 2.5)

---

## Notes for Authors

- **High priority items** are errors that would confuse readers or cause compilation issues
- **Medium priority items** improve clarity and maintainability
- **Low priority items** are minor enhancements
- Issue 6.1 (protein group definition) may require input from the lead author about intended meaning
- Consider adding equation labels proactively to all numbered equations for future referencing
- The document currently has 80 equations - ensure all important ones have labels

---

**Review completed:** 2025-10-30
**Reviewer notes:** This is a well-structured mathematical document with excellent use of formalism. The issues identified are primarily housekeeping items (references, typos) rather than fundamental problems with the methodology.