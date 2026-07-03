---
title: "Phosphosite Localization Benchmark"
subtitle: "The experiment, the sample, the spike-ins, and how we score 'site localization'"
author: "Pioneer.jl — prepared for N. Wamsley"
date: "2026-07-03"
geometry: margin=2.3cm
fontsize: 11pt
colorlinks: true
linkcolor: RoyalBlue
urlcolor: RoyalBlue
---

# 1. The question we are testing

Two distinct abilities, often conflated:

1. **Identification** — did we detect that *this peptide* is phosphorylated (somewhere)?
2. **Localization** — did we place the phosphate on the *correct residue* (which S, T, or Y)?

A search engine can be excellent at (1) and poor at (2). Phospho-DIA makes (2) especially hard:
two "positional isomers" of the same peptide (e.g. phospho on S3 vs. T8) have the **same precursor
m/z** and share **most** of their fragment ions — they differ only in the handful of
"site-determining" fragments that span the modified residue. If those specific ions are weak or
missing in the spectrum, the engine cannot tell the two apart.

This document explains the benchmark we used to measure both, on real data, on our target
instrument class (Orbitrap Astral).

# 2. The sample

We use the public dataset from the paper that defines the de-facto DIA phospho-localization
benchmark:

> **"Fast and deep phosphoproteome analysis with the Orbitrap Astral mass spectrometer"**,
> Nature Communications 2024, DOI **10.1038/s41467-024-51274-0** (Coon lab).
> Data: MassIVE **MSV000093613** (the ProteomeXchange accession PXD049028 cited in some summaries
> is mis-attributed — it is a different HAP1 study).

The specific runs we searched are the **synthetic-standard dilution series**:

`20240306_NML_15min_2Th_250ngYeast_{10000,2500,625,156,39}amol_Stds0{1,2,3}.raw`

Each run is a **15-minute, 2-Th-window Orbitrap Astral DIA** acquisition of a two-component sample:

- **Background:** 250 ng of a **yeast** tryptic **phosphopeptide** sample (phospho-enriched). This
  provides realistic proteome complexity — thousands of co-eluting, co-isolated phosphopeptides
  that compete during signal deconvolution.
- **Spike-in:** **225 synthetic phosphopeptide standards** at a known amount, titrated across a
  **5-point, 4-fold dilution series** (10,000 -> 2,500 -> 625 -> 156 -> 39 amol on-column), in
  **triplicate**. The dilution tests how localization degrades as the analyte approaches the limit
  of detection. We ran the **10,000 amol triplicate** (highest load = best case) for this report.

The yeast standards were **not** re-acquired; they are real Astral DIA files we downloaded
(~7 GB each) and converted to Pioneer's Arrow format with your **PioneerConverter** (Thermo
`.raw` -> `.arrow`, no mzML step, ~57 s/file).

# 3. The spike-in sequences (the ground truth)

The 225 standards come from five commercial kits — **JPT SpikeMix PTM-Kit 52 & 54** and
**Sigma MS PhosphoMix 1, 2, 3** (the "JPTstySig123" in the file names). Their **exact sequences
and exact phospho-site positions are published** in the paper's **Supplementary Data 1**
("Reference Sheet — Peptide Level"). That known-site list is what makes this a *ground-truth*
benchmark: for every spike-in we know the right answer.

Representative entries (site in Spectronaut notation; `C##` = fixed carbamidomethyl, not phospho):

| Modified sequence | Phospho site | # phospho |
|---|---|---|
| `YNTYAYVGLTEGPSPGDFR` | Y6 | 1 |
| `SVQEIQATFFYFTPNK` | Y11 | 1 |
| `VFLDCCNYITELR` | Y8 (C5,C6 = carbamidomethyl) | 1 |
| `LVDQNIFSFYLSR` | Y10 | 1 |

**Multiplicity:** 210 are **mono**-phospho, 11 **di**, 3 **tri**, 1 **tetra**. The mono set is the
clean core of the localization test.

**Built-in positional isomers.** Some standards appear as *pairs with the same sequence but
different true sites* (e.g. `ADEPSSEESDLEIDK`, `LASEYLTPEEMVTFK` each occur twice). These are the
hardest, most direct localization challenge — the engine must distinguish two forms that differ
only in which residue carries the phosphate.

# 4. How this lets us assess "site localization"

The logic is:

1. **Put every possibility in the library.** We build the search library with **variable phospho
   on *every* S/T/Y** (`max_var_mods=2`, phospho-only), so for each standard the library contains
   *all* positional isomers — phospho@S3, phospho@T8, phospho@Y11, ... — as separate precursors with
   the same precursor m/z but different predicted fragment patterns (Prosit-PTM). Because the wrong
   isomers are present, the search *can* get it wrong — which is exactly what we want to measure.
2. **Let the search choose.** Each isomer is scored against the DIA spectra; the deconvolution +
   scoring picks a winner (the best-scoring isomer that passes FDR).
3. **Compare to truth.** For each identified standard we take its best-scoring isomer and compare
   its phospho position(s) to the known site(s) from Supplementary Data 1.

This yields two numbers (the **success criteria**):

- **Identification recovery** = (standards identified with ≥1 passing isomer) / (standards in the
  library). "Did we find it at all?"
- **False-localization rate (FLR)** = (identified standards whose best isomer is on the *wrong*
  residue) / (identified standards). "When we found it, did we place the phosphate correctly?"
  For the positional-isomer pairs, localizing to *either* real site counts correct.

**Reference bar:** the paper reports **1–5 % localization error** on this same standard set (more on
that below), so "good" localization means FLR approaching a few percent.

# 5. What the paper used, and how well it did

The paper is a **Spectronaut**-based study (not our engine):

- **Search engine:** **Spectronaut v17.6 / v18.6** — library-based for the HEK293T method
  development, and **directDIA (library-free)** for the mouse atlas and the entrapment/phospho-
  proline validation experiments. They **cross-checked with Proteome Discoverer 3.1 + CHIMERYS**.
- **Localization scoring:** Spectronaut's **`PTM.SiteProbability`** — a per-site localization
  *probability*. They **filter to sites with probability ≥ 0.75** and report FLR on that confident
  subset. They independently validated the FLR using the **phospho-proline decoy** trick (allow
  phospho on proline, which can't be phosphorylated; any proline hit is a false localization).
- **Reported success:** *"depending on localization probability cutoff, the error rate ranges
  between 1 and 5 percent"* on the 225-standard dilution series, while retaining the majority of
  correct precursors even at stringent cutoffs.
- **Overall depth (context):** ~29,190 phosphosites from the Spectronaut DIA method; a single
  30-min HEK run gave 12,327 phosphopeptides / 9,537 sites; the mouse atlas reached 81,120 sites.

**Two things to note about their 1–5 %:** it is (a) produced by a **dedicated localization
probability model**, and (b) reported **after filtering** at a probability cutoff (≥0.75) — i.e. it
is the FLR of the *confidently-localized* subset, trading some sensitivity for localization
confidence.

# 6. Our result (Pioneer), and the honest comparison

Searching the 10,000-amol run with a Pioneer library that includes the standards (Prosit-PTM
fragments + Chronologer RT, 1 missed cleavage, all positional isomers):

| Metric | Pioneer | Paper (Spectronaut) |
|---|---|---|
| Standards **identified** | **198 / 212 = 93.4 %** | very high (deep phospho-ID) |
| **FLR** (localization error) | **23.2 %** (best-isomer, no cutoff) | **1–5 %** (`PTM.SiteProbability` ≥ 0.75) |

**Identification is excellent** — Pioneer + Prosit-PTM recovers 93 % of the ground-truth standards.
**Localization is the gap**, but the comparison is **not yet apples-to-apples**, and the reason *is*
the finding:

- The paper's low FLR relies on a **localization probability score** and a **stringent cutoff** that
  *discards* low-confidence site calls. Pioneer, today, has **neither** — our 23 % is measured on
  **every** identified standard, "best isomer wins," with **no localization-confidence filter**.
- So a large part of the 23 % vs 1–5 % gap is the **absence of a dedicated site-localization scorer
  + probability cutoff**, not a failure of identification. Adding one (Phase 3) — a site-determining-
  ion / AScore-style probability plus a cutoff — is the principled path to the field's 1–5 %.

We also ran a quick side-experiment adding **L1/L2 regularization to the deconvolution solver**: L2
(ridge) at λ≈1 nudged FLR to ~19.6 %, but across a dose–response the effect was noisy — a cheap knob,
not a substitute for a real localization model.

# 7. Bottom line

- **Sample:** 225 synthetic phosphopeptides (known sequence + known site) spiked into a yeast
  phosphopeptide background, 15-min 2-Th Astral DIA, a 5-point dilution series (we used the 10k-amol
  triplicate).
- **Ground truth:** the published site list (Supplementary Data 1) — we know the correct residue for
  every spike-in.
- **Assessment:** build a library with all positional isomers, search, and compare the chosen site
  to truth -> **ID recovery** and **false-localization rate**.
- **Success bar:** the paper's Spectronaut pipeline achieves **1–5 % FLR** (with a localization
  probability score, filtered at ≥0.75).
- **Where we stand:** Pioneer **identifies** 93 % of standards but sits at **23 % FLR** because it
  has **no localization scorer/cutoff yet** — precisely the Phase-3 build this harness now lets us
  measure against.
