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

- **Identification recovery** = (standards identified with >=1 passing isomer) / (standards in the
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
  *probability* in [0,1]. They **filter to sites with probability >= 0.75** and report FLR on that
  confident subset. They independently validated the FLR using the **phospho-proline decoy** trick
  (allow phospho on proline, which can't be phosphorylated; any proline hit is a false localization).

**What "cutoff" means.** Every localized site carries a confidence score. A *cutoff* is a threshold
on it: report only sites scoring >= the cutoff, set the rest aside as ambiguous. So "FLR X% @ 0.75"
means *among sites we were >= 0.75-confident about, X% were on the wrong residue*. The cutoff trades
coverage for accuracy, so **every FLR number is meaningless without its cutoff** — raising the
cutoff keeps fewer, more-trustworthy sites (lower FLR, lower coverage).

- **Reported success — the exact numbers (Source Data, Fig. 2I; the 225 standards):**

  | `PTM.SiteProbability` cutoff | 0.75 | 0.80 | 0.85 | 0.90 | 0.95 | 0.99 | 0.999 |
  |---|---|---|---|---|---|---|---|
  | **FLR %** (Astral + Spectronaut) | 3.95 | 3.37 | 2.91 | 2.51 | 1.80 | 0.99 | 0.38 |

  Field context (Suppl. Fig. 3A, FLR @ 0.75): **5.8 %** Bekker-Jensen *reported*, **4.6 %** their
  re-analysis, **3.95 %** this Astral study.
- **Coverage is NOT published as a fraction.** Suppl. Fig. 3B plots only an *absolute count* of
  "average correct precursors" vs cutoff (196.8 at 0.75 → 158.8 at 0.999), with **no denominator**
  (total identified before filtering) and no point below 0.75. Deriving total-passing from the FLR
  gives ~205 precursors at 0.75. So there is **no coverage % to quote** — the paper's coverage claim
  is only the qualitative "retain the majority of correct precursors."
- **Overall depth (context):** ~29,190 phosphosites from the Spectronaut DIA method; a single
  30-min HEK run gave 12,327 phosphopeptides / 9,537 sites; the mouse atlas reached 81,120 sites.

Source: Nat Commun 2024, s41467-024-51274-0, Source Data (MOESM8) sheets `Fig2I`, `SFig3A`, `SFig3B`.

# 6. Our result (Pioneer), and the honest comparison

Searching the 10,000-amol run with a Pioneer library that includes the standards (Prosit-PTM
fragments + Chronologer RT, 1 missed cleavage, all positional isomers). Localization = **select the
isomer with the highest apex deconvolution weight, then use the winner's per-scan sibling
weight-fraction (`iso_weight_fraction_at_scan`) as the localization confidence** (Idea-1 Phase A).

| localization rule | FLR | coverage |
|---|---|---|
| best q-value (no confidence) | 20.7 % | 100 % |
| **max apex weight** (no cutoff) | **6.6 %** | 100 % |
| max weight + **fraction cutoff @ 0.75** | **2.5 %** | 81.3 % |

- **Identification:** 198 / 212 standards = **93 %** — excellent, unchanged.
- **Localization** now has both a decision rule (max weight) and a **filterable confidence**
  (sibling weight-fraction). The full FLR-vs-cutoff curve: 6.6 % @ 0 (all), 3.4 % @ 0.60,
  2.9 % @ 0.70, **2.5 % @ 0.75**, ~2.7 % beyond (small-N floor).

**Head-to-head at the 0.75 operating point (same standards benchmark):**

| | FLR @ 0.75 | coverage @ 0.75 |
|---|---|---|
| Paper (Astral + Spectronaut) | 3.95 % | ~205 precursors (no fraction published) |
| **Pioneer (max wt + fraction cut)** | **2.5 %** | 81.3 % (161/198 sequences) |

Pioneer's FLR sits **inside** the field's 1–5 % band and is a touch lower than the paper's 3.95 % at
0.75. **Two caveats keep this honest:** (a) the two "0.75" cutoffs are **not the same stringency** —
Spectronaut's is a calibrated *probability*, ours an uncalibrated *weight fraction*, so equal cutoff
values != equal operating points; and (b) units differ — their coverage is per-precursor, ours
per-sequence (198), and we credit either real site for the isomer-pair standards. The defensible
claim is *matched FLR range at a comparable-but-not-identical operating point*, not "we beat
Spectronaut."

Note the earlier plan gap is closed: the old "23 % FLR, no cutoff" number reflected best-*q-value*
selection with no confidence filter. The finding was that (i) **max apex weight**, not q-value, is
the right decision rule (23 %→6.6 %), and (ii) the **sibling weight-fraction of the winner** is a
good confidence for the cutoff (6.6 %→2.5 % @ 0.75) — but a bad *selection* rule on its own
(confounded by isomers that find a lonely scan and score fraction=1 trivially). A light **L1/L2
regularization** sweep on the deconvolution is being evaluated for whether ridge sharpens the
fraction/cutoff curve.

# 7. Bottom line

- **Sample:** 225 synthetic phosphopeptides (known sequence + known site) spiked into a yeast
  phosphopeptide background, 15-min 2-Th Astral DIA, a 5-point dilution series (we used the 10k-amol
  triplicate).
- **Ground truth:** the published site list (Supplementary Data 1) — we know the correct residue for
  every spike-in.
- **Assessment:** build a library with all positional isomers, search, and compare the chosen site
  to truth -> **ID recovery** and **false-localization rate**.
- **Success bar:** the paper's Spectronaut pipeline achieves **3.95 % FLR @ 0.75** (range 1–5 %
  across cutoffs), with a localization probability score.
- **Where we stand:** Pioneer **identifies** 93 % of standards and, with max-weight selection + the
  sibling weight-fraction as a confidence cutoff, reaches **2.5 % FLR @ 0.75 (81 % coverage)** —
  inside the field's band. Caveat: the 0.75 cutoffs are comparable but not identical stringency
  (calibrated probability vs uncalibrated fraction).
