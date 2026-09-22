# timsTOF quantification: the 2D chromatogram and how to integrate it (2026-09-21/22)

Branch `feat/tdfs-reader`. Scripts in `docs/bruker/` (`chrom2d.jl`, `hye_quant.jl`, the `plot_*.jl` set); they run
offline on the weight dumps that `PIONEER_CHROM_DUMP_DIR` writes, so integration strategies can be compared
without re-searching. Working data on Nathan's machine: `~/BrukerTims/pride_hye/`.

## 1. What a timsTOF chromatogram is

With slice data, a precursor's deconvolved weights form a **grid**: retention-time cycles by ion-mobility slices
(stride 8 scans), about 10 cycles × 9–16 slices per precursor. Peaks are compact: 3–4 cycles wide, 2–4 slices
tall, 27–30% of the weight in the apex cycle, 73–75% within ±1 slice of the apex, 97–98% within ±3. Some ridges
tilt (the mobility apex drifts a slice across the RT peak). Pioneer's 1D integrator collapses the mobility axis
first and integrates the trace, which is what the numbers in §4 say is the problem.

## 2. Benchmark

PXD070049 hybrid-proteome runs, Condition A (human/yeast/E. coli 65/30/5 w/w) vs B (65/15/20): expected
log2(A/B) = 0 / +1 / −2. Library: the three-proteome Pioneer library with Koina `AlphaPept_ccs_generic` mobility
added (`lib/add_im_to_lib.jl`), 10.5 M precursors. Twelve files per gradient (6 A + 6 B) in **one** search, MBR
off, converted at IM σ5 / m/z σ3 / stride 8 / top-1500. Raw files are on RIS `NTW/PXD070049_Bruker/` (~1 min per
file to copy; PRIDE takes ~25).

## 3. Strategies compared

All start from the precursor's grid. The smoothed ones solve a weighted 2D Whittaker-Henderson on the unit square
(second divided differences on both axes, λ = 1e-5, the 1D convention).

| name | boundary | baseline | values summed |
|---|---|---|---|
| raw sum | whole window | none | raw |
| **raw 5×5 core** | ±2 cycles × ±2 slices around the PSM-seeded apex | median raw cell outside the core | raw |
| smoothed 5×5 / 7×7 core | same | median smoothed cell outside | smoothed |
| footprint ≥ f·apex | connected cells above a fraction of the apex | median raw cell outside | raw |
| 2D v3 | region grown from the apex by descent (≤ 1.15× parent) | plane fitted on the region's rim | smoothed, 2D trapezoid |
| Pioneer peak_area | current 1D integrator | | |

## 4. Results (12-file joint searches, 6 replicates per condition)

**15 min 50 ng**, 54,371 precursors quantified in all twelve runs:

| strategy | human CV A / B | yeast log2 | E. coli log2 | E. coli within ±0.5 |
|---|---|---|---|---|
| raw 5×5 core − flat baseline | 12.5% / 7.8% | 0.89 | −1.71 | 73% |
| smoothed 5×5 core − flat base | 12.7% / 7.8% | 0.90 | −1.73 | 75% |
| smoothed 7×7 core − flat base | 12.6% / 7.7% | 0.88 | −1.69 | 71% |
| 2D v3 | 16.2% / 12.1% | 0.90 | −1.73 | 71% |
| **Pioneer peak_area (current)** | **84% / 85%** | 0.84 | −1.51 | 28% |
| DIA-NN 2.6.0 (6 reps, MBR, normalised) | 10.5% / 9.8% | 0.866 | −1.576 | 53% |

**5 min 50 ng**, 22,574 precursors: raw 5×5 core 13.6% / 10.9%, yeast 0.83, E. coli −1.62; Pioneer peak_area
106% / 106%. Search time, one warm session: 3,007 s for the twelve 15-min files, 1,460 s for the twelve 5-min
files (DIA-NN: 460 and 231 min on an 18-thread workstation, with MBR).

**Conclusions.**

1. The current 1D integrator is the quantification problem: 84–106% replicate CVs and only a fifth of the
   precursors get an area at all. Any grid-based strategy is a large improvement.
2. Precision and accuracy trade off through *boundary placement*, not smoothing. Smoothing over a fixed block
   changes nothing (WH nearly conserves the sum over a region); what costs variance is re-estimating the apex,
   the region and the baseline per replicate, which is why 2D v3 is the most accurate and the least precise.
3. Recommended: **the raw 5×5 core around the PSM-seeded apex, minus a flat baseline** (the median cell outside
   it). Best combined score, no smoothing needed beyond locating the apex.
4. The residual compression (yeast 0.89 rather than 1.00) survives every strategy, so it is upstream of
   integration, in the deconvolution weights themselves.
5. Next: share the region across replicates (one template per precursor) or fit a parametric surface, which
   should keep v3's accuracy without its placement variance.

## 5. Extraction windows (fixed along the way)

Measured from the dumps: the RT window (10 cycles for a 3–4 cycle peak) was always adequate; the mobility window
was not, and the cause was **centring, not width**.

| configuration | apex on IM edge | weight in outer IM ring | edge cell > ½ apex |
|---|---|---|---|
| ±32 scans + library-line gate | 4.5% | 6.5% | 14.0% |
| ±64 scans + library-line gate | 4.8% | 4.8% | 13.8% |
| **±64 scans, observed mobility only** | **1.9%** | **1.6%** | **6.5%** |

The library gate is centred on the *predicted* mobility, the empirical window on the *observed* one. When a
prediction is off — 5% of PSMs are beyond 2σ — the intersection is lopsided: it collects empty slices on one side
and cuts the peak on the other. Worked example (precursor 8812954, 2+): apex at IM scan 381, the line puts its
library 1/K0 at 355, so the gate admitted 314–395 and collection stopped at 389 with the weight still at 78% of
apex. A PSM exists by this stage, so the observed mobility is strictly better information; the library gate now
applies only in the fragment index. The ion's mobility is highly reproducible (SD of the weighted centroid across
12 runs: median 1.1 scans) while the window centre (the best PSM's scan) varies by 3.3 — centring on the weighted
centroid would be better still.

## 6. Notes and open bugs

- `precursors_long.arrow` came out **0 bytes** on both 12-file searches while `precursors_long.tsv` of the same
  table is complete (1,005,001 rows for the 15-min set) and every other output wrote normally. Silent, so any
  Arrow consumer gets nothing. Not yet diagnosed; suspicion is a record-batch/offset limit on the largest long
  table we have produced.
- One replicate (15 min, Condition A REP2) is a mild outlier: it disagrees with all five other runs (median
  |log2| 0.22–0.32 vs 0.08–0.12 among clean pairs) and carries ~12% less summed precursor signal. Condition A's
  CV drops from 11.6% to 7.2% without it. Worth keeping in mind when reading single-condition CVs.
- Six replicates per condition is the honest basis; three left individual injections too much leverage.
