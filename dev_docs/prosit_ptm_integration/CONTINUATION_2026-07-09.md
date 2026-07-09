# Continuation memo — 2026-07-09

**Purpose:** enable a fresh session to pick up the site-determining-ion
null-model work exactly where this one left off. Read this + the two PDFs
committed alongside it and you have the full picture.

## Where we are (one paragraph)

Started from the AScore paper (Beausoleil *et al.*, *Nat. Biotechnol.* 2006).
Its localization null is a cumulative binomial with $p = i/100$ — implicit
uniform-fragmentation assumption. Question: is that assumption defensible in
a Prosit-Pioneer setting? Built an offline exploration
(`prototype/explore_matched_counts.jl`) that, for the 12-file Coon dilution
search, enumerates every pair (target PSM, competing target positional
localization) and records per pair how many b/y site-determining ions match
the observed scan on the target side and on the competitor side, with two
definitions of "available" (all theoretical b/y; Prosit library-retained
subset). Ran it on a 10k-PSM sample → 156,715 pair rows. Findings + framing
in `matched_ions_empirical_exploration.pdf`. User correctly pushed back that
the reported "competitor rate" is a **mixture** — coincidental density matches
plus real co-eluting-isomer signal — not a clean null. Doc was corrected. A
design proposal for a defensible null model followed
(`localization_null_model_design.pdf`), initially with `p_chance` computed
analytically from per-scan peak density. User then proposed replacing
analytical `p_chance` with a **sequence-shuffle-based `p_chance`** matched on
the same scan. That's the current recommendation. **No null-model code has
been written yet.**

## Key files added this session

- `dev_docs/prosit_ptm_integration/prototype/explore_matched_counts.jl` —
  target-vs-target site-det ion match counter. Reads `precursors_long.arrow`
  from a search result dir; joins to a clean library by
  `(seq, sorted phospho positions, sorted CAM positions, charge)`; opens
  per-file MS arrow to match at ±20 ppm. Writes per-pair Arrow and prints
  stratified summaries.
- `dev_docs/prosit_ptm_integration/matched_ions_empirical_exploration.md/.pdf` —
  the exploration write-up (corrected re: co-elution mixture).
- `dev_docs/prosit_ptm_integration/localization_null_model_design.md/.pdf` —
  Option A (sequence decoys) vs Option B (`p_real` regression + `p_chance`
  density); LLR construction; where sequence decoys still calibrate. Note
  that this doc's `p_chance` section describes the ANALYTICAL null; the
  active recommendation is now the SHUFFLE null (see "Open decisions" below).

## Data locations

- **Search results:** `/Users/n.t.wamsley/prosit_phospho_test/dilution_curve2/combined/results/`
    - `precursors_long.arrow` — 969,006 target-passing PSMs across 12 files
    - `explore_matched_counts_v1.arrow` — 156,715 per-pair rows produced today
- **Analysis library (clean, no positional-isomer decoys):**
  `/Users/n.t.wamsley/prosit_phospho_test/yeast_plus_standards_1mc.poin/`
    - 2,653,366 targets, 2,418,808 with $C > 1$, no `is_loc_decoy`
- **Search-time library (has shadow decoys, only for `precursor_idx` →
  `(seq, mods, chg)` lookup):**
  `/Users/n.t.wamsley/prosit_phospho_test/yps_decoy2.poin/`
- **Per-file MS arrow:** `/Users/n.t.wamsley/prosit_phospho_test/dilution_curve2/combined/msdata/`
  (12 files; short-form file names in results, prefixed on disk — script
  auto-resolves)

## Key empirical findings to keep in mind

- On the theoretical distinguishing set: target-side match rate 0.28–0.36,
  competitor-side 0.14–0.23.
- On the Prosit library-retained subset: **target 0.69–0.78, competitor
  0.25–0.29**. Much cleaner separation; fewer trials per pair (median 4–8
  vs 25+); this is the right subset to build the probability model on.
- Peak-count covariate is large: competitor rate spans **3% (sparse
  spectra) to 20% (dense)** — a single global `p` cannot represent both.
- Signal:null gap narrows with $C$ (2.1× at $C=2$; 1.4× at $C \ge 21$).
- Fraction distribution is right-skewed, mode 0.2–0.3 — not a single-`p`
  binomial.

## Framing correction we made mid-session

The "competitor rate" reported above is NOT a clean null. Positional isomers
of the same base peptide co-elute in DIA (same mass, near-identical RT), so
their b/y m/z lines are frequently co-fragmented in the same isolation
window. The competitor rate mixes (a) coincidental density matches with (b)
real co-eluting-isomer signal. It is an **upper bound** on the true null
$p_0$. Any null-model construction has to either (i) exclude co-eluting
pairs, or (ii) model the co-elution explicitly. The LLR construction in the
design doc handles (ii) automatically because it ratios likelihoods — in a
true co-elution scan the LLR goes to zero, correctly signalling ambiguity.

## Open decisions (what to think through before writing more code)

1. **Which shuffle scheme for `p_chance`?**
    - **Residue shuffle:** permute interior residues (preserve terminal K/R),
      recompute b/y ladders (including PTMs), match on the same scan.
      Broader null, robust.
    - **Phospho-position shuffle:** keep target residues, re-draw phospho
      positions from STY (or from non-STY for a forced wrong-loc). Tighter
      match to the actual wrong-loc scenario because the m/z shifts are
      exactly ±(phospho mass / z) at the same b/y indices a real competitor
      would produce.
    - Recommendation from the last exchange: **run both on a modest sample
      first**; if per-scan `p_chance` values agree, use the simpler residue
      shuffle; if they disagree, use phospho-position shuffle as the honest
      answer for our specific problem.
    - The phospho-position shuffle is essentially the shadow-site decoy
      construction already in `src/Routines/BuildSpecLib/fragments/isomer_decoys.jl`.
      Reuse that generator (analysis-only, not written to library).

2. **Whether to also fit `p_real` regression.**
    - The LLR needs a `p_real` per ion (probability of match given the
      target's loc is right).
    - Options: (a) fit a regression on all confident-target b/y ions
      (features listed in the design PDF §"`p_real`"); (b) use empirical
      per-Prosit-intensity-decile match rate as a coarse `p_real`; (c) start
      even simpler — treat `p_real` as a constant on the library-retained
      subset (approx 0.75 from the empirical exploration) and only vary
      `p_chance` per scan.
    - (c) is a legitimate v0 that lets us test the LLR framework before
      investing in fitting.

3. **How to treat the mixture.**
    - Even with a good `p_chance`, some pairs will show competitor match
      rate $\gg$ `p_chance`. The LLR handles this by going to zero
      (ambiguous); but for **reporting** we may want to expose an
      "ambiguity flag" per PSM rather than a single localization score.
    - Aligns with the meta-peptidoform pivot (2026-07-08 in memory): report
      ambiguity sets rather than pick a winner. LLR + ambiguity flag is a
      natural marriage.

## Immediate next steps (roughly in order)

1. **Write shuffle-based `p_chance` estimator.** Adds a per-pair column
   `p_chance_shuffle` to the existing Arrow. Implementation:
    - For each pair (target `t`, competitor `c`) at scan `s`:
        - Generate 5–20 phospho-position shuffles of `t`'s peptidoform
          (retain STY constraint; skip shuffles that reproduce `t`'s or
          `c`'s actual phospho set).
        - For each shuffle, compute the site-determining ion set (relative
          to `t`), map to m/z, match against scan `s`, and record match
          fraction.
        - Average across shuffles → `p_chance_shuffle`.
    - Also compute analytical `p_chance` (per-scan peak-density fraction)
      as a sanity check; expect them to differ mostly in dense spectra.
2. **Split existing pairs into "chance regime" vs "co-elution regime"** by
   (observed competitor rate − `p_chance_shuffle`). Report the split
   stratified by peak count and by $C$.
3. **Decide whether to fit `p_real`** based on whether the split reveals a
   clean chance regime worth modeling.
4. **Build LLR scorer,** initially with constant-`p_real` (v0), then with
   fitted-`p_real` (v1) if warranted.
5. **Validate against sequence decoys as an out-of-sample calibration check
   only** (per the null-model design PDF): decoy site-determining ion match
   rate on their own scans should track `p_chance`, not `p_real`.

## Related prior work still relevant

- `dev_docs/prosit_ptm_integration/prototype/extract_sds_features.jl` — the
  S-ratio (intensity-based) feature the LGBM ablation validated. LLR is
  the successor of this feature; don't rebuild what's here.
- `dev_docs/prosit_ptm_integration/prototype/loc_score_core.jl` +
  `test_loc_score_core.jl` — the `S(A,B)` scoring core + tests. Same
  distinguishing-ion machinery; reuse if we move LLR into a package.
- `dev_docs/prosit_ptm_integration/prototype/PROTOTYPE_STATUS.md` — S-score
  validation status prior to this exploration. Notes the phase-2 gate
  (standards calibration) that has not been done yet; LLR would take the
  same phase-2 gate.
- `src/Routines/BuildSpecLib/fragments/isomer_decoys.jl` — shadow-site decoy
  generator. Reuse for the phospho-position shuffle.

## Running the exploration script (for reference)

```
julia --project=Pioneer.jl dev_docs/prosit_ptm_integration/prototype/explore_matched_counts.jl \
  /Users/n.t.wamsley/prosit_phospho_test/dilution_curve2/combined/results \
  /Users/n.t.wamsley/prosit_phospho_test/yeast_plus_standards_1mc.poin \
  /Users/n.t.wamsley/prosit_phospho_test/yps_decoy2.poin \
  /Users/n.t.wamsley/prosit_phospho_test/dilution_curve2/combined/msdata \
  10000 \
  /Users/n.t.wamsley/prosit_phospho_test/dilution_curve2/combined/results/explore_matched_counts_v1
```

Args: `<results_dir> <analysis_lib> <search_lib> <arrow_dir> [n_sample]
[out_prefix]`.

## Things to explicitly NOT do

- Do not use the observed competitor rate as a null $p_0$ — it's a mixture
  upper bound; the framing correction above is load-bearing.
- Do not build sequence-decoy-based null estimation as the primary path —
  the conditioning bias (decoy PSMs enrich for scans with accidental
  support) is real. Sequence decoys are for calibration only.
- Do not re-search the data with the clean library before doing (1)–(3)
  above; the target-space join is exact and re-searching gains nothing.
