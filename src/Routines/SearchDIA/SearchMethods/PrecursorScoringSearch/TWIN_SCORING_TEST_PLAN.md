# Twin-Scoring Test Plan

Goal: measure whether two-round twin_score/delta_irt (cross-run consistency) helps or
harms, using ground-truth where we have it. All A/Bs are on the **develop lineage** so the
only variable is the two-round feature (develop's native FDR is already the simultaneous
`(:global_qval AND :qval)` filter, so twin_score expresses without any env hacks).

**Branches** (see memory `reference_pyr_twin_branches`):
- **baseline** = `develop`, config MBR off (`global.match_between_runs=false`)
- **twin** = `feat/twin-scoring-on-develop`, MBR off + env `TWO_ROUND=1`

## Test 1 — KEAP1 (ground-truth FALSE-TRANSFER test) — DO FIRST

KEAP1 (UniProt Q14145) is knocked out: baseline shows **16 KEAP1 precursors in WT (3 reps),
0 in KO (3 reps)**; KEAP1 protein group present in WT, absent in KO. So KO is a known
true-negative for KEAP1.

- **Question:** does twin-scoring **transfer KEAP1 WT→KO**? A KO run's nearest run is a WT run
  (KO/WT are compositionally near-identical — presence Jaccard sep 0.009), so twin_score for a
  low-scoring KEAP1 candidate in a KO run would be HIGH (from WT, where KEAP1 is confident) →
  round-2 could boost it → KEAP1 appears in KO = false transfer.
- **Metric:** count KEAP1 (Q14145) precursors + protein group passing in each KO run, baseline
  vs twin. Baseline = 0. Any KO KEAP1 IDs under twin = a demonstrable false transfer.
- **Why it matters:** this is the adversarial case for single-twin — conditions that are
  compositionally indistinguishable but have a few real differences. Tests whether twin-scoring
  is safe when clustering CANNOT separate the conditions (it can't here: silhouette 0.044).
- Data: `RegressionTestsLite/KEAP1KO/ms_data` (6 H292 files), human lib, both runs already have
  the develop baseline (`keap1_run`). Just need the twin run.

## Test 2 — EWZ (designed BENEFIT test)

The false-yeast dataset: twin-scoring's intended job is to cut false cross-species transfers.

- **Question:** on develop's native AND FDR, does twin-scoring reduce false yeast at similar IDs?
- **Metric:** total target instances + unique precursors, and false-yeast instances/rate
  (species==YEAST in GO113 human-only runs), baseline vs twin. Prior (v2 branch, env AND):
  twin recovered +94k IDs at +327 fy. Re-measure cleanly on develop.
- Data: `EWZ_local/ms_data` (40 files), `yeast_human.poin`.

## Test 3 — Pyr (COMPLETENESS test, no false set)

Homogeneous patient cohort (clustering ≈1). No labeled false set, so measure ID/completeness.

- **Question:** does twin-scoring help or hurt on near-homogeneous single-species data?
- **Metric:** total IDs, unique precursors, data completeness (instances / unique / n_runs),
  baseline vs twin. Watch for spurious inflation.
- Data: `Pyr/ms_data` (40 files, downloading), human lib.

## Test 4 — CLUSTER-CONSENSUS variant (after 1-3)

Replace single-twin with cluster-consensus (SVD+KMeans on stringent-presence → within-cluster
LOO mean s1; silhouette K=1 floor ~0.2; singleton = self-twin). Re-run Tests 1-3.

- KEAP1: does clustering prevent the WT→KO transfer? (Likely NOT — clustering can't separate
  KO/WT at silhouette 0.044, so it merges them → still references across. This is the honest
  limit: a 16/130k composition difference is below any clustering's resolution. See `CLUSTER_
  CONSENSUS_CLUSTERING.md` §KEAP1.)
- EWZ: does cluster-consensus beat single-twin? (offline suggested ~2x, but leakage-corrected.)
- **PREREQUISITE — clustering record must be intact before starting Test 4:**
  `CLUSTER_CONSENSUS_CLUSTERING.md` (recipe, method sweep, Pyr/KEAP1 cases), memory
  `project_two_round_scoring` (findings) + `reference_pyr_twin_branches` (branches), and plots
  in `RegressionTestsLite/DIANN_vs_Pioneer_PEP/` (svd_kmeans_silhouette, presence_jaccard_
  clustering, cluster_recovery_viz, cluster_quality_viz, olsen3p_clustering, pyr_subset_clustering).

## Run recipe (each A/B)

```bash
# baseline
git checkout develop
julia --project=. --threads 10 --gcthreads 5,1 -e 'using Pioneer; SearchDIA("<cfg>.json")'
# twin
git checkout feat/twin-scoring-on-develop
TWO_ROUND=1 julia --project=. --threads 10 --gcthreads 5,1 -e 'using Pioneer; SearchDIA("<cfg>.json")'
```
Config: MBR off, delete_temp=false (keep main_search dump for diagnostics), q≤1%.
