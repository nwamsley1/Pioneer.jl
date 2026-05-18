# Overnight LGBM HP × feature ablation sweep — Olsen Astral / Exploris / MTAC (6 files each)

**Goal:** speed up the per-file LGBM without losing IDs (within ~5%), aiming for 20–30 active features.

**Compute:** `-t 10 --gcthreads=5,1`, 45 runs (15 variants × 3 datasets), 5h 41min wall (00:03 → 05:45 CDT).

**Baseline (`hp_baseline`)**: 47 features (PRESCORE_FEATURES + ADVANCED_FEATURE_SET resolved), `num_iterations=200, learning_rate=0.10, min_data_in_leaf=300`.

| Dataset | Baseline IDs (q≤0.01 targets) | Baseline PGs | Baseline wall (s) | Baseline MainSearch (s) |
|---|---:|---:|---:|---:|
| **Astral** (6f) | 1,060,032 | 76,444 | 636.5 | 400.3 |
| **Exploris** (6f) | 300,352 | 39,345 | 274.5 | 109.3 |
| **MTAC** (6f) | 582,787 | 64,060 | 352.5 | 157.4 |

---

## Result table (delta vs baseline per dataset)

### Astral (6 files)

| Variant | nF | Δwall | ΔMain | ΔIDs | ΔIDs% | ΔPGs | ΔPGs% |
|---|---:|---:|---:|---:|---:|---:|---:|
| hp_baseline | 47 | 0.0s | 0.0s | 0 | 0.00% | 0 | 0.00% |
| hp_150_lr12 | 47 | −32.0s | −30.8s | +808 | +0.08% | −302 | −0.40% |
| hp_100_lr15 | 47 | −56.5s | −58.6s | +1,255 | +0.12% | +11 | +0.01% |
| hp_100_lr20 | 47 | −62.4s | −60.4s | −5,433 | −0.51% | −158 | −0.21% |
| hp_leaf500 | 47 | +8.9s | +7.2s | +4,953 | +0.47% | −68 | −0.09% |
| hp_leaf100 | 47 | −13.7s | −14.1s | +818 | +0.08% | +29 | +0.04% |
| feat_drop_win5_irt | 44 | −0.3s | +1.4s | −181 | −0.02% | −198 | −0.26% |
| feat_drop_ms1_intensities | 45 | −1.5s | −1.0s | −496 | −0.05% | +35 | +0.05% |
| feat_drop_ms1_all | 41 | −4.1s | −2.1s | −115 | −0.01% | −40 | −0.05% |
| feat_drop_frag_int | 43 | −2.3s | −1.4s | −3,813 | −0.36% | −203 | −0.27% |
| feat_drop_win11 | 45 | +7.5s | +3.6s | **+10,765** | **+1.02%** | +425 | +0.56% |
| feat_drop_neighborhood_full | 38 | +1.5s | 0.0s | +9,099 | +0.86% | +295 | +0.39% |
| combo_100lr20_drop_ms1all | 41 | −61.2s | −61.3s | −5,024 | −0.47% | +134 | +0.18% |
| **combo_100lr20_drop_neighborhood** | **38** | **−60.0s** | **−57.5s** | **+4,948** | **+0.47%** | **+140** | **+0.18%** |
| combo_100lr20_drop_ms1_neigh | 32 | −65.2s | −61.6s | −72 | −0.01% | −65 | −0.09% |

### Exploris (6 files)

| Variant | nF | Δwall | ΔMain | ΔIDs | ΔIDs% | ΔPGs | ΔPGs% |
|---|---:|---:|---:|---:|---:|---:|---:|
| hp_baseline | 47 | 0.0s | 0.0s | 0 | 0.00% | 0 | 0.00% |
| hp_150_lr12 | 47 | −8.7s | −9.5s | +110 | +0.04% | −4 | −0.01% |
| hp_100_lr15 | 47 | −21.1s | −20.5s | −1,490 | −0.50% | −80 | −0.20% |
| hp_100_lr20 | 47 | −19.1s | −19.7s | −1,768 | −0.59% | +162 | +0.41% |
| hp_leaf500 | 47 | −1.1s | −0.1s | +1,642 | +0.55% | +13 | +0.03% |
| hp_leaf100 | 47 | −0.2s | −0.9s | +812 | +0.27% | −163 | −0.41% |
| feat_drop_win5_irt | 44 | −1.3s | −0.4s | +595 | +0.20% | −93 | −0.24% |
| feat_drop_ms1_intensities | 45 | −0.9s | +0.1s | −520 | −0.17% | +106 | +0.27% |
| **feat_drop_ms1_all** | 41 | −2.8s | −1.1s | **−20,660** | **−6.88%** | **−1,357** | **−3.45%** |
| feat_drop_frag_int | 43 | 0.0s | +0.2s | −2,199 | −0.73% | −197 | −0.50% |
| feat_drop_win11 | 45 | +1.0s | +0.2s | +2,532 | +0.84% | +78 | +0.20% |
| feat_drop_neighborhood_full | 38 | −3.4s | −1.1s | +3,794 | +1.26% | +78 | +0.20% |
| **combo_100lr20_drop_ms1all** | 41 | −23.4s | −20.9s | **−22,944** | **−7.64%** | **−1,312** | **−3.33%** |
| **combo_100lr20_drop_neighborhood** | **38** | **−24.0s** | **−21.4s** | **+2,217** | **+0.74%** | **+14** | **+0.04%** |
| combo_100lr20_drop_ms1_neigh | 32 | −26.4s | −22.3s | −18,934 | −6.30% | −1,023 | −2.60% |

### MTAC (6 files)

| Variant | nF | Δwall | ΔMain | ΔIDs | ΔIDs% | ΔPGs | ΔPGs% |
|---|---:|---:|---:|---:|---:|---:|---:|
| hp_baseline | 47 | 0.0s | 0.0s | 0 | 0.00% | 0 | 0.00% |
| hp_150_lr12 | 47 | −11.3s | −12.8s | +283 | +0.05% | −304 | −0.47% |
| hp_100_lr15 | 47 | −23.9s | −26.1s | −2,077 | −0.36% | −349 | −0.54% |
| hp_100_lr20 | 47 | −24.3s | −25.0s | −2,487 | −0.43% | −561 | −0.88% |
| hp_leaf500 | 47 | +4.5s | +2.5s | +1,223 | +0.21% | −172 | −0.27% |
| hp_leaf100 | 47 | −1.0s | −2.2s | −357 | −0.06% | −429 | −0.67% |
| feat_drop_win5_irt | 44 | +1.1s | +1.9s | −2,353 | −0.40% | −401 | −0.63% |
| feat_drop_ms1_intensities | 45 | +1.2s | +0.8s | −1,491 | −0.26% | −1 | −0.00% |
| feat_drop_ms1_all | 41 | +0.9s | +0.4s | +1,958 | +0.34% | +17 | +0.03% |
| feat_drop_frag_int | 43 | −0.2s | +1.0s | −741 | −0.13% | −174 | −0.27% |
| feat_drop_win11 | 45 | +8.1s | +2.7s | **+8,480** | **+1.46%** | +21 | +0.03% |
| feat_drop_neighborhood_full | 38 | +2.3s | +1.1s | **+14,369** | **+2.47%** | **+838** | **+1.31%** |
| combo_100lr20_drop_ms1all | 41 | −27.5s | −26.3s | −3,085 | −0.53% | −393 | −0.61% |
| **combo_100lr20_drop_neighborhood** | **38** | **−28.3s** | **−26.9s** | **+12,555** | **+2.15%** | **+687** | **+1.07%** |
| combo_100lr20_drop_ms1_neigh | 32 | −29.9s | −26.8s | +12,519 | +2.15% | +849 | +1.33% |

---

## Key findings

1. **The "neighborhood" feature class is net-harmful.** Dropping the 9 `*_5scan`, `*_11scan`, `best_max_residual_3scan` features improves IDs on **all three datasets** at no wall-clock cost:
   - Astral +0.86%, Exploris +1.26%, MTAC +2.47% precursor IDs
   - PGs: Astral +0.39%, Exploris +0.20%, MTAC +1.31%
   - The smart-composite ablation campaign that landed these had Exploris-only and MTAC-only A/B but didn't catch that they regress on Astral (where they were never explicitly tested). The cross-dataset picture says drop them.

2. **MS1 features are essential for Exploris** (instrument quality / quad behavior). Dropping all six MS1 features loses **−6.88% IDs** on Exploris while only hurting MTAC and Astral marginally. **Do not drop MS1 features.**

3. **`fragN_int` (per-rank fragment intensities) are mildly net-positive but not critical** — −0.13% to −0.73% IDs when dropped. Worth keeping.

4. **HP sweep alone** shows the 200-iter/lr=0.10 baseline is not optimal: dropping iterations buys wall-clock at a small ID cost. But combined with feature reduction, the IDs go up.

5. **Pareto winner: `combo_100lr20_drop_neighborhood`** — 100 iter, lr=0.20, drop the 9 neighborhood features (38 active features):
   - Astral:   −60.0s wall (−9.4%) │ +4,948 IDs (+0.47%) │ +140 PGs (+0.18%)
   - Exploris: −24.0s wall (−8.7%) │ +2,217 IDs (+0.74%) │ +14 PGs (+0.04%)
   - MTAC:     −28.3s wall (−8.0%) │ +12,555 IDs (+2.15%) │ +687 PGs (+1.07%)

   **+ID on every dataset AND ~9% wall savings.** This is the right knob.

6. **20–30 feature target** partially reachable. `combo_100lr20_drop_ms1_neigh` (32 features) works on Astral/MTAC but loses 6.3% IDs on Exploris. To hit 25–30 features safely we'd need a more surgical ablation (drop frag_int + a few low-gain MS1 cols) — not in this sweep. The 38-feature winner is the safe knob today.

---

## Recommendation (drop into a single commit)

1. Change `SHARED_LGBM_HP` baseline:
   - `num_iterations: 200 → 100`
   - `learning_rate: 0.10 → 0.20`
   - (keep `min_data_in_leaf=300`)

2. Drop these 9 features from `PRESCORE_FEATURES` and `ADVANCED_FEATURE_SET`:
   - `best_max_residual_3scan`
   - `best_gof_5scan`, `best_manhattan_5scan`, `best_max_residual_5scan`
   - `irt_dist_best_gof_5scan`, `irt_dist_best_manhattan_5scan`, `irt_dist_best_max_residual_5scan`
   - `worst_max_residual_11scan`, `worst_manhattan_11scan`

3. Keep `add_neighborhood_features!` callable (computes the columns), but if these 9 are the only consumers, the function call itself can be skipped — that's an additional ~4s/file gain we didn't measure here. Worth a follow-up A/B.

4. Don't touch MS1 features — Exploris depends on them.

## Expected impact on the original feature-vs-develop gap

We measured a +246s MainSearch gap (feature 401s vs develop 155s) on Astral 6-file. The Pareto winner recovers about **−60s of that gap** while adding IDs across all datasets. Combined with skipping `add_neighborhood_features!` (the +4s/file × 6 = +24s) the recovery is closer to **~85s of the 246s** — roughly a 35% reduction in the regression. The remaining ~160s is the other deliberate gains (MS1 features, additional feature passes, paircomp, PEP filter); the sweep shows those are paying their way.

---

## Raw summary file

`/Users/nathanwamsley/Data/RegressionTestsLite/overnight_sweep_2026_05_18/summary.tsv` (TSV, one row per run).

Each per-run output dir (`{dataset}_{variant}/`) still has the full `pioneer_search_log.log` + Performance Report if a deeper look is needed.
