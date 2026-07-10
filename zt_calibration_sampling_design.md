# ZT calibration scan sampling — design notes

Branch `feat/sciex-zt-scanning`. How parameter-tuning / BitVec-calibration draw
MS2 scans, why ZT scanning breaks the default, and the fix.

## How Pioneer normally samples calibration scans

`get_ms2_scan_priority_order(spectra)` (src/structs/MassSpecData/FilteredMassSpecData.jl):

1. `compute_rt_bins` — split the run into `SCAN_PRIORITY_N_RT_BINS` (=15) RT bins.
2. `sort_scans_by_peak_density` — within each RT bin, sort scans by **TIC descending**.
3. `create_priority_order` — **round-robin** one scan from each RT bin: highest-TIC
   scan of bin1, bin2, …, bin15, then 2nd-highest of each, …

So the default is **both TIC-dependent (within-bin sort) AND RT-stratified
(round-robin across 15 bins)**. The tuning loop then walks this order accumulating
scans until it reaches the PSM target (`TUNING_MIN_SAMPLES`), keeping RT coverage
while preferring high-signal scans.

Used by BOTH ParameterTuningSearch.jl:583 and BitVecCalibration.jl:311.

## Why ZT scanning breaks it

Adjacent Q1 bins in the SAME cycle (scan indices ±1..±6) share physical ion
current (bin step S=1.022 m/z, effective window W≈7 m/z ⇒ a precursor spans ~±6
bins). So their **TICs are correlated** → the top-TIC scans within an RT bin are
clusters of adjacent bins. Top-TIC sampling therefore drew highly redundant
neighboring bins → IDs concentrated in RT "blips".

Measured redundancy (fraction of sampled scans within ±6 of another sampled scan):
- Collection N=26,721 (6.6% of MS2): **55.7%** (analytic 1−(1−f)¹² = 55.9%)
- Wide scout N=40,667 (10%): **71.7%** (analytic 71.9%)

Even a full random shuffle can't drop this much, because to reach the PSM floor at
ZT's low yield (~0.05 PSMs/scan) we must sample 6–10% of all scans, and at a ±6
correlation width that density alone forces neighbor collisions.

## Current mitigation (committed) and its cost

ZT-gated full `shuffle!` of scan_priority in ParameterTuningSearch.jl:
```julia
scan_priority = get_ms2_scan_priority_order(spectra)
if _zt
    shuffle!(MersenneTwister(1_800_017 + total_ms2), scan_priority)
end
```
Fixes the *clustering* (18,854→20,318 IDs) but **destroys both the TIC priority and
the RT stratification** — it's uniform random. For ZT's low yield we actually WANT
high-TIC (identifiable) scans, so the shuffle throws away useful signal.

## Proposed fix: gap-dedup (one scan per metascan neighborhood)

Greedy walk with a "blocked" BitVector — O(n·gap), ~50 KB, one pass:
```julia
# scan_priority : priority order (Int scan indices)
# n_scans       : total spectra; gap : exclusion half-width (metascan reach)
function dedup_metascan_order(scan_priority::Vector{Int}, n_scans::Int, gap::Int)
    blocked = falses(n_scans)
    out = similar(scan_priority, 0)
    @inbounds for s in scan_priority
        blocked[s] && continue          # already inside a claimed metascan
        push!(out, s)
        lo = max(1, s - gap); hi = min(n_scans, s + gap)
        for j in lo:hi
            blocked[j] = true           # claim ±gap neighborhood
        end
    end
    return out
end
```
Each accepted scan claims its ±gap neighborhood so no later scan inside it is
accepted. Greedy-in-priority-order keeps the highest-priority (highest-TIC) scan of
each neighborhood.

**Key refinement:** run this on the *original* TIC+RT-stratified order — NOT the
shuffled one. Then we get independent metascans WHILE preserving TIC priority and RT
stratification, i.e. it can **replace** the shuffle rather than stack on it.

**gap choice:**
- `gap=6` (=k): matches the "within ±6" redundancy criterion; metascans 7–12 apart
  still partially overlap.
- `gap=12` (=2k): metascans fully DISJOINT (`[s−6,s+6]` ∩ `[s′−6,s′+6]` = ∅) — true
  independence. Preferred. Pool = 405k/13 ≈ 31k independent metascans, ≫ PSM budget.

Cycles are ~497 scans apart, so a ±12 scan-index block never crosses a cycle
boundary — no explicit cycle bookkeeping needed (negligible edge effect).

Plug-in point (ParameterTuningSearch.jl, ZT-gated), replacing the shuffle:
```julia
scan_priority = get_ms2_scan_priority_order(spectra)
if _zt
    scan_priority = dedup_metascan_order(collect(Int, scan_priority), length(spectra), 12)
end
```
Start with tuning only; consider BitVecCalibration.jl:311 after measuring.

## ±3 tuning-tolerance experiment result (worse)

Widening the tuning precursor-m/z tolerance ±1→±3 (commit 6fc6c6a03):
- IDs q≤.01: **20,318 (±1) → 18,783 (±3)  = −7.6%**
- NCE precursors: 1260 → 363
- collection: 26,721 scans (±1) → **6,915 scans** (±3) to hit the same 1300 PSMs

Mechanism: wider tolerance → ~more candidate precursors per scan → PSM target hit
in ~4× fewer scans → calibration trained on a much narrower RT/precursor slice AND
on more ambiguous matches → worse mass-error/NCE calibration → fewer final IDs.
⇒ Revert to ±1. ±6 likely worse still (same ambiguity mechanism), though its
rationale (match full metascan width) differs — cheap to test if wanted.

---

# ZT chromatogram integration — IMPLEMENTED (Option A, theoretical triangle)

Problem: `IntegrateChromatogramsSearch` was NOT ZT-aware. `build_chromatograms`
(utils.jl:707) extracted per-(precursor,scan) deconv weights using the FITTED quad
(narrow ~1 m/z) → (a) dropped off-center metascan bins (kept ~1–3 of 13),
(b) emitted multiple points per cycle at ~identical RT → WH-smoothing / trapezoid
(`dt=t[i+1]−t[i]`) mishandled them. Quant would be both under-counted and mis-integrated.

Fix (env-gated on `PIONEER_ZT_METASCAN_K`, byte-identical when k=0):
1. **Wide extraction quad** (utils.jl build_chromatograms): for ZT use
   `SquareQuadModel(k·S)` (S = `_zt_bin_step`, median isolation width) so each
   precursor is enumerated across all 2k+1 bins → per-bin deconv weights. Mirrors
   main-search `qtm_deconv`.
2. **`collapse_chromatograms_to_metascans`** (utils.jl, new): group a precursor's
   rows by (precursor_idx, cycle) via `getCycleIdx`; representative bin = min
   |centerMz − prec_mz|; metascan intensity = Σ wⱼ·tⱼ / Σ tⱼ, tⱼ = 1 − |j|/(k+1)
   (bin offset j = scan_idx − center_scan, Σtⱼ = k+1). Emit one point/cycle at the
   center-scan rt/scan_idx. Keeping center scan_idx = the main-search-collapse apex
   scan preserves `find_nearest_scan` apex mapping.
3. **process_file!** (IntegrateChromatogramsSearch.jl): call the collapse after
   extract_chromatograms; for ZT set `precursor_fraction_transmitted = 1` (metascan
   already integrates the full window; skip `get_isotopes_captured!` so the fitted
   quad can't zero the point). Existing WH→2nd-deriv-boundary→baseline→trapezoid
   then integrates the clean one-point-per-cycle chromatogram UNCHANGED.

CombineTraces mode (config `trace_mode="combined"`) sorts by [:precursor_idx,:rt] —
schema sufficient; separate-trace isotope path not used. Precompiles clean.
Efficiency deferred: wide quad = 13× candidate columns per scan (same blowup as
main search); collapse does its own O(n log n) sort. Option B (inline accumulate)
later. **NOT yet run end-to-end** — full pipeline has never completed on ZT.
