# Plan: Fix OOM vs In-Memory Scoring Discrepancy

## Problem

OOM mode produces ~3% fewer protein groups than in-memory mode, even with nearly all PSMs used for training.

## Root Cause Analysis

### Verified: Pairing Algorithm is IDENTICAL ✓
- Both use `assign_pair_ids()` with seed `PAIRING_RANDOM_SEED = 1844`
- Both use grouping keys: `[:irt_bin_idx, :cv_fold, :isotopes_captured]`
- Same algorithm, deterministic results

### Critical Difference 1: Q-Value Computation Timing

**In-Memory** (`percolatorSortOf.jl`):
- Line 374: Computes q-values **every iteration** after train predictions
- Line 875: Additional q-value computation in `update_mbr_features!`

**OOM Phase A** (`percolatorSortOf.jl`):
- Line 544: Computes q-values **only when `itr >= mbr_start_iter - 1`** (late iterations only)

**Impact**: The in-memory model sees q-value-informed feature updates throughout training; OOM model only sees them late.

### Critical Difference 2: Transfer Candidate Logic

**In-Memory** (lines 414-420):
```julia
qvals_prev = get_qvalues!(nonMBR_estimates, target)
pass_mask = (qvals_prev .<= max_q_value_lightgbm_rescore)
psms[!, :MBR_transfer_candidate] = .!pass_mask .& (psms.MBR_max_pair_prob .>= prob_thresh)
```
- Transfer = **NOT passing q-value** AND MBR score >= threshold

**OOM** (lines 1154-1159 in `update_mbr_probs_oom!`):
```julia
psms[!, :MBR_transfer_candidate] = (q > qval_thresh) .& (psms.MBR_max_pair_prob .>= trace_prob_thresh)
```
- Transfer = q-value **exceeds threshold** AND MBR score >= threshold

**These should be equivalent but need verification that thresholds match.**

---

## Proposed Fixes

### Fix 1: Align Q-Value Computation Timing

**File**: `src/utils/ML/percolatorSortOf.jl`

In the OOM Phase A training loop (around line 544), compute q-values every iteration like in-memory:

**Current** (line 544):
```julia
if itr >= mbr_start_iter - 1 && match_between_runs
```

**Change to**:
```julia
# Compute q-values every iteration (like in-memory)
qvals_test = get_qvalues!(prob_test, targets_test)
if itr >= mbr_start_iter - 1 && match_between_runs
```

Or better: ensure the q-value computation in Phase A matches the timing in the in-memory loop.

### Fix 2: Verify Transfer Candidate Threshold Consistency

**Files**:
- `percolatorSortOf.jl` lines 414-420 (in-memory)
- `percolatorSortOf.jl` lines 1154-1159 (`update_mbr_probs_oom!`)

Verify that:
1. `max_q_value_lightgbm_rescore` (in-memory) equals `qval_thresh` (OOM)
2. `prob_thresh` (in-memory) equals `trace_prob_thresh` (OOM)
3. The logical conditions are equivalent

### Fix 3: Match MBR Feature Computation Flow

The in-memory mode computes MBR features during training (line 389), while OOM computes them in Phase B.

Consider whether Phase A should also update MBR features during each iteration to match the in-memory learning dynamics.

---

## Implementation Steps

1. **Read the q-value computation code** in both modes to confirm the timing difference
2. **Read the transfer candidate code** to verify threshold variable names and values
3. **Add q-value computation to OOM Phase A** to match in-memory timing
4. **Verify threshold consistency** between modes
5. **Test** with near-full training data to confirm results converge

---

## Files to Modify

| File | Lines | Change |
|------|-------|--------|
| `percolatorSortOf.jl` | ~544 | Add q-value computation every iteration |
| `percolatorSortOf.jl` | 1154-1159 | Verify transfer candidate logic matches |

---

## Verification

1. Run OOM with 290mb threshold (triggers in-memory) - baseline
2. Run OOM with 275mb threshold (OOM mode with most PSMs)
3. Compare protein group counts - should be within <1% after fixes
4. Compare MBR transfer candidate counts between modes
