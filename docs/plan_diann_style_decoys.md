# Plan: DIA-NN-Style Decoy Generation in Pioneer

## Motivation

Pioneer currently generates decoys by shuffling peptide sequences and then sending both targets and decoys through deep learning models (Chronologer for RT, Altimeter for spectra). This produces decoys with **independently predicted, realistic spectra and retention times** — challenging negative controls.

DIA-NN takes a different approach: it generates decoys by **mutating terminal amino acids** and shifting fragment m/z values accordingly, while **copying fragment intensities and retention times directly from the target**. The hypothesis is that DIA-NN's decoys are easier to distinguish from targets, leading to more IDs at the same FDR.

This plan adds DIA-NN-style decoy generation as an option in `BuildSpecLib`.

## How DIA-NN Generates Decoys (Default Method)

### 1. Sequence Mutation (NOT shuffle/reverse)

DIA-NN mutates exactly **2 amino acid positions** — one near each terminus — using a fixed lookup table (`diann.cpp:1370`):

```
Original:  G  A  V  L  I  F  M  P  W  S  C  T  Y  H  K  R  Q  E  N  D
Mutated:   L  L  L  V  V  L  L  L  L  T  S  S  S  S  L  N  D  Q  E  Q
```

**Position selection logic** (from `generate_decoy()`, `diann.cpp:3621-3636`):

Given amino acid sequence `aas` (0-indexed) with `m = length - 2`:

**N-terminal mutation** — picks ONE position with this priority:
1. Position 1 (2nd AA) — preferred, avoids actual N-terminal AA
2. Position min(2, m+1) — falls back inward
3. Position 0 (1st AA) — only if positions 1 and 2 are modified
4. Position 1 — final fallback

**C-terminal mutation** — picks ONE position with this priority:
1. Position m (2nd-to-last AA) — preferred
2. Position max(0, m-1) — falls back inward
3. Position m+1 (last AA) — only if m and m-1 are modified
4. Position m — final fallback

The `mod_state` check (`!mod_state[pos]`) ensures it avoids positions with modifications, stepping inward or outward to find an unmodified residue.

Each mutation produces a deterministic mass shift: `shift = mass(mutated_AA) - mass(original_AA)`.

### 2. Fragment m/z Shifting

Rather than regenerating fragments, DIA-NN applies a **uniform shift** to all fragments of each type (`diann.cpp:3638-3642`):

```
N_shift = mass(mutated_AA) - mass(original_AA)  at the chosen N-terminal position
C_shift = mass(mutated_AA) - mass(original_AA)  at the chosen C-terminal position

For each fragment:
    if Y-ion:  mz += C_shift / fragment_charge
    else:      mz += N_shift / fragment_charge    (b-ions and everything else)
```

All Y-ions shift by the **same amount** (C-terminal mutation divided by fragment charge), all b-ions shift by the **same amount** (N-terminal mutation divided by fragment charge). The fragment m/z values are not shuffled or randomized — they are deterministically shifted. The precursor m/z shifts by `(N_shift + C_shift) / precursor_charge`.

### 3. Fragment Intensities: Copied from Target

`decoy = target` (`diann.cpp:3593`) copies all fragment intensities. They are **never modified**. The decoy spectrum has the identical intensity pattern as the target, just with shifted m/z positions.

### 4. Retention Time: Copied from Target

Decoy RT = Target RT. The `decoy = target` copy at line 3593 copies iRT/sRT. No explicit RT modification anywhere in `generate_decoy()`.

### 5. DIA-NN also has `--reverse-decoys` mode

Reverses middle of sequence (keeps termini fixed), then **regenerates all fragments from scratch** using `generate_fragments()` on the reversed sequence (`diann.cpp:3643-3653`). This is a completely different approach from the default mutation method. We won't implement this variant initially.

### 6. Source Reference

All of this is in `diann.cpp`:
- Mutation table: line 1370
- `generate_decoy()`: lines 3590-3654
- `get_aas()`: lines 1546-1560 (strips sequence to bare amino acids)
- `get_mod_state()`: lines 1562-1575 (returns bool vector of modified positions)

## How Pioneer Currently Generates Decoys

1. **Shuffle sequence** (keep C-terminal AA fixed)
2. **Predict RT** via Chronologer for the shuffled sequence
3. **Predict spectra** via Altimeter for the shuffled sequence
4. Decoys get fully realistic, independently predicted spectra and RTs

## Why DIA-NN's Approach May Give More IDs

Pioneer's decoys are "too good" as negative controls:
- A shuffled sequence sent through Altimeter produces a **realistic spectrum** that looks like a real peptide
- Chronologer gives it a **plausible RT** for that sequence
- The ML scoring model (LightGBM) has to work harder to distinguish true PSMs from these realistic decoys
- This is conservative: fewer false positives, but also fewer true positives at the same FDR threshold

DIA-NN's decoys are structurally simpler:
- Same intensities as the target → the spectral "shape" is identical, just shifted in m/z
- Same RT → no RT discrimination possible between target and decoy
- The only difference is fragment m/z positions
- This means the decoy can only match if the m/z-shifted fragments happen to align with peaks in the spectrum — a much more constrained event

## Implementation Plan

### Step 1: Add `decoy_method` Option to BuildSpecLib Config

Add a new option `"diann_mutation"` alongside existing `"shuffle"` and `"reverse"`:

```json
"fasta_digest_params": {
    "add_decoys": true,
    "decoy_method": "diann_mutation"  // or "shuffle" (default), "reverse"
}
```

**File**: `src/Routines/BuildSpecLib/utils/buildParamDefaults.jl`

### Step 2: Implement the Amino Acid Mutation Table

**File**: `src/Routines/BuildSpecLib/fasta/fasta_utils.jl`

Add the DIA-NN mutation mapping as a constant:

```julia
const DIANN_MUTATION_TABLE = Dict{Char, Char}(
    'G' => 'L', 'A' => 'L', 'V' => 'L', 'L' => 'V', 'I' => 'V',
    'F' => 'L', 'M' => 'L', 'P' => 'L', 'W' => 'L', 'S' => 'T',
    'C' => 'S', 'T' => 'S', 'Y' => 'S', 'H' => 'S', 'K' => 'L',
    'R' => 'N', 'Q' => 'D', 'E' => 'Q', 'N' => 'E', 'D' => 'Q'
)
```

Add amino acid monoisotopic masses (if not already defined):

```julia
const AA_MASS = Dict{Char, Float64}(
    'G' => 57.02146, 'A' => 71.03711, 'V' => 99.06841, 'L' => 113.08406,
    'I' => 113.08406, 'F' => 147.06841, 'M' => 131.04049, 'P' => 97.05276,
    'W' => 186.07931, 'S' => 87.03203, 'C' => 103.00919, 'T' => 101.04768,
    'Y' => 163.06333, 'H' => 137.05891, 'K' => 128.09496, 'R' => 156.10111,
    'Q' => 128.05858, 'E' => 129.04259, 'N' => 114.04293, 'D' => 115.02694
)
```

### Step 3: Create `mutate_diann_style()` Function

**File**: `src/Routines/BuildSpecLib/fasta/fasta_utils.jl`

The position selection must match DIA-NN's priority logic exactly (1-indexed for Julia):

```julia
function mutate_diann_style(sequence::String, structural_mods)
    chars = collect(sequence)
    n = length(chars)
    m = n - 2  # DIA-NN's 'm' = aas.size() - 2 (0-indexed), = n - 2 in 1-indexed

    mod_at = position -> has_mod_at_position(structural_mods, position)

    # N-terminal position selection (DIA-NN priority: pos 1 → min(2,m+1) → 0 → 1, 0-indexed)
    # In 1-indexed Julia: pos 2 → min(3, m+2) → 1 → 2
    n_pos = if !mod_at(2)
        2
    elseif !mod_at(min(3, m + 2))
        min(3, m + 2)
    elseif !mod_at(1)
        1
    else
        2  # fallback
    end

    # C-terminal position selection (DIA-NN priority: m → max(0,m-1) → m+1 → m, 0-indexed)
    # In 1-indexed Julia: m+1 → max(1, m) → m+2 → m+1
    c_pos = if !mod_at(m + 1)
        m + 1
    elseif !mod_at(max(1, m))
        max(1, m)
    elseif !mod_at(m + 2)
        m + 2
    else
        m + 1  # fallback
    end

    # Compute mass shifts
    n_shift = 0.0
    if n_pos > 0 && haskey(DIANN_MUTATION_TABLE, chars[n_pos])
        orig = chars[n_pos]
        mutated = DIANN_MUTATION_TABLE[orig]
        n_shift = AA_MASS[mutated] - AA_MASS[orig]
        chars[n_pos] = mutated
    end

    c_shift = 0.0
    if c_pos > 0 && c_pos != n_pos && haskey(DIANN_MUTATION_TABLE, chars[c_pos])
        orig = chars[c_pos]
        mutated = DIANN_MUTATION_TABLE[orig]
        c_shift = AA_MASS[mutated] - AA_MASS[orig]
        chars[c_pos] = mutated
    end

    return String(chars), n_shift, c_shift
end
```

### Step 4: Implement DIA-NN Decoy Library Generation

This is the critical step. With the `"diann_mutation"` method, the workflow changes:

**Current flow (shuffle):**
```
targets → shuffle sequences → predict RT for ALL → predict spectra for ALL → build library
```

**New flow (diann_mutation):**
```
targets → predict RT for TARGETS ONLY → predict spectra for TARGETS ONLY →
mutate sequences → copy RT from target → shift fragment m/z → build library
```

This requires restructuring the pipeline so that decoys are created **after** spectral prediction, not before.

#### Option A: Post-prediction decoy injection

After `predict_fragments()` returns target spectra, create decoy entries by:

1. For each target precursor:
   - Mutate the sequence → get `n_shift`, `c_shift`
   - Create a decoy `FastaEntry` with the mutated sequence
   - Copy the target's RT as the decoy's RT
   - Copy the target's fragment list, shifting m/z:
     - Y-ions: `mz += c_shift / frag_charge`
     - B-ions: `mz += n_shift / frag_charge`
   - Keep fragment intensities unchanged
   - Compute new precursor m/z: `target_prec_mz + (n_shift + c_shift) / prec_charge`

2. Insert decoy entries into the library alongside targets

**File to modify**: `src/Routines/BuildSpecLib.jl` (main pipeline) and
`src/Routines/BuildSpecLib/build/build_poin_lib.jl` (library construction)

#### Option B: Pre-prediction with skip

Keep the current pipeline structure but:
1. Create mutated decoy sequences as usual in `add_decoy_sequences_grouped()`
2. Send only targets to Chronologer and Altimeter (skip decoys)
3. After prediction, copy target RT/spectra to paired decoys with m/z shifts

Option A is cleaner because it avoids creating unnecessary FastaEntry objects for decoys during the prediction phase.

### Step 5: Fragment m/z Shifting Implementation

**File**: New function in `src/Routines/BuildSpecLib/build/build_poin_lib.jl` or `fragments/`

```julia
function create_diann_decoy_fragments(
    target_fragments::Vector{DetailedFrag{Float32}},
    target_frag_range::UnitRange,
    n_shift::Float64,
    c_shift::Float64,
    decoy_prec_id::UInt32
)
    decoy_frags = similar(target_fragments[target_frag_range])
    for (i, frag) in enumerate(target_fragments[target_frag_range])
        mz_shift = if frag.is_y
            Float32(c_shift / frag.frag_charge)
        elseif frag.is_b
            Float32(n_shift / frag.frag_charge)
        else
            Float32(0)  # other ion types: no shift (or handle as needed)
        end

        decoy_frags[i] = DetailedFrag{Float32}(
            decoy_prec_id,
            frag.mz + mz_shift,
            frag.intensity,      # UNCHANGED from target
            frag.ion_type,
            frag.is_y,
            frag.is_b,
            frag.is_p,
            frag.is_isotope,
            frag.frag_charge,
            frag.ion_position,
            frag.prec_charge,
            frag.rank,
            frag.sulfur_count
        )
    end
    return decoy_frags
end
```

### Step 6: Target-Decoy Pairing

The existing pairing system (`pair_decoys.jl`) pairs by `(base_target_id, precursor_charge)`. For DIA-NN-style decoys, this pairing is even more natural since each decoy is derived directly from its target.

The pairing code should work as-is, but we need to ensure `base_target_id` is correctly set during the post-prediction decoy injection.

### Step 7: Precursor m/z Update

The decoy precursor m/z must also be shifted:

```julia
decoy_prec_mz = target_prec_mz + Float32((n_shift + c_shift) / prec_charge)
```

This happens when creating the decoy precursor entry in the precursors table.

## Key Design Decisions

### Where in the pipeline to inject decoys

**Recommended: After fragment prediction, before `buildPionLib()`**

The current pipeline in `BuildSpecLib()`:
1. `prepare_chronologer_input()` — creates FastaEntry objects, currently adds decoys here
2. `predict_retention_times()` — Chronologer
3. `parse_chronologer_output()` — loads RT, filters by m/z range
4. `predict_fragments()` — Altimeter/UniSpec via Koina
5. `parse_koina_fragments()` — loads fragments
6. `buildPionLib()` or `buildPionLibSpline()` — builds indexes

For `diann_mutation` mode:
- Step 1: Skip decoy creation (only create targets + entrapment)
- Steps 2-5: Only predict for targets
- **New Step 5.5**: Create decoy entries by mutating targets, copying RT, shifting fragment m/z
- Step 6: Build library with both targets and decoys

### What to do about entrapment sequences

Entrapment sequences are separate from decoys and should continue to be predicted independently (they serve a different purpose — estimating entrapment FDR). Only decoy generation changes.

### How to handle modifications

The DIA-NN mutation only changes unmodified terminal positions. If both terminal AAs have modifications, it skips to the next position inward. Our implementation should do the same — `has_mod_at_position()` checks whether a `structural_mods` entry exists at that residue index.

### What about p-ions, isotope fragments, internal fragments?

For non-b/y ions:
- **p-ions (precursor isotopes)**: Shift by `(n_shift + c_shift) / prec_charge` (same as precursor)
- **Internal fragments**: Not currently used; skip
- **Isotope fragments**: Shift same as parent b/y ion

## Files to Modify

| File | Change |
|------|--------|
| `src/Routines/BuildSpecLib/utils/buildParamDefaults.jl` | Add `"diann_mutation"` as valid `decoy_method` |
| `src/Routines/BuildSpecLib/fasta/fasta_utils.jl` | Add mutation table, `mutate_diann_style()` |
| `src/Routines/BuildSpecLib.jl` | Conditional pipeline: skip decoy prediction when `diann_mutation` |
| `src/Routines/BuildSpecLib/build/build_poin_lib.jl` | Post-prediction decoy injection, fragment m/z shifting |
| `src/Routines/BuildSpecLib/chronologer/pair_decoys.jl` | May need minor adjustments for new pairing flow |

## Testing Plan

1. **Build library with `diann_mutation` decoys**: Verify decoy count matches target count, m/z shifts are correct
2. **Spot-check fragment m/z**: For a known peptide, verify Y-ion shifts match expected C-terminal mass difference
3. **Verify RT copying**: Decoy RT should exactly match paired target RT
4. **SearchDIA comparison on OlsenAstral**:
   - Library A: Current (shuffle + predict)
   - Library B: DIA-NN mutation (shift + copy)
   - Compare unique targets at 1% FDR
5. **SearchDIA comparison on OlsenEclipse**: Same comparison on wider-window data

## Expected Outcome

If the hypothesis is correct:
- DIA-NN-style decoys → higher scores for true targets relative to decoys → more IDs at same FDR
- The effect may be larger on Eclipse (wide isolation, more chimeric spectra) than Astral (narrow isolation)
- The FDR estimate may be slightly optimistic compared to shuffle decoys (a tradeoff to be aware of)

## Risk Assessment

**FDR calibration concern**: If DIA-NN-style decoys are "too easy," the FDR estimate may be anti-conservative (actual FDR higher than reported). This is the fundamental tradeoff — easier decoys give more IDs but potentially less accurate FDR control. The entrapment analysis can help validate this: if entrapment FDR matches reported FDR, the decoys are well-calibrated.

**Mitigation**: Run both decoy methods on the same data and compare entrapment FDR vs reported FDR to assess calibration quality.
