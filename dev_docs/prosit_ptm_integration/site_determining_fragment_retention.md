---
title: "Retaining Site-Determining Fragments in the Library"
subtitle: "Worked on a real peptide: why 'top unique fragment' fails, and boundary-crossing coverage"
author: "Pioneer.jl -- prepared for N. Wamsley"
date: "2026-07-03"
geometry: margin=2.2cm
fontsize: 11pt
colorlinks: true
linkcolor: RoyalBlue
---

# 1. The problem

For localization (and for the isomer-decoy scheme) to work, the library must physically contain the
**site-determining ions** -- the fragments whose mass says *where* the phosphate sits. Today we keep
the **top-10 fragments by intensity per peptidoform**, and those are often all *shared* ions that
carry the phosphate no matter which isomer is true, so two isomers end up with identical stored
fragments and become indistinguishable. This note works one concrete peptide to show exactly which
ions matter and gives an algorithm to retain them.

# 2. The example peptide

$$\textbf{VLSAGSPESIK}\qquad (\text{length } L = 11)$$

| pos | 1 | 2 | **3** | 4 | 5 | **6** | 7 | 8 | **9** | 10 | 11 |
|-----|---|---|---|---|---|---|---|---|---|----|----|
| aa  | V | L | **S** | A | G | **S** | P | E | **S** | I  | K  |

The phospho-acceptor residues (S/T/Y) are at **positions 3, 6, 9** (all serine). With one phosphate,
there are three positional isomers: **mod@3, mod@6, mod@9** -- same sequence, same precursor m/z,
differing only in which serine carries the +79.97 Da.

## Which fragment carries the phosphate?

A fragment "carries" the phosphate (its m/z is +79.97) iff it **contains the modified residue**:

- b-ion `b_i` = residues `1..i` -> contains site `k` iff `k <= i`
- y-ion `y_j` = residues `(L-j+1)..L` = `(12-j)..11` -> contains site `k` iff `k >= 12-j`

Below, **P** = that ion carries the phosphate for that isomer; **.** = it does not.

**b-ions** (N-terminal):

| ion | residues | mod@3 | mod@6 | mod@9 |
|-----|----------|:---:|:---:|:---:|
| b2  | VL       | .   | .   | .   |
| **b3**  | VLS  | **P** | . | . |
| **b4**  | VLSA | **P** | . | . |
| **b5**  | VLSAG| **P** | . | . |
| **b6**  | VLSAGS | P | **P** | . |
| **b7**  | VLSAGSP| P | **P** | . |
| **b8**  | VLSAGSPE | P | **P** | . |
| b9  | ...S9    | P   | P   | P   |
| b10 | ...      | P   | P   | P   |

**y-ions** (C-terminal):

| ion | residues | mod@3 | mod@6 | mod@9 |
|-----|----------|:---:|:---:|:---:|
| y2  | IK       | .   | .   | .   |
| **y3**  | SIK  | . | . | **P** |
| **y4**  | ESIK | . | . | **P** |
| **y5**  | PESIK| . | . | **P** |
| **y6**  | SPESIK | . | **P** | P |
| **y7**  | GSPESIK| . | **P** | P |
| **y8**  | AGSPESIK | . | **P** | P |
| y9  | SAGSPESIK| P   | P   | P   |
| y10 | ...      | P   | P   | P   |

`b9,b10,y9,y10` carry P for **all three** isomers (they span every site) -- useless for localization.
`b2,y2` carry P for none. The action is in the middle rows.

# 3. Why "keep each peptidoform's top *unique* fragment" fails

Read the columns. A fragment is **unique** to an isomer if only that isomer's column has a P (or
only it has a `.`) in that row.

- **mod@3** is easy: `b3,b4,b5` are **P for mod@3 only** -> unique. (Also `y6,y7,y8` are `.` for
  mod@3 only.) So mod@3 has a unique, retainable ion.
- **mod@9** is easy too: `y3,y4,y5` are **P for mod@9 only**; `b6,b7,b8` are `.` for mod@9 only.
- **mod@6 (the middle isomer) has NO unique ion:**
  - `b3,b4,b5`: mod@6 is `.` -- but so is **mod@9**. Shared.
  - `b6,b7,b8`: mod@6 is `P` -- but so is **mod@3**. Shared.
  - `y3,y4,y5`: mod@6 is `.` -- but so is **mod@3**. Shared.
  - `y6,y7,y8`: mod@6 is `P` -- but so is **mod@9**. Shared.

  Every single row where mod@6 could stand out, it is tied with a neighbour. **There is no fragment
  at an m/z unique to mod@6.** "Keep the top unique fragment per peptidoform" would keep *nothing*
  for it, and mod@6 stays invisible.

The fix in intuition: mod@6 is still perfectly identifiable -- it is the only isomer that is
**`.` at b3-b5 AND `P` at b6-b8** (unmodified early, modified late). No single ion says that, but a
**pair** does. The deconvolution doesn't need a unique ion per column; it needs the columns to be
**distinct**, which a combination provides.

# 4. The algorithm that works: boundary-crossing coverage

Reframe: instead of "unique fragment," retain an observable ion **crossing every boundary between
consecutive candidate sites.** Two isomers differ exactly on the ions whose cut lies *between* their
sites, so covering every boundary separates every pair -- including the hard adjacent ones.

For VLSAGSPESIK, the candidate sites are `S = {3, 6, 9}`, giving two boundaries:

- **Boundary (3,6)** -- ions that flip between a "site <=3" and a "site 6" placement:
  crossing b-ions `b3,b4,b5` and mirror y-ions `y6,y7,y8`. Keep the **single most intense** of
  `{b3,b4,b5,y6,y7,y8}` -- say **b4**.
- **Boundary (6,9)** -- crossing b-ions `b6,b7,b8` and mirror y-ions `y3,y4,y5`. Keep the most intense
  of `{b6,b7,b8,y3,y4,y5}` -- say **b7**.

So we force-keep `D = {b4, b7}` (2 ions) for **every** peptidoform. Read off the pattern of the three
isomers across just these two ions:

| isomer | b4 | b7 |
|--------|:--:|:--:|
| mod@3  | P  | P  |
| mod@6  | .  | P  |
| mod@9  | .  | .  |

All three patterns are **distinct** -- mod@6 = `(., P)` stands apart even though it has no unique ion.
Two retained ions make all three isomer columns linearly independent, so the deconvolution can
separate them. Cost: ~`|S|` extra ions on top of the top-10 (here 2).

# 5. Why the candidate-site set must include the DECOY positions

The isomer-decoy plan moves a phosphate onto a **non-S/T/Y neighbour** of a real site. Suppose we add
a decoy for mod@6 on the adjacent **proline P7** (the classic phospho-proline decoy). Is `decoy@7`
distinguishable from the real `mod@6` using our `D = {b4, b7}`?

| ion | real mod@6 | decoy@7 | same? |
|-----|:---:|:---:|:---:|
| b4  | .   | .   | same |
| b7  | P   | P (7 <= 7) | same |

**No** -- `D = {b4, b7}` cannot tell the real site 6 from a decoy on 7. The ions that *do* separate
them are `b6` (mod@6 = P, decoy@7 = `.`) and `y5` (mod@6 = `.`, decoy@7 = P) -- and neither is in `D`,
because the boundary between 6 and 7 never existed in `S = {3,6,9}`.

**Fix: put the decoy candidate positions into `S`.** With `S = {3, 6, 7, 9}` there is now a
**boundary (6,7)**, whose crossing ions are `b6` and `y5`; keeping `b6` separates real@6 from decoy@7.
This is the crucial point: *the site set for fragment retention must be `(real S/T/Y) union
(decoy-candidate neighbours)`* -- otherwise the decoys we deliberately add are indistinguishable from
the reals, and the whole scheme is inert.

# 6. Optional refinement: minimum test-cover

Boundary coverage is simple and always sufficient, but may keep a redundant ion when one high-intensity
ion happens to cover two boundaries. If library size matters, solve the exact **minimum discriminating
set**:

1. Build a matrix `M[ion, pair]` = 1 if that ion separates that pair of peptidoforms (differs in P/.).
2. Greedily add the **highest-intensity** ion covering the most still-uncovered pairs, until every
   pair is covered.

With only a handful of sites the problem is tiny (exact or greedy is instant), and it yields the
smallest, most-intense retained set. Boundary coverage is the guaranteed-correct special case; use the
set-cover only if the ion budget is tight.

# 7. Implementation

One change, in `filter_fragments!(::InstrumentAgnosticModel)` (BuildSpecLib
`fragments/fragment_predict.jl:530`, the per-precursor top-N cap at ~line 608). We already hold the
**full** predicted b/y intensities there (Prosit returns the whole series; we are the ones
truncating), so no re-prediction is needed:

1. From the precursor sequence, compute `S = (S/T/Y positions) union (decoy-candidate neighbours)`.
2. Compute the boundary-crossing ion set `D` (Section 4), picking the max-intensity crossing ion per
   boundary.
3. Keep `top-N union D` instead of `top-N`.

Because a decoy is an m/z-permutation of a real isomer's *retained* fragments, retaining the reals'
boundary ions gives every decoy its shifted site-determining ions **for free** -- so this is a genuine
prerequisite for `localization_decoy_library_plan.md`, not a separate effort. Cost is bounded
(top-10 -> ~13-15 ions for 3-5 sites) and gated to peptides that contain S/T/Y.
