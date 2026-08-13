/** Modifications each Koina prediction model actually accepts.
 *
 *  Not guesswork: Koina encodes every supported PTM in
 *  `Prosit_Preprocess_ac_gain/1/model.py` as `dict_ptm_atom_count_gain`, keyed
 *  `<residue>_<unimodId>` and valued as atom counts over {H,C,N,O,P,S}. The
 *  2024 and 2025 PTM models share that preprocessing step, so their support is
 *  identical. Verified against the live endpoint in both directions: all 60
 *  table pairs are accepted by both models, and 66 probes using a listed
 *  modification on an *unlisted* residue were all rejected — the residue
 *  scoping is real, which is why `sites` is per modification.
 *
 *  Excluded from the table on purpose:
 *    - UNIMOD 5634 — Koina accepts it, but it is absent from unimod.obo, so
 *      there is no mass delta to write into the config.
 *    - UNIMOD 411 (Phenylisocyanate) and the other N-terminal entries — Koina
 *      takes those only as `[UNIMOD:n]-PEPTIDE`, and Pioneer's `pattern` is a
 *      plain Regex matched per residue (`eachmatch` in add_mods), with no
 *      N-terminal concept to express them.
 *    - Koina's own broken rows: R_267 has an empty atom string, and
 *      12118/19903/129317 are not valid UNIMOD accessions.
 *
 *  Altimeter and Prosit 2020 HCD were probed the same way and take only the
 *  two defaults.
 */

export interface KoinaMod {
  unimod: number
  /** UNIMOD's own name, shown as the row description. */
  label: string
  /** Monoisotopic delta from unimod.obo. */
  mass: number
  /** Residues this model accepts the modification on, conventional one first. */
  sites: string[]
}

// Every entry is both accepted by the model and a real UNIMOD record: the two
// are separate constraints and a modification has to satisfy both. Masses and
// sites were checked against UNIMOD's own unimod.obo (1561 terms) on
// 2026-08-13 -- all 29 monoisotopic masses agreed to better than 5e-7, which
// matters because Pioneer takes the mass straight from the config rather than
// looking it up (`chronologer_prep.jl` reads params["fixed_mods"]["mass"]).
//
// Two sites were removed as not being UNIMOD sites at all: Phospho on P, which
// is not a modification that exists, and Gly on C, which was additionally
// isobaric with the carbamidomethyl already pinned to every cysteine. Every
// remaining site is a subset of what UNIMOD lists, as it should be -- the
// model's support is narrower than the chemistry.
//
// `Glutarylation` keeps its correct spelling; UNIMOD:1848's own record name
// reads "Gluratylation". Only the accession is written to a config, so the
// label never leaves the GUI.
const BASE_MODS: KoinaMod[] = [
  { unimod: 4, label: 'Carbamidomethyl', mass: 57.021464, sites: ['C'] },
  { unimod: 35, label: 'Oxidation', mass: 15.994915, sites: ['M'] },
]

const PTM_MODS: KoinaMod[] = [
  { unimod: 1, label: 'Acetyl', mass: 42.010565, sites: ['K'] },
  { unimod: 4, label: 'Carbamidomethyl', mass: 57.021464, sites: ['C', 'K'] },
  { unimod: 7, label: 'Deamidated', mass: 0.984016, sites: ['N', 'Q', 'R'] },
  { unimod: 21, label: 'Phospho', mass: 79.966331, sites: ['S', 'T', 'Y', 'H'] },
  { unimod: 27, label: 'Glu->pyro-Glu', mass: -18.010565, sites: ['E'] },
  { unimod: 28, label: 'Gln->pyro-Glu', mass: -17.026549, sites: ['Q'] },
  { unimod: 34, label: 'Methyl', mass: 14.01565, sites: ['K', 'R', 'D', 'E', 'C', 'H', 'I', 'L', 'N', 'Q'] },
  { unimod: 35, label: 'Oxidation', mass: 15.994915, sites: ['M', 'W', 'C', 'H', 'K', 'P'] },
  { unimod: 36, label: 'Dimethyl', mass: 28.0313, sites: ['K', 'R'] },
  { unimod: 37, label: 'Trimethyl', mass: 42.04695, sites: ['K'] },
  { unimod: 43, label: 'HexNAc', mass: 203.079373, sites: ['S', 'T'] },
  { unimod: 56, label: 'Acetyl:2H(3)', mass: 45.029395, sites: ['K'] },
  { unimod: 58, label: 'Propionyl', mass: 56.026215, sites: ['K'] },
  { unimod: 59, label: 'Propionyl:13C(3)', mass: 59.036279, sites: ['K'] },
  { unimod: 121, label: 'GG', mass: 114.042927, sites: ['K'] },
  { unimod: 214, label: 'iTRAQ4plex', mass: 144.102063, sites: ['K'] },
  { unimod: 312, label: 'Cysteinyl', mass: 119.004099, sites: ['C'] },
  { unimod: 535, label: 'LRGG', mass: 383.228103, sites: ['K'] },
  { unimod: 730, label: 'iTRAQ8plex', mass: 304.20536, sites: ['K'] },
  { unimod: 737, label: 'TMT6plex', mass: 229.162932, sites: ['K'] },
  { unimod: 1263, label: 'Gly', mass: 57.021464, sites: ['K'] },
  { unimod: 1289, label: 'Butyryl', mass: 70.041865, sites: ['K'] },
  { unimod: 1293, label: 'QTGG', mass: 343.149184, sites: ['K'] },
  { unimod: 1848, label: 'Glutarylation', mass: 114.031694, sites: ['K'] },
  { unimod: 1990, label: 'CIGG', mass: 330.136176, sites: ['K'] },
  { unimod: 2016, label: 'TMTpro', mass: 304.207146, sites: ['K'] },
  { unimod: 2062, label: 'DBIA', mass: 296.184841, sites: ['C'] },
]

/** Keyed by PredictionModel.id. */
export const KOINA_MODS: Record<string, KoinaMod[]> = {
  altimeter: BASE_MODS,
  prosit_2020_hcd: BASE_MODS,
  prosit_2024_ptm: PTM_MODS,
  prosit_2025_40ptm: PTM_MODS,
}

export function modsForModel(modelId: string): KoinaMod[] {
  return KOINA_MODS[modelId] ?? BASE_MODS
}

/** "Unimod:35" -> 35. Tolerates the "UNIMOD:35" and bare "35" spellings that
 *  turn up in hand-written configs. */
export function unimodId(name: string): number | null {
  const m = /(\d+)\s*$/.exec(name.trim())
  return m ? Number(m[1]) : null
}

export function findMod(modelId: string, name: string): KoinaMod | null {
  const id = unimodId(name)
  return id === null ? null : (modsForModel(modelId).find((m) => m.unimod === id) ?? null)
}

/** The `pattern` regex for a set of residues: one residue is a literal, several
 *  need a character class. `'STY'` would match the literal tripeptide. */
export function sitePattern(sites: string[]): string {
  return sites.length === 1 ? sites[0] : `[${sites.join('')}]`
}

/** The site choices offered for one modification: each residue on its own, plus
 *  an all-residues entry when there is more than one. */
export function siteOptions(mod: KoinaMod): { value: string; label: string }[] {
  const one = mod.sites.map((s) => ({ value: s, label: s }))
  if (mod.sites.length < 2) return one
  const all = [...mod.sites].sort()
  return [...one, { value: sitePattern(all), label: all.join('') }]
}

/** A ModEntry for a modification, on its conventional site. Structurally a
 *  `ModEntry`; not typed as one so this module stays free of imports and
 *  `types.ts` can build its defaults from it without a cycle. */
export function modEntry(
  modelId: string,
  unimod: number,
): { pattern: string; label: string; name: string; mass: string } {
  const def = modsForModel(modelId).find((m) => m.unimod === unimod)
  if (!def) throw new Error(`${modelId} has no UNIMOD ${unimod}`)
  return { pattern: def.sites[0], label: def.label, name: `Unimod:${def.unimod}`, mass: String(def.mass) }
}

/** True when `pattern` only names residues this model allows for `mod`. */
export function siteAllowed(mod: KoinaMod, pattern: string): boolean {
  const residues = pattern.replace(/[[\]]/g, '')
  return residues.length > 0 && [...residues].every((c) => mod.sites.includes(c))
}

/** UNIMOD 4 on C, Carbamidomethyl.
 *
 *  Every Koina model currently supported is trained on alkylated cysteine, so
 *  the C site is not a choice: always fixed, never variable. A library built
 *  without it does not match what the model predicts, and the mismatch is
 *  silent -- the search simply finds less.
 *
 *  The pin is on the *site*, not the modification. Models whose UNIMOD 4 also
 *  covers K leave K entirely free: fixed, variable or absent, as the user
 *  likes. Multi-site mods are one row with a combined pattern, so "fixed on C
 *  and K" is the single fixed row carrying `[CK]` rather than a second row.
 *
 *  Kept here rather than in the form so a future model which does not want it
 *  opts out in one place.
 */
export const REQUIRED_FIXED_UNIMOD = 4
export const REQUIRED_FIXED_SITE = 'C'

/** Models that require it. Currently all of them; named explicitly so adding a
 *  model is a decision rather than an inheritance. */
const REQUIRES_FIXED_ALKYLATION = new Set(Object.keys(KOINA_MODS))

export function requiresFixedAlkylation(modelId: string): boolean {
  return REQUIRES_FIXED_ALKYLATION.has(modelId)
}

/** The residues a pattern names, `[CK]` -> `['C','K']`. */
function residuesOf(pattern: string): string[] {
  return [...pattern.replace(/[[\]]/g, '')]
}

function isAlkylation(name: string): boolean {
  return unimodId(name) === REQUIRED_FIXED_UNIMOD
}

/** True for the one fixed row the model pins: UNIMOD 4 covering C. It may also
 *  cover K -- that is still the pinned row, and still not removable. */
export function isRequiredFixedMod(modelId: string, name: string, pattern: string): boolean {
  return (
    requiresFixedAlkylation(modelId) &&
    isAlkylation(name) &&
    residuesOf(pattern).includes(REQUIRED_FIXED_SITE)
  )
}

/** Residues already taken by the fixed modifications.
 *
 *  A fixed modification applies to every matching residue, so nothing variable
 *  can sit on one: the residue would be modified twice. Pioneer rejects such a
 *  config outright (`check_mod_site_conflicts`), so the point of computing this
 *  is to keep the choice from being offered in the first place. */
export function occupiedResidues(fixed: { name: string; pattern: string }[]): Set<string> {
  const taken = new Set<string>()
  for (const m of fixed) for (const r of residuesOf(m.pattern)) taken.add(r)
  return taken
}

/** Site choices permitted for one row.
 *
 *  A fixed alkylation row must keep C. A variable row may not name any residue
 *  a fixed modification already holds — which covers C automatically, since the
 *  alkylation row is always present, so that case needs no rule of its own. */
export function allowedSiteValues(
  modelId: string,
  kind: 'fixed' | 'variable',
  name: string,
  values: string[],
  occupied: Set<string> = new Set(),
): string[] {
  if (kind === 'fixed') {
    if (requiresFixedAlkylation(modelId) && isAlkylation(name)) {
      return values.filter((v) => residuesOf(v).includes(REQUIRED_FIXED_SITE))
    }
    return values
  }
  return values.filter((v) => !residuesOf(v).some((r) => occupied.has(r)))
}

/** The site a newly added row should start on: C for a fixed alkylation row,
 *  and the first site that is not C for a variable one. */
export function initialSite(
  modelId: string,
  kind: 'fixed' | 'variable',
  mod: KoinaMod,
  occupied: Set<string> = new Set(),
): string | null {
  if (kind === 'fixed') {
    if (requiresFixedAlkylation(modelId) && mod.unimod === REQUIRED_FIXED_UNIMOD) {
      return REQUIRED_FIXED_SITE
    }
    return mod.sites[0] ?? null
  }
  const free = mod.sites.filter((site) => !occupied.has(site))
  // null means every site this model allows is already held by a fixed
  // modification, so there is no variable form of it left to add.
  return free.length ? free[0] : null
}

/** Force the rule onto a pair of mod lists: exactly one fixed row covering C,
 *  and no variable row that names C.
 *
 *  Applied wherever the lists can change from outside the mod editor -- a model
 *  switch, a loaded config -- so the invariant cannot be dodged by a route that
 *  bypasses the UI. A fixed row already covering C keeps its pattern, so `[CK]`
 *  survives; one that covers only K gains C rather than being replaced, since
 *  dropping it would silently discard a modification the user chose.
 */
export function enforceRequiredMods(
  modelId: string,
  fixed: { pattern: string; label: string; name: string; mass: string }[],
  variable: { pattern: string; label: string; name: string; mass: string }[],
) {
  if (!requiresFixedAlkylation(modelId)) return { fixed, variable }

  const existing = fixed.filter((m) => isAlkylation(m.name))
  const rest = fixed.filter((m) => !isAlkylation(m.name))
  let required: { pattern: string; label: string; name: string; mass: string }
  if (existing.length === 0) {
    required = modEntry(modelId, REQUIRED_FIXED_UNIMOD)
  } else {
    const sites = new Set(existing.flatMap((m) => residuesOf(m.pattern)))
    sites.add(REQUIRED_FIXED_SITE)
    required = { ...existing[0], pattern: sitePattern([...sites].sort()) }
  }

  const taken = occupiedResidues([required, ...rest])
  return {
    fixed: [required, ...rest],
    // Every residue held by a fixed modification is stripped, not just C: a
    // loaded or hand-written config can put oxidation on a cysteine that
    // carbamidomethyl already occupies, which Pioneer rejects outright.
    variable: variable.flatMap((m) => {
      const free = residuesOf(m.pattern).filter((r) => !taken.has(r))
      // Nothing survives, so the row described only occupied residues.
      return free.length ? [{ ...m, pattern: sitePattern(free) }] : []
    }),
  }
}
