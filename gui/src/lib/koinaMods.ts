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

const BASE_MODS: KoinaMod[] = [
  { unimod: 4, label: 'Carbamidomethyl', mass: 57.021464, sites: ['C'] },
  { unimod: 35, label: 'Oxidation', mass: 15.994915, sites: ['M'] },
]

const PTM_MODS: KoinaMod[] = [
  { unimod: 1, label: 'Acetyl', mass: 42.010565, sites: ['K'] },
  { unimod: 4, label: 'Carbamidomethyl', mass: 57.021464, sites: ['C', 'K'] },
  { unimod: 7, label: 'Deamidated', mass: 0.984016, sites: ['N', 'Q', 'R'] },
  { unimod: 21, label: 'Phospho', mass: 79.966331, sites: ['S', 'T', 'Y', 'H', 'P'] },
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
  { unimod: 1263, label: 'Gly', mass: 57.021464, sites: ['C', 'K'] },
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
