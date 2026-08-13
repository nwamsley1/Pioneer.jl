/** Cleavage rules offered on the BuildSpecLib page.
 *
 *  Pioneer digests with a regex whose *first* element is the P1 residue — the
 *  peptide ends where the match starts — and whose remainder is context the
 *  match consumes. Enzymes that cut N-terminal to a residue are therefore
 *  written with a lookahead instead. Every pattern here is checked against a
 *  hand-derived digest in `test/UnitTests/test_cleavage_presets.jl`; change one
 *  and that test should be the thing that objects.
 *
 *  The `_` in each excluded set matches what
 *  `assets/example_config/defaultBuildLibParams.json` ships, so a preset
 *  round-trips through a config unchanged. It never occurs in a protein
 *  sequence.
 *
 *  Specificities are taken from ExPASy PeptideCutter's enzyme table
 *  (https://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html), checked
 *  against Comet's `search_enzyme_number` definitions for what search engines
 *  conventionally do. Where the two differ it is noted on the entry.
 *
 *  Three simplifications, all in the direction of cleaving more rather than
 *  less: PeptideCutter blocks tryptic cleavage for a few P2/P1' combinations
 *  (K-W, R-M, and some involving D, C and H), blocks chymotryptic cleavage
 *  after W when M follows, and blocks it after H when D, M or W follows. None
 *  of those are modelled here, and no search engine models them either.
 *
 *  CNBr is deliberately absent. It cleaves after M, but converts that
 *  methionine to homoserine lactone, and nothing in Pioneer models the residue
 *  change -- every peptide's C-terminal mass would be wrong. Anyone who needs
 *  it can write the rule by hand and accept that.
 */
export interface Enzyme {
  id: string
  label: string
  /** What it does, in the words a method section would use. */
  rule: string
  pattern: string
}

export const ENZYMES: Enzyme[] = [
  { id: 'trypsin_p', label: 'Trypsin/P', rule: 'after K or R, including before proline',
    pattern: '[KR][^_|$]' },
  { id: 'trypsin', label: 'Trypsin', rule: 'after K or R, but not before proline',
    pattern: '[KR][^P_|$]' },
  { id: 'lysc', label: 'Lys-C', rule: 'after K, but not before proline',
    pattern: 'K[^P_|$]' },
  { id: 'lysc_p', label: 'Lys-C/P', rule: 'after K, including before proline',
    pattern: 'K[^_|$]' },
  { id: 'lysn', label: 'Lys-N', rule: 'before K',
    pattern: '[^K](?=K)' },
  // PeptideCutter says P1' has minimal effect on Arg-C, i.e. proline does not
  // block it; Comet blocks proline. The search-engine convention is followed
  // here, since that is what a library is compared against.
  { id: 'argc', label: 'Arg-C', rule: 'after R, but not before proline',
    pattern: 'R[^P_|$]' },
  { id: 'aspn', label: 'Asp-N', rule: 'before D',
    pattern: '[^D](?=[D])' },
  { id: 'aspn_glu', label: 'Asp-N + Glu', rule: 'before D or E',
    pattern: '[^DE](?=[DE])' },
  { id: 'gluc_bicarb', label: 'Glu-C (ammonium bicarbonate)', rule: 'after E, but not before proline',
    pattern: 'E[^P_|$]' },
  { id: 'gluc_phos', label: 'Glu-C (phosphate)', rule: 'after D or E, but not before proline',
    pattern: '[DE][^P_|$]' },
  { id: 'chymotrypsin', label: 'Chymotrypsin (high specificity)',
    rule: 'after F, W or Y, but not before proline',
    pattern: '[FWY][^P_|$]' },
  // PeptideCutter's low-specificity set is F, Y, W, L, M *and H*. Comet's
  // Chymotrypsin is FWYL, omitting both M and H. The fuller set is used, since
  // the point of this entry is the permissive one.
  { id: 'chymotrypsin_broad', label: 'Chymotrypsin (low specificity)',
    rule: 'after F, W, Y, L, M or H, but not before proline',
    pattern: '[FWYLMH][^P_|$]' },
]

/** Selected when nothing matches: the pattern is then edited directly. */
export const CUSTOM_ENZYME = 'custom'

/** Pioneer's own default, so a library built without touching this field is the
 *  one the CLI would have built. Note it is Trypsin/P rather than Trypsin. */
export const DEFAULT_CLEAVAGE = ENZYMES[0].pattern

export function enzymeByPattern(pattern: string): Enzyme | null {
  return ENZYMES.find((e) => e.pattern === pattern.trim()) ?? null
}

export function enzymeById(id: string): Enzyme | null {
  return ENZYMES.find((e) => e.id === id) ?? null
}

/** A sequence carrying one motif for every rule above, used for the preview. */
export const SAMPLE_SEQUENCE = 'AKGRSTKDEKPARQGGSEAADFPGWSYLMTFD'

/**
 * Digest `sequence` the way `digest_sequence` does, for previewing a rule.
 *
 * A peptide ends at the *start* of each match, and matching continues from one
 * past that start rather than past the whole match, which is the effect of
 * Julia's `overlap = true` — without it the context a match consumes would hide
 * the site immediately after it.
 *
 * This is a preview, not the digest: Pioneer has the last word, and the point
 * here is to show what a hand-written rule does before a library is built on
 * it.
 */
export function previewDigest(sequence: string, pattern: string): string[] | null {
  let re: RegExp
  try {
    re = new RegExp(pattern, 'g')
  } catch {
    return null
  }

  const peptides: string[] = []
  let previous = 0
  for (let i = 0; i < sequence.length; i++) {
    re.lastIndex = i
    const m = re.exec(sequence)
    if (!m || m.index !== i) continue
    // Zero-width matches would otherwise cut at every position.
    if (m[0].length === 0) continue
    const end = i + 1
    if (end > previous) peptides.push(sequence.slice(previous, end))
    previous = end
  }
  if (previous < sequence.length) peptides.push(sequence.slice(previous))
  return peptides
}
