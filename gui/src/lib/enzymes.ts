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

export type DigestSpecificity = 'full' | 'semi' | 'semi-n' | 'semi-c'

/** The digestion settings the preview depends on. */
export interface DigestSettings {
  pattern: string
  specificity: DigestSpecificity
  minLength: number
  maxLength: number
  missedCleavages: number
}

/**
 * The cleavage sites a rule finds: index `i` is true when a peptide may end at
 * residue `i`.
 *
 * A site sits at the *start* of each match, and matching continues from one past
 * that start rather than past the whole match, which is the effect of Julia's
 * `overlap = true` — without it the context a match consumes would hide the site
 * immediately after it. Returns null when the pattern will not compile.
 */
function cleavageMask(sequence: string, pattern: string): boolean[] | null {
  let re: RegExp
  try {
    re = new RegExp(pattern, 'g')
  } catch {
    return null
  }
  const mask = new Array<boolean>(sequence.length).fill(false)
  for (let i = 0; i < sequence.length; i++) {
    re.lastIndex = i
    const m = re.exec(sequence)
    // Zero-width matches count. They were skipped here, on the grounds that an
    // empty pattern would otherwise cut at every position -- true, but so does
    // it in Julia, and the price was that a bare lookahead like `(?=D)` found
    // no sites at all while `digest_sequence` found six. A rule that really
    // does cut everywhere is caught by the "cleaves almost everywhere" warning
    // rather than by silently finding nothing.
    if (m && m.index === i) mask[i] = true
  }
  return mask
}

/**
 * The positions where `pattern` finds a cleavage site, or null when it will not
 * compile.
 *
 * Separate from the digest because judging the *rule* must not depend on the
 * length and specificity filters: a sound rule can still yield no peptides once
 * a narrow length range is applied, and that is a different complaint.
 */
export function cleavageSites(sequence: string, pattern: string): number[] | null {
  const mask = cleavageMask(sequence, pattern)
  if (!mask) return null
  const sites: number[] = []
  for (let i = 0; i < mask.length; i++) if (mask[i]) sites.push(i)
  return sites
}

/**
 * Digest `sequence` the way `digest_sequence` does, for previewing the settings.
 *
 * A port of `fasta_digest.jl`, including the length and missed-cleavage filters,
 * so the preview lists the peptides a library would actually carry rather than
 * the fragments the enzyme cuts out of the sequence. That distinction matters
 * most under semi specificity, where the peptide count runs an order of
 * magnitude above the fully-specific one and is otherwise invisible until a
 * build takes all night.
 *
 * Protein termini count as enzymatic, as they do in Pioneer. Still a preview:
 * Pioneer has the last word, and the point is to show what a rule does before a
 * library is built on it.
 *
 * Checked against `digest_sequence` over every preset and several hand-written
 * rules: 30 of 31 pattern/specificity pairs agree peptide-for-peptide. The
 * exception is a *bare* zero-width rule such as `(?=D)` under full specificity,
 * where Julia's separate `_digest_fully_specific_sequence` path emits the final
 * peptide twice; this returns it once. Showing a duplicate in a preview would be
 * the worse of the two behaviours, so the divergence is deliberate. No preset is
 * affected — they all consume a residue before their lookahead.
 */
export function previewDigest(
  sequence: string,
  settings: DigestSettings,
): string[] | null {
  const mask = cleavageMask(sequence, settings.pattern)
  if (!mask) return null

  const n = sequence.length
  const { specificity, minLength, maxLength, missedCleavages } = settings
  if (n === 0 || minLength < 1) return []

  // Prefix counts make the missed-cleavage test O(1) per candidate.
  const prefix = new Array<number>(n).fill(0)
  let running = 0
  for (let i = 0; i < n; i++) {
    if (mask[i]) running++
    prefix[i] = running
  }

  // Every position a peptide may end at. The final residue always qualifies:
  // a protein's C terminus counts as enzymatic.
  const specificEnds: number[] = []
  for (let i = 0; i < n; i++) if (mask[i]) specificEnds.push(i)
  if (specificEnds[specificEnds.length - 1] !== n - 1) specificEnds.push(n - 1)

  const startEnzymatic = (start: number) => start === 0 || mask[start - 1]
  const internalCleavages = (start: number, end: number) =>
    end <= start ? 0 : prefix[end - 1] - (start > 0 ? prefix[start - 1] : 0)

  const peptides: string[] = []
  for (let start = 0; start < n; start++) {
    const minEnd = start + minLength - 1
    if (minEnd > n - 1) continue
    const maxEnd = Math.min(start + maxLength - 1, n - 1)
    const startIsEnzymatic = startEnzymatic(start)

    // semi-c allows an arbitrary C terminus but requires an enzymatic N
    // terminus; general semi does the same for enzymatic starts.
    if (specificity === 'semi-c' || (specificity === 'semi' && startIsEnzymatic)) {
      if (!startIsEnzymatic) continue
      for (let end = minEnd; end <= maxEnd; end++) {
        if (internalCleavages(start, end) <= missedCleavages) {
          peptides.push(sequence.slice(start, end + 1))
        }
      }
      continue
    }

    // full requires both termini enzymatic; semi-n and semi with a
    // non-enzymatic start require an enzymatic C terminus. All three therefore
    // end only at a cleavage site.
    if (specificity === 'full' && !startIsEnzymatic) continue
    for (const end of specificEnds) {
      if (end < minEnd) continue
      if (end > maxEnd) break
      if (internalCleavages(start, end) <= missedCleavages) {
        peptides.push(sequence.slice(start, end + 1))
      }
    }
  }
  return peptides
}
