/** FASTA header presets.
 *
 *  Ported from `presetDefs` / `detectPreset` / `deriveName` / `matchPreset` in
 *  the design. The UniProt regexes are byte-identical to the
 *  ones Pioneer's own `GetBuildLibParams` template emits, so the default row
 *  produces exactly the config Pioneer would have written itself.
 */
import { KOINA_MODS, unimodId } from './koinaMods'
import type { FastaEntry, HeaderPresetId, HeaderRegex } from './types'

export const HEADER_PRESETS: Record<HeaderPresetId, { label: string } & HeaderRegex> = {
  uniprot: {
    label: 'UniProt',
    accessions: '^\\w+\\|(\\w+(?:-\\d+)?)\\|',
    genes: ' GN=(\\S+)',
    proteins: '^\\w+\\|(?:\\w+(?:-\\d+)?)\\|[^ ]+ (.*?) [^ ]+=',
    organisms: ' OS=([^ ]+.*?) [^ ]+=',
  },
  // The design's Ensembl and RefSeq patterns were written against invented
  // headers and do not match real files — see the note below. These are
  // verified against Ensembl release-current peptide FASTAs (yeast R64-1-1 and
  // human GRCh38) and an NCBI RefSeq protein FASTA (E. coli K-12 ASM584v2).
  // An empty pattern is meaningful: Pioneer maps it to `nothing`, i.e. "this
  // format carries no such field", which is better than a wrong capture.
  ensembl: {
    label: 'Ensembl',
    accessions: '^(\\S+)',
    // `gene:` is always present; `gene_symbol:` is friendlier but exists only
    // in some species (human has it, yeast does not), and it appears after
    // `gene:` so a single pattern cannot prefer it.
    genes: ' gene:(\\S+)',
    proteins: ' description:(.*?)(?: \\[Source:|$)',
    organisms: '',
  },
  refseq: {
    label: 'RefSeq',
    accessions: '^(\\S+)',
    genes: '',
    proteins: '^\\S+ (.*?) \\[',
    organisms: '\\[([^\\]]+)\\]',
  },
  custom: { label: 'Custom', accessions: '', genes: '', proteins: '', organisms: '' },
}

const regexOf = (id: HeaderPresetId): HeaderRegex => {
  const p = HEADER_PRESETS[id]
  return { accessions: p.accessions, genes: p.genes, proteins: p.proteins, organisms: p.organisms }
}

/** Guess the header format from the filename. */
export function detectPreset(path: string): HeaderPresetId {
  const f = String(path || '').toLowerCase()
  if (/\.pep\.all|ensembl|ensp\d|grch3|grcm3|\.cdna\./.test(f)) return 'ensembl'
  if (/refseq|gcf_|gca_|_protein\.faa|\.faa(\.gz)?$/.test(f)) return 'refseq'
  return 'uniprot'
}

/** Guess the species tag from the filename. */
export function deriveName(path: string): string {
  const f = String(path || '').toLowerCase()
  if (/human|sapiens|9606/.test(f)) return 'HUMAN'
  if (/mouse|musculus|10090|grcm/.test(f)) return 'MOUSE'
  if (/yeast|cerevisiae|4932/.test(f)) return 'YEAST'
  if (/coli|83333/.test(f)) return 'ECOLI'
  const base = (String(path || '').split(/[\\/]/).pop() || '').replace(
    /\.(fasta|fa|faa|fna|fas)(\.gz)?$/i,
    '',
  )
  return base ? base.toUpperCase().slice(0, 18) : 'SPECIES'
}

/** Which preset a set of regexes corresponds to, or 'custom'. */
export function matchPreset(regex: HeaderRegex): HeaderPresetId {
  for (const id of ['uniprot', 'ensembl', 'refseq'] as const) {
    const p = HEADER_PRESETS[id]
    if (
      p.accessions === regex.accessions &&
      p.genes === regex.genes &&
      p.proteins === regex.proteins &&
      p.organisms === regex.organisms
    ) {
      return id
    }
  }
  return 'custom'
}

export function makeFastaRow(path: string): FastaEntry {
  const presetId = detectPreset(path)
  return {
    path,
    name: deriveName(path),
    presetId,
    regex: regexOf(presetId),
    auto: true,
    showRegex: false,
  }
}

export { regexOf as presetRegex }

/** The human label for a Unimod accession, when we know it.
 *
 *  Looks across every model's catalogue rather than the selected model's, so a
 *  config naming a modification the current model cannot predict still shows
 *  its name — ModTable needs that to say what it is rejecting.
 */
export function unimodLabel(name: string): string {
  const id = unimodId(name)
  if (id === null) return ''
  for (const mods of Object.values(KOINA_MODS)) {
    const hit = mods.find((m) => m.unimod === id)
    if (hit) return hit.label
  }
  return ''
}

/** One of `library_info`'s modification entries, with the accession replaced by
 *  the modification's name.
 *
 *  The Rust side formats these as `"<name> on <pattern>"`, where `<name>` is
 *  whatever the library's config recorded — in practice `"Unimod:4"`. That is
 *  what the config stores but not what anyone calls the modification, so
 *  `"Unimod:4 on C"` reads as `"Carbamidomethyl on C"`.
 *
 *  An accession we do not recognise is left exactly as it came: opaque beats
 *  wrong, and the number is still something the reader can look up.
 */
export function unimodDisplay(entry: string): string {
  const sep = entry.indexOf(' on ')
  const label = unimodLabel(sep === -1 ? entry : entry.slice(0, sep))
  if (!label) return entry
  return sep === -1 ? label : label + entry.slice(sep)
}
