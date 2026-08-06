/** FASTA header presets and modification presets.
 *
 *  Ported from `presetDefs` / `detectPreset` / `deriveName` / `matchPreset` /
 *  `modPresetDefs` in the design. The UniProt regexes are byte-identical to the
 *  ones Pioneer's own `GetBuildLibParams` template emits, so the default row
 *  produces exactly the config Pioneer would have written itself.
 */
import type { FastaEntry, HeaderPresetId, HeaderRegex, ModEntry } from './types'

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

export const MOD_PRESETS: Record<string, ModEntry & { menuLabel: string }> = {
  carbamidomethyl: {
    pattern: 'C',
    label: 'Carbamidomethyl',
    name: 'Unimod:4',
    mass: '57.021464',
    menuLabel: 'Carbamidomethyl (C) · +57.02',
  },
  oxidation: {
    pattern: 'M',
    label: 'Oxidation',
    name: 'Unimod:35',
    mass: '15.994915',
    menuLabel: 'Oxidation (M) · +15.99',
  },
  acetylk: {
    pattern: 'K',
    label: 'Acetyl',
    name: 'Unimod:1',
    mass: '42.010565',
    menuLabel: 'Acetyl (K) · +42.01',
  },
  phospho: {
    pattern: 'STY',
    label: 'Phospho',
    name: 'Unimod:21',
    mass: '79.966331',
    menuLabel: 'Phospho (STY) · +79.97',
  },
  deamidation: {
    pattern: 'NQ',
    label: 'Deamidation',
    name: 'Unimod:7',
    mass: '0.984016',
    menuLabel: 'Deamidation (NQ) · +0.98',
  },
  tmt6: {
    pattern: 'K',
    label: 'TMT6plex',
    name: 'Unimod:737',
    mass: '229.162932',
    menuLabel: 'TMT6plex (K) · +229.16',
  },
  custom: { pattern: '', label: '', name: '', mass: '', menuLabel: 'Custom…' },
}

/** The human label for a Unimod accession, when we know it. */
export function unimodLabel(name: string): string {
  for (const k of Object.keys(MOD_PRESETS)) {
    const d = MOD_PRESETS[k]
    if (d.name && d.name === name) return d.label
  }
  return ''
}
