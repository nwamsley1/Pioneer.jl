/** Reading the published library names.
 *
 *  Catalog libraries are named
 *  `{Model}_{Genus}_{species}_{UPID...}_up{YYYY}_{RR}_dl{YYYY-MM-DD}.poin`,
 *  which is unambiguous but long. This splits a name into the columns the
 *  download table shows. A name that does not follow the convention parses to
 *  `null`, and the table falls back to showing it whole.
 */

export interface ParsedLibraryName {
  /** Prediction model, e.g. "Altimeter". */
  model: string
  /** Name a user would search for, e.g. "Human". */
  common: string
  /** Latin binomial, or '' for a mixed-proteome library. */
  latin: string
  /** UniProt proteome IDs, e.g. ["UP000005640"]. */
  proteomes: string[]
  /** UniProt release, e.g. "2026_03". */
  release: string
  /** Date the FASTA was downloaded, YYYY-MM-DD. */
  downloaded: string
}

const COMMON_NAMES: Record<string, string> = {
  'Homo sapiens': 'Human',
  'Mus musculus': 'Mouse',
  'Rattus norvegicus': 'Rat',
  'Saccharomyces cerevisiae': 'Yeast',
  'Escherichia coli': 'E. coli',
  'Arabidopsis thaliana': 'Arabidopsis',
  'Drosophila melanogaster': 'Fruit fly',
  'Caenorhabditis elegans': 'C. elegans',
}

/** Mixed-proteome libraries, keyed by the token that follows the model. */
const MIXES: Record<string, string> = {
  HYE: 'Human + Yeast + E. coli',
}

const NAME_RE = /^([A-Za-z0-9]+)_(.+)_up(\d{4}_\d{2})_dl(\d{4}-\d{2}-\d{2})$/
const PROTEOME_RE = /^UP\d{9}$/

export function parseLibraryName(name: string): ParsedLibraryName | null {
  const m = NAME_RE.exec(name.replace(/\.poin$/, ''))
  if (!m) return null
  const [, model, middle, release, downloaded] = m
  const tokens = middle.split('_')
  const proteomes = tokens.filter((t) => PROTEOME_RE.test(t))
  const words = tokens.filter((t) => !PROTEOME_RE.test(t))
  if (proteomes.length === 0 || words.length === 0) return null

  const mix = MIXES[words[0]]
  if (mix) return { model, common: mix, latin: '', proteomes, release, downloaded }

  const latin = words.join(' ')
  return { model, common: COMMON_NAMES[latin] ?? latin, latin, proteomes, release, downloaded }
}
