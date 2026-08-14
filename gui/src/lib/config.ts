/** Translating between the SearchDIA form and Pioneer's params JSON.
 *
 *  Ported from `buildJsonBase` / `configToState` / `computeExtras` in
 *  "Pioneer GUI - 1 Console.dc.html". The round-trip matters: a user can load
 *  a config.json from an earlier run, tweak two fields in the form, and run it
 *  again — and any keys the form does not model must survive that trip. They
 *  are held aside in `extraConfig` and merged back on the way out.
 */
import { DEFAULT_CLEAVAGE } from './enzymes'
import { makeFastaRow, matchPreset, presetRegex, unimodLabel } from './fasta'
import type { BuildParams, ConvertParams, DownloadParams, ModEntry, SearchParams } from './types'

export type Json = Record<string, unknown>

const num = (v: string, fallback: number): number => {
  const n = parseFloat(v)
  return Number.isNaN(n) ? fallback : n
}

/** Placeholder shown in the JSON preview when a path is not yet set. */
const disp = (v: string, fallback: string): string => (v && v.length ? v : fallback)

/** The JSON paths the SearchDIA form owns. Anything else in a loaded config is
 *  "extra" and is preserved verbatim. */
export const SEARCH_OWNED_PATHS = [
  'paths.ms_data',
  'paths.library',
  'paths.results',
  'output.write_csv',
  'output.write_decoys',
  'output.delete_temp',
  'global.q_value_threshold',
  'global.match_between_runs',
  'search.n_isotopes',
  'acquisition.nce',
  'optimization.chromatogram_integration.trace_mode',
  'proteinScoring.min_peptides',
  'maxLFQ.run_to_run_normalization',
  'logging.debug_console_level',
]

export function buildSearchJsonBase(s: SearchParams): Json {
  return {
    paths: {
      ms_data: disp(s.msData, '/path/to/ms/data'),
      library: disp(s.library, '/path/to/library.poin'),
      results: disp(s.results, '/path/to/results'),
    },
    output: {
      write_csv: s.writeCsv,
      write_decoys: s.writeDecoys,
      delete_temp: s.deleteTemp,
    },
    global: {
      q_value_threshold: num(s.qValue, 0.01),
      // Emitted explicitly even though Pioneer's fallback matches the default.
      // The config is the record of what a run was told to do, and MBR changes
      // the result enough that "it was left unset" is not a useful answer later.
      match_between_runs: s.matchBetweenRuns,
    },
    search: { n_isotopes: num(s.nIsotopes, 2) },
    acquisition: { nce: num(s.nce, 26) },
    // Always combined. The GUI no longer offers a choice, and the value it
    // used to send for the other option -- "separated" -- was not one Pioneer
    // accepts. It expects "combined" or "separate", so choosing it failed the
    // run at IntegrateChromatogramsSearch.
    optimization: { chromatogram_integration: { trace_mode: 'combined' } },
    proteinScoring: { min_peptides: num(s.minPeptides, 1) },
    maxLFQ: { run_to_run_normalization: s.runToRunNorm },
    // Console verbosity only. The debug *file* level is left alone: Pioneer
    // defaults it to 1 because the log is useful after the fact whether or not
    // the console was noisy at the time.
    logging: { debug_console_level: s.debugLogging ? 1 : 0 },
  }
}

export function buildSearchJson(s: SearchParams, extra: Json | null): Json {
  const base = buildSearchJsonBase(s)
  return extra ? deepMerge(base, extra) : base
}

const isObj = (x: unknown): x is Json =>
  typeof x === 'object' && x !== null && !Array.isArray(x)

export function deepMerge(base: Json, extra: Json): Json {
  const out: Json = { ...base }
  for (const k of Object.keys(extra)) {
    const e = extra[k]
    const b = out[k]
    out[k] = isObj(e) && isObj(b) ? deepMerge(b, e) : e
  }
  return out
}

function delPath(obj: unknown, parts: string[]): void {
  if (!isObj(obj)) return
  if (parts.length === 1) {
    delete obj[parts[0]]
    return
  }
  delPath(obj[parts[0]], parts.slice(1))
}

function pruneEmpty(obj: Json): Json {
  for (const k of Object.keys(obj)) {
    const v = obj[k]
    if (isObj(v)) {
      pruneEmpty(v)
      if (Object.keys(v).length === 0) delete obj[k]
    }
  }
  return obj
}

/** Strip the keys the form owns; whatever remains is carried through untouched. */
export function computeExtras(obj: unknown, ownedPaths: string[]): Json | null {
  if (!isObj(obj)) return null
  const clone = JSON.parse(JSON.stringify(obj)) as Json
  ownedPaths.forEach((p) => delPath(clone, p.split('.')))
  pruneEmpty(clone)
  return Object.keys(clone).length ? clone : null
}

/** Dotted leaf paths of the preserved extras, for the "custom keys" hint. */
export function extraLeafPaths(obj: Json | null, prefix = ''): string[] {
  if (!obj) return []
  let out: string[] = []
  for (const k of Object.keys(obj)) {
    if (/^comment/i.test(k)) continue
    const path = prefix ? `${prefix}.${k}` : k
    const v = obj[k]
    if (isObj(v)) out = out.concat(extraLeafPaths(v, path))
    else out.push(path)
  }
  return out
}

// ---------------------------------------------------------------------------
// BuildSpecLib
// ---------------------------------------------------------------------------

/** The JSON paths the BuildSpecLib form owns. Everything else in a loaded
 *  config — `nce_params`, `library_params`, `isotope_mod_groups`, the digestion
 *  keys the form does not surface — is preserved as extras. */
export const BUILD_OWNED_PATHS = [
  'library_path',
  'library_params.prediction_model',
  'library_params.auto_detect_frag_bounds',
  'library_params.frag_mz_min',
  'library_params.frag_mz_max',
  'library_params.prec_mz_min',
  'library_params.prec_mz_max',
  'library_params.frag_bounds',
  'calibration_raw_file',
  'fasta_paths',
  'fasta_names',
  'fasta_header_regex_accessions',
  'fasta_header_regex_genes',
  'fasta_header_regex_proteins',
  'fasta_header_regex_organisms',
  'fasta_digest_params.min_length',
  'fasta_digest_params.max_length',
  'fasta_digest_params.min_charge',
  'fasta_digest_params.max_charge',
  'fasta_digest_params.missed_cleavages',
  'fasta_digest_params.specificity',
  'fasta_digest_params.max_var_mods',
  'fasta_digest_params.add_decoys',
  'variable_mods',
  'fixed_mods',
  'include_contaminants',
  'predict_fragments',
  'logging',
]

/** Pioneer stores modifications column-wise: three parallel arrays. */
function modsToJson(mods: ModEntry[]): Json {
  return {
    pattern: mods.map((m) => m.pattern),
    mass: mods.map((m) => num(m.mass, 0)),
    name: mods.map((m) => m.name),
  }
}

function modsFromJson(mm: unknown): ModEntry[] | null {
  if (!isObj(mm) || !Array.isArray(mm.pattern)) return null
  const mass = Array.isArray(mm.mass) ? mm.mass : []
  const name = Array.isArray(mm.name) ? mm.name : []
  return mm.pattern.map((p, i) => {
    const nm = name[i] == null ? '' : String(name[i])
    return {
      pattern: p == null ? '' : String(p),
      label: unimodLabel(nm),
      name: nm,
      mass: mass[i] == null ? '' : String(mass[i]),
    }
  })
}

/** The optional `library_params.frag_bounds` key, or nothing at all.
 *
 *  Only a sloped ceiling produces a key. 'constant' is the absence of one,
 *  which is exactly how Pioneer reads a config written before the key existed.
 */
function fragBoundsJson(s: BuildParams): Json {
  if (s.fragBoundsRule === 'constant') return {}
  if (s.fragBoundsRule === 'custom') {
    return {
      frag_bounds: {
        // The floor is flat on every Thermo method measured, so only the
        // ceiling is offered; low is emitted explicitly so the config records
        // the whole rule rather than half of it.
        low: { slope: 0, intercept: 0 },
        high: {
          slope: num(s.fragCeilingSlope, 2),
          intercept: num(s.fragCeilingIntercept, 10),
        },
      },
    }
  }
  return { frag_bounds: s.fragBoundsRule }
}

export function buildLibJsonBase(s: BuildParams): Json {
  // Keep the JSON preview meaningful before any FASTA has been added.
  const files = s.fastaFiles.length ? s.fastaFiles : [makeFastaRow('/path/to/fasta/file.fasta')]
  return {
    library_path: disp(s.libPath, '/path/to/output/my_library'),
    library_params: {
      prediction_model: s.predictionModel,
      auto_detect_frag_bounds: s.autoDetectFragBounds,
      frag_mz_min: num(s.fragMzMin, 150),
      frag_mz_max: num(s.fragMzMax, 2020),
      prec_mz_min: num(s.precMzMin, 390),
      prec_mz_max: num(s.precMzMax, 1010),
      // Omitted for 'constant': absence is what Pioneer reads as flat bounds,
      // and emitting the name would put a key into every config that did not
      // have one for no change in behaviour.
      ...fragBoundsJson(s),
    },
    // Pioneer's own template carries this key with an empty default, so emit it
    // either way rather than omitting it when unset.
    calibration_raw_file: s.calibrationFile.trim(),
    fasta_paths: files.map((f) => f.path || '/path/to/fasta/file.fasta'),
    fasta_names: files.map((f) => f.name || 'SPECIES'),
    fasta_header_regex_accessions: files.map((f) => f.regex.accessions),
    fasta_header_regex_genes: files.map((f) => f.regex.genes),
    fasta_header_regex_proteins: files.map((f) => f.regex.proteins),
    fasta_header_regex_organisms: files.map((f) => f.regex.organisms),
    fasta_digest_params: {
      min_length: num(s.minLen, 7),
      max_length: num(s.maxLen, 40),
      min_charge: num(s.minCharge, 2),
      max_charge: num(s.maxCharge, 3),
      missed_cleavages: num(s.missedCleav, 1),
      // Written explicitly rather than left to Pioneer's default, so the
      // config records the rule the library was actually built with.
      cleavage_regex: s.cleavageRegex.trim() || DEFAULT_CLEAVAGE,
      specificity: s.digestSpecificity,
      max_var_mods: num(s.maxVarMods, 1),
      add_decoys: s.addDecoys,
    },
    variable_mods: modsToJson(s.variableMods),
    fixed_mods: modsToJson(s.fixedMods),
    include_contaminants: s.includeContaminants,
    predict_fragments: s.predictFragments,
    logging: { debug_console_level: s.debugLogging ? 1 : 0 },
  }
}

export function buildLibJson(s: BuildParams, extra: Json | null): Json {
  const base = buildLibJsonBase(s)
  return extra ? deepMerge(base, extra) : base
}

/** Does this object look like a BuildSpecLib config? */
export function isBuildConfig(obj: unknown): boolean {
  return (
    isObj(obj) &&
    ('fasta_paths' in obj || 'library_path' in obj || 'fasta_digest_params' in obj)
  )
}

export function buildConfigToState(obj: unknown): Partial<BuildParams> | null {
  if (!isBuildConfig(obj) || !isObj(obj)) return null
  const set: Partial<BuildParams> = {}
  const str = (x: unknown): string | undefined =>
    x === null || x === undefined ? undefined : String(x)

  if (str(obj.library_path) !== undefined) set.libPath = str(obj.library_path)
  if (str(obj.calibration_raw_file) !== undefined) {
    set.calibrationFile = str(obj.calibration_raw_file)
  }

  const lp = isObj(obj.library_params) ? obj.library_params : {}
  if ('auto_detect_frag_bounds' in lp) {
    set.autoDetectFragBounds = !!lp.auto_detect_frag_bounds
  }
  if (lp.frag_mz_min != null) set.fragMzMin = String(lp.frag_mz_min)
  if (lp.frag_mz_max != null) set.fragMzMax = String(lp.frag_mz_max)
  if (lp.prec_mz_min != null) set.precMzMin = String(lp.prec_mz_min)
  if (lp.prec_mz_max != null) set.precMzMax = String(lp.prec_mz_max)
  // Absent means flat bounds, which is what 'constant' represents here.
  const fb = lp.frag_bounds
  if (fb === undefined || fb === null) {
    set.fragBoundsRule = 'constant'
  } else if (typeof fb === 'string') {
    const name = fb.trim().toLowerCase()
    set.fragBoundsRule =
      name === 'thermo_auto' || name === 'thermo_auto_documented' ? name : 'constant'
  } else if (isObj(fb) && isObj(fb.high)) {
    // Explicit coefficients round-trip as Custom, even when they happen to
    // match a preset: the config said coefficients, so the form should too.
    set.fragBoundsRule = 'custom'
    if (fb.high.slope != null) set.fragCeilingSlope = String(fb.high.slope)
    if (fb.high.intercept != null) set.fragCeilingIntercept = String(fb.high.intercept)
  }

  const paths = Array.isArray(obj.fasta_paths) ? obj.fasta_paths : []
  const names = Array.isArray(obj.fasta_names) ? obj.fasta_names : []
  const ra = Array.isArray(obj.fasta_header_regex_accessions) ? obj.fasta_header_regex_accessions : []
  const rg = Array.isArray(obj.fasta_header_regex_genes) ? obj.fasta_header_regex_genes : []
  const rp = Array.isArray(obj.fasta_header_regex_proteins) ? obj.fasta_header_regex_proteins : []
  const ro = Array.isArray(obj.fasta_header_regex_organisms) ? obj.fasta_header_regex_organisms : []

  if (paths.length) {
    set.fastaFiles = paths.map((p, i) => {
      const path = String(p)
      const regex = {
        accessions: ra[i] == null ? '' : String(ra[i]),
        genes: rg[i] == null ? '' : String(rg[i]),
        proteins: rp[i] == null ? '' : String(rp[i]),
        organisms: ro[i] == null ? '' : String(ro[i]),
      }
      const hasRegex = !!(regex.accessions || regex.genes || regex.proteins || regex.organisms)
      const final = hasRegex ? regex : presetRegex('uniprot')
      return {
        path,
        name: names[i] == null ? makeFastaRow(path).name : String(names[i]),
        presetId: matchPreset(final),
        regex: final,
        auto: false,
        showRegex: false,
      }
    })
  }

  const d = isObj(obj.fasta_digest_params) ? obj.fasta_digest_params : {}
  if (str(d.min_length) !== undefined) set.minLen = str(d.min_length)
  if (str(d.max_length) !== undefined) set.maxLen = str(d.max_length)
  if (str(d.min_charge) !== undefined) set.minCharge = str(d.min_charge)
  if (str(d.max_charge) !== undefined) set.maxCharge = str(d.max_charge)
  if (str(d.missed_cleavages) !== undefined) set.missedCleav = str(d.missed_cleavages)
  if (typeof d.cleavage_regex === 'string' && d.cleavage_regex.trim()) {
    set.cleavageRegex = d.cleavage_regex.trim()
  }
  const specificity = str(d.specificity)
    ?.trim()
    .toLowerCase()
    .replace('_', '-')
  if (specificity && ['full', 'semi', 'semi-n', 'semi-c'].includes(specificity)) {
    set.digestSpecificity = specificity as BuildParams['digestSpecificity']
  }
  if (str(d.max_var_mods) !== undefined) set.maxVarMods = str(d.max_var_mods)
  if ('add_decoys' in d) set.addDecoys = !!d.add_decoys

  if ('include_contaminants' in obj) set.includeContaminants = !!obj.include_contaminants
  if ('predict_fragments' in obj) set.predictFragments = !!obj.predict_fragments

  const vmods = modsFromJson(obj.variable_mods)
  if (vmods) set.variableMods = vmods
  const fmods = modsFromJson(obj.fixed_mods)
  if (fmods) set.fixedMods = fmods

  return set
}

/** Does this object look like a SearchDIA config? */
export function isSearchConfig(obj: unknown): boolean {
  return (
    isObj(obj) &&
    ('paths' in obj || 'global' in obj || 'maxLFQ' in obj || 'proteinScoring' in obj)
  )
}

/** Map a parsed config.json onto form fields. Returns null if unrecognized. */
export function searchConfigToState(obj: unknown): Partial<SearchParams> | null {
  if (!isSearchConfig(obj) || !isObj(obj)) return null
  const set: Partial<SearchParams> = {}
  const str = (x: unknown): string | undefined =>
    x === null || x === undefined ? undefined : String(x)

  const p = isObj(obj.paths) ? obj.paths : {}
  if (str(p.ms_data) !== undefined) set.msData = str(p.ms_data)
  if (str(p.library) !== undefined) set.library = str(p.library)
  if (str(p.results) !== undefined) set.results = str(p.results)

  const o = isObj(obj.output) ? obj.output : {}
  if ('write_csv' in o) set.writeCsv = !!o.write_csv
  if ('write_decoys' in o) set.writeDecoys = !!o.write_decoys
  if ('delete_temp' in o) set.deleteTemp = !!o.delete_temp

  if (isObj(obj.global) && obj.global.q_value_threshold != null) {
    set.qValue = String(obj.global.q_value_threshold)
  }
  // Guarded on presence, like every other optional key here: a config written
  // before this field existed leaves the toggle as it stands rather than being
  // read as off, which would be the wrong reading — Pioneer's fallback is on.
  if (isObj(obj.global) && 'match_between_runs' in obj.global) {
    set.matchBetweenRuns = !!obj.global.match_between_runs
  }
  if (isObj(obj.search) && obj.search.n_isotopes != null) {
    set.nIsotopes = String(obj.search.n_isotopes)
  }
  if (isObj(obj.acquisition) && obj.acquisition.nce != null) {
    set.nce = String(obj.acquisition.nce)
  }
  if (isObj(obj.proteinScoring) && obj.proteinScoring.min_peptides != null) {
    set.minPeptides = String(obj.proteinScoring.min_peptides)
  }
  if (isObj(obj.logging) && 'debug_console_level' in obj.logging) {
    set.debugLogging = Number(obj.logging.debug_console_level) > 0
  }
  if (isObj(obj.maxLFQ) && 'run_to_run_normalization' in obj.maxLFQ) {
    set.runToRunNorm = !!obj.maxLFQ.run_to_run_normalization
  }
  return set
}

// ---------------------------------------------------------------------------
// ConvertRAW
// ---------------------------------------------------------------------------

/** Build PioneerConverter's argv.
 *
 *  Unlike the two Julia commands this has no params file — the converter takes
 *  a positional RAW path plus flags:
 *
 *    PioneerConverter RAW_PATH [-o DIR] [--skip-existing]
 *                     [-n CONCURRENT] [-t PER_FILE] [-b BATCH] [--scan-chunk-size N]
 *
 *  Flags equal to the converter's own defaults are still emitted, so the logged
 *  command line is an exact, re-runnable record of what was executed.
 */
export function buildConvertArgs(s: ConvertParams): string[] {
  const args: string[] = [s.input.trim()]
  if (s.outputDir.trim()) args.push('--output-dir', s.outputDir.trim())
  if (s.skipExisting) args.push('--skip-existing')
  // One file at a time, split across `threads` scan readers. The converter can
  // work on several files concurrently, but exposing both knobs meant the two
  // multiplied and it was easy to oversubscribe the machine without noticing.
  // Pinned to 1 explicitly rather than left to the converter's own default.
  args.push('--concurrent-files', '1')
  args.push('--threads-per-file', String(Math.max(1, parseInt(s.threadsPerFile, 10) || 1)))
  args.push('--batch-size', s.batchSize.trim())
  args.push('--scan-chunk-size', s.scanChunkSize.trim())
  return args
}

/** Argv for DownloadSpecLib.
 *
 *  `--dest` is required by the binary and by validateDownloadRun, so it is
 *  always present here; there is deliberately no default destination. */
export function buildDownloadArgs(s: DownloadParams): string[] {
  const args: string[] = [s.selected.trim(), '--dest', s.dest.trim()]
  if (s.force) args.push('--force')
  return args
}

/** The command line as a user would type it, for the preview panel. */
export function downloadCommandLine(s: DownloadParams): string {
  const quote = (a: string) => (/[\s"']/.test(a) ? JSON.stringify(a) : a)
  return ['DownloadSpecLib', ...buildDownloadArgs(s).map(quote)].join(' ')
}

/** The command line as a user would type it, for the preview panel. */
export function convertCommandLine(s: ConvertParams): string {
  const quote = (a: string) => (/[\s"']/.test(a) ? JSON.stringify(a) : a)
  return ['PioneerConverter', ...buildConvertArgs(s).map(quote)].join(' ')
}
