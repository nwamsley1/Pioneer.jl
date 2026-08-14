import { DEFAULT_CLEAVAGE } from './enzymes'
import { modEntry } from './koinaMods'

export type CommandId = 'searchdia' | 'buildspeclib' | 'downloadspeclib' | 'convertraw'

/** One library offered by the Hugging Face repository, as reported by
 *  `DownloadSpecLib --list --json`. Mirrors LibraryEntry in catalog.jl — the
 *  contract between the two halves of the feature. */
export interface RemoteLibrary {
  name: string
  title: string
  model: string
  description: string
  recommended_for: string
  total_bytes: number
  size_human: string
  n_files: number
  details: Record<string, string>
}

/** Everything the DownloadSpecLib page owns. */
export interface DownloadParams {
  /** Directory the library is written into. Required; there is no default. */
  dest: string
  /** Selected library name, '' when nothing is picked yet. */
  selected: string
  /** Replace an existing directory at the destination. */
  force: boolean
}

export const DOWNLOAD_DEFAULTS: DownloadParams = {
  dest: '',
  selected: '',
  force: false,
}

/** Every field the SearchDIA form owns. Kept as strings where the design keeps
 *  strings, so a half-typed number ("0.0") survives a re-render intact. */
export interface SearchParams {
  msData: string
  library: string
  results: string
  writeCsv: boolean
  writeDecoys: boolean
  deleteTemp: boolean
  qValue: string
  nIsotopes: string
  nce: string
  minPeptides: string
  runToRunNorm: boolean
  /** global.match_between_runs. Pioneer falls back to true when the key is
   *  absent, and the GUI did not emit it before this field existed — so on is
   *  what every GUI run has done to date, and is the default here. */
  matchBetweenRuns: boolean
  /** logging.debug_console_level: 0 off, 1 on. Pioneer's own default is 0. */
  debugLogging: boolean
}

export const SEARCH_DEFAULTS: SearchParams = {
  msData: '',
  library: '',
  results: '',
  writeCsv: true,
  writeDecoys: false,
  deleteTemp: true,
  qValue: '0.01',
  nIsotopes: '2',
  nce: '26',
  minPeptides: '1',
  runToRunNorm: false,
  matchBetweenRuns: true,
  debugLogging: false,
}

/** One FASTA input row in the BuildSpecLib form. */
export interface FastaEntry {
  path: string
  /** Short species tag written to `fasta_names`, e.g. HUMAN. */
  name: string
  presetId: HeaderPresetId
  regex: HeaderRegex
  /** True while the preset is still being inferred from the filename; a manual
   *  choice pins it and stops re-detection. */
  auto: boolean
  showRegex: boolean
}

export type HeaderPresetId = 'uniprot' | 'ensembl' | 'refseq' | 'custom'

export interface HeaderRegex {
  accessions: string
  genes: string
  proteins: string
  organisms: string
}

/** One fixed or variable modification. Masses stay strings so a half-typed
 *  value survives a re-render. */
export interface ModEntry {
  pattern: string
  label: string
  name: string
  mass: string
}

/** A fragment-prediction model Pioneer can drive. Mirrors MODEL_CONFIGS in
 *  src/Pioneer.jl; `peptideLength` is that entry's `peptide_length`, with null
 *  meaning unconstrained.
 *
 *  Pioneer clamps the digest to this range and warns
 *  (clamp_digest_length_to_model), so a mismatch is never silent — but the
 *  warning only appears once the build is already running. Surfacing the range
 *  here lets the user see it while they are still choosing. */
export interface PredictionModel {
  id: string
  label: string
  note: string
  peptideLength: { min: number; max: number } | null
}

export const PREDICTION_MODELS: PredictionModel[] = [
  {
    id: 'altimeter',
    label: 'Altimeter',
    note: 'Spline coefficients across collision energies; instrument-aware.',
    peptideLength: null,
  },
  {
    id: 'prosit_2020_hcd',
    label: 'Prosit 2020 HCD',
    note: 'Fixed collision energy, no PTM support beyond the defaults.',
    peptideLength: { min: 7, max: 30 },
  },
  {
    id: 'prosit_2024_ptm',
    label: 'Prosit 2024 PTM',
    note: 'Fixed collision energy, HCD, with PTM support.',
    peptideLength: { min: 7, max: 30 },
  },
  {
    id: 'prosit_2025_40ptm',
    label: 'Prosit 2025 40-PTM',
    note: 'As 2024 PTM, with a wider (40) modification vocabulary.',
    peptideLength: { min: 7, max: 30 },
  },
]

/** True for the Prosit families.
 *
 *  Searching a Prosit-predicted library that carries variable modifications is
 *  experimental: Pioneer does not score site-localization confidence, so a
 *  modified residue is placed but the placement is not evidenced. Altimeter
 *  libraries are unaffected -- the caveat is about the combination, which is
 *  why it is asked of the library rather than of the search.
 */
export function isPrositModel(id: string): boolean {
  return id.startsWith('prosit')
}

export function predictionModelById(id: string): PredictionModel {
  return PREDICTION_MODELS.find((m) => m.id === id) ?? PREDICTION_MODELS[0]
}

export interface BuildParams {
  fastaFiles: FastaEntry[]
  libPath: string
  /** Key into PREDICTION_MODELS; emitted as `library_params.prediction_model`. */
  predictionModel: string
  /** Optional MS data file used to auto-detect fragment and precursor m/z
   *  bounds. Without it Pioneer falls back to fixed defaults. */
  calibrationFile: string
  /** library_params.auto_detect_frag_bounds. Off means the four m/z bounds
   *  below are used as given, rather than read from the reference file. */
  autoDetectFragBounds: boolean
  fragMzMin: string
  fragMzMax: string
  precMzMin: string
  precMzMax: string
  /** library_params.frag_bounds. 'constant' emits no key at all, which is the
   *  behaviour of every config written before the key existed. */
  fragBoundsRule: 'constant' | 'thermo_auto_documented' | 'thermo_auto' | 'custom'
  /** Only read when fragBoundsRule is 'custom'. */
  fragCeilingSlope: string
  fragCeilingIntercept: string
  minLen: string
  maxLen: string
  minCharge: string
  maxCharge: string
  missedCleav: string
  /** The digestion rule, as the regex Pioneer takes. Presets in `enzymes.ts`
   *  set it; "Custom" lets it be typed. Stored as the pattern rather than a
   *  preset id so a config round-trips even when it matches nothing. */
  cleavageRegex: string
  /** How many termini must obey that rule. Orthogonal to it: the enzyme says
   *  where cleavage may occur, this says how much of the peptide has to
   *  respect it. */
  digestSpecificity: 'full' | 'semi' | 'semi-n' | 'semi-c'
  maxVarMods: string
  addDecoys: boolean
  includeContaminants: boolean
  predictFragments: boolean
  variableMods: ModEntry[]
  fixedMods: ModEntry[]
  /** logging.debug_console_level: 0 off, 1 on. */
  debugLogging: boolean
}

/** Defaults match Pioneer's own simplified template from GetBuildLibParams. */
export const BUILD_DEFAULTS: BuildParams = {
  fastaFiles: [],
  libPath: '',
  predictionModel: 'altimeter',
  calibrationFile: '',
  // Mirrors assets/example_config/defaultBuildLibParams.json, so an untouched
  // form emits what Pioneer would have defaulted to anyway.
  autoDetectFragBounds: true,
  fragMzMin: '150',
  fragMzMax: '2020',
  precMzMin: '390',
  precMzMax: '1010',
  fragBoundsRule: 'constant',
  // Thermo's documented rule, shown when Custom is first chosen so the fields
  // start somewhere real rather than at zero.
  fragCeilingSlope: '2.0',
  fragCeilingIntercept: '10',
  minLen: '7',
  maxLen: '40',
  minCharge: '2',
  maxCharge: '3',
  missedCleav: '1',
  cleavageRegex: DEFAULT_CLEAVAGE,
  digestSpecificity: 'full',
  maxVarMods: '1',
  addDecoys: true,
  includeContaminants: true,
  predictFragments: true,
  // Built from the catalogue so the starting rows are exactly what the picker
  // would add; hand-written copies drifted (Oxidation was 15.99491).
  variableMods: [modEntry('altimeter', 35)],
  fixedMods: [modEntry('altimeter', 4)],
  debugLogging: false,
}

/** ConvertRAW is a .NET program driven entirely by CLI flags — there is no
 *  params JSON for it. Defaults are PioneerConverter's own. */
export interface ConvertParams {
  /** A single .raw file, or a folder of them. The binary accepts one path. */
  inputMode: 'file' | 'folder'
  input: string
  /** Blank means the converter's default of <input_dir>/arrow_out. */
  outputDir: string
  skipExisting: boolean
  /** Scan-reader threads within the single file being converted.
   *
   *  PioneerConverter parallelises on two levels and the knobs multiply, so
   *  files-at-a-time stays pinned at 1 (see buildConvertArgs) and this is the
   *  only one exposed. It is deliberately not the sidebar thread count: that
   *  drives JULIA_NUM_THREADS, and the converter is a .NET program that never
   *  reads it. */
  threadsPerFile: string
  batchSize: string
  scanChunkSize: string
}

export const CONVERT_DEFAULTS: ConvertParams = {
  inputMode: 'folder',
  input: '',
  outputDir: '',
  skipExisting: false,
  threadsPerFile: '3',
  batchSize: '10000',
  scanChunkSize: '128',
}

/** Mirrors the Rust `runner::Invocation` enum. */
export type Invocation =
  | { kind: 'paramsFile'; json: string }
  | { kind: 'args'; args: string[] }

/** `interrupted` is not reported by a run — it is inferred at startup for a
 *  row still marked queued or running, which can only mean the app went away
 *  before it finished. */
export type JobStatus =
  | 'queued'
  | 'running'
  | 'done'
  | 'failed'
  | 'cancelled'
  | 'interrupted'

/** `app` marks lines the GUI itself wrote (the command line, the params path,
 *  the exit message). Terminal cursor control only ever rewrites the stream it
 *  came from, so these are never overwritten. */
export type LogStream = 'stdout' | 'stderr' | 'app'

export interface LogLine {
  text: string
  stream: LogStream
  /** Ended in a carriage return: the next line from this stream replaces it. */
  transient: boolean
}

/** The form state a run was launched from, kept so the queue and history can put
 *  it back on screen.
 *
 *  Stored as the form's own params rather than the serialized params file:
 *  ConvertRAW's `paramsJson` is a display command line and cannot be parsed back,
 *  and even for the Julia commands round-tripping through the Pioneer JSON would
 *  lose anything the form models but the config does not. */
export type JobSnapshot =
  | { cmd: 'searchdia'; search: SearchParams }
  | { cmd: 'buildspeclib'; build: BuildParams }
  | { cmd: 'downloadspeclib'; download: DownloadParams }
  | { cmd: 'convertraw'; convert: ConvertParams }

export interface Job {
  id: string
  /** Its number in the run history: the next one at creation, kept for life.
   *  Unlike a queue position this never changes, so deleting run 3 leaves runs
   *  4 and 5 as 4 and 5. Zero for entries restored from before numbering. */
  runNo: number
  cmd: CommandId
  /** Generated adjective-noun name, e.g. "brisk-otter". */
  title: string
  snapshot: JobSnapshot
  /** The output path this run writes to — shown next to the title. */
  target: string
  threads: number
  status: JobStatus
  logLines: LogLine[]
  failMsg: string
  /** Serialized params, held so the job can be started when it reaches the
   *  front of the queue rather than at enqueue time. */
  paramsJson: string
  /** How to invoke: a params file for the Julia commands, argv for the .NET
   *  converter, which has no params file. */
  invocation: Invocation
  /** Paths needed to open this run's output in the viewer. */
  viewerPaths: { results: string; msData: string; library: string } | null
  /** Where the params file was written, once the job actually started. */
  paramsPath: string
  /** Unix seconds when the run reached a final status. 0 while it is still
   *  going, or for history written before this was recorded. */
  finishedAt: number
}

/** What the Rust `inspect_path` command reports. */
export interface PathInfo {
  exists: boolean
  is_dir: boolean
  is_file: boolean
  extension: string
  entry_count: number
  ms_file_count: number
  raw_count: number
  mzml_count: number
  arrow_count: number
  has_config_json: boolean
  is_pion_library: boolean
  error: string | null
}

export const EMPTY_PATH_INFO: PathInfo = {
  exists: false,
  is_dir: false,
  is_file: false,
  extension: '',
  entry_count: 0,
  ms_file_count: 0,
  raw_count: 0,
  mzml_count: 0,
  arrow_count: 0,
  has_config_json: false,
  is_pion_library: false,
  error: null,
}

export interface PioneerInfo {
  home: string
  source: string
  executables: string[]
  has_wrapper: boolean
}
