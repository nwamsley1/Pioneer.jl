export type CommandId = 'searchdia' | 'buildspeclib' | 'convertraw'

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
  traceMode: 'combined' | 'separated'
  minPeptides: string
  runToRunNorm: boolean
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
  traceMode: 'combined',
  minPeptides: '1',
  runToRunNorm: false,
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
  minLen: string
  maxLen: string
  minCharge: string
  maxCharge: string
  missedCleav: string
  maxVarMods: string
  addDecoys: boolean
  includeContaminants: boolean
  predictFragments: boolean
  variableMods: ModEntry[]
  fixedMods: ModEntry[]
}

/** Defaults match Pioneer's own simplified template from GetBuildLibParams. */
export const BUILD_DEFAULTS: BuildParams = {
  fastaFiles: [],
  libPath: '',
  predictionModel: 'altimeter',
  calibrationFile: '',
  minLen: '7',
  maxLen: '40',
  minCharge: '2',
  maxCharge: '3',
  missedCleav: '1',
  maxVarMods: '1',
  addDecoys: true,
  includeContaminants: true,
  predictFragments: true,
  variableMods: [{ pattern: 'M', label: 'Oxidation', name: 'Unimod:35', mass: '15.99491' }],
  fixedMods: [{ pattern: 'C', label: 'Carbamidomethyl', name: 'Unimod:4', mass: '57.021464' }],
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
  batchSize: string
  scanChunkSize: string
}

export const CONVERT_DEFAULTS: ConvertParams = {
  inputMode: 'folder',
  input: '',
  outputDir: '',
  skipExisting: false,
  batchSize: '10000',
  scanChunkSize: '128',
}

/** Mirrors the Rust `runner::Invocation` enum. */
export type Invocation =
  | { kind: 'paramsFile'; json: string }
  | { kind: 'args'; args: string[] }

export type JobStatus = 'queued' | 'running' | 'done' | 'failed' | 'cancelled'

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

export interface Job {
  id: string
  cmd: CommandId
  title: string
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
