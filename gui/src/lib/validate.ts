/** Field validation.
 *
 *  The numeric rules are ported verbatim from `numSpecs()` / `numError()` in
 *  the design. The path rules are NOT: the prototype only checked that a string
 *  contained a slash and compared against a hard-coded list of three folders.
 *  These consult the real filesystem through the `inspect_path` command.
 */
import type {
  BuildParams,
  ConvertParams,
  DownloadParams,
  ModEntry,
  PathInfo,
} from './types'

export interface NumSpec {
  label: string
  min: number | null
  max: number | null
  step: number
  int: boolean
  /** Shown on an info dot beside the label, for fields whose name cannot carry
   *  the whole story. */
  info?: string
}

export const NUM_SPECS: Record<string, NumSpec> = {
  qValue: { label: 'q-value threshold', min: 0, max: 1, step: 0.005, int: false },
  nIsotopes: {
    label: 'Fragment isotopes',
    min: 1,
    max: 3,
    step: 1,
    int: true,
    info: 'How many isotope peaks of each fragment are considered when matching. Fragment, not precursor: the parameter is search.n_isotopes, flattened from search.fragment_settings. More costs time; 2 suits most data.',
  },
  nce: {
    label: 'Initial NCE',
    min: 0,
    max: 100,
    step: 1,
    int: false,
    info: 'The starting guess for normalized collision energy. Pioneer refines it during the parameter tuning search, so it does not have to be exact — it only has to be close enough for tuning to converge.',
  },
  minPeptides: {
    label: 'Min peptides',
    min: 1,
    max: null,
    step: 1,
    int: true,
    info: 'How many distinct peptides a protein group needs before it is reported. Counted per run, not across the experiment: a protein seen by two peptides in one file and one in another is kept for the first and dropped from the second at a threshold of 2.',
  },
  fragMzMin: { label: 'Fragment m/z min', min: 0, max: null, step: 10, int: false },
  fragMzMax: { label: 'Fragment m/z max', min: 0, max: null, step: 10, int: false },
  precMzMin: { label: 'Precursor m/z min', min: 0, max: null, step: 10, int: false },
  precMzMax: { label: 'Precursor m/z max', min: 0, max: null, step: 10, int: false },
  fragCeilingSlope: { label: 'Slope', min: 0, max: null, step: 0.01, int: false },
  fragCeilingIntercept: { label: 'Intercept', min: -1000, max: null, step: 1, int: false },
  minLen: { label: 'Min length', min: 7, max: 40, step: 1, int: true },
  maxLen: { label: 'Max length', min: 7, max: 40, step: 1, int: true },
  minCharge: { label: 'Min charge', min: 1, max: 4, step: 1, int: true },
  maxCharge: { label: 'Max charge', min: 1, max: 4, step: 1, int: true },
  missedCleav: { label: 'Missed cleav.', min: 0, max: 9, step: 1, int: true },
  maxVarMods: { label: 'Max var. mods', min: 0, max: 5, step: 1, int: true },
  // PioneerConverter's own defaults are 2 / 3 / 10000 / 128.
  threadsPerFile: { label: 'Threads per file', min: 1, max: null, step: 1, int: true },
  batchSize: { label: 'Batch size (scans)', min: 1, max: null, step: 1000, int: true },
  scanChunkSize: { label: 'Scan chunk size', min: 1, max: null, step: 16, int: true },
  // convertMzML's own default is 2.
  concurrentFiles: {
    label: 'Files at a time',
    min: 1,
    max: null,
    step: 1,
    int: true,
    info: 'How many .mzML files are converted at the same time. Each one is read and written whole, so this is bounded by disk and memory rather than by cores — the default of 2 is a deliberately conservative starting point.',
  },
}

/** The numeric fields PioneerConverter reads. */
export const CONVERT_NUM_KEYS = [
  'threadsPerFile',
  'batchSize',
  'scanChunkSize',
] as const

/** The numeric fields convertMzML reads. Deliberately disjoint from
 *  CONVERT_NUM_KEYS: the two converters share no tuning parameter, and
 *  validating a field the running binary ignores would block a valid run. */
export const MZML_NUM_KEYS = ['concurrentFiles'] as const

export const DIGEST_KEYS = [
  'minLen',
  'maxLen',
  'minCharge',
  'maxCharge',
  'missedCleav',
  'maxVarMods',
] as const

export function numError(key: string, raw: string): string {
  const sp = NUM_SPECS[key]
  if (!sp) return ''
  const v = (raw ?? '').trim()
  if (v === '') return 'required'
  if (/[^0-9.eE+-]/.test(v) || !Number.isFinite(Number(v))) return 'must be a number'
  const n = Number(v)
  if (sp.int && !Number.isInteger(n)) return 'whole number'
  if (sp.min != null && n < sp.min) return `min ${sp.min}`
  if (sp.max != null && n > sp.max) return `max ${sp.max}`
  return ''
}

/** A field-level note. `error` blocks the run; `warn` does not. */
export interface Note {
  level: 'error' | 'warn' | ''
  msg: string
}

const NONE: Note = { level: '', msg: '' }

const MS_EXTENSIONS = ['raw', 'mzml', 'arrow']

export function msDataNote(value: string, info: PathInfo): Note {
  if (!value.trim()) return NONE
  if (info.error) return { level: 'error', msg: info.error }
  if (!info.exists) return { level: 'error', msg: 'This path does not exist.' }
  if (info.is_file) {
    return MS_EXTENSIONS.includes(info.extension)
      ? NONE
      : {
          level: 'error',
          msg: 'Unsupported type — Pioneer reads .raw, .mzML or .arrow (a folder is fine too).',
        }
  }
  if (info.ms_file_count === 0) {
    return { level: 'error', msg: 'No .raw, .mzML or .arrow files in this folder.' }
  }
  if (info.raw_count > 0 && info.arrow_count === 0) {
    // SearchDIA reads Arrow; .raw still needs converting first.
    return {
      level: 'warn',
      msg: `${info.raw_count} .raw file${info.raw_count > 1 ? 's' : ''} found and no .arrow — convert them with ConvertRAW first.`,
    }
  }
  return NONE
}

/**
 * `libraryNote`, except that the path a queued build or download is about to
 * produce is reported as pending rather than missing.
 *
 * Without this the field contradicts the button: `validateSearchRun` lets the
 * search be queued behind the job that makes its library, while the note under
 * the field still calls the path an error.
 */
export function pendingLibraryNote(
  value: string,
  info: PathInfo,
  pendingLibrary: string | null,
): Note {
  const note = libraryNote(value, info)
  if (note.level !== 'error') return note
  if (pendingLibrary === null || value.trim() !== pendingLibrary) return note
  if (info.exists) return note
  return { level: 'warn', msg: 'Being built by a queued run — this search will wait for it.' }
}

/** A spectral library is a *directory* of Arrow/JLD2/JLS tables, not a single
 *  file, so the design's extension test is not sufficient: an empty folder
 *  named `lib.poin` passes it and then fails deep inside Pioneer with
 *  "Fragment file not found".
 *
 *  The contents are also the more reliable signal than the name. The canonical
 *  extension is `.poin` (BuildSpecLib appends it; the design's `.pion` was a
 *  typo), but libraries named `.pion` exist in the wild — so the marker files
 *  decide, and an unexpected extension is only a warning. */
export function libraryNote(value: string, info: PathInfo): Note {
  if (!value.trim()) return NONE
  if (info.error) return { level: 'error', msg: info.error }
  if (!info.exists) return { level: 'error', msg: 'This path does not exist.' }
  if (info.is_file) {
    return { level: 'error', msg: 'A spectral library is a folder, not a file.' }
  }
  if (!info.is_pion_library) {
    return {
      level: 'error',
      msg: 'This folder is not a spectral library — no precursors_table or detailed_fragments inside.',
    }
  }
  if (info.extension !== 'poin') {
    return {
      level: 'warn',
      msg: `Spectral libraries are normally named .poin — this one ends in .${info.extension || '(none)'}.`,
    }
  }
  return NONE
}

export function resultsNote(value: string, info: PathInfo): Note {
  if (!value.trim()) return NONE
  if (info.error) return { level: 'error', msg: info.error }
  if (!info.exists) return NONE // created on demand
  if (info.is_file) {
    return { level: 'error', msg: 'A file exists at this path — choose a folder.' }
  }
  if (info.entry_count > 0) {
    return {
      level: 'warn',
      msg: "This folder already exists and isn't empty — existing files may be overwritten.",
    }
  }
  return NONE
}

const FASTA_EXTENSIONS = ['fasta', 'fa', 'faa', 'fna', 'fas']

export function fastaNote(value: string, info: PathInfo): Note {
  if (!value.trim()) return NONE
  if (info.error) return { level: 'error', msg: info.error }
  if (!info.exists) return { level: 'error', msg: 'This file does not exist.' }
  if (info.is_dir) {
    return { level: 'error', msg: 'Choose a FASTA file, not a folder.' }
  }
  // `extension_of` in paths.rs looks through a trailing `.gz`, so a
  // `proteins.fasta.gz` reports "fasta" and passes. It did not always: the
  // backend returned "gz" and every compressed FASTA was rejected here, even
  // though the picker offered them and Pioneer reads them.
  if (info.extension && !FASTA_EXTENSIONS.includes(info.extension)) {
    return { level: 'error', msg: 'Expected a FASTA file (.fasta, .fa, .faa…).' }
  }
  return NONE
}

/** The calibration file is optional, but Pioneer warns when it is missing and
 *  falls back to fixed m/z bounds, so an absent one is surfaced as a standing
 *  warning rather than left silent. */
export function calibrationNote(value: string, info: PathInfo): Note {
  if (!value.trim()) {
    return {
      level: 'warn',
      msg: 'No reference MS file — m/z bounds fall back to defaults (fragment 150–2020, precursor 390–1010) instead of being detected from your data.',
    }
  }
  if (info.error) return { level: 'error', msg: info.error }
  if (!info.exists) return { level: 'error', msg: 'This file does not exist.' }
  if (info.is_dir) return { level: 'error', msg: 'Choose a single MS data file, not a folder.' }
  if (info.extension && !MS_EXTENSIONS.includes(info.extension)) {
    return { level: 'error', msg: 'Expected a .raw, .mzML or .arrow file.' }
  }
  return NONE
}

/** BuildSpecLib appends `.poin` to `library_path` if it is not already there,
 *  so the directory it will actually create is what we warn about. */
export function libraryTargetPath(libPath: string): string {
  const v = libPath.trim()
  if (!v) return ''
  return v.endsWith('.poin') ? v : `${v}.poin`
}

/** Where a downloaded library will land.
 *
 *  DownloadSpecLib writes `<dest>/<library name>`, and the name already carries
 *  its .poin/.pion extension, so nothing is appended here. */
export function downloadTargetPath(dest: string, library: string): string {
  const d = dest.trim()
  const l = library.trim()
  if (!d || !l) return ''
  return d.endsWith('/') ? `${d}${l}` : `${d}/${l}`
}

/** The file name without its directory or extension: what a per-file search
 *  names its run and its results folder after.
 *
 *  Splits on both separators rather than the platform's own, because a path can
 *  reach the console from a drag-and-drop or a stored run rather than from this
 *  machine's picker. */
export function fileStem(path: string): string {
  const name = path.trim().split(/[\\/]/).pop() ?? ''
  const dot = name.lastIndexOf('.')
  return dot > 0 ? name.slice(0, dot) : name
}

/** Join a directory and one child name, tolerating a trailing separator on the
 *  directory. Same shape as downloadTargetPath, kept separate because that one
 *  is about a specific pair of fields and this one is not. */
export function joinPath(dir: string, child: string): string {
  const d = dir.trim()
  if (!d) return child
  return /[\\/]$/.test(d) ? `${d}${child}` : `${d}/${child}`
}

export function libPathNote(value: string, targetInfo: PathInfo): Note {
  if (!value.trim()) return NONE
  if (!/[\\/]/.test(value)) return { level: 'error', msg: 'Enter a full path.' }
  if (targetInfo.is_file) {
    return { level: 'error', msg: 'A file already exists at this path.' }
  }
  if (targetInfo.exists && targetInfo.entry_count > 0) {
    return {
      level: 'warn',
      msg: `${libraryTargetPath(value)} already exists — it will be overwritten.`,
    }
  }
  return NONE
}

/** The first blocking problem, or null when the run can proceed. */
export interface RunBlock {
  key: string
  msg: string
}

/**
 * @param pendingLibrary Path a queued or running BuildSpecLib/DownloadSpecLib
 *   job will produce, or null when no such job is waiting.
 */
export function validateSearchRun(
  values: {
    msDataMode?: 'folder' | 'files'
    msDataFiles?: string[]
    msData: string
    library: string
    results: string
    qValue: string
    nIsotopes: string
    nce: string
    minPeptides: string
  },
  notes: { msData: Note; library: Note; results: Note },
  pendingLibrary: string | null = null,
): RunBlock | null {
  // In file-list mode the folder field is unused and its note is meaningless:
  // what has to be there is at least one file. Everything below is shared --
  // each fanned-out run gets the same library, the same thresholds, and a
  // results folder under the same root.
  if (values.msDataMode === 'files') {
    if (!(values.msDataFiles ?? []).length) {
      return { key: 'msData', msg: 'Add at least one file to search.' }
    }
  } else {
    if (!values.msData.trim()) return { key: 'msData', msg: 'Set the MS data path before running.' }
    if (notes.msData.level === 'error') return { key: 'msData', msg: notes.msData.msg }
  }

  // A library that is still being predicted is not a missing library. With a
  // build or download ahead of it in the queue, the search that consumes its
  // output can be queued behind it -- the folder does not exist while the
  // search sits in the queue, and does by the time it starts. Only that one
  // path earns the exemption: any other path that does not exist is still the
  // mistake it always was.
  const awaitsPending =
    pendingLibrary !== null && values.library.trim() === pendingLibrary

  if (!awaitsPending) {
    if (!values.library.trim()) {
      return { key: 'library', msg: 'Set the spectral library path before running.' }
    }
    if (notes.library.level === 'error') return { key: 'library', msg: notes.library.msg }
  }

  if (!values.results.trim()) {
    return { key: 'results', msg: 'Set the results path before running.' }
  }
  if (notes.results.level === 'error') return { key: 'results', msg: notes.results.msg }

  for (const key of ['qValue', 'nIsotopes', 'nce', 'minPeptides'] as const) {
    const err = numError(key, values[key])
    if (err) return { key, msg: `${NUM_SPECS[key].label}: ${err}.` }
  }
  return null
}

export function validateBuildRun(
  p: BuildParams,
  fastaNotes: Note[],
  libNote: Note,
  calibNote: Note,
): RunBlock | null {
  // Absent is fine (warned about inline); a path that is set but wrong is not.
  if (calibNote.level === 'error') return { key: 'calibrationFile', msg: calibNote.msg }
  if (!p.fastaFiles.length) {
    return { key: 'fastaAdd', msg: 'Add at least one FASTA file before building.' }
  }
  for (let i = 0; i < p.fastaFiles.length; i++) {
    if (!p.fastaFiles[i].path.trim()) {
      return { key: 'fastaAdd', msg: 'Every FASTA file needs a path.' }
    }
    if (fastaNotes[i]?.level === 'error') {
      return { key: 'fastaAdd', msg: fastaNotes[i].msg }
    }
  }
  if (!p.libPath.trim()) {
    return { key: 'libPath', msg: 'Set the output library path before building.' }
  }
  if (libNote.level === 'error') return { key: 'libPath', msg: libNote.msg }

  for (const key of DIGEST_KEYS) {
    const err = numError(key, p[key])
    if (err) return { key, msg: `${NUM_SPECS[key].label}: ${err}.` }
  }
  if (Number(p.minLen) > Number(p.maxLen)) {
    return { key: 'minLen', msg: 'Min length cannot exceed max length.' }
  }
  if (Number(p.minCharge) > Number(p.maxCharge)) {
    return { key: 'minCharge', msg: 'Min charge cannot exceed max charge.' }
  }
  // Pioneer indexes mods by three parallel arrays, so a blank row would produce
  // a mod with an empty pattern rather than being ignored.
  for (const [kind, mods] of [
    ['fixed', p.fixedMods],
    ['variable', p.variableMods],
  ] as const) {
    for (const m of mods) {
      if (!m.pattern.trim() || !m.name.trim() || !m.mass.trim()) {
        return { key: `${kind}Mods`, msg: `Every ${kind} modification needs a site, Unimod name and mass.` }
      }
      if (!Number.isFinite(Number(m.mass))) {
        return { key: `${kind}Mods`, msg: `"${m.mass}" is not a valid ${kind} modification mass.` }
      }
    }
  }
  const conflict = modSiteConflict(p.fixedMods, p.variableMods)
  if (conflict) return { key: 'variableMods', msg: conflict }
  // Only when they are actually used: with auto-detection on these come from
  // the reference file and whatever is in the fields is ignored.
  if (!p.autoDetectFragBounds) {
    for (const key of ['fragMzMin', 'fragMzMax', 'precMzMin', 'precMzMax'] as const) {
      const err = numError(key, p[key])
      if (err) return { key, msg: `${NUM_SPECS[key].label}: ${err}.` }
    }
    if (Number(p.fragMzMin) >= Number(p.fragMzMax)) {
      return { key: 'fragMzMin', msg: 'Fragment m/z min must be below max.' }
    }
    if (Number(p.precMzMin) >= Number(p.precMzMax)) {
      return { key: 'precMzMin', msg: 'Precursor m/z min must be below max.' }
    }
    if (p.fragBoundsRule === 'custom') {
      for (const key of ['fragCeilingSlope', 'fragCeilingIntercept'] as const) {
        const err = numError(key, p[key])
        if (err) return { key, msg: `Fragment ceiling ${NUM_SPECS[key].label}: ${err}.` }
      }
    }
  }
  return null
}

/** The residues a modification pattern applies to.
 *
 *  Port of `mod_pattern_residues` (check_params.jl): test the pattern against
 *  each amino acid on its own, so a plain `"C"` and a class `"[ST]"` both
 *  resolve, while a multi-residue rule matches nothing and is left to another
 *  check. An unparseable pattern yields no residues, likewise.
 */
export function modPatternResidues(pattern: string): Set<string> {
  const residues = new Set<string>()
  let re: RegExp
  try {
    re = new RegExp(pattern)
  } catch {
    return residues
  }
  for (const aa of 'ACDEFGHIKLMNPQRSTVWY') {
    if (re.test(aa)) residues.add(aa)
  }
  return residues
}

/** Whether a variable modification shares a residue with a fixed one.
 *
 *  Port of `check_mod_site_conflicts`. A fixed modification occupies *every*
 *  matching residue, so a variable one on the same residue describes a peptide
 *  that cannot exist — `fillVarModStrings!` stacks them and the library carries
 *  impossible masses without erroring. Pioneer rejects the config outright, so
 *  catching it here turns a failed build into an inline message.
 *
 *  Any shared residue counts, not just the same modification twice: a fixed
 *  carbamidomethyl on C conflicts with a variable oxidation on C.
 */
/** Residues claimed by a fixed modification and a variable one at once.
 *
 *  The same rule `modSiteConflict` reports, as a set, so a table can mark the
 *  rows involved while they are being edited rather than waiting for Run. Both
 *  sides are implicated: removing either resolves it.
 */
export function conflictingResidues(fixed: ModEntry[], variable: ModEntry[]): Set<string> {
  const shared = new Set<string>()
  const fixedResidues = fixed.flatMap((f) => [...modPatternResidues(f.pattern)])
  for (const v of variable) {
    const vres = modPatternResidues(v.pattern)
    for (const r of fixedResidues) if (vres.has(r)) shared.add(r)
  }
  return shared
}

export function modSiteConflict(fixed: ModEntry[], variable: ModEntry[]): string | null {
  for (const v of variable) {
    const vres = modPatternResidues(v.pattern)
    if (!vres.size) continue
    for (const f of fixed) {
      const shared = [...modPatternResidues(f.pattern)].filter((r) => vres.has(r))
      if (!shared.length) continue
      const vname = v.label || v.name || 'that modification'
      const fname = f.label || f.name || 'a fixed modification'
      return (
        `${vname} is variable on ${shared.sort().join(', ')}, but ${fname} is already ` +
        `fixed there. A fixed modification takes every matching residue, so the ` +
        `variable one would land on top of it. Remove one, or narrow a site.`
      )
    }
  }
  return null
}

// ---------------------------------------------------------------------------
// ConvertRAW
// ---------------------------------------------------------------------------

/** What the Input field is saying about itself.
 *
 *  Folder mode: the format picker, not the path, decides which converter runs,
 *  so a `.raw` folder handed to the mzML converter is caught here rather than
 *  by the binary, whose failure would arrive seconds later as a stack trace.
 *
 *  File-list mode: the format is read off each name instead, so what has to be
 *  checked is that every file is one a converter can actually read.
 */
export function convertInputNote(p: ConvertParams, info: PathInfo): Note {
  if (p.inputMode === 'files') {
    const files = p.inputFiles
    if (files.length === 0) return NONE
    const bad = files.filter((f) => !/\.(raw|mzml)$/i.test(f.trim()))
    if (bad.length) {
      // Named, not counted: with a list on screen, "2 files are not supported"
      // still leaves you hunting for which.
      const names = bad.slice(0, 3).map((f) => f.trim().split(/[\\/]/).pop())
      return {
        level: 'error',
        msg: `Not a .raw or .mzML file: ${names.join(', ')}${bad.length > 3 ? `, and ${bad.length - 3} more` : ''}.`,
      }
    }
    // Two files of the same name cannot share the one folder each converter is
    // handed, so they are refused here rather than by paths::stage_files once
    // the run has already been queued.
    const seen = new Set<string>()
    for (const f of files) {
      const name = (f.trim().split(/[\\/]/).pop() ?? '').toLowerCase()
      if (seen.has(name)) {
        return {
          level: 'error',
          msg: `Two of these files are named ${name}. Rename one, or convert them separately.`,
        }
      }
      seen.add(name)
    }
    return NONE
  }

  if (!p.input.trim()) return NONE
  if (info.error) return { level: 'error', msg: info.error }
  if (!info.exists) return { level: 'error', msg: 'This path does not exist.' }

  const mzml = p.format === 'mzml'
  const label = mzml ? '.mzML' : '.raw'
  const count = mzml ? info.mzml_count : info.raw_count
  const other = mzml ? info.raw_count : info.mzml_count
  const otherLabel = mzml ? '.raw' : '.mzML'

  if (info.is_file) return { level: 'error', msg: 'That is a file — switch to Files mode.' }
  if (count === 0) {
    return {
      level: 'error',
      msg: other
        ? `No ${label} files in this folder, but ${other} ${otherLabel} file${other > 1 ? 's' : ''} — switch the format to ${otherLabel}.`
        : `No ${label} files in this folder.`,
    }
  }
  return {
    level: '',
    msg: `${count} ${label} file${count > 1 ? 's' : ''} to convert.`,
  }
}

/** Takes the whole `convert` object rather than just the path, because whether
 *  existing output is at risk depends on "Skip existing" as much as on what is
 *  in the folder. */
export function convertOutputNote(p: ConvertParams, info: PathInfo): Note {
  if (!p.outputDir.trim()) return NONE // converter defaults to <input_dir>/arrow_out
  if (info.is_file) return { level: 'error', msg: 'A file exists at this path — choose a folder.' }
  if (info.exists && info.arrow_count > 0) {
    // Nothing to warn about once the files are being left alone: the warning
    // existed only to offer this setting, and telling someone to enable what
    // they have already enabled reads as the app not knowing its own state.
    if (p.skipExisting) {
      return {
        level: '',
        msg: `This folder already holds ${info.arrow_count} .arrow file${info.arrow_count > 1 ? 's' : ''}. They will be left alone.`,
      }
    }
    return {
      level: 'warn',
      msg: `This folder already holds ${info.arrow_count} .arrow file${info.arrow_count > 1 ? 's' : ''} — they may be overwritten. Enable "Skip existing" to keep them.`,
    }
  }
  return NONE
}

/** The two parallelism knobs multiply, so individually reasonable values can
 *  still ask for more threads than the machine has. The Julia thread picker is
 *  clamped at its control; this product cannot be, so it is enforced here. */
/** Blocks a download that cannot be started.
 *
 *  A destination is required rather than defaulted: a 3 GiB download landing
 *  somewhere the user did not choose is worse than being asked. */
export function validateDownloadRun(p: DownloadParams, targetExists: boolean): RunBlock | null {
  if (!p.selected.trim()) return { key: 'downloadSelected', msg: 'Choose a library to download.' }
  if (!p.dest.trim()) return { key: 'downloadDest', msg: 'Choose where the library should be saved.' }
  // Caught here rather than by the binary: it refuses to overwrite, and that
  // refusal otherwise arrives only after the job has been queued and started.
  if (targetExists && !p.force) {
    return {
      key: 'downloadDest',
      msg: `${p.selected} already exists in that folder. Turn on “Replace an existing copy”, or choose another folder.`,
    }
  }
  return null
}

export function validateConvertRun(
  p: ConvertParams,
  inputNote: Note,
  outputNote: Note,
): RunBlock | null {
  if (p.inputMode === 'files') {
    if (p.inputFiles.length === 0) {
      return { key: 'convertInput', msg: 'Add at least one file to convert.' }
    }
    // Each converter is handed a staging folder under the system temp
    // directory, so its own <input_dir>/arrow_out default would write there --
    // somewhere nobody would think to look. Required rather than guessed at:
    // the files can come from several folders and none of them is the obvious
    // answer.
    if (!p.outputDir.trim()) {
      return { key: 'convertOutput', msg: 'Choose an output folder for the converted files.' }
    }
  } else if (!p.input.trim()) {
    return { key: 'convertInput', msg: 'Choose the folder to convert.' }
  }
  if (inputNote.level === 'error') return { key: 'convertInput', msg: inputNote.msg }
  if (outputNote.level === 'error') return { key: 'convertOutput', msg: outputNote.msg }
  // Only the fields the converters that will actually run read. A list can hold
  // both formats, so it is checked against both; a stale batch size left over
  // from a RAW run must not block an mzML-only conversion that ignores it.
  const keys =
    p.inputMode === 'files'
      ? [
          ...(p.inputFiles.some((f) => /\.raw$/i.test(f.trim())) ? CONVERT_NUM_KEYS : []),
          ...(p.inputFiles.some((f) => /\.mzml$/i.test(f.trim())) ? MZML_NUM_KEYS : []),
        ]
      : p.format === 'mzml'
        ? MZML_NUM_KEYS
        : CONVERT_NUM_KEYS
  for (const key of keys) {
    const err = numError(key, p[key])
    if (err) return { key, msg: `${NUM_SPECS[key].label}: ${err}.` }
  }
  return null
}
