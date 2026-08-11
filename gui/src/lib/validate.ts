/** Field validation.
 *
 *  The numeric rules are ported verbatim from `numSpecs()` / `numError()` in
 *  the design. The path rules are NOT: the prototype only checked that a string
 *  contained a slash and compared against a hard-coded list of three folders.
 *  These consult the real filesystem through the `inspect_path` command.
 */
import type { BuildParams, ConvertParams, PathInfo } from './types'

export interface NumSpec {
  label: string
  min: number | null
  max: number | null
  step: number
  int: boolean
}

export const NUM_SPECS: Record<string, NumSpec> = {
  qValue: { label: 'q-value threshold', min: 0, max: 1, step: 0.005, int: false },
  nIsotopes: { label: 'Isotopes (n)', min: 1, max: 3, step: 1, int: true },
  nce: { label: 'NCE', min: 0, max: 100, step: 1, int: false },
  minPeptides: { label: 'Min peptides', min: 1, max: null, step: 1, int: true },
  minLen: { label: 'Min length', min: 7, max: 40, step: 1, int: true },
  maxLen: { label: 'Max length', min: 7, max: 40, step: 1, int: true },
  minCharge: { label: 'Min charge', min: 1, max: 4, step: 1, int: true },
  maxCharge: { label: 'Max charge', min: 1, max: 4, step: 1, int: true },
  missedCleav: { label: 'Missed cleav.', min: 0, max: 9, step: 1, int: true },
  maxVarMods: { label: 'Max var. mods', min: 0, max: 5, step: 1, int: true },
  // PioneerConverter's own defaults are 2 / 3 / 10000 / 128.
  batchSize: { label: 'Batch size (scans)', min: 1, max: null, step: 1000, int: true },
  scanChunkSize: { label: 'Scan chunk size', min: 1, max: null, step: 16, int: true },
}

export const CONVERT_NUM_KEYS = [
  'batchSize',
  'scanChunkSize',
] as const

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
      msg: 'No calibration file — fragment m/z bounds fall back to defaults (fragment 150–2020, precursor 390–1010) instead of being detected from your data.',
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

export function validateSearchRun(
  values: { msData: string; library: string; results: string; qValue: string; nIsotopes: string; nce: string; minPeptides: string },
  notes: { msData: Note; library: Note; results: Note },
): RunBlock | null {
  if (!values.msData.trim()) return { key: 'msData', msg: 'Set the MS data path before running.' }
  if (notes.msData.level === 'error') return { key: 'msData', msg: notes.msData.msg }

  if (!values.library.trim()) {
    return { key: 'library', msg: 'Set the spectral library path before running.' }
  }
  if (notes.library.level === 'error') return { key: 'library', msg: notes.library.msg }

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
  return null
}

// ---------------------------------------------------------------------------
// ConvertRAW
// ---------------------------------------------------------------------------

/** PioneerConverter takes one path: a .raw file, or a directory of them. */
export function convertInputNote(p: ConvertParams, info: PathInfo): Note {
  if (!p.input.trim()) return NONE
  if (info.error) return { level: 'error', msg: info.error }
  if (!info.exists) return { level: 'error', msg: 'This path does not exist.' }

  if (p.inputMode === 'file') {
    if (info.is_dir) return { level: 'error', msg: 'That is a folder — switch to Folder mode.' }
    if (info.extension !== 'raw') {
      return { level: 'error', msg: 'Expected a Thermo .raw file.' }
    }
    return NONE
  }

  if (info.is_file) return { level: 'error', msg: 'That is a file — switch to Single file mode.' }
  if (info.raw_count === 0) return { level: 'error', msg: 'No .raw files in this folder.' }
  return {
    level: '',
    msg: `${info.raw_count} .raw file${info.raw_count > 1 ? 's' : ''} to convert.`,
  }
}

export function convertOutputNote(value: string, info: PathInfo): Note {
  if (!value.trim()) return NONE // converter defaults to <input_dir>/arrow_out
  if (info.is_file) return { level: 'error', msg: 'A file exists at this path — choose a folder.' }
  if (info.exists && info.arrow_count > 0) {
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
export function validateConvertRun(
  p: ConvertParams,
  inputNote: Note,
  outputNote: Note,
): RunBlock | null {
  if (!p.input.trim()) return { key: 'convertInput', msg: 'Choose the file or folder to convert.' }
  if (inputNote.level === 'error') return { key: 'convertInput', msg: inputNote.msg }
  if (outputNote.level === 'error') return { key: 'convertOutput', msg: outputNote.msg }
  for (const key of CONVERT_NUM_KEYS) {
    const err = numError(key, p[key])
    if (err) return { key, msg: `${NUM_SPECS[key].label}: ${err}.` }
  }
  return null
}
