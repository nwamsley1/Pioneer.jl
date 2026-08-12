/** Thin typed wrappers over the Rust commands and events. */
import { invoke } from '@tauri-apps/api/core'
import { listen, type UnlistenFn } from '@tauri-apps/api/event'
import { open, save } from '@tauri-apps/plugin-dialog'
import { homeDir } from '@tauri-apps/api/path'

import {
  EMPTY_PATH_INFO,
  type CommandId,
  type Invocation,
  type PathInfo,
  type PioneerInfo,
} from './types'

export interface LineEvent {
  job_id: string
  line: string
  stream: 'stdout' | 'stderr'
  /** The line ended in a carriage return: the writer means to overwrite it
   *  (a progress bar redrawing), so replace the previous line rather than
   *  appending. */
  transient: boolean
  /** Committed lines this one replaces, from an `ESC[<n>A` cursor-up. */
  overwrite: number
}

export interface ExitEvent {
  job_id: string
  code: number | null
  success: boolean
  cancelled: boolean
  message: string
}

export const pioneerInfo = (): Promise<PioneerInfo> => invoke('pioneer_info')

export const cpuCount = (): Promise<number> => invoke('cpu_count')

export async function inspectPath(path: string): Promise<PathInfo> {
  if (!path.trim()) return EMPTY_PATH_INFO
  try {
    return await invoke<PathInfo>('inspect_path', { path })
  } catch {
    return EMPTY_PATH_INFO
  }
}

export const readConfig = (path: string): Promise<string> => invoke('read_config', { path })

export interface Started {
  params_path: string
  /** The environment the backend actually set, for the log header. */
  env_summary: string
}

export const startJob = (
  jobId: string,
  command: CommandId,
  invocation: Invocation,
  threads: number,
): Promise<Started> => invoke('start_job', { jobId, command, invocation, threads })

export const cancelJob = (jobId: string): Promise<void> => invoke('cancel_job', { jobId })

export const onJobLine = (cb: (e: LineEvent) => void): Promise<UnlistenFn> =>
  listen<LineEvent>('job-line', (e) => cb(e.payload))

export const onJobExit = (cb: (e: ExitEvent) => void): Promise<UnlistenFn> =>
  listen<ExitEvent>('job-exit', (e) => cb(e.payload))

/** Where pickers open.
 *
 *  Without a `defaultPath` the native dialog opens wherever the OS decides,
 *  which in a packaged build is the install directory — never a useful place
 *  to start.
 *
 *  Resolution order:
 *    1. where the last picker landed *this session*
 *    2. the configured default directory, if one has been set
 *    3. the home directory
 *
 *  The session part is deliberately not persisted: each launch should start
 *  from the default rather than wherever you happened to finish last time.
 *  Within a session it still tracks you, because MS data, library and results
 *  for one experiment normally live near each other.
 */
const DEFAULT_DIR_KEY = 'pioneerConsole.defaultDir'

/** Session-scoped, so it resets on launch. */
let sessionDir: string | undefined

let homeDirCache: string | undefined

/** Resolve and cache the home directory. Called once at startup. */
export async function initHomeDir(): Promise<void> {
  try {
    homeDirCache = await homeDir()
  } catch {
    /* leave undefined; the dialog falls back to the OS default */
  }
}

/** The configured default browse directory, or '' when unset. */
export function defaultDir(): string {
  try {
    return localStorage.getItem(DEFAULT_DIR_KEY) || ''
  } catch {
    return ''
  }
}

export function setDefaultDir(dir: string): void {
  try {
    if (dir) localStorage.setItem(DEFAULT_DIR_KEY, dir)
    else localStorage.removeItem(DEFAULT_DIR_KEY)
  } catch {
    /* private mode — the setting just will not stick */
  }
}

function lastDir(): string | undefined {
  return sessionDir || defaultDir() || homeDirCache
}

/** Let the user choose the directory every picker starts from. */
export async function pickDefaultDir(): Promise<string | null> {
  const picked = await open({
    directory: true,
    multiple: false,
    title: 'Choose the folder Pioneer starts browsing from',
    defaultPath: defaultDir() || homeDirCache,
  })
  if (typeof picked !== 'string') return null
  setDefaultDir(picked)
  return picked
}

/** Remember where a pick landed. Files remember their parent folder. */
function rememberDir(picked: string, isDirectory: boolean): void {
  const sep = picked.includes('\\') && !picked.includes('/') ? '\\' : '/'
  const dir = isDirectory ? picked : picked.slice(0, picked.lastIndexOf(sep))
  if (dir) sessionDir = dir
}

/** The spectral library picked last, remembered across sessions.
 *
 *  Unlike the general browse location this one is persisted and field-specific.
 *  A library is reused across many searches while data and results folders
 *  change every time, so it is the one path worth defaulting to by name rather
 *  than by neighbourhood.
 */
const LAST_LIBRARY_KEY = 'pioneerConsole.lastLibrary'

export function lastLibrary(): string {
  try {
    return localStorage.getItem(LAST_LIBRARY_KEY) || ''
  } catch {
    return ''
  }
}

/** Library picker. Starts at the library used last if it is still there, and
 *  falls back to the ordinary browse location if it has been moved or deleted. */
export async function pickLibrary(title: string): Promise<string | null> {
  const previous = lastLibrary()
  const start = previous && (await inspectPath(previous)).exists ? previous : lastDir()
  const picked = await open({ directory: true, multiple: false, title, defaultPath: start })
  if (typeof picked !== 'string') return null
  rememberDir(picked, true)
  try {
    localStorage.setItem(LAST_LIBRARY_KEY, picked)
  } catch {
    /* private mode — it just will not be remembered */
  }
  return picked
}

/** Native folder picker. Returns null when the user cancels. */
export async function pickFolder(title: string): Promise<string | null> {
  const picked = await open({ directory: true, multiple: false, title, defaultPath: lastDir() })
  if (typeof picked !== 'string') return null
  rememberDir(picked, true)
  return picked
}

/** Where to write a new library.
 *
 *  A save dialog, not a folder picker. `pickFolder` can only return a
 *  directory that already exists, so naming a new library in the place you
 *  want it was impossible -- you had to pick the parent and then edit the path
 *  by hand.
 *
 *  Returns the path with `.poin` appended when the user did not type it. The
 *  extension is not left to the dialog's filter: on macOS the filter is a
 *  suggestion the user can override, and Pioneer expects the suffix.
 */
export async function pickLibraryTarget(title: string): Promise<string | null> {
  const picked = await save({
    title,
    defaultPath: lastDir(),
    filters: [{ name: 'Pioneer library', extensions: ['poin'] }],
  })
  if (typeof picked !== 'string' || !picked) return null
  rememberDir(picked, false)
  return picked.toLowerCase().endsWith('.poin') ? picked : `${picked}.poin`
}

/** FASTA picker. Returns [] when cancelled. */
export async function pickFastaFiles(multiple = true): Promise<string[]> {
  const picked = await open({
    directory: false,
    multiple,
    title: multiple ? 'Choose FASTA files' : 'Choose a FASTA file',
    filters: [{ name: 'FASTA', extensions: ['fasta', 'fa', 'faa', 'fna', 'fas', 'gz'] }],
    defaultPath: lastDir(),
  })
  if (!picked) return []
  const paths = Array.isArray(picked) ? picked : [picked]
  if (paths.length) rememberDir(paths[0], false)
  return paths
}

/** Native file picker restricted to `extensions`. */
export async function pickFile(
  title: string,
  name: string,
  extensions: string[],
): Promise<string | null> {
  const picked = await open({
    directory: false,
    multiple: false,
    title,
    filters: [{ name, extensions }],
    defaultPath: lastDir(),
  })
  if (typeof picked !== 'string') return null
  rememberDir(picked, false)
  return picked
}
