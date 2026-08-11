/** Thin typed wrappers over the Rust commands and events. */
import { invoke } from '@tauri-apps/api/core'
import { listen, type UnlistenFn } from '@tauri-apps/api/event'
import { open } from '@tauri-apps/plugin-dialog'

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

/** Where the last picker landed, so the next one opens there.
 *
 *  Without a `defaultPath` the native dialog opens wherever the OS decides,
 *  which in a packaged build is the install directory — never a useful place to
 *  start. Remembering one location across every picker matches how these are
 *  actually used: MS data, library and results for a given experiment normally
 *  live near each other.
 *
 *  Kept in localStorage rather than component state so it survives a restart,
 *  same as the form parameters.
 */
const LAST_DIR_KEY = 'pioneerConsole.lastDir'

function lastDir(): string | undefined {
  try {
    return localStorage.getItem(LAST_DIR_KEY) || undefined
  } catch {
    return undefined
  }
}

/** Remember where a pick landed. Files remember their parent folder. */
function rememberDir(picked: string, isDirectory: boolean): void {
  const sep = picked.includes('\\') && !picked.includes('/') ? '\\' : '/'
  const dir = isDirectory ? picked : picked.slice(0, picked.lastIndexOf(sep))
  if (!dir) return
  try {
    localStorage.setItem(LAST_DIR_KEY, dir)
  } catch {
    /* private mode / quota — the dialog just opens at the OS default */
  }
}

/** Native folder picker. Returns null when the user cancels. */
export async function pickFolder(title: string): Promise<string | null> {
  const picked = await open({ directory: true, multiple: false, title, defaultPath: lastDir() })
  if (typeof picked !== 'string') return null
  rememberDir(picked, true)
  return picked
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
