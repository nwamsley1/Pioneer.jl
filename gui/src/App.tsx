import { useCallback, useEffect, useMemo, useRef, useState } from 'react'

import { BuildSpecLibForm } from './components/BuildSpecLibForm'
import { ConfirmDialog } from './components/ConfirmDialog'
import { ConvertRawForm } from './components/ConvertRawForm'
import { DownloadSpecLibForm } from './components/DownloadSpecLibForm'
import { JobNameField } from './components/JobNameField'
import { JsonModal } from './components/JsonModal'
import { LoadConfigModal } from './components/LoadConfigModal'
import { LogDrawer } from './components/LogDrawer'
import { SearchDiaForm } from './components/SearchDiaForm'
import { TopBar } from './components/TopBar'
import { Sidebar } from './components/Sidebar'
import * as backend from './lib/backend'
import {
  buildConvertArgs,
  downloadCommandLine,
  buildDownloadArgs,
  buildLibJson,
  buildSearchJson,
  convertCommandLine,
  buildConfigToState,
  computeExtras,
  extraLeafPaths,
  isBuildConfig,
  isSearchConfig,
  searchConfigToState,
  BUILD_OWNED_PATHS,
  SEARCH_OWNED_PATHS,
  type Json,
} from './lib/config'
import { makeFastaRow, presetRegex } from './lib/fasta'
import { moveQueuedJob } from './lib/queue'
import { libraryOf, recentLibraries } from './lib/recent'
import {
  enforceRequiredMods,
  findMod,
  initialSite,
  modsForModel,
  occupiedResidues,
  siteAllowed,
  unimodId,
} from './lib/koinaMods'
import { SAMPLE_SEQUENCE, cleavageSites, enzymeById } from './lib/enzymes'
import { listenForDrops } from './lib/dragdrop'
import { resolveRunName } from './lib/names'
import { TITLEBAR_H } from './lib/styles'
import { applyTheme, loadTheme, type ThemeId } from './lib/theme'
import {
  BUILD_DEFAULTS,
  CONVERT_DEFAULTS,
  EMPTY_PATH_INFO,
  SEARCH_DEFAULTS,
  DOWNLOAD_DEFAULTS,
  type BuildParams,
  type CommandId,
  type ConvertParams,
  type FastaEntry,
  type HeaderPresetId,
  type Job,
  type JobSnapshot,
  type JobStatus,
  type LogLine,
  type ModEntry,
  type PathInfo,
  type DownloadParams,
  type RemoteLibrary,
  type SearchParams,
  predictionModelById,
} from './lib/types'
import {
  calibrationNote,
  convertInputNote,
  convertOutputNote,
  fastaNote,
  libPathNote,
  pendingLibraryNote,
  libraryTargetPath,
  msDataNote,
  resultsNote,
  validateBuildRun,
  validateConvertRun,
  type Note,
  validateDownloadRun,
  downloadTargetPath,
  validateSearchRun,
} from './lib/validate'

const PERSIST_KEY = 'pioneerConsole.v2'

/** Finished runs, kept across sessions so their parameters can be recalled.
 *
 *  Separate key from the form draft so a corrupt or oversized history cannot
 *  take the working parameters down with it.
 *
 *  Parameters only — no log output. A single SearchDIA log is megabytes, and
 *  localStorage is a low-single-digit-megabyte budget shared with the draft; a
 *  handful of runs would evict everything. Restored entries therefore show an
 *  empty log, which is why the drawer says so rather than looking broken.
 */
const HISTORY_KEY = 'pioneerConsole.history.v1'

/** The next history number to hand out.
 *
 *  History numbers are an identity, not a position: run #1000 stays #1000 when
 *  earlier runs are deleted, and the counter never resets. It is therefore kept
 *  on its own rather than derived from the stored history, which would rewind
 *  the moment anything was removed.
 */
const RUN_COUNTER_KEY = 'pioneerConsole.runCounter.v1'

function nextRunNo(): number {
  let n = 1
  try {
    n = Math.max(0, Number(localStorage.getItem(RUN_COUNTER_KEY)) || 0) + 1
    localStorage.setItem(RUN_COUNTER_KEY, String(n))
  } catch {
    /* private mode — numbering restarts, which beats refusing to run */
  }
  return n
}

/** What survives a restart. A subset of Job: the identity, the outcome, and the
 *  form state needed to put the run back on screen. */
interface PersistedRun {
  id: string
  /** Its number in the history. Stable for the life of the run. */
  runNo: number
  cmd: CommandId
  title: string
  target: string
  threads: number
  status: JobStatus
  snapshot: JobSnapshot
}

/** A stored row as a Job the sidebar can render. */
function runToJob(r: backend.StoredRun): Job | null {
  let snapshot: JobSnapshot
  try {
    snapshot = JSON.parse(r.snapshot) as JobSnapshot
  } catch {
    return null
  }
  if (!snapshot?.cmd) return null
  return {
    id: r.id,
    runNo: r.run_no,
    cmd: r.cmd as CommandId,
    title: r.title,
    target: r.target,
    threads: r.threads,
    status: r.status as JobStatus,
    snapshot,
    logLines: [],
    failMsg: '',
    paramsJson: '',
    invocation: { kind: 'paramsFile' as const, json: '' },
    viewerPaths:
      snapshot.cmd === 'searchdia'
        ? {
            results: snapshot.search.results,
            msData: snapshot.search.msData,
            library: snapshot.search.library,
          }
        : null,
    paramsPath: '',
    finishedAt: r.finished_at,
  }
}

function jobToRun(j: Job): backend.StoredRun {
  return {
    id: j.id,
    run_no: j.runNo,
    cmd: j.cmd,
    title: j.title,
    target: j.target,
    threads: j.threads,
    status: j.status,
    snapshot: JSON.stringify(j.snapshot),
    finished_at: j.finishedAt,
  }
}

/** The old localStorage history, read only to hand it to the store once. */
function loadHistory(): Job[] {
  try {
    const raw = localStorage.getItem(HISTORY_KEY)
    if (!raw) return []
    const saved = JSON.parse(raw) as PersistedRun[]
    if (!Array.isArray(saved)) return []
    return saved
      .filter((r) => r && r.id && r.snapshot && r.snapshot.cmd)
      .map((r) => ({
        ...r,
        runNo: r.runNo ?? 0,
        finishedAt: 0,
        logLines: [],
        failMsg: '',
        paramsJson: '',
        // Never used: Run always builds a fresh invocation from the form. Present
        // only because Job requires it.
        invocation: { kind: 'paramsFile' as const, json: '' },
        viewerPaths:
          r.snapshot.cmd === 'searchdia'
            ? {
                results: r.snapshot.search.results,
                msData: r.snapshot.search.msData,
                library: r.snapshot.search.library,
              }
            : null,
        paramsPath: '',
      }))
  } catch {
    return []
  }
}

const TITLES: Record<CommandId, string> = {
  searchdia: 'SearchDIA',
  buildspeclib: 'BuildSpecLib',
  downloadspeclib: 'DownloadSpecLib',
  convertraw: 'ConvertRAW',
}

const SUBTITLES: Record<CommandId, string> = {
  searchdia: 'Find & quantify proteins',
  buildspeclib: 'Predict a spectral library',
  downloadspeclib: 'Download a prebuilt spectral library',
  convertraw: 'Convert .raw files to Arrow',
}

/** Append one line of process output, emulating the bits of a terminal that
 *  Julia's progress bars rely on.
 *
 *  ProgressBars.jl commits each redraw with `\n` and then emits `ESC[1A` to put
 *  the cursor back over it. A plain append turns one progress bar into dozens
 *  of lines, so a cursor-up (or a pending carriage-return redraw) removes the
 *  lines it was meant to overwrite. Only lines from the same stream are
 *  considered — output interleaved from the other pipe must not be eaten.
 */
function applyLine(lines: LogLine[], ev: backend.LineEvent): LogLine[] {
  const out = lines.slice()
  let toDrop = ev.overwrite
  for (let i = out.length - 1; i >= 0; i--) {
    if (out[i].stream !== ev.stream) continue
    if (toDrop > 0) {
      out.splice(i, 1)
      toDrop--
      continue
    }
    if (out[i].transient) out.splice(i, 1)
    break
  }
  out.push({ text: ev.line, stream: ev.stream, transient: ev.transient })
  return out
}


export default function App() {
  const [command, setCommand] = useState<CommandId>('searchdia')
  const [search, setSearch] = useState<SearchParams>(SEARCH_DEFAULTS)
  const [build, setBuild] = useState<BuildParams>(BUILD_DEFAULTS)
  const [convert, setConvert] = useState<ConvertParams>(CONVERT_DEFAULTS)
  const [download, setDownload] = useState<DownloadParams>(DOWNLOAD_DEFAULTS)
  /** The remote catalog, fetched by DownloadSpecLib rather than by the webview,
   *  which the CSP forbids from reaching any external origin. */
  const [libraries, setLibraries] = useState<RemoteLibrary[]>([])
  const [catalogLoading, setCatalogLoading] = useState(false)
  const [catalogError, setCatalogError] = useState('')
  /** Keys of a loaded config that the form does not model, per command. */
  const [extras, setExtras] = useState<Record<string, Json | null>>({})
  /** Optional name for the next run. Empty means "generate one". Shared across
   *  the command tabs: it names the run you are about to start, not the form. */
  const [jobName, setJobName] = useState('')
  const [advancedOpen, setAdvancedOpen] = useState(false)
  const [navCollapsed, setNavCollapsed] = useState(false)
  const [runError, setRunError] = useState('')
  const [modNote, setModNote] = useState({ fixed: '', variable: '' })

  const [jobs, setJobs] = useState<Job[]>([])

  /** Shown in the sidebar footer. `pioneer` versions the CLI tools and the
   *  converter together — they ship as one distribution — and `app` is this
   *  window. Either can be blank: an older distribution has no VERSION file. */
  const [versions, setVersions] = useState({ app: '', pioneer: '' })
  const [uninstall, setUninstall] = useState<backend.UninstallInfo | null>(null)
  const [uninstallOpen, setUninstallOpen] = useState(false)
  const [uninstalling, setUninstalling] = useState(false)
  const [uninstallError, setUninstallError] = useState('')
  const [viewJobId, setViewJobId] = useState<string | null>(null)
  const [drawerOpen, setDrawerOpen] = useState(false)
  const [drawerHeight, setDrawerHeight] = useState(300)
  const [confirmCancel, setConfirmCancel] = useState(false)

  // 0 means "not yet determined" — a non-zero value here would make the
  // derive-from-core-count branch below unreachable.
  const [threads, setThreads] = useState(0)
  const [maxThreads, setMaxThreads] = useState(8)
  const [pioneerError, setPioneerError] = useState('')

  const [jsonOpen, setJsonOpen] = useState(false)
  const [loadOpen, setLoadOpen] = useState(false)
  const [jobConfirm, setJobConfirm] = useState<{
    id: string
    kind: 'cancel' | 'delete'
    title: string
    /** Queued but never started, so there is no process to stop. */
    queued: boolean
  } | null>(null)
  const [overwriteOpen, setOverwriteOpen] = useState(false)

  const [pathInfos, setPathInfos] = useState<Record<string, PathInfo>>({})
  /** What the selected .poin records about itself, or null when there is none. */
  const [libInfo, setLibInfo] = useState<backend.LibraryInfo | null>(null)
  const [theme, setTheme] = useState<ThemeId>(loadTheme)

  /** Set while the form is showing a past run's parameters rather than your own
   *  draft. The drafts are stashed here on entry and put back when you click a
   *  workflow tab, so inspecting a run — or tweaking it and running the tweak —
   *  never costs you the work in progress. */
  const [inspectingJobId, setInspectingJobId] = useState<string | null>(null)
  const stashedDraft = useRef<{
    command: CommandId
    search: SearchParams
    build: BuildParams
    convert: ConvertParams
  } | null>(null)

  const jobSeq = useRef(0)
  const startedIds = useRef(new Set<string>())
  /** Job ids restart at 1 each launch, but the params file each job writes is
   *  kept so a run can be reproduced. Without a per-launch prefix a new
   *  session's job1 would overwrite the previous session's. */
  const sessionId = useRef(Date.now().toString(36))

  const isSearch = command === 'searchdia'
  const isConvert = command === 'convertraw'
  const isDownload = command === 'downloadspeclib'
  /** Driven by argv, not a params JSON. These two have no config file to show
   *  or edit, so the preview panel shows the command line instead. */
  const isArgv = isConvert || isDownload

  /** Read the catalog. Called when the page is first opened and by Refresh,
   *  rather than on a timer: the repository changes when a library is
   *  published, which is rare and not something to poll for. */
  const refreshCatalog = useCallback(() => {
    setCatalogLoading(true)
    setCatalogError('')
    backend
      .listSpecLibs()
      .then((text) => {
        const parsed = JSON.parse(text) as { libraries: RemoteLibrary[] }
        setLibraries(parsed.libraries ?? [])
      })
      .catch((e) => setCatalogError(String(e)))
      .finally(() => setCatalogLoading(false))
  }, [])

  useEffect(() => {
    if (isDownload && libraries.length === 0 && !catalogLoading && !catalogError) {
      refreshCatalog()
    }
  }, [isDownload, libraries.length, catalogLoading, catalogError, refreshCatalog])

  // ---- startup -----------------------------------------------------------

  useEffect(() => {
    try {
      const raw = localStorage.getItem(PERSIST_KEY)
      if (raw) {
        const saved = JSON.parse(raw)
        if (saved.search) setSearch((p) => ({ ...p, ...saved.search }))
        if (saved.build) setBuild((p) => ({ ...p, ...saved.build }))
        if (saved.convert) setConvert((p) => ({ ...p, ...saved.convert }))
        if (saved.command) setCommand(saved.command)
        if (typeof saved.threads === 'number') setThreads(saved.threads)
        if (typeof saved.jobName === 'string') setJobName(saved.jobName)
      }
    } catch {
      /* corrupt or absent — fall back to defaults */
    }

    backend
      .cpuCount()
      .then((n) => {
        setMaxThreads(n)
        // Default to all cores but one, leaving the machine usable. A value
        // already restored from disk wins, clamped to what this machine has.
        setThreads((t) => (t > 0 ? Math.min(t, n) : Math.max(1, n - 1)))
      })
      .catch(() => setThreads((t) => (t > 0 ? t : 1)))

    // Resolve at startup purely to surface a failure: pioneerError blocks the
    // Run button with the resolver's own message ("SearchDIA not found in ...,
    // looked for bin/SearchDIA and the pioneer wrapper"), which is more use than
    // the detail this used to print at the bottom of every page.
    backend
      .pioneerInfo()
      .then((info) => {
        setPioneerError('')
        setVersions((v) => ({ ...v, pioneer: info.version ?? '' }))
      })
      .catch((e) => setPioneerError(String(e)))
    backend
      .appVersion()
      .then((app) => setVersions((v) => ({ ...v, app })))
      .catch(() => {})
    backend
      .uninstallInfo()
      .then(setUninstall)
      .catch(() => {})
  }, [])

  useEffect(() => {
    // Not while inspecting a past run: the form is showing that run's params,
    // and saving them here would overwrite the draft this session is meant to
    // give back when you click a workflow tab.
    if (inspectingJobId) return
    try {
      localStorage.setItem(PERSIST_KEY, JSON.stringify({ command, search, build, convert, threads, jobName }))
    } catch {
      /* private mode / quota — persistence is a convenience, not a requirement */
    }
  }, [command, search, build, convert, threads, inspectingJobId])

  // A restored run has no log, so the drawer offers its output folder instead.
  // Whether that folder still exists is not knowable from the stored history --
  // it may have been moved, deleted, or be on a volume that is not mounted -- so
  // stat it rather than claiming.
  useEffect(() => {
    if (viewJobId) {
      const j = jobs.find((x) => x.id === viewJobId)
      if (j && j.logLines.length === 0 && j.target && j.target !== '—') {
        inspect('jobTarget', j.target)
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [viewJobId])

  /** Put a past run's parameters on screen, switching to its workflow. */
  const inspectJob = useCallback((job: Job) => {
    setInspectingJobId((current) => {
      // Stash only on the way in. Clicking from one run straight to another must
      // not overwrite the draft with the first run's params.
      if (!current) {
        stashedDraft.current = { command, search, build, convert }
      }
      return job.id
    })
    setCommand(job.cmd)
    if (job.snapshot.cmd === 'searchdia') setSearch(job.snapshot.search)
    else if (job.snapshot.cmd === 'buildspeclib') {
      // Stored runs from before digestion specificity was added do not carry
      // the field. Merge defaults so recalling one still renders and emits a
      // full-specific build. The same applies to the cleavage rule.
      setBuild({ ...BUILD_DEFAULTS, ...job.snapshot.build })
    } else if (job.snapshot.cmd === 'downloadspeclib') setDownload(job.snapshot.download)
    else setConvert(job.snapshot.convert)
    setRunError('')
  }, [command, search, build, convert])

  /** Leave inspection and restore the draft, if there is one to restore. */
  const restoreDraft = useCallback(() => {
    const d = stashedDraft.current
    if (!d) return
    setSearch(d.search)
    setBuild(d.build)
    setConvert(d.convert)
    stashedDraft.current = null
    setInspectingJobId(null)
  }, [])

  /** `${id}:${status}` pairs already written, so a re-render does not rewrite
   *  the history — but a status *change* does write, which is what lets an
   *  interrupted run be recognised later. */
  const savedRuns = useRef(new Set<string>())

  useEffect(() => {
    for (const j of jobs) {
      // Queued and running rows are written too, not just finished ones. That
      // is how interruption is detected: if the app goes away, the row stays
      // pending on disk and startup can see it never finished. Waiting for a
      // close event would miss a crash or a force quit.
      // runNo is part of the key because reordering the queue reassigns it,
      // and a status that has not changed would otherwise suppress the write.
      const key = `${j.id}:${j.status}:${j.runNo}`
      if (savedRuns.current.has(key)) continue
      savedRuns.current.add(key)
      backend.historySave(jobToRun(j)).catch(() => {
        savedRuns.current.delete(key)
      })
    }
  }, [jobs])

  /** Read the store, importing the old localStorage history the first time. */
  const refreshHistory = useCallback(async () => {
    try {
      if (await backend.historyNeedsImport()) {
        const legacy = loadHistory()
        const counter = Number(localStorage.getItem(RUN_COUNTER_KEY)) || 0
        await backend.historyImport(legacy.map(jobToRun), counter)
      }
      const rows = await backend.historyLoad()
      setJobs((prev) => {
        // A job this session started is already in memory, and that copy is the
        // better one: the store has no log output, and for a run still going it
        // does not yet know the outcome. So the in-memory job always wins, and
        // rows are only used for runs this session has never seen.
        //
        // Every status is written to the store now, so a live run *is* on disk.
        // Merging by id rather than concatenating is what keeps it from
        // appearing twice, and from being read back as interrupted while it is
        // still going.
        const live = new Map(prev.map((j) => [j.id, j]))
        const merged = rows
          .map(runToJob)
          .filter((j): j is Job => j !== null)
          .map((j) => {
            const mine = live.get(j.id)
            if (mine) return mine
            // Not this session's, and still pending on disk: the app that
            // started it went away. Say so rather than re-queuing — starting a
            // multi-hour search because someone opened the app is not a
            // reasonable thing to do unasked.
            return j.status === 'queued' || j.status === 'running'
              ? { ...j, status: 'interrupted' as JobStatus }
              : j
          })
        merged.forEach((j) => savedRuns.current.add(`${j.id}:${j.status}`))
        // Persist only the reinterpretations, so they are not redone next launch.
        const onDisk = new Map(rows.map((r) => [r.id, r.status]))
        merged
          .filter((j) => j.status === 'interrupted' && onDisk.get(j.id) !== 'interrupted')
          .forEach((j) => void backend.historySave(jobToRun(j)).catch(() => undefined))
        // Anything in memory the store has not caught up with yet.
        const seen = new Set(rows.map((r) => r.id))
        return [...merged, ...prev.filter((j) => !seen.has(j.id))]
      })
    } catch {
      /* store unavailable — leave whatever is on screen alone */
    }
  }, [])

  useEffect(() => {
    void refreshHistory()
  }, [refreshHistory])

  // Re-read when the window regains focus. Redundant while the app is
  // single-instance, but the database is a file: anything else that writes it
  // shows up on the next focus rather than needing a restart.
  useEffect(() => {
    const onFocus = () => void refreshHistory()
    window.addEventListener('focus', onFocus)
    return () => window.removeEventListener('focus', onFocus)
  }, [refreshHistory])

  useEffect(() => {
    applyTheme(theme)
  }, [theme])

  // Read the library's own config whenever the path changes. Debounced with the
  // same idea as path validation: the field is typed into, and each keystroke
  // would otherwise be a filesystem read.
  useEffect(() => {
    const path = search.library.trim()
    if (!path) {
      setLibInfo(null)
      return
    }
    let cancelled = false
    const t = setTimeout(() => {
      backend
        .libraryInfo(path)
        .then((i) => !cancelled && setLibInfo(i))
        .catch(() => !cancelled && setLibInfo(null))
    }, 250)
    return () => {
      cancelled = true
      clearTimeout(t)
    }
  }, [search.library])

  // Resolve home once, so pickers have somewhere sensible to start before the
  // user has picked anything or configured a default.
  useEffect(() => {
    void backend.initHomeDir()
  }, [])

  useEffect(() => {
    let stop: (() => void) | undefined
    let cancelled = false
    listenForDrops((key, kind, paths) => {
      void (async () => {
        const path = paths[0]
        // Trust the drop for the path but not for what it is: a folder named
        // lib.poin and a file named lib.poin are indistinguishable by string.
        const info = await backend.inspectPath(path)
        const isDir = info.is_dir
        if (kind === 'dir' && !isDir) {
          setRunError('That is a file — this field takes a folder.')
          return
        }
        if (kind === 'file' && isDir) {
          setRunError('That is a folder — this field takes a single file.')
          return
        }
        setRunError('')
        if (key === 'fastaAdd') {
          // Every dropped path, not just the first: picking several FASTAs at
          // once is the normal case.
          setBuild((b) => ({ ...b, fastaFiles: [...b.fastaFiles, ...paths.map(makeFastaRow)] }))
          return
        }
        if (key === 'convertInput') {
          setConvert((c) => ({ ...c, input: path, inputMode: isDir ? 'folder' : 'file' }))
          return
        }
        // `data-key` doubles as the scroll target for validation errors, and
        // ConvertRAW's output field is keyed `convertOutput` while the state it
        // writes is `outputDir`. Every other droppable field happens to use the
        // same name for both, so a plain onParam(key, ...) silently wrote a
        // field nothing reads and the drop appeared to do nothing at all.
        const DROP_KEY_TO_FIELD: Record<string, string> = { convertOutput: 'outputDir' }
        onParam(DROP_KEY_TO_FIELD[key] ?? key, path)
      })()
    }).then((un) => {
      if (cancelled) un()
      else stop = un
    })
    return () => {
      cancelled = true
      stop?.()
    }
    // onParam closes over the current command, so the listener is rebound when
    // the tab changes; otherwise a drop would land in the previous form.
  }, [command])

  // ---- live path validation ---------------------------------------------

  /** How many runs have reached a terminal state. Only its changing matters. */
  const finishedCount = jobs.filter(
    (j) => j.status === 'done' || j.status === 'failed' || j.status === 'cancelled',
  ).length

  const inspect = useCallback((key: string, value: string) => {
    backend.inspectPath(value).then((info) => {
      setPathInfos((prev) => ({ ...prev, [key]: info }))
    })
  }, [])

  const fastaPaths = build.fastaFiles.map((f) => f.path).join('\0')

  // Debounced so typing a long path doesn't stat on every keystroke.
  useEffect(() => {
    const t = setTimeout(() => {
      inspect('msData', search.msData)
      inspect('library', search.library)
      inspect('results', search.results)
      inspect('libTarget', libraryTargetPath(build.libPath))
      inspect('calibrationFile', build.calibrationFile)
      inspect('convertInput', convert.input)
      inspect('convertOutput', convert.outputDir)
      inspect('downloadTarget', downloadTargetPath(download.dest, download.selected))
      build.fastaFiles.forEach((f) => inspect(`fasta:${f.path}`, f.path))
    }, 250)
    return () => clearTimeout(t)
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [
    search.msData,
    search.library,
    search.results,
    build.libPath,
    build.calibrationFile,
    convert.input,
    convert.outputDir,
    download.dest,
    download.selected,
    fastaPaths,
    inspect,
    // Not a path, but a reason to re-stat them all: a run that has just
    // finished may have created the folder one of these fields points at, and
    // nothing else would prompt another look.
    finishedCount,
  ])

  const info = (k: string): PathInfo => pathInfos[k] ?? EMPTY_PATH_INFO

  /** The library a queued or running build/download is about to produce, so a
   *  search that consumes it can be queued behind it rather than waiting for
   *  the folder to appear. Most recently queued wins when there are several. */
  const pendingLibrary = useMemo(() => {
    for (let i = jobs.length - 1; i >= 0; i--) {
      const j = jobs[i]
      if (j.status !== 'queued' && j.status !== 'running') continue
      if (j.snapshot.cmd !== 'buildspeclib' && j.snapshot.cmd !== 'downloadspeclib') continue
      const path = libraryOf(j)
      if (path) return path
    }
    return null
  }, [jobs])

  // Fill an empty library field from that job. The Run handler already does
  // this at the moment a build is queued; this also covers arriving at
  // SearchDIA with the field cleared, or a build queued before it was.
  useEffect(() => {
    if (!pendingLibrary) return
    setSearch((p) => (p.library.trim() ? p : { ...p, library: pendingLibrary }))
  }, [pendingLibrary])

  const searchNotes = useMemo(
    () => ({
      msData: msDataNote(search.msData, info('msData')),
      library: pendingLibraryNote(search.library, info('library'), pendingLibrary),
      results: resultsNote(search.results, info('results')),
    }),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [search.msData, search.library, search.results, pathInfos, pendingLibrary],
  )

  const fastaNotes = useMemo(
    () => build.fastaFiles.map((f) => fastaNote(f.path, info(`fasta:${f.path}`))),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [fastaPaths, pathInfos],
  )

  const calibNote = useMemo(
    () => calibrationNote(build.calibrationFile, info('calibrationFile')),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [build.calibrationFile, pathInfos],
  )

  const convertNotes = useMemo(
    () => ({
      input: convertInputNote(convert, info('convertInput')),
      output: convertOutputNote(convert.outputDir, info('convertOutput')),
    }),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [
      convert.input,
      convert.inputMode,
      convert.outputDir,
      pathInfos,
    ],
  )

  const libNote = useMemo(
    () => libPathNote(build.libPath, info('libTarget')),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [build.libPath, pathInfos],
  )

  // ---- job events --------------------------------------------------------

  useEffect(() => {
    const unlisteners: Array<() => void> = []
    // The listeners are registered asynchronously. If this effect is torn down
    // before they resolve (React StrictMode remounts it immediately), the
    // unlisten handles would land after cleanup had already run and every log
    // line would be delivered twice.
    let disposed = false
    const keep = (u: () => void) => (disposed ? u() : unlisteners.push(u))

    backend
      .onJobLine((ev) => {
        setJobs((prev) =>
          prev.map((j) => (j.id === ev.job_id ? { ...j, logLines: applyLine(j.logLines, ev) } : j)),
        )
      })
      .then(keep)

    backend
      .onJobExit(({ job_id, success, cancelled, message }) => {
        setJobs((prev) =>
          prev.map((j) => {
            if (j.id !== job_id) return j
            const status: JobStatus = cancelled ? 'cancelled' : success ? 'done' : 'failed'
            return {
              ...j,
              status,
              finishedAt: Math.floor(Date.now() / 1000),
              failMsg: success ? '' : message,
              logLines: message
                ? [...j.logLines, { text: message, stream: 'app' as const, transient: false }]
                : j.logLines,
            }
          }),
        )
      })
      .then(keep)

    return () => {
      disposed = true
      unlisteners.forEach((u) => u())
    }
  }, [])

  // ---- the queue: one job at a time --------------------------------------

  useEffect(() => {
    if (jobs.some((j) => j.status === 'running')) return
    const next = jobs.find((j) => j.status === 'queued' && !startedIds.current.has(j.id))
    if (!next) return

    startedIds.current.add(next.id)
    const appLine = (text: string): LogLine => ({ text, stream: 'app', transient: false })

    setJobs((prev) =>
      prev.map((j) =>
        j.id === next.id
          ? {
              ...j,
              status: 'running',
              logLines: [...j.logLines, appLine(`$ ${j.title}`)],
            }
          : j,
      ),
    )

    backend
      .startJob(next.id, next.cmd, next.invocation, next.threads)
      .then(({ params_path, env_summary }) => {
        setJobs((prev) =>
          prev.map((j) =>
            j.id === next.id
              ? {
                  ...j,
                  paramsPath: params_path,
                  logLines: [
                    ...j.logLines,
                    ...(env_summary ? [appLine(env_summary)] : []),
                    appLine(`params: ${params_path}`),
                  ],
                }
              : j,
          ),
        )
      })
      .catch((e) => {
        setJobs((prev) =>
          prev.map((j) =>
            j.id === next.id
              ? {
                  ...j,
                  status: 'failed',
                  failMsg: String(e),
                  logLines: [...j.logLines, appLine(`ERROR: ${e}`)],
                }
              : j,
          ),
        )
      })
  }, [jobs])

  // ---- shared form handlers ----------------------------------------------

  const onParam = (key: string, value: string) => {
    if (isSearch) setSearch((p) => ({ ...p, [key]: value }))
    else if (isConvert) setConvert((p) => ({ ...p, [key]: value }))
    else if (key === 'predictionModel') switchModel(value)
    else setBuild((p) => ({ ...p, [key]: value }))
    setRunError('')
  }

  /** Modification support is per model, so switching models can strand entries
   *  the new model would have Koina reject. Keep what we can — a modification
   *  the new model knows on a site it does not is moved to its first allowed
   *  site — and drop only what it has no entry for at all. */
  const switchModel = (modelId: string) => {
    const label = predictionModelById(modelId).label
    const retarget = (mods: ModEntry[]): { kept: ModEntry[]; note: string } => {
      const kept: ModEntry[] = []
      const moved: string[] = []
      const dropped: string[] = []
      for (const m of mods) {
        const def = findMod(modelId, m.name)
        if (!def) {
          dropped.push(m.label || m.name)
        } else if (siteAllowed(def, m.pattern)) {
          kept.push(m)
        } else {
          moved.push(`${def.label} to ${def.sites[0]}`)
          kept.push({ ...m, pattern: def.sites[0] })
        }
      }
      const parts: string[] = []
      if (dropped.length) parts.push(`removed ${dropped.join(', ')}`)
      if (moved.length) parts.push(`moved ${moved.join(', ')}`)
      return { kept, note: parts.length ? `${label}: ${parts.join('; ')}.` : '' }
    }
    const f = retarget(build.fixedMods)
    const v = retarget(build.variableMods)
    setModNote({ fixed: f.note, variable: v.note })
    const req = enforceRequiredMods(modelId, f.kept, v.kept)
    setBuild((p) => ({
      ...p,
      predictionModel: modelId,
      fixedMods: req.fixed,
      variableMods: req.variable,
    }))
  }

  const onToggle = (key: string) => {
    if (isSearch) setSearch((p) => ({ ...p, [key]: !p[key as keyof SearchParams] }))
    else if (isConvert) setConvert((p) => ({ ...p, [key]: !p[key as keyof ConvertParams] }))
    else setBuild((p) => ({ ...p, [key]: !p[key as keyof BuildParams] }))
  }

  const browseConvertInput = async () => {
    const picked =
      convert.inputMode === 'file'
        ? await backend.pickFile('Choose a .raw file', 'Thermo RAW', ['raw'])
        : await backend.pickFolder('Choose a folder of .raw files')
    if (picked) onParam('input', picked)
  }

  const browseConvertOutput = async () => {
    const picked = await backend.pickFolder('Choose the output folder')
    if (picked) onParam('outputDir', picked)
  }

  const browseDownloadDest = async () => {
    const picked = await backend.pickFolder('Choose where to save the library')
    if (picked) setDownload((p) => ({ ...p, dest: picked }))
  }

  /** What is already at `<dest>/<library>`. Checked before the run rather than
   *  after: the binary refuses to overwrite, and finding that out only once the
   *  job has started wastes the trip and reads like a failure. */
  const downloadTargetExists = info('downloadTarget').exists

  /** The destination note. Deliberately shows the size that is about to land:
   *  a 3 GiB fetch should say so before it starts, not once it is running. */
  const downloadDestNote: Note = !download.dest.trim()
    ? { level: '', msg: 'Required — the library is written into this folder.' }
    : downloadTargetExists && download.selected
      ? {
          level: download.force ? 'warn' : 'error',
          msg: download.force
            ? `${download.selected} already exists here and will be replaced.`
            : `${download.selected} already exists here. Turn on “Replace an existing copy” to overwrite it, or choose another folder.`,
        }
      : (() => {
          const lib = libraries.find((l) => l.name === download.selected)
          return lib
            ? { level: '', msg: `${lib.size_human} will be written to ${download.dest}` }
            : { level: '', msg: `Libraries are written to ${download.dest}` }
        })()

  const onBrowseSearch = async (key: 'msData' | 'library' | 'results') => {
    const titles = {
      msData: 'Choose the MS data folder',
      library: 'Choose the spectral library folder (.poin)',
      results: 'Choose the results folder',
    }
    // All three are folders — a .poin library is a directory of tables. The
    // library gets its own picker so it can start at the one used last.
    const picked =
      key === 'library'
        ? await backend.pickLibrary(titles.library)
        : await backend.pickFolder(titles[key])
    if (picked) onParam(key, picked)
  }

  // ---- BuildSpecLib handlers ---------------------------------------------

  const patchFasta = (idx: number, fn: (f: FastaEntry) => FastaEntry) => {
    setBuild((p) => ({ ...p, fastaFiles: p.fastaFiles.map((f, i) => (i === idx ? fn(f) : f)) }))
    setRunError('')
  }

  const addFasta = async () => {
    const picked = await backend.pickFastaFiles()
    setBuild((p) => ({
      ...p,
      // Cancelling still adds a blank row, so the button always does something
      // and a path can be typed by hand.
      fastaFiles: [...p.fastaFiles, ...(picked.length ? picked.map(makeFastaRow) : [makeFastaRow('')])],
    }))
    setRunError('')
  }

  const browseFasta = async (idx: number) => {
    const picked = await backend.pickFastaFiles(false)
    if (!picked.length) return
    patchFasta(idx, (f) => {
      const path = picked[0]
      if (!f.auto) return { ...f, path }
      // While auto, re-derive the species tag and header preset from the name.
      const fresh = makeFastaRow(path)
      return { ...f, path, name: fresh.name, presetId: fresh.presetId, regex: fresh.regex }
    })
  }

  const onFastaField = (idx: number, field: keyof FastaEntry, value: string) => {
    patchFasta(idx, (f) => {
      const next = { ...f, [field]: value } as FastaEntry
      if (field === 'path' && f.auto) {
        // The design re-derived only the preset and regex here, leaving the
        // species tag at its default. That was invisible in the prototype
        // (Browse filled the path) but writes fasta_names: ["SPECIES"] as soon
        // as a path is typed. While `auto` holds, everything is derived.
        const fresh = makeFastaRow(value)
        next.presetId = fresh.presetId
        next.regex = fresh.regex
        next.name = fresh.name
      }
      // Naming the species by hand is an explicit choice; stop deriving.
      if (field === 'name') next.auto = false
      return next
    })
  }

  const onFastaPreset = (idx: number, preset: HeaderPresetId) => {
    patchFasta(idx, (f) => ({
      ...f,
      presetId: preset,
      auto: false,
      regex: preset === 'custom' ? f.regex : presetRegex(preset),
      showRegex: preset === 'custom' ? true : f.showRegex,
    }))
  }

  const onFastaRegex = (idx: number, key: keyof FastaEntry['regex'], value: string) => {
    patchFasta(idx, (f) => ({
      ...f,
      auto: false,
      presetId: 'custom',
      regex: { ...f.regex, [key]: value },
    }))
  }

  const toggleFastaRegex = (idx: number) => patchFasta(idx, (f) => ({ ...f, showRegex: !f.showRegex }))

  const removeFasta = (idx: number) => {
    setBuild((p) => ({ ...p, fastaFiles: p.fastaFiles.filter((_, i) => i !== idx) }))
    setRunError('')
  }

  const browseCalibration = async () => {
    const picked = await backend.pickFile('Choose one run from this experiment', 'MS data', [
      'arrow',
      'mzML',
      'mzml',
      'raw',
    ])
    if (picked) onParam('calibrationFile', picked)
  }

  const browseLibPath = async () => {
    const picked = await backend.pickLibraryTarget('Name the library to build')
    if (picked) onParam('libPath', picked)
  }

  const modKey = (kind: 'fixed' | 'variable') => (kind === 'fixed' ? 'fixedMods' : 'variableMods')

  const onModField = (kind: 'fixed' | 'variable', idx: number, field: keyof ModEntry, value: string) => {
    const key = modKey(kind)
    let displaced: string[] = []
    setBuild((p) => {
      const next = { ...p, [key]: p[key].map((m, i) => (i === idx ? { ...m, [field]: value } : m)) }
      if (kind !== 'fixed') return next
      // Widening a fixed mod onto a residue a variable one uses is the one way
      // left to build a conflict, since the variable side only ever offers free
      // residues. The variable rows yield, and are named rather than vanishing.
      const req = enforceRequiredMods(next.predictionModel, next.fixedMods, next.variableMods)
      displaced = next.variableMods
        .filter((m) => !req.variable.some((k) => k.name === m.name && k.pattern === m.pattern))
        .map((m) => `${m.label || m.name}${m.pattern ? ` on ${m.pattern}` : ''}`)
      return { ...next, fixedMods: req.fixed, variableMods: req.variable }
    })
    setModNote((n) => ({
      ...n,
      [kind]: '',
      ...(displaced.length
        ? { variable: `Fixed modifications now cover ${displaced.join(', ')} — removed from variable.` }
        : null),
    }))
    setRunError('')
  }

  const removeMod = (kind: 'fixed' | 'variable', idx: number) => {
    const key = modKey(kind)
    setBuild((p) => ({ ...p, [key]: p[key].filter((_, i) => i !== idx) }))
    setModNote((n) => ({ ...n, [kind]: '' }))
    setRunError('')
  }

  /** Choosing a preset writes its pattern; choosing Custom keeps whatever is
   *  there, so switching to Custom to make a small edit does not wipe the rule
   *  you were editing.
   *
   *  That second half used to be the whole handler, and did nothing visible:
   *  the picker reads its value back from the pattern, so keeping a preset's
   *  pattern meant the picker snapped straight back to that preset and the text
   *  field never appeared. The choice is now recorded rather than inferred. */
  const onEnzyme = (id: string) => {
    const preset = enzymeById(id)
    setBuild((p) =>
      preset
        ? { ...p, cleavageRegex: preset.pattern, customEnzyme: false }
        : { ...p, customEnzyme: true },
    )
    setRunError('')
  }

  /** Whether the typed rule is usable, judged by trying it. A pattern that
   *  compiles can still be wrong, which is what the preview beside it is for;
   *  this catches only what can be decided. */
  const cleavageNote: Note = (() => {
    const pattern = build.cleavageRegex.trim()
    if (!pattern) return { level: 'error', msg: 'Enter a cleavage rule, or choose an enzyme.' }
    // Cut sites, not peptides: whether the rule is sound is a separate question
    // from whether the length limits leave anything, and the preview beside it
    // now answers the second one.
    const sites = cleavageSites(SAMPLE_SEQUENCE, pattern)
    if (sites === null) return { level: 'error', msg: 'Not a valid regular expression.' }
    if (sites.length === 0) {
      return { level: 'warn', msg: 'This rule finds no cleavage site in the sample sequence.' }
    }
    if (sites.length > SAMPLE_SEQUENCE.length / 2) {
      return { level: 'warn', msg: 'This rule cleaves almost everywhere — check it is what you meant.' }
    }
    return { level: '', msg: '' }
  })()

  const addMod = (kind: 'fixed' | 'variable', unimod: string) => {
    if (!unimod) return
    const key = modKey(kind)
    setBuild((p) => {
      const def = modsForModel(p.predictionModel).find((d) => d.unimod === Number(unimod))
      if (!def) return p
      // C is reserved for the pinned fixed row, so a variable alkylation row
      // starts on the next site the model allows.
      const site = initialSite(p.predictionModel, kind, def, occupiedResidues(p.fixedMods))
      if (site === null) return p
      if (p[key].some((m) => unimodId(m.name) === def.unimod)) {
        setModNote((n) => ({ ...n, [kind]: `${def.label} is already in the list.` }))
        return p
      }
      setModNote((n) => ({ ...n, [kind]: '' }))
      return {
        ...p,
        [key]: [
          ...p[key],
          {
            pattern: site,
            label: def.label,
            name: `Unimod:${def.unimod}`,
            mass: String(def.mass),
          },
        ],
      }
    })
    setRunError('')
  }

  // ---- running -----------------------------------------------------------

  const currentExtras = extras[command] ?? null

  const jsonText = useMemo(
    () =>
      isDownload
        ? downloadCommandLine(download)
        : isConvert
        ? // Clamped: while the threads field is being cleared it holds 0, and
          // the preview should not show `--threads-per-file 0` as if that were
          // the command we would run. Run refuses in that state anyway.
          convertCommandLine(convert)
        : JSON.stringify(
            isSearch ? buildSearchJson(search, currentExtras) : buildLibJson(build, currentExtras),
            null,
            2,
          ),
    // `threads` belongs here: the ConvertRAW preview renders --threads-per-file
    // from it, so without it the command line goes stale when the stepper moves.
    [isConvert, isDownload, isSearch, search, build, convert, download, currentExtras, threads],
  )

  const overwriteNote = isConvert
    ? convertNotes.output
    // With force on, replacing an existing library is a warning rather than a
    // block, so it routes into the same confirm dialog the other commands use.
    // 3 GiB is too much to overwrite on a single click.
    : isDownload
      ? downloadDestNote
      : isSearch
        ? searchNotes.results
        : libNote

  /** Names already spoken for, across this session and persisted history. */
  const takenTitles = useMemo(() => new Set(jobs.map((j) => j.title)), [jobs])
  // Straight off the run history, which already records each run's parameters.
  const recentLibraryPaths = useMemo(() => recentLibraries(jobs), [jobs])

  const trimmedJobName = jobName.trim()
  /** Only computed for a typed name: for an empty one the answer is random, and
   *  a preview that reshuffles on every keystroke is noise. */
  const resolvedJobName = trimmedJobName ? resolveRunName(trimmedJobName, takenTitles) : ''

  const enqueue = async () => {
    jobSeq.current += 1
    const id = `${sessionId.current}-job${jobSeq.current}`
    // From the store, so it keeps climbing across restarts and cannot
    // diverge from the rows it numbers. Falls back to the local counter if
    // the store is unavailable.
    const runNo = await backend.historyNextRunNo().catch(() => nextRunNo())
    const job: Job = {
      id,
      runNo,
      finishedAt: 0,
      cmd: command,
      title: resolveRunName(jobName, jobs.map((j) => j.title)),
      snapshot: isConvert
        ? { cmd: 'convertraw' as const, convert }
        : isDownload
          ? { cmd: 'downloadspeclib' as const, download }
          : isSearch
            ? { cmd: 'searchdia' as const, search }
            : { cmd: 'buildspeclib' as const, build },
      target:
        (isConvert
          ? convert.outputDir || convert.input
          : isDownload
            ? download.dest
            : isSearch
              ? search.results
              : libraryTargetPath(build.libPath)) || '—',
      threads: Math.max(1, threads),
      status: 'queued',
      logLines: [],
      failMsg: '',
      paramsJson: jsonText,
      invocation: isConvert
        ? { kind: 'args' as const, args: buildConvertArgs(convert) }
        : isDownload
          ? { kind: 'args' as const, args: buildDownloadArgs(download) }
          : { kind: 'paramsFile' as const, json: jsonText },
      viewerPaths: isSearch
        ? { results: search.results, msData: search.msData, library: search.library }
        : null,
      paramsPath: '',
    }
    setJobs((prev) => [...prev, job])
    // If this run came from tweaking a past one, follow the new job rather than
    // staying pointed at the old one. The stashed draft is deliberately kept, so
    // a workflow tab still gives back the work in progress.
    setInspectingJobId((current) => (current ? id : null))
    setViewJobId(id)
    setDrawerOpen(true)
    setRunError('')
  }

  const run = (skipOverwriteCheck = false) => {
    if (pioneerError) {
      setRunError(pioneerError)
      return
    }
    // Not on ConvertRAW: the picker is hidden there and the value is unused,
    // so blocking on it would be an error about an invisible control.
    if (!isConvert && threads < 1) {
      setRunError(`Enter a thread count between 1 and ${maxThreads}.`)
      const el = document.querySelector('[data-key="threads"]')
      if (el) (el as HTMLElement).focus?.()
      return
    }
    const block = isConvert
      ? validateConvertRun(convert, convertNotes.input, convertNotes.output)
      : isDownload
        ? validateDownloadRun(download, downloadTargetExists)
        : isSearch
          ? validateSearchRun(search, searchNotes, pendingLibrary)
          : validateBuildRun(build, fastaNotes, libNote, calibNote)
    if (block) {
      setRunError(block.msg)
      const el = document.querySelector(`[data-key="${block.key}"]`)
      if (el) {
        el.scrollIntoView({ behavior: 'smooth', block: 'center' })
        setTimeout(() => (el as HTMLElement).focus?.(), 350)
      }
      return
    }
    if (!skipOverwriteCheck && overwriteNote.level === 'warn') {
      setOverwriteOpen(true)
      return
    }
    // Point SearchDIA at what this run is about to produce. Done here rather
    // than on success so the fields are ready while the job runs -- the usual
    // shape of a session is convert, then build or download, then search, and
    // filling them in afterwards means going back to look. The path does not
    // exist yet, so the field will show "does not exist" until the job lands;
    // the effect below re-checks every path when one finishes, which clears it
    // without the user touching anything.
    if (isConvert) {
      setSearch((p) => ({ ...p, msData: convert.outputDir || convert.input }))
    } else if (isDownload) {
      const produced = downloadTargetPath(download.dest, download.selected)
      if (produced) setSearch((p) => ({ ...p, library: produced }))
    } else if (!isSearch) {
      const produced = libraryTargetPath(build.libPath)
      if (produced) setSearch((p) => ({ ...p, library: produced }))
    }
    void enqueue()
  }

  /** Apply an edited or loaded config. Returns an error message, or null. */
  const applyConfig = (draft: string, switchCommand: boolean): string | null => {
    let obj: unknown
    try {
      obj = JSON.parse(draft)
    } catch (e) {
      return (e as Error).message
    }

    if (isSearchConfig(obj)) {
      const patch = searchConfigToState(obj)
      if (!patch) return 'Unrecognized configuration.'
      if (!switchCommand && !isSearch) return 'That is a SearchDIA config, not a BuildSpecLib one.'
      setSearch((p) => ({ ...p, ...patch }))
      setExtras((e) => ({ ...e, searchdia: computeExtras(obj, SEARCH_OWNED_PATHS) }))
      if (switchCommand) setCommand('searchdia')
      setRunError('')
      return null
    }

    if (isBuildConfig(obj)) {
      const patch = buildConfigToState(obj)
      if (!patch) return 'Unrecognized configuration.'
      if (!switchCommand && isSearch) return 'That is a BuildSpecLib config, not a SearchDIA one.'
      // A loaded config is outside the mod editor entirely, and older or
      // hand-written ones carry carbamidomethyl as variable, or not at all.
      setBuild((p) => {
        const next = { ...p, ...patch }
        const req = enforceRequiredMods(next.predictionModel, next.fixedMods, next.variableMods)
        return { ...next, fixedMods: req.fixed, variableMods: req.variable }
      })
      setExtras((e) => ({ ...e, buildspeclib: computeExtras(obj, BUILD_OWNED_PATHS) }))
      if (switchCommand) setCommand('buildspeclib')
      setRunError('')
      return null
    }

    return 'Unrecognized configuration — expected a SearchDIA or BuildSpecLib config.'
  }

  const applyLoad = (draft: string): string | null => {
    const err = applyConfig(draft, true)
    if (err) return err
    setLoadOpen(false)
    setAdvancedOpen(false)
    return null
  }

  // ---- drawer drag -------------------------------------------------------

  const startDragDrawer = (e: React.MouseEvent) => {
    e.preventDefault()
    const startY = e.clientY
    const startH = drawerHeight
    const onMove = (ev: MouseEvent) => {
      const maxH = window.innerHeight - 140
      setDrawerHeight(Math.max(140, Math.min(maxH, startH - (ev.clientY - startY))))
    }
    const onUp = () => {
      window.removeEventListener('mousemove', onMove)
      window.removeEventListener('mouseup', onUp)
      window.removeEventListener('blur', onUp)
      document.body.style.userSelect = ''
    }
    document.body.style.userSelect = 'none'
    window.addEventListener('mousemove', onMove)
    window.addEventListener('mouseup', onUp)
    window.addEventListener('blur', onUp)
  }

  // ---- keyboard ----------------------------------------------------------

  const anyModalOpen = jsonOpen || loadOpen || overwriteOpen || uninstallOpen || !!jobConfirm

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if ((e.metaKey || e.ctrlKey) && !e.shiftKey && !e.altKey && ['1', '2', '3', '4'].includes(e.key)) {
        if (anyModalOpen) return
        e.preventDefault()
        setCommand(
          e.key === '1'
            ? 'convertraw'
            : e.key === '2'
              ? 'buildspeclib'
              : e.key === '3'
                ? 'downloadspeclib'
                : 'searchdia',
        )
        return
      }
      if (e.key !== 'Enter' || e.shiftKey || anyModalOpen) return
      const t = e.target as HTMLElement | null
      if (t && (t.tagName === 'TEXTAREA' || t.tagName === 'BUTTON' || t.tagName === 'SELECT')) return
      if (t && t.tagName === 'INPUT') {
        e.preventDefault()
        t.blur()
      }
      run()
    }
    window.addEventListener('keydown', onKey)
    return () => window.removeEventListener('keydown', onKey)
  })

  // ---- derived -----------------------------------------------------------

  const viewJob = jobs.find((j) => j.id === viewJobId) ?? null
  const statusText = viewJob
    ? viewJob.status === 'running'
      ? `Running ${viewJob.title}…`
      : {
          queued: 'Queued',
          done: 'Completed',
          failed: 'Failed',
          cancelled: 'Cancelled',
          interrupted: 'Interrupted',
          running: 'Running',
        }[viewJob.status]
    : 'Ready'
  const extraKeys = extraLeafPaths(currentExtras)
  const modKeyLabel = /Mac|iPhone|iPad/.test(navigator.platform || navigator.userAgent) ? '⌘' : 'Ctrl '
  const runLabel = `Run ${TITLES[command]}`

  return (
    <div
      style={{
        display: 'flex',
        height: '100vh',
        width: '100vw',
        background: '#F4F6F8',
        color: '#1A2230',
        overflow: 'hidden',
      }}
    >
      <Sidebar
        collapsed={navCollapsed}
        versions={versions}
        uninstall={uninstall}
        uninstallBlocked={jobs.some((j) => j.status === 'queued' || j.status === 'running')}
        onUninstall={() => {
          setUninstallError('')
          setUninstallOpen(true)
        }}
        selected={command}
        jobs={jobs}
        viewJobId={viewJobId}
        modKey={modKeyLabel}
        theme={theme}
        onTheme={setTheme}
        onSelect={(id) => {
          // Clicking a workflow tab always means "back to my own parameters".
          restoreDraft()
          setCommand(id)
          setRunError('')
        }}
        onToggleCollapsed={() => setNavCollapsed((c) => !c)}
        onViewJob={(id) => {
          const job = jobs.find((j) => j.id === id)
          if (job) inspectJob(job)
          if (drawerOpen && viewJobId === id) setDrawerOpen(false)
          else {
            setViewJobId(id)
            setDrawerOpen(true)
          }
        }}
        onJobAction={(id, kind) => {
          const job = jobs.find((j) => j.id === id)
          setJobConfirm({
            id,
            kind,
            title: job ? job.title : '',
            queued: job?.status === 'queued',
          })
        }}
        onReorderQueued={(dragId, dropId) =>
          setJobs((prev) => moveQueuedJob(prev, dragId, dropId))
        }
        onRenameJob={(id, title) => {
          const wanted = title.trim()
          if (!wanted) return
          setJobs((prev) =>
            prev.map((j) => {
              if (j.id !== id) return j
              // Resolved against every other run, this one excluded -- keeping
              // its own name in `taken` would turn a no-op edit into "name-2".
              // Same rule as naming a new run, so the two cannot disagree about
              // what counts as taken.
              const resolved = resolveRunName(
                wanted,
                prev.filter((o) => o.id !== id).map((o) => o.title),
              )
              return resolved === j.title ? j : { ...j, title: resolved }
            }),
          )
        }}
      />

      <main style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0, position: 'relative' }}>
        <header
          data-tauri-drag-region
          style={{
            // 62 of header plus the overlay title bar's strip, so the title and
            // Run button keep the same clearance they had under the OS bar
            // rather than riding up into it.
            height: 62 + TITLEBAR_H,
            flex: 'none',
            borderBottom: '1px solid #E5E9ED',
            background: '#fff',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            padding: `${TITLEBAR_H}px 28px 0`,
          }}
        >
          <div style={{ minWidth: 0 }}>
            <div style={{ fontSize: 16, fontWeight: 700 }}>{TITLES[command]}</div>
            <div style={{ fontSize: 12.5, color: '#667085', marginTop: 1 }}>{SUBTITLES[command]}</div>
          </div>
          <TopBar
            // ConvertRAW is a .NET program: it never reads JULIA_NUM_THREADS,
            // and its own thread count lives in the form. A picker here would
            // be a control that changes nothing.
            showThreads={!isConvert}
            threads={threads}
            maxThreads={maxThreads}
            runLabel={runLabel}
            previewLabel={isArgv ? 'View the command line' : 'Edit the .json directly'}
            onThreads={setThreads}
            onEditJson={() => setJsonOpen(true)}
            onRun={() => run()}
          />
        </header>

        {/* minHeight: 0 is load-bearing. A flex item defaults to min-height:auto,
            which refuses to shrink below its content, so without this the drawer
            would push the form off the bottom instead of taking space from it. */}
        <div style={{ flex: 1, minHeight: 0, overflowY: 'auto', padding: '24px 44px 40px' }}>
          <div style={{ maxWidth: 740, margin: '0 auto' }}>
            {(runError || pioneerError) && (
              <div
                style={{
                  marginBottom: 22,
                  padding: '13px 15px',
                  borderRadius: 11,
                  background: '#FEF2F2',
                  border: '1px solid #FECACA',
                  display: 'flex',
                  alignItems: 'center',
                  gap: 11,
                  // Never intercept a drop. The banner is not interactive, but
                  // it appears above the fields and would otherwise be what
                  // elementFromPoint returns for a release aimed at the field
                  // underneath it.
                  pointerEvents: 'none',
                }}
              >
                <svg width="18" height="18" viewBox="0 0 24 24" fill="none" style={{ flex: 'none' }}>
                  <path d="M12 8v5M12 16.5v.5" stroke="#DC2626" strokeWidth="2" strokeLinecap="round" />
                  <circle cx="12" cy="12" r="9" stroke="#DC2626" strokeWidth="1.8" />
                </svg>
                <span style={{ fontSize: 13, color: '#B42318', fontWeight: 500 }}>
                  {runError || pioneerError}
                </span>
              </div>
            )}

            {/* SearchDIA carries this inside its Essentials card, beside the
                paths the run will be remembered alongside. The other commands
                have no equivalent card to host it, so for them it keeps one of
                its own. */}
            {!isSearch && (
              <section
                style={{
                  background: '#fff',
                  border: '1px solid #E7EAEE',
                  borderRadius: 13,
                  padding: '18px 20px',
                  marginBottom: 14,
                }}
              >
                <JobNameField
                  value={jobName}
                  resolved={resolvedJobName}
                  onChange={setJobName}
                />
              </section>
            )}

            {isConvert ? (
              <ConvertRawForm
                params={convert}
                inputNote={convertNotes.input}
                outputNote={convertNotes.output}
                advancedOpen={advancedOpen}
                onParam={onParam}
                onToggle={onToggle}
                onBrowseInput={browseConvertInput}
                onBrowseOutput={browseConvertOutput}
                onToggleAdvanced={() => setAdvancedOpen((o) => !o)}
              />
            ) : isDownload ? (
              <DownloadSpecLibForm
                params={download}
                libraries={libraries}
                loading={catalogLoading}
                error={catalogError}
                destNote={downloadDestNote}
                onParam={(key, value) => setDownload((p) => ({ ...p, [key]: value }))}
                onToggleForce={() => setDownload((p) => ({ ...p, force: !p.force }))}
                onSelect={(name) => setDownload((p) => ({ ...p, selected: name }))}
                onBrowseDest={browseDownloadDest}
                onRefresh={refreshCatalog}
              />
            ) : isSearch ? (
              <SearchDiaForm
                params={search}
                notes={searchNotes}
                libInfo={libInfo}
                recentLibraries={recentLibraryPaths}
                jobName={jobName}
                resolvedJobName={resolvedJobName}
                onParam={onParam}
                onToggle={onToggle}
                onBrowse={onBrowseSearch}
                onJobName={setJobName}
                onOpenLoad={() => setLoadOpen(true)}
                onGoToBuild={() => setCommand('buildspeclib')}
              />
            ) : (
              <BuildSpecLibForm
                params={build}
                fastaNotes={fastaNotes}
                libNote={libNote}
                calibNote={calibNote}
                modNote={modNote}
                onParam={onParam}
                onToggle={onToggle}
                onOpenLoad={() => setLoadOpen(true)}
                onAddFasta={addFasta}
                onBrowseFasta={browseFasta}
                onFastaField={onFastaField}
                onFastaPreset={onFastaPreset}
                onFastaRegex={onFastaRegex}
                onToggleFastaRegex={toggleFastaRegex}
                onRemoveFasta={removeFasta}
                onBrowseLibPath={browseLibPath}
                onBrowseCalibration={browseCalibration}
                onModField={onModField}
                onRemoveMod={removeMod}
                onAddMod={addMod}
                onEnzyme={onEnzyme}
                cleavageNote={cleavageNote}
              />
            )}

          </div>
        </div>

        {drawerOpen && (
          <LogDrawer
            job={viewJob}
            targetInfo={pathInfos['jobTarget'] ?? null}
            height={drawerHeight}
            statusText={statusText}
            confirmCancel={confirmCancel}
            onStartDrag={startDragDrawer}
            onAskCancel={() => setConfirmCancel(true)}
            onKeepRunning={() => setConfirmCancel(false)}
            onConfirmCancel={() => {
              setConfirmCancel(false)
              if (viewJobId) backend.cancelJob(viewJobId).catch(() => undefined)
            }}
            onOpenTarget={() => {
              const job = jobs.find((j) => j.id === viewJobId)
              if (!job?.target) return
              // The folder can go between the last stat and the click, so the
              // failure is reported rather than swallowed.
              backend.openFolder(job.target).catch((e) => setRunError(String(e)))
            }}
            onClose={() => setDrawerOpen(false)}
          />
        )}

        {jsonOpen && (
          <JsonModal
            fileName={isArgv ? 'command line' : `${command}.json`}
            editable={!isArgv}
            text={jsonText}
            hasExtras={!!currentExtras}
            extraKeys={extraKeys}
            onClose={() => setJsonOpen(false)}
            onApply={(draft) => applyConfig(draft, false)}
            onRevert={() => setExtras((e) => ({ ...e, [command]: null }))}
          />
        )}

        {loadOpen && <LoadConfigModal onClose={() => setLoadOpen(false)} onApply={applyLoad} />}

        {jobConfirm && (
          <ConfirmDialog
            tone="danger"
            title={
              jobConfirm.kind === 'delete'
                ? 'Delete this run?'
                : jobConfirm.queued
                  ? 'Take this run out of the queue?'
                  : 'Stop this run?'
            }
            body={
              jobConfirm.kind === 'delete'
                ? 'This removes the run and its log from the history. This can’t be undone.'
                : jobConfirm.queued
                  ? 'It has not started, so nothing is interrupted. It leaves the queue and will not run.'
                  : 'This stops the Pioneer process and lets the next queued job start.'
            }
            detail={jobConfirm.title}
            dismissLabel={jobConfirm.queued ? 'Leave in queue' : 'Keep'}
            confirmLabel={
              jobConfirm.kind === 'delete'
                ? 'Delete'
                : jobConfirm.queued
                  ? 'Remove'
                  : 'Stop run'
            }
            onDismiss={() => setJobConfirm(null)}
            onConfirm={() => {
              const { id, kind } = jobConfirm
              if (kind === 'delete') {
                // Also from the store, or the next read brings it back.
                backend.historyDelete(id).catch(() => undefined)
                setJobs((prev) => {
                  const rest = prev.filter((j) => j.id !== id)
                  if (viewJobId === id) {
                    setViewJobId(rest.length ? rest[rest.length - 1].id : null)
                    if (!rest.length) setDrawerOpen(false)
                  }
                  return rest
                })
              } else if (jobConfirm.queued) {
                // No process to cancel -- cancelJob would find nothing to kill,
                // which is why this used to leave the row sitting in the queue.
                // Marking it cancelled takes it out of the queue the scheduler
                // reads and moves it into the history, keeping its run number:
                // "cancelled before it starts" is a state the numbering already
                // accounts for.
                setJobs((prev) =>
                  prev.map((j) =>
                    j.id === id ? { ...j, status: 'cancelled' as JobStatus } : j,
                  ),
                )
              } else {
                backend.cancelJob(id).catch(() => undefined)
              }
              setJobConfirm(null)
            }}
          />
        )}

        {overwriteOpen && (
          <ConfirmDialog
            tone="warning"
            title={isSearch ? 'Results folder already exists' : 'A library already exists here'}
            body={overwriteNote.msg}
            detail={isSearch ? search.results : libraryTargetPath(build.libPath)}
            dismissLabel="Cancel"
            confirmLabel={isSearch ? 'Overwrite & run' : 'Overwrite & build'}
            onDismiss={() => setOverwriteOpen(false)}
            onConfirm={() => {
              setOverwriteOpen(false)
              run(true)
            }}
          />
        )}

        {uninstallOpen && uninstall?.available && (
          <ConfirmDialog
            tone="danger"
            title={`Uninstall Pioneer ${uninstall.version}?`}
            body="This removes this version’s macOS app and command-line tools. Your settings, run history, libraries, and analysis data are kept."
            detail={`${uninstall.app_path}\n${uninstall.install_root}`}
            dismissLabel="Keep Pioneer"
            confirmLabel={uninstalling ? 'Uninstalling…' : 'Uninstall'}
            pending={uninstalling}
            error={uninstallError}
            onDismiss={() => {
              setUninstallOpen(false)
              setUninstallError('')
            }}
            onConfirm={() => {
              setUninstalling(true)
              setUninstallError('')
              backend
                .uninstallThisVersion()
                .then(() => setUninstallOpen(false))
                .catch((e) => {
                  setUninstallError(String(e))
                  setUninstalling(false)
                })
            }}
          />
        )}
      </main>
    </div>
  )
}
