import { useCallback, useEffect, useMemo, useRef, useState } from 'react'

import { BuildSpecLibForm } from './components/BuildSpecLibForm'
import { ConfirmDialog } from './components/ConfirmDialog'
import { ConvertRawForm } from './components/ConvertRawForm'
import { JsonModal } from './components/JsonModal'
import { LoadConfigModal } from './components/LoadConfigModal'
import { LogDrawer } from './components/LogDrawer'
import { SearchDiaForm } from './components/SearchDiaForm'
import { TopBar } from './components/TopBar'
import { Sidebar } from './components/Sidebar'
import * as backend from './lib/backend'
import {
  buildConvertArgs,
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
import { MOD_PRESETS, makeFastaRow, presetRegex } from './lib/fasta'
import { generateRunName } from './lib/names'
import { applyTheme, loadTheme, type ThemeId } from './lib/theme'
import {
  BUILD_DEFAULTS,
  CONVERT_DEFAULTS,
  EMPTY_PATH_INFO,
  SEARCH_DEFAULTS,
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
  type SearchParams,
} from './lib/types'
import {
  calibrationNote,
  convertInputNote,
  convertOutputNote,
  fastaNote,
  libPathNote,
  libraryNote,
  libraryTargetPath,
  msDataNote,
  resultsNote,
  validateBuildRun,
  validateConvertRun,
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
const HISTORY_LIMIT = 100

/** What survives a restart. A subset of Job: the identity, the outcome, and the
 *  form state needed to put the run back on screen. */
interface PersistedRun {
  id: string
  cmd: CommandId
  title: string
  target: string
  threads: number
  status: JobStatus
  snapshot: JobSnapshot
}

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
  convertraw: 'ConvertRAW',
}

const SUBTITLES: Record<CommandId, string> = {
  searchdia: 'Find & quantify proteins · all settings on one page',
  buildspeclib: 'Predict a spectral library · all settings on one page',
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
  /** Keys of a loaded config that the form does not model, per command. */
  const [extras, setExtras] = useState<Record<string, Json | null>>({})
  const [advancedOpen, setAdvancedOpen] = useState(false)
  const [navCollapsed, setNavCollapsed] = useState(false)
  const [runError, setRunError] = useState('')
  const [modNote, setModNote] = useState({ fixed: '', variable: '' })

  const [jobs, setJobs] = useState<Job[]>(loadHistory)
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
  const [jobConfirm, setJobConfirm] = useState<{ id: string; kind: 'cancel' | 'delete'; title: string } | null>(null)
  const [overwriteOpen, setOverwriteOpen] = useState(false)

  const [pathInfos, setPathInfos] = useState<Record<string, PathInfo>>({})
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
      .then(() => setPioneerError(''))
      .catch((e) => setPioneerError(String(e)))
  }, [])

  useEffect(() => {
    // Not while inspecting a past run: the form is showing that run's params,
    // and saving them here would overwrite the draft this session is meant to
    // give back when you click a workflow tab.
    if (inspectingJobId) return
    try {
      localStorage.setItem(PERSIST_KEY, JSON.stringify({ command, search, build, convert, threads }))
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
    else if (job.snapshot.cmd === 'buildspeclib') setBuild(job.snapshot.build)
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

  useEffect(() => {
    const finished = jobs.filter(
      (j) => j.status === 'done' || j.status === 'failed' || j.status === 'cancelled',
    )
    // Keep the most recent HISTORY_LIMIT. Older entries fall off the front, so a
    // long-running session cannot grow the stored blob without bound.
    const keep = finished.slice(-HISTORY_LIMIT).map<PersistedRun>((j) => ({
      id: j.id,
      cmd: j.cmd,
      title: j.title,
      target: j.target,
      threads: j.threads,
      status: j.status,
      snapshot: j.snapshot,
    }))
    try {
      localStorage.setItem(HISTORY_KEY, JSON.stringify(keep))
    } catch {
      /* quota — history is a convenience, and the draft lives under its own key */
    }
  }, [jobs])

  useEffect(() => {
    applyTheme(theme)
  }, [theme])

  // ---- live path validation ---------------------------------------------

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
    fastaPaths,
    inspect,
  ])

  const info = (k: string): PathInfo => pathInfos[k] ?? EMPTY_PATH_INFO

  const searchNotes = useMemo(
    () => ({
      msData: msDataNote(search.msData, info('msData')),
      library: libraryNote(search.library, info('library')),
      results: resultsNote(search.results, info('results')),
    }),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [search.msData, search.library, search.results, pathInfos],
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
    else setBuild((p) => ({ ...p, [key]: value }))
    setRunError('')
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

  const onBrowseSearch = async (key: 'msData' | 'library' | 'results') => {
    const titles = {
      msData: 'Choose the MS data folder',
      library: 'Choose the spectral library folder (.poin)',
      results: 'Choose the results folder',
    }
    // All three are folders — a .poin library is a directory of tables.
    const picked = await backend.pickFolder(titles[key])
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
    const picked = await backend.pickFile('Choose one MS data file', 'MS data', [
      'arrow',
      'mzML',
      'mzml',
      'raw',
    ])
    if (picked) onParam('calibrationFile', picked)
  }

  const browseLibPath = async () => {
    const picked = await backend.pickFolder('Choose where to write the library')
    if (picked) onParam('libPath', picked)
  }

  const modKey = (kind: 'fixed' | 'variable') => (kind === 'fixed' ? 'fixedMods' : 'variableMods')

  const onModField = (kind: 'fixed' | 'variable', idx: number, field: keyof ModEntry, value: string) => {
    const key = modKey(kind)
    setBuild((p) => ({
      ...p,
      [key]: p[key].map((m, i) => (i === idx ? { ...m, [field]: value } : m)),
    }))
    setModNote((n) => ({ ...n, [kind]: '' }))
    setRunError('')
  }

  const removeMod = (kind: 'fixed' | 'variable', idx: number) => {
    const key = modKey(kind)
    setBuild((p) => ({ ...p, [key]: p[key].filter((_, i) => i !== idx) }))
    setModNote((n) => ({ ...n, [kind]: '' }))
    setRunError('')
  }

  const addMod = (kind: 'fixed' | 'variable', preset: string) => {
    if (!preset) return
    const def = MOD_PRESETS[preset] ?? MOD_PRESETS.custom
    const key = modKey(kind)
    setBuild((p) => {
      const exists = p[key].some(
        (m) => (m.name && m.name === def.name) || (m.pattern === def.pattern && m.label === def.label),
      )
      if (exists) {
        setModNote((n) => ({ ...n, [kind]: `${def.label} is already in the list.` }))
        return p
      }
      setModNote((n) => ({ ...n, [kind]: '' }))
      return {
        ...p,
        [key]: [...p[key], { pattern: def.pattern, label: def.label, name: def.name, mass: def.mass }],
      }
    })
    setRunError('')
  }

  // ---- running -----------------------------------------------------------

  const currentExtras = extras[command] ?? null

  const jsonText = useMemo(
    () =>
      isConvert
        ? convertCommandLine(convert, threads)
        : JSON.stringify(
            isSearch ? buildSearchJson(search, currentExtras) : buildLibJson(build, currentExtras),
            null,
            2,
          ),
    // `threads` belongs here: the ConvertRAW preview renders --threads-per-file
    // from it, so without it the command line goes stale when the stepper moves.
    [isConvert, isSearch, search, build, convert, currentExtras, threads],
  )

  const overwriteNote = isConvert
    ? convertNotes.output
    : isSearch
      ? searchNotes.results
      : libNote

  const enqueue = () => {
    jobSeq.current += 1
    const id = `${sessionId.current}-job${jobSeq.current}`
    const job: Job = {
      id,
      cmd: command,
      title: generateRunName(jobs.map((j) => j.title)),
      snapshot: isConvert
        ? { cmd: 'convertraw' as const, convert }
        : isSearch
          ? { cmd: 'searchdia' as const, search }
          : { cmd: 'buildspeclib' as const, build },
      target:
        (isConvert
          ? convert.outputDir || convert.input
          : isSearch
            ? search.results
            : libraryTargetPath(build.libPath)) || '—',
      threads: Math.max(1, threads),
      status: 'queued',
      logLines: [],
      failMsg: '',
      paramsJson: jsonText,
      invocation: isConvert
        ? { kind: 'args' as const, args: buildConvertArgs(convert, threads) }
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
    const block = isConvert
      ? validateConvertRun(convert, convertNotes.input, convertNotes.output)
      : isSearch
        ? validateSearchRun(search, searchNotes)
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
    enqueue()
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
      setBuild((p) => ({ ...p, ...patch }))
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

  const anyModalOpen = jsonOpen || loadOpen || overwriteOpen || !!jobConfirm

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if ((e.metaKey || e.ctrlKey) && !e.shiftKey && !e.altKey && ['1', '2', '3'].includes(e.key)) {
        if (anyModalOpen) return
        e.preventDefault()
        setCommand(
          e.key === '1' ? 'convertraw' : e.key === '2' ? 'buildspeclib' : 'searchdia',
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
      : { queued: 'Queued', done: 'Completed', failed: 'Failed', cancelled: 'Cancelled' }[viewJob.status]
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
          setJobConfirm({ id, kind, title: job ? job.title : '' })
        }}
      />

      <main style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0, position: 'relative' }}>
        <header
          style={{
            height: 62,
            flex: 'none',
            borderBottom: '1px solid #E5E9ED',
            background: '#fff',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            padding: '0 28px',
          }}
        >
          <div style={{ minWidth: 0 }}>
            <div style={{ fontSize: 16, fontWeight: 700 }}>{TITLES[command]}</div>
            <div style={{ fontSize: 12.5, color: '#667085', marginTop: 1 }}>{SUBTITLES[command]}</div>
          </div>
          <TopBar
            showThreads
            threads={threads}
            maxThreads={maxThreads}
            runLabel={runLabel}
            previewLabel={isConvert ? 'View the command line' : 'Edit the .json directly'}
            onThreads={setThreads}
            onEditJson={() => setJsonOpen(true)}
            onRun={() => run()}
          />
        </header>

        {/* minHeight: 0 is load-bearing. A flex item defaults to min-height:auto,
            which refuses to shrink below its content, so without this the drawer
            would push the form off the bottom instead of taking space from it. */}
        <div style={{ flex: 1, minHeight: 0, overflowY: 'auto', padding: '24px 28px 40px' }}>
          <div style={{ maxWidth: 680, margin: '0 auto' }}>
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
            ) : isSearch ? (
              <SearchDiaForm
                params={search}
                notes={searchNotes}
                advancedOpen={advancedOpen}
                onParam={onParam}
                onToggle={onToggle}
                onBrowse={onBrowseSearch}
                onToggleAdvanced={() => setAdvancedOpen((o) => !o)}
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
            onClose={() => setDrawerOpen(false)}
          />
        )}

        {jsonOpen && (
          <JsonModal
            fileName={isConvert ? 'command line' : `${command}.json`}
            editable={!isConvert}
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
            title={jobConfirm.kind === 'delete' ? 'Delete this run?' : 'Stop this run?'}
            body={
              jobConfirm.kind === 'delete'
                ? 'This removes the run and its log from the queue. This can’t be undone.'
                : 'This stops the Pioneer process and lets the next queued job start.'
            }
            detail={jobConfirm.title}
            dismissLabel="Keep"
            confirmLabel={jobConfirm.kind === 'delete' ? 'Delete' : 'Stop run'}
            onDismiss={() => setJobConfirm(null)}
            onConfirm={() => {
              const { id, kind } = jobConfirm
              if (kind === 'delete') {
                setJobs((prev) => {
                  const rest = prev.filter((j) => j.id !== id)
                  if (viewJobId === id) {
                    setViewJobId(rest.length ? rest[rest.length - 1].id : null)
                    if (!rest.length) setDrawerOpen(false)
                  }
                  return rest
                })
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
      </main>
    </div>
  )
}
