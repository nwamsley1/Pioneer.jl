/** The navy sidebar: brand, workflow picker, live job queue, thread picker
 *  and the Run button. Ported from the `<aside>` block of the design.
 *
 *  The design's "Analysis" section (viewer tabs mounting Pioneer Viewer.dc.html)
 *  is not included yet — that is a separate design file and a later phase.
 */
import { useEffect, useRef, useState } from 'react'

import { defaultDir, pickDefaultDir, setDefaultDir } from '../lib/backend'
import { InfoDot } from './InfoDot'
import { TITLEBAR_H } from '../lib/styles'
import { THEMES, type ThemeId } from '../lib/theme'
import type { CommandId, Job } from '../lib/types'

/** The sidebar's small square icon buttons: config and collapse. Replaces the
 *  full-width "Collapse" pill, which was the heaviest thing in the sidebar for
 *  something pressed once a session. */
const iconBtn: React.CSSProperties = {
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  width: 30,
  height: 30,
  flex: 'none',
  padding: 0,
  border: 'none',
  borderRadius: 8,
  background: 'none',
  color: 'var(--pio-nav-fg-dim)',
  cursor: 'pointer',
}

const DOT_COLORS: Record<Job['status'], string> = {
  queued: '#F59E0B',
  running: 'var(--pio-accent-soft)',
  done: '#10B981',
  failed: '#DC2626',
  cancelled: 'var(--pio-nav-fg-faint)',
  interrupted: '#F59E0B',
}

/** Everything about a run worth matching a search against.
 *
 *  Not just the name and output folder: the input paths are usually what
 *  someone remembers a run by — "the one on the yeast data", "the one that
 *  used the phospho library". Every path in the snapshot is included, so a
 *  match on any part of any of them finds the run.
 */
function searchableText(j: Job): string {
  const parts: string[] = [j.title, CMD_TEXT[j.cmd], j.target, String(j.runNo), j.status]
  const s = j.snapshot
  if (s.cmd === 'searchdia') {
    parts.push(s.search.msData, s.search.library, s.search.results)
  } else if (s.cmd === 'buildspeclib') {
    parts.push(s.build.libPath, s.build.calibrationFile)
    parts.push(...s.build.fastaFiles.map((f) => f.path))
    parts.push(...s.build.fastaFiles.map((f) => f.name))
    // The modifications a library was built with, so "phospho" finds it.
    parts.push(...s.build.variableMods.map((m) => m.label))
    parts.push(...s.build.fixedMods.map((m) => m.label))
    parts.push(s.build.predictionModel)
  } else if (s.cmd === 'downloadspeclib') {
    parts.push(s.download.selected, s.download.dest)
  } else {
    parts.push(s.convert.input, s.convert.outputDir)
  }
  return parts.filter(Boolean).join(' \u0000 ').toLowerCase()
}

/** "12 Aug 2026, 14:32" in the viewer's locale, or '' when unknown.
 *
 *  Runs restored from before completion times were recorded have 0, and so do
 *  runs still going — both should say nothing rather than 1 Jan 1970. */
function finishedText(unixSeconds: number): string {
  if (!unixSeconds) return ''
  return new Date(unixSeconds * 1000).toLocaleString(undefined, {
    day: 'numeric',
    month: 'short',
    year: 'numeric',
    hour: '2-digit',
    minute: '2-digit',
  })
}

/** The row's hover text. Free to be verbose — a title costs nothing until it is
 *  hovered, and collapsed the row is a bare dot with no other explanation. */
function rowTitle(j: Job): string {
  const when = finishedText(j.finishedAt)
  return [
    `${j.runNo ? `Run ${j.runNo} · ` : ''}${j.title}`,
    `${CMD_TEXT[j.cmd]} · ${STATUS_TEXT[j.status]}`,
    when && `Finished ${when}`,
    j.target,
  ]
    .filter(Boolean)
    .join('\n')
}

const CMD_TEXT: Record<CommandId, string> = {
  searchdia: 'SearchDIA',
  buildspeclib: 'BuildSpecLib',
  downloadspeclib: 'DownloadSpecLib',
  convertraw: 'ConvertRAW',
}

const STATUS_TEXT: Record<Job['status'], string> = {
  queued: 'Queued',
  running: 'Running…',
  done: 'Completed',
  failed: 'Failed',
  cancelled: 'Cancelled',
  interrupted: 'Interrupted — click to reload its parameters',
}

/** The indeterminate "barber pole" bar the design uses while a job runs — an
 *  honest choice, since Pioneer's stdout gives no reliable percentage. */
const barberStyle: React.CSSProperties = {
  position: 'absolute',
  inset: 0,
  borderRadius: 2,
  background:
    'linear-gradient(45deg,var(--pio-accent-soft) 25%,var(--pio-accent-softer) 25%,var(--pio-accent-softer) 50%,var(--pio-accent-soft) 50%,var(--pio-accent-soft) 75%,var(--pio-accent-softer) 75%)',
  backgroundSize: '20px 20px',
  animation: 'pio-barber .6s linear infinite',
}

/** The settings popover under the gear: where Browse starts, and the theme.
 *
 *  Both were previously loose in the sidebar — the folder behind a bare icon
 *  with no way to see its current value, the theme as a lone swatch above the
 *  collapse button. Neither is touched often enough to earn standing space.
 */
function SettingsPanel({ theme, onTheme }: { theme: ThemeId; onTheme: (id: ThemeId) => void }) {
  const [dir, setDir] = useState(defaultDir())
  const heading: React.CSSProperties = {
    fontSize: 10,
    fontWeight: 700,
    letterSpacing: '0.07em',
    textTransform: 'uppercase',
    color: 'var(--pio-nav-fg-faint)',
    marginBottom: 7,
  }
  return (
    <div
      role="dialog"
      aria-label="Settings"
      style={{
        position: 'absolute',
        top: 36,
        left: 0,
        zIndex: 40,
        width: 232,
        padding: '13px 14px 14px',
        borderRadius: 11,
        background: 'var(--pio-nav)',
        border: '1px solid var(--pio-nav-hair-strong)',
        boxShadow: '0 12px 32px rgba(8,12,20,0.45)',
      }}
    >
      <div style={heading}>Browse from</div>
      <div
        title={dir || 'Your home folder'}
        style={{
          fontSize: 11.5,
          fontFamily: "'IBM Plex Mono'",
          color: 'var(--pio-nav-fg)',
          overflow: 'hidden',
          textOverflow: 'ellipsis',
          whiteSpace: 'nowrap',
          direction: 'rtl',
          textAlign: 'left',
          marginBottom: 8,
        }}
      >
        {/* rtl so a long path truncates at the front, keeping the leaf visible */}
        {dir || 'Home folder'}
      </div>
      <div style={{ display: 'flex', gap: 6 }}>
        <button
          type="button"
          className="pio-ghost"
          onClick={async () => {
            const picked = await pickDefaultDir()
            if (picked) setDir(picked)
          }}
          style={panelBtn}
        >
          Change…
        </button>
        {dir && (
          <button
            type="button"
            className="pio-ghost"
            onClick={() => {
              setDefaultDir('')
              setDir('')
            }}
            title="Go back to starting at your home folder"
            style={panelBtn}
          >
            Reset
          </button>
        )}
      </div>

      <div style={{ ...heading, marginTop: 15 }}>Theme</div>
      <div style={{ display: 'flex', gap: 7, flexWrap: 'wrap' }}>
        {THEMES.map((t) => {
          const active = t.id === theme
          return (
            <button
              key={t.id}
              type="button"
              onClick={() => onTheme(t.id)}
              title={t.label}
              aria-label={`${t.label} theme`}
              aria-pressed={active}
              style={{
                width: 20,
                height: 20,
                padding: 0,
                borderRadius: '50%',
                cursor: 'pointer',
                // Half surface, half accent: previews both the sidebar colour
                // and what the Run button becomes.
                background: `linear-gradient(135deg, ${t.swatch} 50%, ${t.accent} 50%)`,
                border: active
                  ? '2px solid rgba(255,255,255,0.9)'
                  : '1px solid var(--pio-nav-hair-strong)',
                boxShadow: active ? '0 0 0 2px rgba(0,0,0,0.4)' : 'none',
              }}
            />
          )
        })}
      </div>
    </div>
  )
}

const panelBtn: React.CSSProperties = {
  padding: '5px 10px',
  borderRadius: 7,
  border: '1px solid var(--pio-nav-hair-strong)',
  background: 'none',
  color: 'var(--pio-nav-fg)',
  font: "600 11.5px 'IBM Plex Sans'",
  cursor: 'pointer',
}

interface NavItemProps {
  id: CommandId
  title: string
  subtitle: string
  chip: string
  active: boolean
  collapsed: boolean
  disabled?: boolean
  icon: React.ReactNode
  onClick: (id: CommandId) => void
}

function NavItem({
  id,
  title,
  subtitle,
  chip,
  active,
  collapsed,
  disabled,
  icon,
  onClick,
}: NavItemProps) {
  const base: React.CSSProperties = {
    display: 'flex',
    alignItems: 'center',
    gap: 11,
    width: '100%',
    padding: collapsed ? '11px 0' : '9px 11px',
    borderRadius: 9,
    border: 'none',
    cursor: disabled ? 'not-allowed' : 'pointer',
    textAlign: 'left',
    transition: 'background .12s',
    justifyContent: collapsed ? 'center' : undefined,
    opacity: disabled ? 0.45 : 1,
    ...(active
      ? {
          background: 'var(--pio-accent-wash)',
          color: 'var(--pio-nav-fg)',
          paddingLeft: 14,
          boxShadow: 'inset 4px 0 0 var(--pio-accent-soft), inset 6px 0 0 rgba(10,18,35,0.6)',
        }
      : { background: 'none', color: 'var(--pio-nav-fg)' }),
  }
  return (
    <button
      type="button"
      className={active ? 'pio-nav-active' : 'pio-nav-item'}
      onClick={() => !disabled && onClick(id)}
      style={base}
      title={disabled ? `${title} — not implemented yet` : title}
    >
      {icon}
      {!collapsed && (
        <span
          style={{
            display: 'flex',
            flexDirection: 'column',
            alignItems: 'flex-start',
            lineHeight: 1.25,
          }}
        >
          <span style={{ fontSize: 13.5, fontWeight: 600 }}>{title}</span>
          <span style={{ fontSize: 11, color: 'var(--pio-nav-fg-dim)' }}>{subtitle}</span>
        </span>
      )}
      {!collapsed && (
        <span
          style={{
            marginLeft: 'auto',
            fontSize: 10,
            fontWeight: 600,
            color: 'var(--pio-nav-fg-faint)',
            border: '1px solid var(--pio-nav-hair)',
            borderRadius: 5,
            padding: '2px 5px',
            whiteSpace: 'nowrap',
          }}
        >
          {chip}
        </span>
      )}
    </button>
  )
}

interface Props {
  collapsed: boolean
  selected: CommandId
  jobs: Job[]
  viewJobId: string | null
  modKey: string
  theme: ThemeId
  onTheme: (id: ThemeId) => void
  onSelect: (id: CommandId) => void
  onToggleCollapsed: () => void
  onViewJob: (id: string) => void
  onJobAction: (id: string, kind: 'cancel' | 'delete') => void
  /** Move the queued job `dragId` into the position currently held by `dropId`. */
  onReorderQueued: (dragId: string, dropId: string) => void
  /** Give a run a different name. Collisions are resolved by the caller. */
  onRenameJob: (id: string, title: string) => void
}

export function Sidebar({
  collapsed,
  selected,
  jobs,
  viewJobId,
  modKey,
  theme,
  onTheme,
  onSelect,
  onToggleCollapsed,
  onViewJob,
  onJobAction,
  onReorderQueued,
  onRenameJob,
}: Props) {
  const [themeOpen, setThemeOpen] = useState(false)
  const [query, setQuery] = useState('')
  /** The queued job being dragged, and the one it is currently over.
   *
   *  Pointer events rather than HTML5 drag-and-drop: Tauri handles the OS drop
   *  natively (`dragDropEnabled`, which the file-drop-onto-fields support
   *  depends on), so `dragstart`/`drop` never reach the webview. See
   *  `lib/dragdrop.ts`. */
  const [dragJobId, setDragJobId] = useState<string | null>(null)
  const [dropJobId, setDropJobId] = useState<string | null>(null)
  /** The run whose name is being edited, and the text so far. Held apart from
   *  `jobs` so abandoning an edit leaves the name alone. */
  const [editingJobId, setEditingJobId] = useState<string | null>(null)
  const [editingTitle, setEditingTitle] = useState('')

  /** Runs that were pending at some point since the sidebar was last open.
   *
   *  Collapsed, the sidebar is a strip of status dots, and a run you were
   *  watching should not vanish the instant it finishes — you want to see it go
   *  green. But the whole history would make the strip useless, so only these
   *  stay. Opening the sidebar shows them in the history proper, which is the
   *  acknowledgement: the set is cleared, and collapsing again leaves them
   *  behind.
   *
   *  A ref rather than state: every change that matters also changes `jobs` or
   *  `collapsed`, so a render is already coming. */
  const watched = useRef(new Set<string>())

  useEffect(() => {
    if (!collapsed) {
      watched.current.clear()
      return
    }
    // Anything pending now, plus anything that becomes pending while collapsed.
    // Runs on collapse too, which is what captures the queue as it stood then.
    for (const j of jobs) {
      if (j.status === 'queued' || j.status === 'running') watched.current.add(j.id)
    }
  }, [collapsed, jobs])
  const themeRef = useRef<HTMLDivElement>(null)

  // Close on a click anywhere else, so an accidental open does not leave the row
  // expanded for the rest of the session.
  useEffect(() => {
    if (!themeOpen) return
    const onDown = (e: MouseEvent) => {
      if (!themeRef.current?.contains(e.target as Node)) setThemeOpen(false)
    }
    const onKey = (e: KeyboardEvent) => {
      if (e.key === 'Escape') setThemeOpen(false)
    }
    document.addEventListener('mousedown', onDown)
    document.addEventListener('keydown', onKey)
    return () => {
      document.removeEventListener('mousedown', onDown)
      document.removeEventListener('keydown', onKey)
    }
  }, [themeOpen])

  const sectionStyle: React.CSSProperties | undefined = collapsed
    ? { display: 'none' }
    : {
        padding: '16px 12px 8px',
        fontSize: 10.5,
        letterSpacing: '0.08em',
        textTransform: 'uppercase',
        color: 'var(--pio-nav-fg-faint)',
        fontWeight: 600,
      }

  // Pending runs and finished runs are different things to look at: the queue is
  // what is about to happen, history is what to go back to. Same row renderer for
  // both, so a finished run stays clickable and still loads its parameters.
  const queue = jobs.filter((j) => j.status === 'queued' || j.status === 'running')
  // Newest first. Sorted here rather than relying on the order rows come back
  // in, because a run that finishes during this session is appended to the end
  // of `jobs` and would otherwise land at the bottom of the history even though
  // it is the most recent. Copy before sorting — `jobs` is props.
  //
  // The queue keeps insertion order: there, position is the point.
  const history = jobs
    .filter((j) => j.status !== 'queued' && j.status !== 'running')
    .slice()
    .sort((a, b) => b.runNo - a.runNo)

  /** What the collapsed strip shows: everything pending, plus whatever finished
   *  while it was collapsed — in `jobs` order, which is the order they were
   *  queued in.
   *
   *  One list rather than pending-then-finished, so a run that goes green stays
   *  where it was instead of jumping to the bottom. It only leaves the strip
   *  when the sidebar is opened, which is when it joins the history proper. */
  const collapsedRows = jobs.filter(
    (j) => j.status === 'queued' || j.status === 'running' || watched.current.has(j.id),
  )
  // Every term must appear somewhere, so "yeast phospho" narrows rather than
  // widening. Matched case-insensitively across the whole of searchableText.
  const terms = query.trim().toLowerCase().split(/\s+/).filter(Boolean)
  const shown = terms.length
    ? history.filter((j) => {
        const hay = searchableText(j)
        return terms.every((t) => hay.includes(t))
      })
    : history
  const emptyHint: React.CSSProperties = {
    padding: '2px 10px 8px',
    fontSize: 11.5,
    color: 'var(--pio-nav-fg-faint)',
  }

  /** Whether the queue is worth reordering: more than one job waiting, and room
   *  to show a handle. Collapsed, the strip is status dots and there is nothing
   *  to grab. */
  const reorderable = !collapsed && jobs.filter((j) => j.status === 'queued').length > 1

  const renderRow = (j: Job) => {
            const running = j.status === 'running'
            const pending = running || j.status === 'queued'
            // Only a queued job can move. A running one has already started and
            // a finished one has no position left to change, so neither is a
            // drag source nor a drop target.
            const movable = reorderable && j.status === 'queued'
            const dragging = dragJobId === j.id
            const dropBefore = dropJobId === j.id && dragJobId !== null && dragJobId !== j.id
            return (
              <div
                key={j.id}
                title={rowTitle(j)}
                className="pio-job pio-row-hover"
                data-job-id={j.id}
                data-job-movable={movable ? '1' : undefined}
                style={{
                  display: 'flex',
                  alignItems: 'center',
                  gap: 10,
                  width: '100%',
                  padding: '8px 11px',
                  border: 'none',
                  borderRadius: 9,
                  background: j.id === viewJobId ? 'var(--pio-accent-wash-strong)' : 'none',
                  // The row being carried fades; the one under the pointer grows
                  // a line along the edge the drop will land on.
                  opacity: dragging ? 0.4 : 1,
                  boxShadow: dropBefore ? 'inset 0 2px 0 0 var(--pio-accent-soft)' : undefined,
                }}
              >
                {reorderable && (
                  <span
                    onPointerDown={
                      movable
                        ? (e) => {
                            // preventDefault stops the press from starting a
                            // text selection that would follow the pointer.
                            e.preventDefault()
                            e.currentTarget.setPointerCapture(e.pointerId)
                            setDragJobId(j.id)
                            setDropJobId(j.id)
                          }
                        : undefined
                    }
                    onPointerMove={
                      movable
                        ? (e) => {
                            if (dragJobId !== j.id) return
                            // Capture routes every move here regardless of what
                            // is under the pointer, so the row being hovered has
                            // to be hit-tested rather than read from the event.
                            const under = document
                              .elementFromPoint(e.clientX, e.clientY)
                              ?.closest('[data-job-id]') as HTMLElement | null
                            const id = under?.dataset.jobId
                            if (id && under?.dataset.jobMovable === '1') setDropJobId(id)
                          }
                        : undefined
                    }
                    onPointerUp={
                      movable
                        ? () => {
                            if (dragJobId && dropJobId && dragJobId !== dropJobId) {
                              onReorderQueued(dragJobId, dropJobId)
                            }
                            setDragJobId(null)
                            setDropJobId(null)
                          }
                        : undefined
                    }
                    onPointerCancel={() => {
                      setDragJobId(null)
                      setDropJobId(null)
                    }}
                    title={movable ? 'Drag to reorder' : undefined}
                    aria-label={movable ? `Reorder ${j.title}` : undefined}
                    style={{
                      flex: 'none',
                      display: 'flex',
                      alignItems: 'center',
                      // Held at low opacity rather than revealed on hover: a
                      // handle nobody can see is a feature nobody finds.
                      opacity: movable ? 0.45 : 0,
                      cursor: movable ? (dragging ? 'grabbing' : 'grab') : 'default',
                      color: 'var(--pio-nav-fg-dim)',
                    }}
                  >
                    <svg width="10" height="14" viewBox="0 0 10 14" aria-hidden="true">
                      {[2, 7, 12].map((cy) =>
                        [2, 8].map((cx) => (
                          <circle key={`${cx}-${cy}`} cx={cx} cy={cy} r="1.2" fill="currentColor" />
                        )),
                      )}
                    </svg>
                  </span>
                )}
                <div
                  onClick={() => onViewJob(j.id)}
                  style={{
                    flex: 1,
                    minWidth: 0,
                    display: 'flex',
                    alignItems: 'center',
                    gap: 10,
                    cursor: 'pointer',
                  }}
                >
                  <span
                    style={{
                      width: 9,
                      height: 9,
                      borderRadius: '50%',
                      flex: 'none',
                      background: DOT_COLORS[j.status],
                      animation: running ? 'pio-pulse 1.1s ease-in-out infinite' : undefined,
                    }}
                  />
                  {!collapsed && (
                    <span
                      style={{
                        display: 'flex',
                        flexDirection: 'column',
                        alignItems: 'flex-start',
                        lineHeight: 1.25,
                      }}
                    >
                      {editingJobId === j.id ? (
                        <input
                          autoFocus
                          value={editingTitle}
                          onChange={(e) => setEditingTitle(e.target.value)}
                          // The row opens the run on click; typing in the field
                          // must not also do that.
                          onClick={(e) => e.stopPropagation()}
                          onKeyDown={(e) => {
                            e.stopPropagation()
                            if (e.key === 'Enter') {
                              onRenameJob(j.id, editingTitle)
                              setEditingJobId(null)
                            } else if (e.key === 'Escape') {
                              setEditingJobId(null)
                            }
                          }}
                          // Commit on blur as well as Enter: clicking away is a
                          // reasonable way to mean "done", and losing the edit
                          // there would be the more surprising outcome. Escape
                          // clears the id first, so it does not commit.
                          onBlur={() => {
                            if (editingJobId === j.id) onRenameJob(j.id, editingTitle)
                            setEditingJobId(null)
                          }}
                          style={{
                            width: 120,
                            padding: '1px 5px',
                            border: '1px solid var(--pio-nav-hair-strong)',
                            borderRadius: 5,
                            background: 'rgba(0,0,0,0.25)',
                            color: 'var(--pio-nav-fg)',
                            font: "600 12.5px 'IBM Plex Sans'",
                            outline: 'none',
                          }}
                        />
                      ) : (
                        <span
                          style={{
                            fontSize: 12.5,
                            fontWeight: 600,
                            color: 'var(--pio-nav-fg)',
                            maxWidth: 120,
                            overflow: 'hidden',
                            textOverflow: 'ellipsis',
                            whiteSpace: 'nowrap',
                          }}
                        >
                          {j.title}
                        </span>
                      )}
                      {/* The descriptor is unconditional. It used to be the
                          else-branch of `running`, so the bar replaced it and a
                          running row was the one place the command and status
                          were not written down -- while being the row most
                          worth identifying. The bar now sits under it. */}
                      <span style={{ fontSize: 10.5, color: 'var(--pio-nav-fg-dim)' }}>
                        {CMD_TEXT[j.cmd]} · {STATUS_TEXT[j.status]}
                      </span>
                      {running && (
                        <div
                          style={{
                            position: 'relative',
                            overflow: 'hidden',
                            height: 3,
                            width: 96,
                            maxWidth: '100%',
                            marginTop: 5,
                            borderRadius: 2,
                            background: 'var(--pio-nav-hair)',
                          }}
                        >
                          <div style={barberStyle} />
                        </div>
                      )}
                    </span>
                  )}
                  {!collapsed && (
                    <span
                      style={{
                        marginLeft: 'auto',
                        flex: 'none',
                        fontSize: 11,
                        fontWeight: 700,
                        fontFamily: "'IBM Plex Mono'",
                        color: 'var(--pio-nav-fg-faint)',
                        minWidth: 16,
                        textAlign: 'right',
                      }}
                    >
                      {/* Handed out when the run is queued, so a queued job
                          already shows the number it will carry into the
                          history -- queue three after run 41 and they are 42,
                          43, 44, including if one is cancelled before it
                          starts.

                          While a run is still waiting the number is its place
                          in line: reordering the queue redistributes the
                          numbers so they keep ascending down it, rather than
                          leaving the column reading 44, 42, 43. Once a run
                          starts, its number is fixed for life -- that is what
                          the history is sorted and searched by. */}
                      {j.runNo || ''}
                    </span>
                  )}
                </div>
                {/* Queued, running or finished alike: the name is a label, and
                    the run worth renaming is often the one that already
                    finished. Hidden until the row is hovered, like the other
                    row actions. */}
                {!collapsed && editingJobId !== j.id && (
                  <button
                    type="button"
                    className="pio-jobact pio-iconbtn"
                    onClick={() => {
                      setEditingJobId(j.id)
                      setEditingTitle(j.title)
                    }}
                    title="Rename"
                    style={{
                      flex: 'none',
                      border: 'none',
                      background: 'none',
                      cursor: 'pointer',
                      padding: 3,
                      color: 'var(--pio-nav-fg-dim)',
                      display: 'flex',
                    }}
                  >
                    <svg width="13" height="13" viewBox="0 0 24 24" fill="none">
                      <path
                        d="M4 20h4L19 9a2.1 2.1 0 0 0-3-3L5 17v3Z"
                        stroke="currentColor"
                        strokeWidth="1.8"
                        strokeLinejoin="round"
                      />
                    </svg>
                  </button>
                )}
                {pending && (
                  <button
                    type="button"
                    className="pio-jobact pio-iconbtn"
                    onClick={() => onJobAction(j.id, 'cancel')}
                    title={running ? 'Cancel run' : 'Remove from queue'}
                    style={{
                      flex: 'none',
                      border: 'none',
                      background: 'none',
                      cursor: 'pointer',
                      padding: 3,
                      color: '#F0A8A8',
                      display: 'flex',
                    }}
                  >
                    <svg width="14" height="14" viewBox="0 0 24 24" fill="none">
                      <rect x="6" y="6" width="12" height="12" rx="2" fill="currentColor" />
                    </svg>
                  </button>
                )}
                {!pending && (
                  <button
                    type="button"
                    className="pio-jobact pio-iconbtn"
                    onClick={() => onJobAction(j.id, 'delete')}
                    title="Delete from history"
                    style={{
                      flex: 'none',
                      border: 'none',
                      background: 'none',
                      cursor: 'pointer',
                      padding: 3,
                      color: 'var(--pio-nav-fg-faint)',
                      display: 'flex',
                    }}
                  >
                    <svg width="14" height="14" viewBox="0 0 24 24" fill="none">
                      <path
                        d="M5 7h14M10 7V5h4v2M6 7l1 13h10l1-13"
                        stroke="currentColor"
                        strokeWidth="1.7"
                        strokeLinecap="round"
                        strokeLinejoin="round"
                      />
                    </svg>
                  </button>
                )}
              </div>
            )
  }

  return (
    <aside
      style={{
        width: collapsed ? 66 : 252,
        flex: 'none',
        background: 'var(--pio-nav)',
        color: 'var(--pio-nav-fg)',
        display: 'flex',
        flexDirection: 'column',
        borderRight: '1px solid var(--pio-nav-border)',
        transition: 'width .16s ease',
      }}
    >
      <div
        data-tauri-drag-region
        style={{
          // On macOS the title bar is an overlay, so the webview starts at y=0
          // and the traffic lights float over this block; TITLEBAR_H is their
          // clearance. It is 0 elsewhere, where the real title bar sits above.
          padding: `${TITLEBAR_H + 8}px 10px 8px`,
          display: 'flex',
          alignItems: 'center',
          justifyContent: collapsed ? 'center' : 'space-between',
          gap: 6,
          borderBottom: '1px solid var(--pio-nav-veil)',
        }}
      >
        {!collapsed && (
          <div ref={themeRef} style={{ position: 'relative' }}>
            <button
              type="button"
              className="pio-navicon"
              onClick={() => setThemeOpen((o) => !o)}
              title="Settings"
              aria-label="Settings"
              aria-expanded={themeOpen}
              style={{ ...iconBtn, color: themeOpen ? 'var(--pio-nav-fg)' : 'var(--pio-nav-fg-dim)' }}
            >
              {/* A cog: hub, ring, and eight radial teeth. Two earlier
                  attempts failed as icons -- a circle with spokes reads as a
                  sun, and an eight-pointed star path reads as an explosion.
                  Teeth as separate strokes off a plain ring is unambiguous at
                  16px. */}
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none">
                <circle cx="12" cy="12" r="6.7" stroke="currentColor" strokeWidth="1.6" />
                <circle cx="12" cy="12" r="2.5" stroke="currentColor" strokeWidth="1.6" />
                <path
                  d="M18.7 12.0L20.9 12.0M16.7 16.7L18.3 18.3M12.0 18.7L12.0 20.9M7.3 16.7L5.7 18.3M5.3 12.0L3.1 12.0M7.3 7.3L5.7 5.7M12.0 5.3L12.0 3.1M16.7 7.3L18.3 5.7"
                  stroke="currentColor"
                  strokeWidth="2.1"
                  strokeLinecap="round"
                />
              </svg>
            </button>
            {themeOpen && (
              <SettingsPanel theme={theme} onTheme={onTheme} />
            )}
          </div>
        )}
        <button
          type="button"
          className="pio-navicon"
          onClick={onToggleCollapsed}
          title={collapsed ? 'Expand sidebar' : 'Collapse sidebar'}
          aria-label={collapsed ? 'Expand sidebar' : 'Collapse sidebar'}
          style={iconBtn}
        >
          {/* The panel-with-a-rail icon every app uses for this, rather than a
              bare chevron: it shows what is being toggled, not just a direction.
              The rail sits on the side the sidebar is on. */}
          <svg width="17" height="17" viewBox="0 0 24 24" fill="none">
            <rect
              x="3"
              y="4.5"
              width="18"
              height="15"
              rx="2.6"
              stroke="currentColor"
              strokeWidth="1.7"
            />
            <path d="M9.2 4.5v15" stroke="currentColor" strokeWidth="1.7" />
          </svg>
        </button>
      </div>

      <div style={sectionStyle}>Workflow</div>
      <nav style={{ display: 'flex', flexDirection: 'column', gap: 3, padding: '0 12px' }}>
        <NavItem
          id="convertraw"
          title="ConvertRAW"
          subtitle=".raw → .arrow"
          chip={`${modKey}1`}
          active={selected === 'convertraw'}
          collapsed={collapsed}
          onClick={onSelect}
          icon={
            // A bare arrow said "something moves right", not "a file becomes
            // another file". Two sheets with the second overlapping the first,
            // and the arrow between them, is the conventional convert glyph.
            <svg width="19" height="19" viewBox="0 0 24 24" fill="none" style={{ flex: 'none' }}>
              <path d="M13.4 3.6H5.6v10.2" stroke="currentColor" strokeWidth="1.6" strokeLinecap="round" />
              <path
                d="M9.4 8.8h5.4l3.3 3.3v8.3H9.4z"
                stroke="currentColor"
                strokeWidth="1.6"
                strokeLinejoin="round"
              />
              <path d="M14.8 8.8v3.3h3.3" stroke="currentColor" strokeWidth="1.6" strokeLinejoin="round" />
              <path
                d="M5.2 16.2h6.4M9.4 14l2.2 2.2-2.2 2.2"
                stroke="currentColor"
                strokeWidth="1.6"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          }
        />
        <NavItem
          id="buildspeclib"
          title="BuildSpecLib"
          subtitle="FASTA → library"
          chip={`${modKey}2`}
          active={selected === 'buildspeclib'}
          collapsed={collapsed}
          onClick={onSelect}
          icon={
            <svg width="19" height="19" viewBox="0 0 24 24" fill="none" style={{ flex: 'none' }}>
              <ellipse cx="12" cy="5.5" rx="7" ry="2.8" stroke="currentColor" strokeWidth="1.6" />
              <path
                d="M5 5.5v13c0 1.55 3.13 2.8 7 2.8s7-1.25 7-2.8v-13"
                stroke="currentColor"
                strokeWidth="1.6"
              />
              <path d="M5 12c0 1.55 3.13 2.8 7 2.8s7-1.25 7-2.8" stroke="currentColor" strokeWidth="1.6" />
            </svg>
          }
        />
        <NavItem
          id="downloadspeclib"
          title="DownloadSpecLib"
          subtitle="Fetch a prebuilt library"
          chip={`${modKey}3`}
          active={selected === 'downloadspeclib'}
          collapsed={collapsed}
          onClick={onSelect}
          icon={
            <svg width="19" height="19" viewBox="0 0 24 24" fill="none" style={{ flex: 'none' }}>
              <path
                d="M12 3v12M7.5 10.5 12 15l4.5-4.5"
                stroke="currentColor"
                strokeWidth="1.7"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
              <path d="M4 19h16" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" />
            </svg>
          }
        />
        <NavItem
          id="searchdia"
          title="SearchDIA"
          subtitle="Find & Quantify Proteins"
          chip={`${modKey}4`}
          active={selected === 'searchdia'}
          collapsed={collapsed}
          onClick={onSelect}
          icon={
            <svg width="19" height="19" viewBox="0 0 24 24" fill="none" style={{ flex: 'none' }}>
              <circle cx="11" cy="11" r="6.2" stroke="currentColor" strokeWidth="1.7" />
              <path d="m20 20-3.7-3.7" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" />
            </svg>
          }
        />
      </nav>

      <div
        style={{
          flex: 1,
          minHeight: 0,
          display: 'flex',
          flexDirection: 'column',
          paddingTop: 2,
        }}
      >
        <div
          style={{
            flex: 1,
            minHeight: 0,
            overflowY: 'auto',
            display: 'flex',
            flexDirection: 'column',
            gap: 2,
            padding: '0 12px',
          }}
        >
          <div style={sectionStyle}>Queue</div>
          {queue.length === 0 && !collapsed && <div style={emptyHint}>Nothing queued.</div>}
          {(collapsed ? collapsedRows : queue).map(renderRow)}
          <div
            style={{
              ...sectionStyle,
              marginTop: 12,
              display: collapsed ? 'none' : 'flex',
              alignItems: 'center',
              gap: 7,
            }}
          >
            History
            {/* What is searchable is not guessable from a box that says
                "Search runs…" — the paths and modifications in particular. */}
            <InfoDot
              tone="dark"
              text={
                'Search matches a run\u2019s name, number, command and status, and every path it used \u2014 MS data folder, spectral library, results folder, FASTA files. ' +
                'For a library build it also matches the modifications and the prediction model, so \u201coxidation\u201d or \u201cprosit\u201d finds one. ' +
                'Several words narrow rather than widen: \u201cyeast phospho\u201d matches only runs with both. ' +
                'Individual MS file names are not searchable \u2014 only the folder they are in was recorded.'
              }
            />
          </div>
          {/* Only once there is enough history to be worth searching — a box
              above two runs is clutter. */}
          {!collapsed && history.length > 4 && (
            <div style={{ padding: '0 12px 6px' }}>
              <input
                value={query}
                onChange={(e) => setQuery(e.target.value)}
                placeholder="Search runs…"
                aria-label="Search run history"
                style={{
                  width: '100%',
                  boxSizing: 'border-box',
                  padding: '6px 9px',
                  borderRadius: 7,
                  border: '1px solid var(--pio-nav-hair-strong)',
                  background: 'var(--pio-nav-veil)',
                  color: 'var(--pio-nav-fg)',
                  font: "12px 'IBM Plex Sans'",
                  outline: 'none',
                }}
              />
            </div>
          )}
          {history.length === 0 && !collapsed && (
            <div style={emptyHint}>No finished runs yet.</div>
          )}
          {history.length > 0 && shown.length === 0 && !collapsed && (
            <div style={emptyHint}>No run matches “{query}”.</div>
          )}
          {/* Collapsed, history is deliberately absent: the strip is for what is
              happening now, plus whatever finished while you were not looking. */}
          {!collapsed && shown.map(renderRow)}
        </div>
      </div>



    </aside>
  )
}
