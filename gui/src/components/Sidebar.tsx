/** The navy sidebar: brand, workflow picker, live job queue, thread picker
 *  and the Run button. Ported from the `<aside>` block of the design.
 *
 *  The design's "Analysis" section (viewer tabs mounting Pioneer Viewer.dc.html)
 *  is not included yet — that is a separate design file and a later phase.
 */
import { useEffect, useRef, useState } from 'react'

import { PioneerLogo } from './PioneerLogo'
import { TITLEBAR_H } from '../lib/styles'
import { THEMES, type ThemeId } from '../lib/theme'
import type { CommandId, Job } from '../lib/types'

const DOT_COLORS: Record<Job['status'], string> = {
  queued: '#F59E0B',
  running: 'var(--pio-accent-soft)',
  done: '#10B981',
  failed: '#DC2626',
  cancelled: 'var(--pio-nav-fg-faint)',
}

const CMD_TEXT: Record<CommandId, string> = {
  searchdia: 'SearchDIA',
  buildspeclib: 'BuildSpecLib',
  convertraw: 'ConvertRAW',
}

const STATUS_TEXT: Record<Job['status'], string> = {
  queued: 'Queued',
  running: 'Running…',
  done: 'Completed',
  failed: 'Failed',
  cancelled: 'Cancelled',
}

/** The indeterminate "barber pole" bar the design uses while a job runs — an
 *  honest choice, since Pioneer's stdout gives no reliable percentage. */
const barberStyle: React.CSSProperties = {
  position: 'absolute',
  inset: 0,
  borderRadius: 2,
  background:
    'repeating-linear-gradient(45deg,var(--pio-accent-soft) 0 7px,var(--pio-accent-softer) 7px 14px)',
  backgroundSize: '28px 100%',
  animation: 'pio-barber .6s linear infinite',
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
}: Props) {
  const [themeOpen, setThemeOpen] = useState(false)
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
  const history = jobs.filter((j) => j.status !== 'queued' && j.status !== 'running')
  const emptyHint: React.CSSProperties = {
    padding: '2px 10px 8px',
    fontSize: 11.5,
    color: 'var(--pio-nav-fg-faint)',
  }

  const renderRow = (j: Job, idx: number) => {
            const running = j.status === 'running'
            const pending = running || j.status === 'queued'
            return (
              <div
                key={j.id}
                className="pio-job pio-row-hover"
                style={{
                  display: 'flex',
                  alignItems: 'center',
                  gap: 10,
                  width: '100%',
                  padding: '8px 11px',
                  border: 'none',
                  borderRadius: 9,
                  background: j.id === viewJobId ? 'var(--pio-accent-wash-strong)' : 'none',
                }}
              >
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
                      {running ? (
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
                      ) : (
                        <span style={{ fontSize: 10.5, color: 'var(--pio-nav-fg-dim)' }}>
                          {CMD_TEXT[j.cmd]} · {STATUS_TEXT[j.status]}
                        </span>
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
                      {idx + 1}
                    </span>
                  )}
                </div>
                {pending && (
                  <button
                    type="button"
                    className="pio-jobact pio-iconbtn"
                    onClick={() => onJobAction(j.id, 'cancel')}
                    title="Cancel run"
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
          // On macOS the title bar is an overlay, so the webview starts at
          // y=0 and the traffic lights float over this block; TITLEBAR_H is
          // their clearance. It is 0 elsewhere, where the real title bar
          // already sits above us.
          padding: `${22 + TITLEBAR_H}px ${collapsed ? 0 : 20}px 20px`,
          display: 'flex',
          alignItems: 'center',
          justifyContent: collapsed ? 'center' : 'flex-start',
          gap: 13,
          borderBottom: '1px solid var(--pio-nav-veil)',
        }}
      >
        <PioneerLogo />
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
            <svg width="17" height="17" viewBox="0 0 24 24" fill="none" style={{ flex: 'none' }}>
              <path
                d="M4 12h12M13 7l5 5-5 5"
                stroke="currentColor"
                strokeWidth="1.7"
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
            <svg width="17" height="17" viewBox="0 0 24 24" fill="none" style={{ flex: 'none' }}>
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
          id="searchdia"
          title="SearchDIA"
          subtitle="Find & Quantify Proteins"
          chip={`${modKey}3`}
          active={selected === 'searchdia'}
          collapsed={collapsed}
          onClick={onSelect}
          icon={
            <svg width="17" height="17" viewBox="0 0 24 24" fill="none" style={{ flex: 'none' }}>
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
          {queue.map(renderRow)}
          <div style={{ ...sectionStyle, marginTop: 12 }}>History</div>
          {history.length === 0 && !collapsed && (
            <div style={emptyHint}>No finished runs yet.</div>
          )}
          {history.map(renderRow)}
        </div>
      </div>


      {/* Collapsed to a single swatch by default. Four permanent circles is a lot
          of standing visual weight for something set once and rarely revisited,
          so the alternatives only appear once asked for. */}
      {!collapsed && (
        <div ref={themeRef} style={{ display: 'flex', alignItems: 'center', gap: 6, padding: '12px 14px 4px' }}>
          {(themeOpen ? THEMES : THEMES.filter((t) => t.id === theme)).map((t) => {
            const active = t.id === theme
            return (
              <button
                key={t.id}
                type="button"
                onClick={() => {
                  if (!themeOpen) setThemeOpen(true)
                  else {
                    onTheme(t.id)
                    setThemeOpen(false)
                  }
                }}
                title={themeOpen ? t.label : `Theme: ${t.label} — click to change`}
                aria-label={themeOpen ? `${t.label} theme` : `Theme: ${t.label}. Click to change.`}
                aria-pressed={themeOpen ? active : undefined}
                style={{
                  width: 18,
                  height: 18,
                  padding: 0,
                  borderRadius: '50%',
                  cursor: 'pointer',
                  // Half surface, half accent: previews both the sidebar colour
                  // and what the Run button becomes.
                  background: `linear-gradient(135deg, ${t.swatch} 50%, ${t.accent} 50%)`,
                  border:
                    themeOpen && active
                      ? '2px solid rgba(255,255,255,0.85)'
                      : '1px solid var(--pio-nav-hair-strong)',
                  boxShadow: themeOpen && active ? '0 0 0 2px rgba(0,0,0,0.35)' : 'none',
                  transition: 'border-color .12s',
                }}
              />
            )
          })}
          {/* Always shown. Closed, the swatch is a lone circle with no hint
              that it does anything — it read as decoration. The word is what
              makes it a control; the tooltip only helps once you suspect it. */}
          <span style={{ fontSize: 10.5, color: 'var(--pio-nav-fg-faint)', marginLeft: 2 }}>theme</span>
        </div>
      )}

      <button
        type="button"
        className="pio-navtoggle"
        onClick={onToggleCollapsed}
        title={collapsed ? 'Expand sidebar' : 'Collapse sidebar'}
        style={{
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          gap: 9,
          margin: '8px 12px 16px',
          padding: 12,
          border: 'none',
          borderRadius: 10,
          background: 'var(--pio-nav-veil)',
          color: 'var(--pio-nav-fg)',
          cursor: 'pointer',
          font: "600 13.5px 'IBM Plex Sans'",
          letterSpacing: '0.01em',
        }}
      >
        <svg
          width="17"
          height="17"
          viewBox="0 0 24 24"
          fill="none"
          style={{
            flex: 'none',
            transition: 'transform .16s',
            transform: collapsed ? 'rotate(180deg)' : undefined,
          }}
        >
          <path
            d="M14 7l-5 5 5 5"
            stroke="currentColor"
            strokeWidth="2"
            strokeLinecap="round"
            strokeLinejoin="round"
          />
        </svg>
        {!collapsed && <span>Collapse</span>}
      </button>
    </aside>
  )
}
