/** The dark log drawer pinned to the bottom of the main pane: drag handle,
 *  title + status, indeterminate progress bar, streamed output, failure banner.
 *  Ported from the LOG DRAWER block of the design.
 */
import { useEffect, useRef } from 'react'
import type { Job, JobStatus, PathInfo } from '../lib/types'

const PILL: React.CSSProperties = {
  fontSize: 11.5,
  fontWeight: 600,
  padding: '4px 11px',
  borderRadius: 20,
}

const drawerPill = (status: JobStatus | 'idle'): React.CSSProperties => {
  const byStatus: Record<string, React.CSSProperties> = {
    running: { background: 'rgba(245,158,11,0.18)', color: '#FBBF24' },
    failed: { background: 'rgba(220,38,38,0.2)', color: '#FCA5A5' },
    done: { background: 'rgba(16,185,129,0.18)', color: '#34D399' },
    queued: { background: 'rgba(245,158,11,0.14)', color: '#FCD9A0' },
  }
  return { ...PILL, ...(byStatus[status] ?? { background: 'rgba(255,255,255,0.08)', color: '#9AA4B0' }) }
}

const GHOST_BTN: React.CSSProperties = {
  display: 'flex',
  alignItems: 'center',
  gap: 6,
  padding: '6px 12px',
  borderRadius: 8,
  border: '1px solid rgba(255,255,255,0.18)',
  background: 'none',
  color: '#C8D0DA',
  font: "600 12px 'IBM Plex Sans'",
  cursor: 'pointer',
}

interface Props {
  job: Job | null
  /** Stat of the run's output folder, or null if not checked yet. */
  targetInfo: PathInfo | null
  height: number
  statusText: string
  confirmCancel: boolean
  onStartDrag: (e: React.MouseEvent) => void
  onAskCancel: () => void
  onKeepRunning: () => void
  onConfirmCancel: () => void
  onClose: () => void
}

export function LogDrawer({
  job,
  targetInfo,
  height,
  statusText,
  confirmCancel,
  onStartDrag,
  onAskCancel,
  onKeepRunning,
  onConfirmCancel,
  onClose,
}: Props) {
  const status: JobStatus | 'idle' = job ? job.status : 'idle'
  const preRef = useRef<HTMLPreElement>(null)
  const pinnedRef = useRef(true)

  const lineCount = job ? job.logLines.length : 0

  // Follow the tail as output arrives, but stop following the moment the user
  // scrolls up to read something — otherwise a long search is unreadable.
  useEffect(() => {
    const el = preRef.current
    if (el && pinnedRef.current) el.scrollTop = el.scrollHeight
  }, [lineCount])

  const onScroll = () => {
    const el = preRef.current
    if (!el) return
    pinnedRef.current = el.scrollHeight - el.scrollTop - el.clientHeight < 24
  }

  const barInner: React.CSSProperties | null =
    status === 'running'
      ? {
          position: 'absolute',
          inset: 0,
          background: 'repeating-linear-gradient(45deg,var(--pio-accent-soft) 0 7px,var(--pio-accent-softer) 7px 14px)',
          backgroundSize: '28px 100%',
          animation: 'pio-barber .6s linear infinite',
        }
      : status === 'done'
        ? { position: 'absolute', inset: 0, background: 'linear-gradient(90deg,#10B981,#34D399)' }
        : status === 'failed'
          ? { position: 'absolute', inset: 0, background: '#DC2626' }
          : null

  const title = job ? `${job.title}${job.target ? ` · ${job.target}` : ''}` : 'Console'

  return (
    <div
      style={{
        // A flex sibling of the form pane, not an overlay. It used to be
        // `position: absolute; bottom: 0`, which floated it over the form: the
        // bottom of the parameter panel became unreachable the moment a run
        // started, because the pane still believed it had the full height and
        // its last rows sat underneath the console.
        flex: 'none',
        height,
        minHeight: 0,
        background: 'var(--pio-nav)',
        borderTop: '1px solid var(--pio-nav-border)',
        display: 'flex',
        flexDirection: 'column',
        boxShadow: '0 -8px 30px rgba(15,20,27,0.18)',
      }}
    >
      <div
        onMouseDown={onStartDrag}
        title="Drag to resize"
        className="pio-drag"
        style={{
          height: 14,
          flex: 'none',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          cursor: 'ns-resize',
          background: 'rgba(255,255,255,0.02)',
        }}
      >
        <div style={{ width: 42, height: 4, borderRadius: 2, background: 'rgba(255,255,255,0.24)' }} />
      </div>

      <div
        style={{
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          padding: '12px 18px',
          borderBottom: '1px solid rgba(255,255,255,0.07)',
        }}
      >
        <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
          <span
            style={{
              fontSize: 12,
              fontWeight: 600,
              color: '#C8D0DA',
              maxWidth: 420,
              overflow: 'hidden',
              textOverflow: 'ellipsis',
              whiteSpace: 'nowrap',
            }}
          >
            {title}
          </span>
          <span style={drawerPill(status)}>{statusText}</span>
        </div>
        <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
          {status === 'running' && !confirmCancel && (
            <button
              type="button"
              className="pio-cancel"
              onClick={onAskCancel}
              style={{ ...GHOST_BTN, color: '#F0A8A8', padding: '6px 13px' }}
            >
              <svg width="12" height="12" viewBox="0 0 24 24" fill="none">
                <rect x="6" y="6" width="12" height="12" rx="2" fill="currentColor" />
              </svg>
              Cancel
            </button>
          )}
          {status === 'running' && confirmCancel && (
            <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
              <span style={{ fontSize: 12, color: '#F0A8A8', fontWeight: 600 }}>Stop this run?</span>
              <button
                type="button"
                className="pio-ghost"
                onClick={onKeepRunning}
                style={{ ...GHOST_BTN, padding: '6px 12px' }}
              >
                Keep running
              </button>
              <button
                type="button"
                className="pio-danger"
                onClick={onConfirmCancel}
                style={{
                  padding: '6px 12px',
                  borderRadius: 8,
                  border: 'none',
                  background: '#DC2626',
                  color: '#fff',
                  font: "600 12px 'IBM Plex Sans'",
                  cursor: 'pointer',
                }}
              >
                Stop
              </button>
            </div>
          )}
          <button type="button" className="pio-ghost" onClick={onClose} title="Hide" style={GHOST_BTN}>
            <svg width="14" height="14" viewBox="0 0 24 24" fill="none">
              <path
                d="M6 9l6 6 6-6"
                stroke="currentColor"
                strokeWidth="2"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
            Hide
          </button>
        </div>
      </div>

      <div
        style={{
          position: 'relative',
          overflow: 'hidden',
          height: 5,
          background: 'rgba(255,255,255,0.14)',
          flex: 'none',
        }}
      >
        {barInner && <div style={barInner} />}
      </div>

      <pre
        ref={preRef}
        onScroll={onScroll}
        style={{
          flex: 1,
          margin: 0,
          overflowY: 'auto',
          padding: '14px 18px',
          font: "12.5px/1.65 'IBM Plex Mono'",
          color: '#AEC2E6',
          whiteSpace: 'pre-wrap',
          wordBreak: 'break-word',
        }}
      >
        {job && job.logLines.length === 0 && status !== 'running' && status !== 'queued' ? (
          // A run restored from a previous session: history keeps the parameters
          // so they can be recalled, but not the output, which is far too large
          // for localStorage. Say so rather than showing an empty pane that looks
          // like a bug.
          <span style={{ color: '#6E7E97', fontStyle: 'italic' }}>
            Output is not kept between sessions — this run's parameters were
            restored from history.
            {targetInfo === null
              ? ` Its results were written to ${job.target}.`
              : targetInfo.exists
                ? ` Its results are still at ${job.target}.`
                : ` Its results are no longer at ${job.target} — the folder has been moved or deleted.`}
          </span>
        ) : (
          job ? job.logLines.map((l) => l.text).join('\n') : ''
        )}
        {status === 'running' && (
          <span style={{ color: 'var(--pio-accent)', animation: 'pio-blink 1s step-end infinite' }}>▋</span>
        )}
      </pre>

      {status === 'failed' && (
        <div
          style={{
            flex: 'none',
            display: 'flex',
            alignItems: 'center',
            gap: 9,
            padding: '11px 18px',
            background: 'rgba(220,38,38,0.14)',
            borderTop: '1px solid rgba(220,38,38,0.3)',
          }}
        >
          <svg width="16" height="16" viewBox="0 0 24 24" fill="none" style={{ flex: 'none' }}>
            <circle cx="12" cy="12" r="9" stroke="#FCA5A5" strokeWidth="1.8" />
            <path d="M12 7v6M12 16v.5" stroke="#FCA5A5" strokeWidth="2" strokeLinecap="round" />
          </svg>
          <span style={{ fontSize: 12.5, color: '#FCA5A5', fontWeight: 600 }}>
            Run failed — {job?.failMsg}
          </span>
        </div>
      )}
    </div>
  )
}
