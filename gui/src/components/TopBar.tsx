/** The controls pinned to the header: threads, the .json editor link, run
 *  status, and the primary Run button.
 *
 *  These used to live in two places at once — a card at the bottom of the form
 *  (RunBar) and again in the sidebar — so there was no single obvious Run
 *  button, and the thread count could be set from either. Worse, the form-bottom
 *  copy was only reachable by scrolling, which stopped working entirely once the
 *  log drawer opened. Pinning them to the header makes them reachable at any
 *  scroll position and gives each control exactly one home.
 */
import type React from 'react'

export interface TopBarProps {
  /** ConvertRAW has no thread control: PioneerConverter takes its own flag. */
  showThreads: boolean
  threads: number
  maxThreads: number
  runLabel: string
  /** "Edit the .json directly", or "View the command line" for ConvertRAW. */
  previewLabel: string
  onThreads: (n: number) => void
  onEditJson: () => void
  onRun: () => void
}

const stepBtn: React.CSSProperties = {
  width: 26,
  height: 28,
  border: 'none',
  background: '#F8FAFB',
  cursor: 'pointer',
  font: "700 14px 'IBM Plex Sans'",
  color: '#344054',
  lineHeight: 1,
}

export function TopBar({
  showThreads,
  threads,
  maxThreads,
  runLabel,
  previewLabel,
  onThreads,
  onEditJson,
  onRun,
}: TopBarProps) {
  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: 14, flex: 'none' }}>
      {/* Icon only. The label ("Edit the .json directly") was the widest thing in
          the header and pushed Run off at the narrow default width; the tooltip
          carries the wording instead.

          There was also a status pill here, mirroring the selected job's state.
          It came out: the log drawer shows the identical string and the queue
          shows a status dot per job, so it was a third copy whose idle text
          ("Ready") read as "ready to run" next to the Run button. */}
      <button
        type="button"
        className="pio-link"
        onClick={onEditJson}
        title={previewLabel}
        aria-label={previewLabel}
        style={{
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          width: 30,
          height: 30,
          borderRadius: 8,
          background: 'none',
          border: '1px solid #E5E9ED',
          cursor: 'pointer',
          color: '#667085',
          padding: 0,
          flex: 'none',
        }}
      >
        <svg width="15" height="15" viewBox="0 0 24 24" fill="none">
          <path
            d="m13 5 6 6M4 20l1-4L16 5l3 3L8 19l-4 1Z"
            stroke="currentColor"
            strokeWidth="1.7"
            strokeLinecap="round"
            strokeLinejoin="round"
          />
        </svg>
      </button>

      {showThreads && (
        <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
          <span style={{ fontSize: 12.5, color: '#667085', fontWeight: 600 }}>Threads</span>
          <div
            style={{
              display: 'flex',
              alignItems: 'center',
              border: '1px solid #D7DBE0',
              borderRadius: 8,
              overflow: 'hidden',
              flex: 'none',
            }}
          >
            <button
              type="button"
              className="pio-browse"
              onClick={() => onThreads(Math.max(1, threads - 1))}
              title="Fewer threads"
              style={stepBtn}
            >
              −
            </button>
            <input
              value={threads}
              onChange={(e) => {
                const n = parseInt(e.target.value.replace(/[^0-9]/g, ''), 10)
                if (Number.isFinite(n)) onThreads(Math.max(1, Math.min(maxThreads, n)))
              }}
              title={`${maxThreads} available · max ${maxThreads}`}
              style={{
                width: 40,
                textAlign: 'center',
                border: 'none',
                borderLeft: '1px solid #E5E9ED',
                borderRight: '1px solid #E5E9ED',
                padding: '6px 0',
                font: "600 13px 'IBM Plex Mono'",
                color: '#1A2230',
                outline: 'none',
              }}
            />
            <button
              type="button"
              className="pio-browse"
              onClick={() => onThreads(Math.min(maxThreads, threads + 1))}
              title="More threads"
              style={stepBtn}
            >
              +
            </button>
          </div>
        </div>
      )}

      <button
        type="button"
        className="pio-run"
        onClick={onRun}
        style={{
          display: 'flex',
          alignItems: 'center',
          gap: 8,
          padding: '10px 18px',
          border: 'none',
          borderRadius: 10,
          cursor: 'pointer',
          font: "700 13.5px 'IBM Plex Sans'",
          color: '#FFFFFF',
          background: 'linear-gradient(135deg,#3E62A0,#2E4D7E)',
          boxShadow: '0 4px 12px rgba(46,77,126,0.28)',
          flex: 'none',
        }}
      >
        <svg width="14" height="14" viewBox="0 0 24 24" fill="#FFFFFF">
          <path d="M7 4.5v15l13-7.5z" />
        </svg>
        {runLabel}
      </button>
    </div>
  )
}
