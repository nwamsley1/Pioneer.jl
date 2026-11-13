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
  /** Short text on the button itself; previewLabel is the tooltip. */
  previewText: string
  onThreads: (n: number) => void
  onEditJson: () => void
  onRun: () => void
}

/** One height for every header control: the .json button, the threads
 *  stepper (inner height plus its 1px border) and Run. */
const CONTROL_H = 32

const stepBtn: React.CSSProperties = {
  width: 26,
  height: CONTROL_H - 2,
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
  previewText,
  onThreads,
  onEditJson,
  onRun,
}: TopBarProps) {
  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: 14, flex: 'none' }}>
      {/* A short text label (".json" or "Command"): an icon alone did not say
          what it opened. The full wording ("Edit the .json directly") was the
          widest thing in the header and pushed Run off at the narrow default
          width, so it stays in the tooltip.

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
          height: CONTROL_H,
          boxSizing: 'border-box',
          borderRadius: 8,
          background: 'none',
          border: '1px solid #E5E9ED',
          cursor: 'pointer',
          color: '#667085',
          padding: '0 10px',
          font: "500 12px 'IBM Plex Mono'",
          flex: 'none',
        }}
      >
        {previewText}
      </button>

      {showThreads && (
        <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
          <span style={{ fontSize: 12.5, color: '#667085', fontWeight: 600 }}>Threads</span>
          <div
            style={{
              display: 'flex',
              alignItems: 'center',
              border: `1px solid ${threads === 0 ? '#E5484D' : '#D7DBE0'}`,
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
              data-key="threads"
              // Empty is a state the field must be able to hold. Rejecting a
              // value that does not parse meant backspacing to nothing never
              // reached onThreads, so the old number came straight back and the
              // field could not be cleared to retype. 0 already meant "not
              // determined" and nothing after startup reads it, so empty
              // reuses it; Run refuses while it is set.
              value={threads === 0 ? '' : String(threads)}
              onChange={(e) => {
                const digits = e.target.value.replace(/[^0-9]/g, '')
                if (!digits) return onThreads(0)
                onThreads(Math.max(1, Math.min(maxThreads, parseInt(digits, 10))))
              }}
              title={`${maxThreads} available · max ${maxThreads}`}
              style={{
                width: 40,
                textAlign: 'center',
                border: 'none',
                borderLeft: '1px solid #E5E9ED',
                borderRight: '1px solid #E5E9ED',
                height: CONTROL_H - 2,
                padding: 0,
                font: "600 13px 'IBM Plex Mono'",
                color: threads === 0 ? '#C0392B' : '#1A2230',
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
          gap: 7,
          height: CONTROL_H,
          padding: '0 14px',
          border: 'none',
          borderRadius: 8,
          cursor: 'pointer',
          font: "700 13px 'IBM Plex Sans'",
          color: 'var(--pio-run-fg)',
          background: 'linear-gradient(135deg,var(--pio-run-from),var(--pio-run-to))',
          boxShadow: '0 1px 2px rgba(15,20,27,0.15)',
          flex: 'none',
        }}
      >
        <svg width="12" height="12" viewBox="0 0 24 24" fill="currentColor">
          <path d="M7 4.5v15l13-7.5z" />
        </svg>
        {runLabel}
      </button>
    </div>
  )
}
