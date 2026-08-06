/** The block beneath the form body: threads card, "Edit the .json directly",
 *  and the primary Run button. Ported from lines 411-429 of the design, shared
 *  by SearchDIA and BuildSpecLib.
 */
interface Props {
  /** The Julia thread control is meaningless for the .NET converter. */
  showThreads?: boolean
  /** ConvertRAW has no params file — its preview shows the command line. */
  previewLabel?: string
  threads: number
  maxThreads: number
  runLabel: string
  onThreads: (n: number) => void
  onEditJson: () => void
  onRun: () => void
}

export function RunBar({
  showThreads = true,
  previewLabel = 'Edit the .json directly',
  threads,
  maxThreads,
  runLabel,
  onThreads,
  onEditJson,
  onRun,
}: Props) {
  const stepBtn: React.CSSProperties = {
    padding: '8px 13px',
    border: 'none',
    background: '#F8FAFB',
    cursor: 'pointer',
    font: "700 15px 'IBM Plex Sans'",
    color: '#344054',
  }

  return (
    <div style={{ marginTop: 18 }}>
      {showThreads && (
      <div
        style={{
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          gap: 12,
          padding: '12px 16px',
          border: '1px solid #E7EAEE',
          borderRadius: 11,
          background: '#fff',
          marginBottom: 12,
        }}
      >
        <div>
          <div style={{ fontSize: 13.5, fontWeight: 600, color: '#1A2230' }}>Threads</div>
          <div style={{ fontSize: 12, color: '#98A2B3', marginTop: 1 }}>
            {maxThreads} available · max {maxThreads}
          </div>
        </div>
        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            border: '1px solid #D7DBE0',
            borderRadius: 9,
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
            style={{
              width: 46,
              textAlign: 'center',
              border: 'none',
              borderLeft: '1px solid #E5E9ED',
              borderRight: '1px solid #E5E9ED',
              padding: '8px 0',
              font: "600 14px 'IBM Plex Mono'",
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
        className="pio-link"
        onClick={onEditJson}
        style={{
          margin: '0 auto 12px',
          display: 'flex',
          alignItems: 'center',
          gap: 7,
          background: 'none',
          border: 'none',
          cursor: 'pointer',
          font: "600 12.5px 'IBM Plex Sans'",
          color: '#667085',
        }}
      >
        <svg width="14" height="14" viewBox="0 0 24 24" fill="none">
          <path
            d="m13 5 6 6M4 20l1-4L16 5l3 3L8 19l-4 1Z"
            stroke="currentColor"
            strokeWidth="1.7"
            strokeLinecap="round"
            strokeLinejoin="round"
          />
        </svg>
        {previewLabel}
      </button>

      <button
        type="button"
        className="pio-run"
        onClick={onRun}
        style={{
          width: '100%',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          gap: 10,
          padding: 15,
          border: 'none',
          borderRadius: 12,
          cursor: 'pointer',
          font: "700 15px 'IBM Plex Sans'",
          color: '#FFFFFF',
          background: 'linear-gradient(135deg,#3E62A0,#2E4D7E)',
          boxShadow: '0 6px 18px rgba(46,77,126,0.3)',
        }}
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="#FFFFFF">
          <path d="M7 4.5v15l13-7.5z" />
        </svg>
        {runLabel}
      </button>
    </div>
  )
}
