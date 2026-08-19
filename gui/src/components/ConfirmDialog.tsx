/** The confirm sheets from the design: job cancel/delete and results-overwrite.
 *  Same shell, two tones. */

const OVERLAY: React.CSSProperties = {
  position: 'absolute',
  inset: 0,
  background: 'rgba(15,20,27,0.45)',
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  zIndex: 31,
  padding: 40,
}

const LIGHT_BTN: React.CSSProperties = {
  padding: '9px 18px',
  borderRadius: 9,
  border: '1px solid #D7DBE0',
  background: '#fff',
  color: '#344054',
  font: "600 13px 'IBM Plex Sans'",
  cursor: 'pointer',
}

interface Props {
  title: string
  body: string
  /** Monospaced detail line under the body — the job title, or a path. */
  detail?: string
  dismissLabel: string
  confirmLabel: string
  tone: 'danger' | 'warning'
  pending?: boolean
  error?: string
  onDismiss: () => void
  onConfirm: () => void
}

export function ConfirmDialog({
  title,
  body,
  detail,
  dismissLabel,
  confirmLabel,
  tone,
  pending = false,
  error,
  onDismiss,
  onConfirm,
}: Props) {
  const warning = tone === 'warning'
  return (
    <div onClick={() => !pending && onDismiss()} style={OVERLAY}>
      <div
        onClick={(e) => e.stopPropagation()}
        style={{
          background: '#fff',
          borderRadius: 14,
          width: warning ? 448 : 420,
          maxWidth: '100%',
          overflow: 'hidden',
          boxShadow: '0 20px 60px rgba(0,0,0,0.3)',
        }}
      >
        <div style={{ padding: '22px 22px 18px', display: 'flex', gap: 14 }}>
          {warning && (
            <span
              style={{
                width: 40,
                height: 40,
                borderRadius: '50%',
                background: '#FEF3C7',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                flex: 'none',
              }}
            >
              <svg width="22" height="22" viewBox="0 0 24 24" fill="none">
                <path d="M12 3 2 20h20L12 3Z" stroke="#B45309" strokeWidth="1.8" strokeLinejoin="round" />
                <path d="M12 10v4M12 17v.5" stroke="#B45309" strokeWidth="2.2" strokeLinecap="round" />
              </svg>
            </span>
          )}
          <div style={{ minWidth: 0 }}>
            <div style={{ fontSize: 15.5, fontWeight: 700, color: '#1A2230' }}>{title}</div>
            <div style={{ fontSize: 13, color: '#667085', marginTop: 5, lineHeight: 1.5 }}>{body}</div>
            {detail && (
              <div
                style={
                  warning
                    ? {
                        fontSize: 12,
                        fontFamily: "'IBM Plex Mono'",
                        color: '#92400E',
                        background: '#FFFBEB',
                        border: '1px solid #FDE68A',
                        borderRadius: 8,
                        padding: '8px 10px',
                        marginTop: 11,
                        wordBreak: 'break-all',
                        whiteSpace: 'pre-line',
                      }
                    : {
                        fontSize: 12.5,
                        fontWeight: 600,
                        color: '#1B2A4A',
                        marginTop: 10,
                        fontFamily: "'IBM Plex Mono'",
                      }
                }
              >
                {detail}
              </div>
            )}
            {error && (
              <div
                role="alert"
                style={{
                  fontSize: 12,
                  color: '#B42318',
                  background: '#FEF3F2',
                  border: '1px solid #FECDCA',
                  borderRadius: 8,
                  padding: '8px 10px',
                  marginTop: 11,
                }}
              >
                {error}
              </div>
            )}
          </div>
        </div>
        <div
          style={{
            padding: '14px 22px',
            borderTop: '1px solid #EEF1F4',
            display: 'flex',
            justifyContent: 'flex-end',
            gap: 10,
            background: '#FAFBFC',
          }}
        >
          <button
            type="button"
            className="pio-light-btn"
            onClick={onDismiss}
            disabled={pending}
            style={{ ...LIGHT_BTN, cursor: pending ? 'not-allowed' : 'pointer', opacity: pending ? 0.6 : 1 }}
          >
            {dismissLabel}
          </button>
          <button
            type="button"
            className={warning ? 'pio-bright' : 'pio-danger'}
            onClick={onConfirm}
            disabled={pending}
            style={{
              padding: '9px 18px',
              borderRadius: 9,
              border: 'none',
              background: warning ? '#B45309' : '#DC2626',
              color: '#fff',
              font: "700 13px 'IBM Plex Sans'",
              cursor: pending ? 'wait' : 'pointer',
              opacity: pending ? 0.72 : 1,
            }}
          >
            {confirmLabel}
          </button>
        </div>
      </div>
    </div>
  )
}
