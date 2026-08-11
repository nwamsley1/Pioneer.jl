/** The params-JSON modal: view the generated config, or edit it by hand and
 *  fold the result back into the form. Ported from the JSON MODAL block.
 *
 *  "Apply" maps the edited JSON onto the form fields and stashes any keys the
 *  form does not model, so hand-written settings survive being edited here.
 */
import { useState } from 'react'

interface Props {
  fileName: string
  /** ConvertRAW has no params file to parse back, so its preview is read-only. */
  editable?: boolean
  /** The currently generated JSON, pretty-printed. */
  text: string
  /** True when keys outside the form's model are being carried through. */
  hasExtras: boolean
  extraKeys: string[]
  onClose: () => void
  /** Returns an error message, or null on success. */
  onApply: (draft: string) => string | null
  onRevert: () => void
}

export function JsonModal({
  fileName,
  editable = true,
  text,
  hasExtras,
  extraKeys,
  onClose,
  onApply,
  onRevert,
}: Props) {
  const [editing, setEditing] = useState(false)
  const [draft, setDraft] = useState(text)
  const [error, setError] = useState('')

  const startEdit = () => {
    setDraft(text)
    setError('')
    setEditing(true)
  }

  const onDraft = (v: string) => {
    setDraft(v)
    try {
      JSON.parse(v)
      setError('')
    } catch (e) {
      setError((e as Error).message)
    }
  }

  const apply = () => {
    const err = onApply(draft)
    if (err) {
      setError(err)
      return
    }
    setEditing(false)
    setError('')
  }

  const ghost: React.CSSProperties = {
    padding: '9px 16px',
    borderRadius: 9,
    border: '1px solid rgba(255,255,255,0.14)',
    background: 'none',
    color: '#C8D0DA',
    font: "600 13px 'IBM Plex Sans'",
    cursor: 'pointer',
  }

  return (
    <div
      onClick={onClose}
      style={{
        position: 'absolute',
        inset: 0,
        background: 'rgba(15,20,27,0.45)',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        zIndex: 20,
        padding: 40,
      }}
    >
      <div
        onClick={(e) => e.stopPropagation()}
        style={{
          background: 'var(--pio-nav)',
          borderRadius: 14,
          width: 600,
          maxWidth: '100%',
          maxHeight: '100%',
          display: 'flex',
          flexDirection: 'column',
          overflow: 'hidden',
          boxShadow: '0 20px 60px rgba(0,0,0,0.4)',
        }}
      >
        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            padding: '16px 20px',
            borderBottom: '1px solid rgba(255,255,255,0.08)',
          }}
        >
          <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
            <span
              style={{
                fontSize: 13,
                fontWeight: 600,
                color: '#E7EBF0',
                fontFamily: "'IBM Plex Mono'",
              }}
            >
              {fileName}
            </span>
            <span
              title={hasExtras ? extraKeys.join(', ') : undefined}
              style={{
                fontSize: 10.5,
                fontWeight: 600,
                letterSpacing: '0.04em',
                textTransform: 'uppercase',
                padding: '2px 8px',
                borderRadius: 20,
                background: hasExtras ? 'var(--pio-accent-wash)' : 'rgba(255,255,255,0.07)',
                color: hasExtras ? '#CBD8EE' : '#98A6BC',
              }}
            >
              {hasExtras ? 'custom keys' : 'generated'}
            </span>
          </div>
          <button
            type="button"
            onClick={onClose}
            title="Close"
            style={{
              background: 'none',
              border: 'none',
              color: '#98A6BC',
              cursor: 'pointer',
              padding: 4,
              borderRadius: 6,
              display: 'flex',
            }}
          >
            <svg width="17" height="17" viewBox="0 0 24 24" fill="none">
              <path d="M6 6l12 12M18 6 6 18" stroke="currentColor" strokeWidth="2" strokeLinecap="round" />
            </svg>
          </button>
        </div>

        {!editing && (
          <pre
            style={{
              margin: 0,
              overflow: 'auto',
              padding: '18px 20px',
              font: "12.5px/1.6 'IBM Plex Mono'",
              color: '#C8D0DA',
              maxHeight: '58vh',
            }}
          >
            {text}
          </pre>
        )}
        {editing && (
          <div style={{ padding: '14px 20px 4px' }}>
            <textarea
              className="pio-textarea"
              value={draft}
              onChange={(e) => onDraft(e.target.value)}
              spellCheck={false}
              style={{
                width: '100%',
                height: '46vh',
                resize: 'vertical',
                boxSizing: 'border-box',
                background: '#080B0F',
                color: '#C8D0DA',
                border: '1px solid rgba(255,255,255,0.12)',
                borderRadius: 9,
                padding: 14,
                font: "12.5px/1.6 'IBM Plex Mono'",
                outline: 'none',
              }}
            />
            <div
              style={{
                minHeight: 20,
                padding: '7px 2px 2px',
                font: "12px 'IBM Plex Mono'",
                color: error ? '#F08C8C' : 'transparent',
              }}
            >
              {error ? `⚠  ${error}` : '.'}
            </div>
          </div>
        )}

        <div
          style={{
            padding: '14px 20px',
            borderTop: '1px solid rgba(255,255,255,0.08)',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            gap: 10,
          }}
        >
          {!editing ? (
            <>
              {editable ? (
              <button type="button" className="pio-ghost" onClick={startEdit} style={{ ...ghost, display: 'flex', alignItems: 'center', gap: 7 }}>
                <svg width="14" height="14" viewBox="0 0 24 24" fill="none">
                  <path
                    d="m13 5 6 6M4 20l1-4L16 5l3 3L8 19l-4 1Z"
                    stroke="currentColor"
                    strokeWidth="1.7"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  />
                </svg>
                Edit
              </button>
              ) : (
                <span style={{ fontSize: 12, color: '#98A6BC' }}>
                  Run this yourself to reproduce the conversion.
                </span>
              )}
              <button
                type="button"
                className="pio-ghost-strong"
                onClick={onClose}
                style={{
                  padding: '9px 18px',
                  borderRadius: 9,
                  border: 'none',
                  background: 'rgba(255,255,255,0.12)',
                  color: '#E7EBF0',
                  font: "600 13px 'IBM Plex Sans'",
                  cursor: 'pointer',
                }}
              >
                Close
              </button>
            </>
          ) : (
            <>
              <button
                type="button"
                className="pio-ghost"
                onClick={() => {
                  onRevert()
                  setEditing(false)
                  setError('')
                }}
                style={{ ...ghost, color: '#AEB9C8' }}
              >
                Revert to parameters
              </button>
              <div style={{ display: 'flex', gap: 10 }}>
                <button
                  type="button"
                  className="pio-ghost"
                  onClick={() => {
                    setEditing(false)
                    setError('')
                  }}
                  style={ghost}
                >
                  Cancel
                </button>
                <button
                  type="button"
                  onClick={apply}
                  disabled={!!error}
                  style={{
                    padding: '9px 18px',
                    borderRadius: 9,
                    border: 'none',
                    font: "600 13px 'IBM Plex Sans'",
                    ...(error
                      ? {
                          background: 'var(--pio-accent-wash-strong)',
                          color: 'rgba(255,255,255,0.5)',
                          cursor: 'not-allowed',
                        }
                      : { background: 'var(--pio-accent)', color: '#FFFFFF', cursor: 'pointer' }),
                  }}
                >
                  Apply
                </button>
              </div>
            </>
          )}
        </div>
      </div>
    </div>
  )
}
