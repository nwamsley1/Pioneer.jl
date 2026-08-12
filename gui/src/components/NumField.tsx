/** Port of NumField.dc.html — a labeled numeric input with clamped
 *  up/down steppers and an inline error beside the label. */
import { LABEL_TIGHT } from '../lib/styles'
import { numError, NUM_SPECS } from '../lib/validate'

interface Props {
  fieldKey: string
  value: string
  onChange: (key: string, value: string) => void
}

/** A hoverable ⓘ beside a label, for a field whose name cannot carry the whole
 *  explanation. `title` rather than a custom popover: the native tooltip is
 *  keyboard- and screen-reader-reachable, and this is supplementary detail
 *  rather than something needed to use the field. */
function InfoDot({ text }: { text: string }) {
  return (
    <span
      title={text}
      aria-label={text}
      tabIndex={0}
      style={{
        display: 'inline-flex',
        alignItems: 'center',
        justifyContent: 'center',
        width: 14,
        height: 14,
        flex: 'none',
        borderRadius: '50%',
        border: '1px solid #B6BFC9',
        color: '#8A93A0',
        font: "italic 700 9.5px 'IBM Plex Serif', Georgia, serif",
        cursor: 'help',
        userSelect: 'none',
      }}
    >
      i
    </span>
  )
}

export function NumField({ fieldKey, value, onChange }: Props) {
  const spec = NUM_SPECS[fieldKey]
  const error = numError(fieldKey, value)

  const clamp = (v: number): number => {
    if (spec.min != null && v < spec.min) return spec.min
    if (spec.max != null && v > spec.max) return spec.max
    return v
  }

  const fmt = (v: number): string =>
    spec.int ? String(Math.round(v)) : String(parseFloat(v.toFixed(6)))

  const stepBy = (dir: number) => {
    let v = Number(value)
    if (!Number.isFinite(v)) v = spec.min ?? 0
    onChange(fieldKey, fmt(clamp(v + dir * spec.step)))
  }

  // The two buttons divide the column rather than asserting a height. A fixed
  // height (18 each + the 1px divider = 37) fell short of the input's ~41px
  // row, so the column's unpainted remainder showed as a white sliver under
  // the down arrow and pushed the divider above the field's centre. Splitting
  // by flex tracks the input's line-height and padding automatically.
  const stepBtn: React.CSSProperties = {
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    flex: 1,
    minHeight: 0,
    width: 26,
    border: 'none',
    background: '#F5F7F9',
    cursor: 'pointer',
    padding: 0,
  }

  return (
    <div style={{ display: 'flex', flexDirection: 'column' }}>
      <div
        style={{
          display: 'flex',
          alignItems: 'center',
          gap: 8,
          marginBottom: 6,
          minHeight: 16,
        }}
      >
        <label style={LABEL_TIGHT}>{spec.label}</label>
        {spec.info && <InfoDot text={spec.info} />}
        {error && (
          <span
            style={{
              fontSize: 11,
              color: '#C0392B',
              fontWeight: 600,
              whiteSpace: 'nowrap',
            }}
          >
            ⚠ {error}
          </span>
        )}
      </div>
      <div
        style={{
          display: 'flex',
          alignItems: 'stretch',
          border: `1px solid ${error ? '#E5484D' : '#D7DBE0'}`,
          borderRadius: 8,
          overflow: 'hidden',
          background: '#fff',
          // Sized to what a number needs, not to the column it sits in. These
          // hold four or five characters; stretched to a full column they read
          // as unfinished next to the path fields, which are full-width because
          // a path genuinely can be any length.
          maxWidth: 150,
        }}
      >
        <input
          data-key={fieldKey}
          value={value}
          onChange={(e) => onChange(fieldKey, e.target.value)}
          inputMode="decimal"
          style={{
            flex: 1,
            minWidth: 0,
            width: '100%',
            padding: '9px 11px',
            border: 'none',
            background: 'transparent',
            font: "13px 'IBM Plex Mono'",
            outline: 'none',
            color: '#1A2230',
          }}
        />
        <div
          style={{
            display: 'flex',
            flexDirection: 'column',
            borderLeft: '1px solid #E2E6EA',
            flex: 'none',
            // Matches the buttons, so a sub-pixel rounding remainder cannot
            // show through as white.
            background: '#F5F7F9',
          }}
        >
          <button
            type="button"
            className="pio-numstep"
            onClick={() => stepBy(1)}
            tabIndex={-1}
            style={stepBtn}
          >
            <svg width="9" height="9" viewBox="0 0 24 24" fill="none">
              <path
                d="M6 15l6-6 6 6"
                stroke="#5E6877"
                strokeWidth="2.6"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </button>
          <button
            type="button"
            className="pio-numstep"
            onClick={() => stepBy(-1)}
            tabIndex={-1}
            style={{ ...stepBtn, borderTop: '1px solid #E2E6EA' }}
          >
            <svg width="9" height="9" viewBox="0 0 24 24" fill="none">
              <path
                d="M6 9l6 6 6-6"
                stroke="#5E6877"
                strokeWidth="2.6"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </button>
        </div>
      </div>
    </div>
  )
}
