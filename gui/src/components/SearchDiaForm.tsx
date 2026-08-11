/** The SearchDIA page: Essentials, Confidence & output, and a collapsible
 *  Advanced section. Ported from the `isSearch` branch of the design template.
 *
 *  Note the design's stepper (currentStep / buildSteps / next / back) is not
 *  reproduced: the Console variant carries that state on its Component class
 *  but never renders it — everything is on one page.
 */
import { NumField } from './NumField'
import { Toggle } from './Toggle'
import type { SearchParams } from '../lib/types'
import type { Note } from '../lib/validate'

const CARD: React.CSSProperties = {
  background: '#fff',
  border: '1px solid #E7EAEE',
  borderRadius: 13,
  padding: '18px 20px',
}

const H2: React.CSSProperties = {
  margin: 0,
  fontSize: 12,
  fontWeight: 700,
  letterSpacing: '0.04em',
  textTransform: 'uppercase',
  color: '#1B2A4A',
}

const LABEL: React.CSSProperties = {
  display: 'block',
  fontSize: 12,
  fontWeight: 600,
  color: '#344054',
  marginBottom: 6,
}

const BROWSE: React.CSSProperties = {
  flex: 'none',
  padding: '0 14px',
  border: '1px solid #D7DBE0',
  borderRadius: 9,
  background: '#F8FAFB',
  font: "600 12.5px 'IBM Plex Sans'",
  color: '#344054',
  cursor: 'pointer',
}

const inputStyle = (note: Note): React.CSSProperties => ({
  flex: 1,
  padding: '9px 12px',
  borderRadius: 9,
  font: "13px 'IBM Plex Mono'",
  color: '#1A2230',
  background: '#fff',
  outline: 'none',
  minWidth: 0,
  border: `1px solid ${note.level === 'error' ? '#E5484D' : '#D7DBE0'}`,
})

function PathRow({
  label,
  fieldKey,
  value,
  placeholder,
  note,
  onChange,
  onBrowse,
  children,
}: {
  label: React.ReactNode
  fieldKey: string
  value: string
  placeholder: string
  note: Note
  onChange: (key: string, value: string) => void
  onBrowse: () => void
  children?: React.ReactNode
}) {
  return (
    <div>
      <label style={LABEL}>{label}</label>
      <div style={{ display: 'flex', gap: 8 }}>
        <input
          className="pio-input"
          data-key={fieldKey}
          value={value}
          placeholder={placeholder}
          onChange={(e) => onChange(fieldKey, e.target.value)}
          style={inputStyle(note)}
        />
        <button type="button" className="pio-browse" onClick={onBrowse} style={BROWSE}>
          Browse
        </button>
      </div>
      {note.msg && (
        <div
          style={{
            marginTop: 9,
            fontSize: 12,
            lineHeight: 1.4,
            color: note.level === 'error' ? '#C0392B' : '#B45309',
          }}
        >
          ⚠&nbsp; {note.msg}
        </div>
      )}
      {children}
    </div>
  )
}

function ToggleRow({
  title,
  hint,
  on,
  fieldKey,
  onToggle,
}: {
  title: string
  hint: string
  on: boolean
  fieldKey: string
  onToggle: (key: string) => void
}) {
  return (
    <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', gap: 14 }}>
      <div>
        <div style={{ fontSize: 13, fontWeight: 600, color: '#344054' }}>{title}</div>
        <div style={{ fontSize: 11.5, color: '#98A2B3' }}>{hint}</div>
      </div>
      <Toggle on={on} fieldKey={fieldKey} onClick={() => onToggle(fieldKey)} />
    </div>
  )
}

const LoadPreviousButton = ({ onClick }: { onClick: () => void }) => (
  <button
    type="button"
    className="pio-link"
    onClick={onClick}
    style={{
      background: 'none',
      border: 'none',
      padding: 0,
      cursor: 'pointer',
      font: "500 12px 'IBM Plex Sans'",
      color: '#8A93A0',
      display: 'inline-flex',
      alignItems: 'center',
      gap: 5,
    }}
  >
    <svg width="13" height="13" viewBox="0 0 24 24" fill="none">
      <path
        d="M12 15V4m0 0L8 8m4-4 4 4M5 16v3a1 1 0 0 0 1 1h12a1 1 0 0 0 1-1v-3"
        stroke="currentColor"
        strokeWidth="1.7"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
    </svg>
    Load previous run
  </button>
)

interface Props {
  params: SearchParams
  notes: { msData: Note; library: Note; results: Note }
  advancedOpen: boolean
  onParam: (key: string, value: string) => void
  onToggle: (key: string) => void
  onBrowse: (key: 'msData' | 'library' | 'results') => void
  onToggleAdvanced: () => void
  onOpenLoad: () => void
  onGoToBuild: () => void
}

export function SearchDiaForm({
  params,
  notes,
  advancedOpen,
  onParam,
  onToggle,
  onBrowse,
  onToggleAdvanced,
  onOpenLoad,
  onGoToBuild,
}: Props) {
  const seg = (active: boolean): React.CSSProperties => ({
    flex: 1,
    padding: '7px 0',
    border: 'none',
    borderRadius: 6,
    font: "600 12.5px 'IBM Plex Sans'",
    cursor: 'pointer',
    transition: 'all .12s',
    ...(active
      ? { background: '#fff', color: '#1B2A4A', boxShadow: '0 1px 3px rgba(15,20,27,0.12)' }
      : { background: 'none', color: '#667085' }),
  })

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 14 }}>
      <section style={CARD}>
        <div
          style={{
            display: 'flex',
            alignItems: 'baseline',
            justifyContent: 'space-between',
            gap: 12,
            marginBottom: 14,
          }}
        >
          <h2 style={H2}>Essentials</h2>
          <LoadPreviousButton onClick={onOpenLoad} />
        </div>
        <div style={{ display: 'flex', flexDirection: 'column', gap: 12 }}>
          <PathRow
            label="MS data folder"
            fieldKey="msData"
            value={params.msData}
            placeholder="/path/to/ms/data"
            note={notes.msData}
            onChange={onParam}
            onBrowse={() => onBrowse('msData')}
          />
          <PathRow
            label="Spectral library"
            fieldKey="library"
            value={params.library}
            placeholder="/path/to/library.poin"
            note={notes.library}
            onChange={onParam}
            onBrowse={() => onBrowse('library')}
          >
            <button
              type="button"
              className="pio-link-underline"
              onClick={onGoToBuild}
              style={{
                marginTop: 8,
                background: 'none',
                border: 'none',
                padding: 0,
                cursor: 'pointer',
                font: "500 12px 'IBM Plex Sans'",
                color: 'var(--pio-accent)',
                display: 'inline-flex',
                alignItems: 'center',
                gap: 5,
              }}
            >
              <svg width="13" height="13" viewBox="0 0 24 24" fill="none">
                <path d="M12 5v14M5 12h14" stroke="currentColor" strokeWidth="2" strokeLinecap="round" />
              </svg>
              No library yet? Build one from a FASTA →
            </button>
          </PathRow>
          <PathRow
            label="Results folder"
            fieldKey="results"
            value={params.results}
            placeholder="/path/to/results"
            note={notes.results}
            onChange={onParam}
            onBrowse={() => onBrowse('results')}
          />
        </div>
      </section>

      <section style={CARD}>
        <div
          style={{
            display: 'flex',
            alignItems: 'baseline',
            justifyContent: 'space-between',
            gap: 12,
            marginBottom: 14,
          }}
        >
          <h2 style={H2}>Confidence &amp; output</h2>
        </div>
        <div style={{ display: 'flex', gap: 22, alignItems: 'flex-start', flexWrap: 'wrap' }}>
          <div style={{ width: 190, flex: 'none' }}>
            <NumField fieldKey="qValue" value={params.qValue} onChange={onParam} />
            <div style={{ fontSize: 11.5, color: '#98A2B3', marginTop: 6 }}>1% FDR = 0.01</div>
          </div>
          <div
            style={{
              flex: 1,
              minWidth: 240,
              display: 'flex',
              flexDirection: 'column',
              gap: 11,
            }}
          >
            <ToggleRow
              title="Write CSV report"
              hint="Human-readable tables"
              on={params.writeCsv}
              fieldKey="writeCsv"
              onToggle={onToggle}
            />
            <ToggleRow
              title="Write decoys"
              hint="Include decoy matches"
              on={params.writeDecoys}
              fieldKey="writeDecoys"
              onToggle={onToggle}
            />
            <ToggleRow
              title="Delete temp files"
              hint="Clean up intermediates"
              on={params.deleteTemp}
              fieldKey="deleteTemp"
              onToggle={onToggle}
            />
            <ToggleRow
              title="Run-to-run normalization"
              hint="Retention-time dependent cross-run scaling"
              on={params.runToRunNorm}
              fieldKey="runToRunNorm"
              onToggle={onToggle}
            />
          </div>
        </div>
      </section>

      <section
        style={{
          background: '#fff',
          border: '1px solid #E7EAEE',
          borderRadius: 13,
          overflow: 'hidden',
        }}
      >
        <button
          type="button"
          className="pio-adv"
          onClick={onToggleAdvanced}
          style={{
            width: '100%',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            padding: '15px 20px',
            background: 'none',
            border: 'none',
            cursor: 'pointer',
            font: 'inherit',
            textAlign: 'left',
          }}
        >
          <div>
            <div
              style={{
                fontSize: 12,
                fontWeight: 700,
                letterSpacing: '0.04em',
                textTransform: 'uppercase',
                color: '#1B2A4A',
              }}
            >
              Advanced parameters
            </div>
            <div style={{ fontSize: 11.5, color: '#98A2B3', marginTop: 2 }}>
              Defaults suit most experiments
            </div>
          </div>
          <span
            style={{
              display: 'flex',
              transition: 'transform .18s',
              transform: advancedOpen ? 'rotate(180deg)' : undefined,
            }}
          >
            <svg width="18" height="18" viewBox="0 0 24 24" fill="none">
              <path
                d="m6 9 6 6 6-6"
                stroke="#667085"
                strokeWidth="2"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </span>
        </button>
        {advancedOpen && (
          <div style={{ padding: '4px 20px 20px', borderTop: '1px solid #EEF1F4' }}>
            <div
              style={{
                display: 'grid',
                gridTemplateColumns: '1fr 1fr 1fr',
                gap: 14,
                alignItems: 'start',
                marginTop: 16,
              }}
            >
              <NumField fieldKey="nIsotopes" value={params.nIsotopes} onChange={onParam} />
              <NumField fieldKey="nce" value={params.nce} onChange={onParam} />
              <NumField fieldKey="minPeptides" value={params.minPeptides} onChange={onParam} />
            </div>
            <div style={{ marginTop: 14, maxWidth: 280 }}>
              <label style={{ display: 'block', fontSize: 12, color: '#475467', marginBottom: 6 }}>
                Trace mode
              </label>
              <div style={{ display: 'flex', padding: 3, background: '#EEF1F4', borderRadius: 8 }}>
                <button
                  type="button"
                  onClick={() => onParam('traceMode', 'combined')}
                  style={seg(params.traceMode === 'combined')}
                >
                  combined
                </button>
                <button
                  type="button"
                  onClick={() => onParam('traceMode', 'separated')}
                  style={seg(params.traceMode === 'separated')}
                >
                  separated
                </button>
              </div>
            </div>
          </div>
        )}
      </section>
    </div>
  )
}
