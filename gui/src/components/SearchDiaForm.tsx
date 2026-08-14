/** The SearchDIA page: Essentials, Confidence & output, and a collapsible
 *  Advanced section. Ported from the `isSearch` branch of the design template.
 *
 *  Note the design's stepper (currentStep / buildSteps / next / back) is not
 *  reproduced: the Console variant carries that state on its Component class
 *  but never renders it — everything is on one page.
 */
import { LibrarySummary } from './LibrarySummary'
import { NumField } from './NumField'
import { Toggle } from './Toggle'
import { BROWSE, HINT, LABEL, LABEL_TIGHT } from '../lib/styles'
import type { LibraryInfo } from '../lib/backend'
import { isPrositModel } from '../lib/types'
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
    // The whole row is the drop zone, not just the input: aiming at "the MS
    // data field" can plausibly land on the label, the Browse button or the
    // padding between them.
    <div data-drop={fieldKey}>
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
        <div style={LABEL_TIGHT}>{title}</div>
        <div style={HINT}>{hint}</div>
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
  /** What the selected .poin records about itself; null when none is chosen. */
  libInfo: LibraryInfo | null
  onParam: (key: string, value: string) => void
  onToggle: (key: string) => void
  onBrowse: (key: 'msData' | 'library' | 'results') => void
  onOpenLoad: () => void
  onGoToBuild: () => void
}

export function SearchDiaForm({
  params,
  notes,
  libInfo,
  onParam,
  onToggle,
  onBrowse,
  onOpenLoad,
  onGoToBuild,
}: Props) {
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
            <LibrarySummary info={libInfo} />
            {libInfo &&
              libInfo.is_library &&
              isPrositModel(libInfo.prediction_model) &&
              libInfo.variable_mods.length > 0 && (
                <div
                  style={{
                    marginTop: 10,
                    padding: '10px 12px',
                    borderRadius: 9,
                    background: '#FFFBEB',
                    border: '1px solid #FDE68A',
                    fontSize: 12,
                    lineHeight: 1.5,
                    color: '#92400E',
                  }}
                >
                  <strong style={{ fontWeight: 600 }}>Experimental.</strong> This library
                  is Prosit-predicted and carries variable modifications
                  ({libInfo.variable_mods.join(', ')}). Pioneer does not report
                  site-localization confidence, so a modified residue is placed but the
                  placement is not scored.
                </div>
              )}
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

      {/* Confidence & output and Advanced sit side by side at half width each.
          Basis 300px, not 340: at the window's 980px minimum the sidebar and
          padding leave about 670px, so 340 each would exceed it and wrap to a
          column at the default size. At 300 they fit and then grow to ~330.
          flex-start so the shorter card does not stretch to match the taller. */}
      <div style={{ display: 'flex', gap: 14, alignItems: 'flex-start', flexWrap: 'wrap' }}>
        <section style={{ ...CARD, flex: '1 1 300px', minWidth: 0 }}>
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
        {/* Stacked, not side by side. One field beside four toggles left a large
            empty quadrant under the field, and no amount of alignment fixes a
            column with nothing else to hold. Full width also lets each toggle
            sit at the card edge, which reads as deliberate. */}
        <div style={{ display: 'flex', flexDirection: 'column', gap: 16 }}>
          <div>
            <NumField fieldKey="qValue" value={params.qValue} onChange={onParam} />
            <div style={{ fontSize: 11.5, color: '#98A2B3', marginTop: 6 }}>1% FDR = 0.01</div>
          </div>
          <div
            style={{
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
            <ToggleRow
              title="Debug logging"
              hint="Verbose console output — the log file keeps its own detail either way"
              on={params.debugLogging}
              fieldKey="debugLogging"
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
          flex: '1 1 300px',
          minWidth: 0,
        }}
      >
        {/* Always open. Collapsed, this card was a half-width empty box beside
            Confidence & output, and the four settings inside it are short
            enough that hiding them saved nothing. */}
        <div style={{ padding: '15px 20px 0' }}>
          <h2 style={H2}>Advanced parameters</h2>
          <div style={{ fontSize: 11.5, color: '#98A2B3', marginTop: 2 }}>
            Defaults suit most experiments
          </div>
        </div>
        <div style={{ padding: '4px 20px 20px' }}>
            <div
              style={{
                display: 'grid',
                gridTemplateColumns: '1fr 1fr 1fr',
                gap: 14,
                alignItems: 'stretch',
                marginTop: 16,
              }}
            >
              <NumField fieldKey="nIsotopes" value={params.nIsotopes} onChange={onParam} />
              <NumField fieldKey="nce" value={params.nce} onChange={onParam} />
              <NumField fieldKey="minPeptides" value={params.minPeptides} onChange={onParam} />
            </div>
        </div>
        </section>
      </div>
    </div>
  )
}
