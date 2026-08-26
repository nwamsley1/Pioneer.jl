/** The ConvertRAW page.
 *
 *  Adapted rather than ported: the design shows a file / file-list / folder
 *  toggle and a JSON config, but neither converter takes a params file -- both
 *  take ONE positional path plus flags. So the file-list mode is dropped (the
 *  binaries cannot express it) and the panel shows the real command line
 *  instead of an invented JSON document.
 *
 *  The page drives two programs, not one: Thermo `.raw` goes through
 *  PioneerConverter (.NET), `.mzML` through Pioneer's own `convertMzML`
 *  (Julia). They are one page because they are one task -- get instrument data
 *  into the Arrow files SearchDIA reads -- but their tuning knobs are disjoint,
 *  so the advanced section swaps rather than greying fields out.
 */
import { NumField } from './NumField'
import { Toggle } from './Toggle'
import { BROWSE, LABEL, SEG_TRACK, seg } from '../lib/styles'
import type { ConvertParams } from '../lib/types'
import { type Note } from '../lib/validate'

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

const noteStyle = (note: Note): React.CSSProperties => ({
  marginTop: 9,
  fontSize: 12,
  lineHeight: 1.4,
  color: note.level === 'error' ? '#C0392B' : note.level === 'warn' ? '#B45309' : '#98A2B3',
})

interface Props {
  params: ConvertParams
  /** Logical cores, for the oversubscription warning. */
  inputNote: Note
  outputNote: Note
  advancedOpen: boolean
  onParam: (key: string, value: string) => void
  onToggle: (key: string) => void
  onBrowseInput: () => void
  onBrowseOutput: () => void
  onToggleAdvanced: () => void
}

export function ConvertRawForm({
  params,
  inputNote,
  outputNote,
  advancedOpen,
  onParam,
  onToggle,
  onBrowseInput,
  onBrowseOutput,
  onToggleAdvanced,
}: Props) {
  const isMzml = params.format === 'mzml'


  const defaultOut = params.input.trim()
    ? `${params.inputMode === 'file' ? params.input.trim().replace(/[^\\/]+$/, '').replace(/[\\/]$/, '') : params.input.trim()}/arrow_out`
    : '<input folder>/arrow_out'

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 14 }}>
      <section style={CARD}>
        <div style={{ display: 'flex', alignItems: 'baseline', gap: 12, marginBottom: 14 }}>
          <h2 style={H2}>{isMzml ? 'Convert mzML files' : 'Convert raw files'}</h2>
        </div>
        <p style={{ margin: '-4px 0 14px', fontSize: 12.5, color: '#667085', lineHeight: 1.5 }}>
          {isMzml ? (
            <>
              Convert <code style={{ fontFamily: "'IBM Plex Mono'", fontSize: 12 }}>.mzML</code>{' '}
              files to Arrow, which is what SearchDIA reads. Runs Pioneer&rsquo;s own converter, so
              it works on any platform. The data must be centroided — profile-mode spectra are
              skipped.
            </>
          ) : (
            <>
              Convert Thermo{' '}
              <code style={{ fontFamily: "'IBM Plex Mono'", fontSize: 12 }}>.raw</code> files to
              Arrow, which is what SearchDIA reads.
            </>
          )}
        </p>

        <label style={LABEL}>Format</label>
        <div style={{ ...SEG_TRACK, marginBottom: 14 }}>
          <button
            type="button"
            onClick={() => onParam('format', 'raw')}
            style={seg(!isMzml)}
          >
            Thermo .raw
          </button>
          <button
            type="button"
            onClick={() => onParam('format', 'mzml')}
            style={seg(isMzml)}
          >
            .mzML
          </button>
        </div>

        <label style={LABEL}>Input</label>
        <div style={{ ...SEG_TRACK, marginBottom: 12 }}>
          <button
            type="button"
            onClick={() => onParam('inputMode', 'file')}
            style={seg(params.inputMode === 'file')}
          >
            Single file
          </button>
          <button
            type="button"
            onClick={() => onParam('inputMode', 'folder')}
            style={seg(params.inputMode === 'folder')}
          >
            Folder
          </button>
        </div>

        <div data-drop="convertInput" style={{ display: 'flex', gap: 8 }}>
          <input
            className="pio-input"
            data-key="convertInput"
            value={params.input}
            onChange={(e) => onParam('input', e.target.value)}
            placeholder={
              params.inputMode === 'file'
                ? isMzml
                  ? '/path/to/file.mzML'
                  : '/path/to/file.raw'
                : isMzml
                  ? '/path/to/mzml/folder'
                  : '/path/to/raw/folder'
            }
            style={{
              flex: 1,
              padding: '9px 12px',
              borderRadius: 9,
              font: "13px 'IBM Plex Mono'",
              outline: 'none',
              minWidth: 0,
              border: `1px solid ${inputNote.level === 'error' ? '#E5484D' : '#D7DBE0'}`,
            }}
          />
          <button type="button" className="pio-browse" onClick={onBrowseInput} style={BROWSE}>
            Browse
          </button>
        </div>
        {inputNote.msg && (
          <div style={noteStyle(inputNote)}>
            {inputNote.level ? '⚠  ' : ''}
            {inputNote.msg}
          </div>
        )}
      </section>

      <section style={CARD}>
        <div style={{ display: 'flex', alignItems: 'baseline', gap: 12, marginBottom: 14 }}>
          <h2 style={H2}>Output</h2>
        </div>
        <label style={LABEL}>
          Output folder{' '}
          <span style={{ fontWeight: 400, color: '#98A2B3' }}>— defaults to the input folder</span>
        </label>
        <div data-drop="convertOutput" style={{ display: 'flex', gap: 8 }}>
          <input
            className="pio-input"
            data-key="convertOutput"
            value={params.outputDir}
            onChange={(e) => onParam('outputDir', e.target.value)}
            placeholder={defaultOut}
            style={{
              flex: 1,
              padding: '9px 12px',
              borderRadius: 9,
              font: "13px 'IBM Plex Mono'",
              outline: 'none',
              minWidth: 0,
              border: `1px solid ${outputNote.level === 'error' ? '#E5484D' : '#D7DBE0'}`,
            }}
          />
          <button type="button" className="pio-browse" onClick={onBrowseOutput} style={BROWSE}>
            Browse
          </button>
        </div>
        {outputNote.msg && <div style={noteStyle(outputNote)}>⚠&nbsp; {outputNote.msg}</div>}

        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            gap: 14,
            marginTop: 16,
          }}
        >
          <div>
            <div style={{ fontSize: 13, fontWeight: 600, color: '#344054' }}>Skip existing</div>
            <div style={{ fontSize: 11.5, color: '#98A2B3' }}>
              Leave files alone when their output already looks complete
            </div>
          </div>
          <Toggle
            on={params.skipExisting}
            fieldKey="skipExisting"
            onClick={() => onToggle('skipExisting')}
          />
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
            {/* Disjoint sets, not a shared set with some fields disabled: none
                of PioneerConverter's knobs has a counterpart in convertMzML,
                so showing them greyed out would only suggest otherwise. */}
            {isMzml ? (
              <>
                <div
                  style={{
                    display: 'grid',
                    gridTemplateColumns: '1fr 1fr',
                    gap: 14,
                    alignItems: 'start',
                    marginTop: 16,
                  }}
                >
                  <NumField
                    fieldKey="concurrentFiles"
                    value={params.concurrentFiles}
                    onChange={onParam}
                  />
                </div>
                <div
                  style={{
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'space-between',
                    gap: 14,
                    marginTop: 18,
                  }}
                >
                  <div>
                    <div style={{ fontSize: 13, fontWeight: 600, color: '#344054' }}>
                      Include scan headers
                    </div>
                    <div style={{ fontSize: 11.5, color: '#98A2B3' }}>
                      Carry each scan&rsquo;s header into the Arrow file. SearchDIA does not read
                      them and they roughly double the output — leave off unless something else
                      needs them.
                    </div>
                  </div>
                  <Toggle
                    on={params.includeScanHeader}
                    fieldKey="includeScanHeader"
                    onClick={() => onToggle('includeScanHeader')}
                  />
                </div>
              </>
            ) : (
              <div
                style={{
                  display: 'grid',
                  gridTemplateColumns: '1fr 1fr',
                  gap: 14,
                  alignItems: 'start',
                  marginTop: 16,
                }}
              >
                <NumField
                  fieldKey="threadsPerFile"
                  value={params.threadsPerFile}
                  onChange={onParam}
                />
                <NumField fieldKey="batchSize" value={params.batchSize} onChange={onParam} />
                <NumField
                  fieldKey="scanChunkSize"
                  value={params.scanChunkSize}
                  onChange={onParam}
                />
              </div>
            )}
          </div>
        )}
      </section>
    </div>
  )
}
