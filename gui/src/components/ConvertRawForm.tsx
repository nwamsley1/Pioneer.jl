/** The ConvertRAW page.
 *
 *  Adapted rather than ported: the design shows a file / file-list / folder
 *  toggle and a JSON config, but PioneerConverter is a .NET program that takes
 *  ONE positional RAW path plus flags and has no params file at all. So the
 *  file-list mode is dropped (the binary cannot express it) and the panel shows
 *  the real command line instead of an invented JSON document.
 */
import { NumField } from './NumField'
import { Toggle } from './Toggle'
import type { ConvertParams } from '../lib/types'
import { convertTotalThreads, type Note } from '../lib/validate'

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

const noteStyle = (note: Note): React.CSSProperties => ({
  marginTop: 9,
  fontSize: 12,
  lineHeight: 1.4,
  color: note.level === 'error' ? '#C0392B' : note.level === 'warn' ? '#B45309' : '#98A2B3',
})

interface Props {
  params: ConvertParams
  /** Logical cores, for the oversubscription warning. */
  maxThreads: number
  inputNote: Note
  outputNote: Note
  threadsNote: Note
  advancedOpen: boolean
  onParam: (key: string, value: string) => void
  onToggle: (key: string) => void
  onBrowseInput: () => void
  onBrowseOutput: () => void
  onToggleAdvanced: () => void
}

export function ConvertRawForm({
  params,
  maxThreads,
  inputNote,
  outputNote,
  threadsNote,
  advancedOpen,
  onParam,
  onToggle,
  onBrowseInput,
  onBrowseOutput,
  onToggleAdvanced,
}: Props) {
  const seg = (active: boolean): React.CSSProperties => ({
    padding: '7px 14px',
    border: 'none',
    borderRadius: 8,
    font: "600 12.5px 'IBM Plex Sans'",
    cursor: 'pointer',
    transition: 'all .12s',
    ...(active
      ? { background: '#fff', color: '#1B2A4A', boxShadow: '0 1px 3px rgba(15,20,27,0.12)' }
      : { background: 'none', color: '#667085' }),
  })

  // The two knobs multiply: N files in flight, each with M scan-reader threads.
  // Shared with the validator so the readout and the block can never disagree.
  const total = convertTotalThreads(params)

  const defaultOut = params.input.trim()
    ? `${params.inputMode === 'file' ? params.input.trim().replace(/[^\\/]+$/, '').replace(/[\\/]$/, '') : params.input.trim()}/arrow_out`
    : '<input folder>/arrow_out'

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 14 }}>
      <section style={CARD}>
        <div style={{ display: 'flex', alignItems: 'baseline', gap: 12, marginBottom: 14 }}>
          <h2 style={H2}>Convert raw files</h2>
        </div>
        <p style={{ margin: '-4px 0 14px', fontSize: 12.5, color: '#667085', lineHeight: 1.5 }}>
          Convert Thermo{' '}
          <code style={{ fontFamily: "'IBM Plex Mono'", fontSize: 12 }}>.raw</code> files to Arrow,
          which is what SearchDIA reads.
        </p>

        <label style={LABEL}>Input</label>
        <div
          style={{
            display: 'inline-flex',
            padding: 3,
            background: '#EEF1F4',
            borderRadius: 10,
            marginBottom: 12,
          }}
        >
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

        <div style={{ display: 'flex', gap: 8 }}>
          <input
            className="pio-input"
            data-key="convertInput"
            value={params.input}
            onChange={(e) => onParam('input', e.target.value)}
            placeholder={
              params.inputMode === 'file' ? '/path/to/file.raw' : '/path/to/raw/folder'
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
        <div style={{ display: 'flex', gap: 8 }}>
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

      <section style={CARD}>
        <div style={{ display: 'flex', alignItems: 'baseline', gap: 12, marginBottom: 4 }}>
          <h2 style={H2}>Parallelism</h2>
        </div>
        <p style={{ margin: '0 0 14px', fontSize: 12, color: '#98A2B3', lineHeight: 1.5 }}>
          The converter works on several files at once, and splits each file across its own scan
          readers — so these multiply. The Julia thread setting does not apply here; this is a .NET
          program.
        </p>
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 14 }}>
          <NumField fieldKey="concurrentFiles" value={params.concurrentFiles} onChange={onParam} />
          <NumField fieldKey="threadsPerFile" value={params.threadsPerFile} onChange={onParam} />
        </div>
        <div style={{ marginTop: 10, fontSize: 12, color: '#475467' }}>
          ≈ <strong style={{ fontFamily: "'IBM Plex Mono'" }}>{total}</strong> threads in total ·{' '}
          {maxThreads} available
          {params.inputMode === 'file' && parseInt(params.concurrentFiles, 10) > 1 && (
            <span style={{ color: '#B45309' }}>
              {' '}
              · converting a single file, so only threads-per-file matters
            </span>
          )}
        </div>
        {threadsNote.msg && (
          <div style={{ marginTop: 9, fontSize: 12, lineHeight: 1.4, color: '#C0392B' }}>
            ⚠&nbsp; {threadsNote.msg}
          </div>
        )}
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
                gridTemplateColumns: '1fr 1fr',
                gap: 14,
                alignItems: 'start',
                marginTop: 16,
              }}
            >
              <NumField fieldKey="batchSize" value={params.batchSize} onChange={onParam} />
              <NumField fieldKey="scanChunkSize" value={params.scanChunkSize} onChange={onParam} />
            </div>
          </div>
        )}
      </section>
    </div>
  )
}
