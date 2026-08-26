/** The SearchDIA page: Files and name, then Common parameters beside Output,
 *  then Advanced parameters across the full width.
 *
 *  The first card holds what is particular to this run — the two inputs, the
 *  destination, and what to call it. The three below it hold settings that
 *  mostly stay put between runs, split by how often you touch them.
 *
 *  Note the design's stepper (currentStep / buildSteps / next / back) is not
 *  reproduced: the Console variant carries that state on its Component class
 *  but never renders it — everything is on one page.
 */
import { JobNameField } from './JobNameField'
import { LibrarySummary } from './LibrarySummary'
import { NumField } from './NumField'
import { RecentLibraries } from './RecentLibraries'
import { Toggle } from './Toggle'
import { unimodDisplay } from '../lib/fasta'
import { BROWSE, HINT, LABEL, LABEL_TIGHT, SEG_TRACK, seg } from '../lib/styles'
import type { LibraryInfo } from '../lib/backend'
import { isPrositModel } from '../lib/types'
import type { SearchParams } from '../lib/types'
import { fileStem } from '../lib/validate'
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
  action,
  children,
}: {
  label: React.ReactNode
  fieldKey: string
  value: string
  placeholder: string
  note: Note
  onChange: (key: string, value: string) => void
  onBrowse: () => void
  /** An extra control beside Browse, e.g. a recent-paths picker. */
  action?: React.ReactNode
  children?: React.ReactNode
}) {
  return (
    // The whole row is the drop zone, not just the input: aiming at "the MS
    // data field" can plausibly land on the label, the Browse button or the
    // padding between them.
    <div data-drop={fieldKey}>
      {label ? <label style={LABEL}>{label}</label> : null}
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
        {action}
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

/** The chosen-files list.
 *
 *  Each row is one search: the file on the left, the results folder it will
 *  write to on the right. Showing the destination per row rather than once
 *  above the list is the point of the mode -- what makes these separate runs is
 *  that they do not share an output tree, and a list that only showed file
 *  names would look like the folder mode with extra steps.
 */
function MsFileList({
  params,
  onAdd,
  onRemove,
}: {
  params: SearchParams
  onAdd: () => void
  onRemove: (index: number) => void
}) {
  const files = params.msDataFiles
  const root = params.results.trim()
  return (
    <div>
      {files.length > 0 && (
        <div
          style={{
            border: '1px solid #E3E8EC',
            borderRadius: 9,
            overflow: 'hidden',
            marginBottom: 9,
          }}
        >
          {files.map((f, i) => (
            <div
              key={`${f}:${i}`}
              style={{
                display: 'flex',
                alignItems: 'center',
                gap: 10,
                padding: '8px 10px',
                borderTop: i ? '1px solid #EEF1F4' : undefined,
              }}
            >
              <div style={{ flex: 1, minWidth: 0 }}>
                <div
                  title={f}
                  style={{
                    font: "12.5px 'IBM Plex Mono'",
                    color: '#1A2230',
                    whiteSpace: 'nowrap',
                    overflow: 'hidden',
                    textOverflow: 'ellipsis',
                    direction: 'rtl',
                    textAlign: 'left',
                  }}
                >
                  {f}
                </div>
                <div style={{ ...HINT, marginTop: 2 }}>
                  {root ? `\u2192 ${root}/${fileStem(f)}` : '\u2192 set a results folder below'}
                </div>
              </div>
              <button
                type="button"
                onClick={() => onRemove(i)}
                title="Remove from the list"
                style={{
                  flex: 'none',
                  border: 'none',
                  background: 'none',
                  cursor: 'pointer',
                  padding: 4,
                  lineHeight: 0,
                  color: '#98A2B3',
                }}
              >
                <svg width="14" height="14" viewBox="0 0 24 24" fill="none">
                  <path
                    d="M6 6l12 12M18 6L6 18"
                    stroke="currentColor"
                    strokeWidth="2"
                    strokeLinecap="round"
                  />
                </svg>
              </button>
            </div>
          ))}
        </div>
      )}
      <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
        <button type="button" className="pio-browse" onClick={onAdd} style={BROWSE}>
          {files.length ? 'Add more files' : 'Choose files'}
        </button>
        <span style={HINT}>
          {files.length
            ? `${files.length} file${files.length > 1 ? 's' : ''}, searched separately \u2014 ${files.length} run${files.length > 1 ? 's' : ''} queued.`
            : 'Each file is searched on its own, into its own results folder.'}
        </span>
      </div>
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
  /** Libraries used by earlier SearchDIA runs, most recent first. */
  recentLibraries: string[]
  jobName: string
  /** The name the run will get once collisions are resolved; empty when unset. */
  resolvedJobName: string
  onParam: (key: string, value: string) => void
  onToggle: (key: string) => void
  onBrowse: (key: 'msData' | 'library' | 'results') => void
  /** Add files to the list, via the multi-select picker. */
  onAddMsFiles: () => void
  /** Drop one file from the list, by index. */
  onRemoveMsFile: (index: number) => void
  onJobName: (value: string) => void
  onOpenLoad: () => void
  onGoToBuild: () => void
}

export function SearchDiaForm({
  params,
  notes,
  libInfo,
  recentLibraries,
  jobName,
  resolvedJobName,
  onParam,
  onToggle,
  onBrowse,
  onAddMsFiles,
  onRemoveMsFile,
  onJobName,
  onOpenLoad,
  onGoToBuild,
}: Props) {
  const byFiles = params.msDataMode === 'files'
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
          <h2 style={H2}>Files and name</h2>
          <LoadPreviousButton onClick={onOpenLoad} />
        </div>
        <div style={{ display: 'flex', flexDirection: 'column', gap: 12 }}>
          <div data-key="msData">
            <label style={LABEL}>MS data</label>
            <div style={{ ...SEG_TRACK, marginBottom: 10 }}>
              <button type="button" onClick={() => onParam('msDataMode', 'folder')} style={seg(!byFiles)}>
                One folder
              </button>
              <button type="button" onClick={() => onParam('msDataMode', 'files')} style={seg(byFiles)}>
                Chosen files
              </button>
            </div>
            {byFiles ? <MsFileList params={params} onAdd={onAddMsFiles} onRemove={onRemoveMsFile} /> : (
              <PathRow
                label=""
                fieldKey="msData"
                value={params.msData}
                placeholder="/path/to/ms/data"
                note={notes.msData}
                onChange={onParam}
                onBrowse={() => onBrowse('msData')}
              />
            )}
          </div>
          <PathRow
            label="Spectral library"
            fieldKey="library"
            value={params.library}
            placeholder="/path/to/library.poin"
            note={notes.library}
            onChange={onParam}
            onBrowse={() => onBrowse('library')}
            action={
              <RecentLibraries
                paths={recentLibraries}
                onPick={(p) => onParam('library', p)}
              />
            }
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
                  ({libInfo.variable_mods.map(unimodDisplay).join(', ')}). Pioneer does not report
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
          <JobNameField value={jobName} resolved={resolvedJobName} onChange={onJobName} />
        </div>
      </section>

      {/* Common and Output share a row; Advanced spans the full width beneath.
          The two in the row are near enough the same height -- the q-value field
          and its hint against four toggle rows -- that neither leaves the other
          with a long empty tail. Advanced is one row of three fields, so full
          width suits it better than half, and it is the section you open least.

          Basis 300px, not 340: at the window's 980px minimum the sidebar and
          padding leave about 670px, so 340 each would exceed it and wrap to a
          single column at the default size. At 300 they fit and then grow to
          ~330.

          `stretch` so both cards take the height of the taller one. Their
          contents differ by only a row or so, and a pair of boxes that end at
          slightly different heights reads as a mistake rather than as a
          consequence of what is in them. When the row wraps at narrow widths
          each card is alone on its line, so stretch has nothing to match and
          they return to their natural heights. */}
      <div style={{ display: 'flex', gap: 14, alignItems: 'stretch', flexWrap: 'wrap' }}>
        <section style={{ ...CARD, flex: '1 1 300px', minWidth: 0 }}>
          <h2 style={{ ...H2, marginBottom: 14 }}>Common parameters</h2>
          {/* Stacked, not side by side. One field beside the toggles left a large
              empty quadrant under the field, and no amount of alignment fixes a
              column with nothing else to hold. Full width also lets each toggle
              sit at the card edge, which reads as deliberate. */}
          <div style={{ display: 'flex', flexDirection: 'column', gap: 16 }}>
            <div>
              <NumField fieldKey="qValue" value={params.qValue} onChange={onParam} />
              <div style={{ fontSize: 11.5, color: '#98A2B3', marginTop: 6 }}>1% FDR = 0.01</div>
            </div>
            <div style={{ display: 'flex', flexDirection: 'column', gap: 11 }}>
              <ToggleRow
                title="Run-to-run normalization"
                hint="Retention-time dependent cross-run intensity scaling"
                on={params.runToRunNorm}
                fieldKey="runToRunNorm"
                onToggle={onToggle}
              />
              <ToggleRow
                title="MBR"
                hint="Match between runs"
                on={params.matchBetweenRuns}
                fieldKey="matchBetweenRuns"
                onToggle={onToggle}
              />
            </div>
          </div>
        </section>

        <section style={{ ...CARD, flex: '1 1 300px', minWidth: 0 }}>
          <h2 style={{ ...H2, marginBottom: 14 }}>Output</h2>
          <div style={{ display: 'flex', flexDirection: 'column', gap: 11 }}>
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
              title="Debug logging"
              hint="Verbose console output"
              on={params.debugLogging}
              fieldKey="debugLogging"
              onToggle={onToggle}
            />
          </div>
        </section>
      </div>

      <section style={CARD}>
        <h2 style={{ ...H2, marginBottom: 16 }}>Advanced parameters</h2>
        <div
          style={{
            display: 'grid',
            gridTemplateColumns: '1fr 1fr 1fr',
            gap: 14,
            alignItems: 'stretch',
          }}
        >
          <NumField fieldKey="nIsotopes" value={params.nIsotopes} onChange={onParam} />
          <NumField fieldKey="nce" value={params.nce} onChange={onParam} />
          <NumField fieldKey="minPeptides" value={params.minPeptides} onChange={onParam} />
        </div>
      </section>
    </div>
  )
}
