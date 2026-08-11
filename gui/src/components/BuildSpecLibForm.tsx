/** The BuildSpecLib page: FASTA input, library output, digestion,
 *  modifications and options. Ported from the `isBuild` branch of the design.
 */
import { NumField } from './NumField'
import { Toggle } from './Toggle'
import { HEADER_PRESETS, MOD_PRESETS } from '../lib/fasta'
import { PREDICTION_MODELS, predictionModelById } from '../lib/types'
import type { BuildParams, FastaEntry, HeaderPresetId, ModEntry } from '../lib/types'
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

const SMALL_INPUT: React.CSSProperties = {
  width: '100%',
  padding: '8px 10px',
  border: '1px solid #D7DBE0',
  borderRadius: 8,
  font: "12.5px 'IBM Plex Sans'",
  outline: 'none',
  boxSizing: 'border-box',
}

const MOD_HEADERS = (
  <div style={{ display: 'flex', gap: 8, alignItems: 'center', marginBottom: 6 }}>
    <span style={{ width: 52, flex: 'none', fontSize: 10.5, fontWeight: 600, color: '#98A2B3', textAlign: 'center' }}>
      Site
    </span>
    <span style={{ flex: 1, minWidth: 0, fontSize: 10.5, fontWeight: 600, color: '#98A2B3' }}>
      Description
    </span>
    <span style={{ width: 120, flex: 'none', fontSize: 10.5, fontWeight: 600, color: '#98A2B3' }}>
      Unimod
    </span>
    <span style={{ width: 96, flex: 'none', fontSize: 10.5, fontWeight: 600, color: '#98A2B3', textAlign: 'right' }}>
      Δ mass (Da)
    </span>
    <span style={{ width: 21, flex: 'none' }} />
  </div>
)

const removeBtn = (onClick: () => void, title: string, width?: number): React.ReactNode => (
  <button
    type="button"
    className="pio-remove"
    onClick={onClick}
    title={title}
    style={{
      flex: 'none',
      width,
      background: 'none',
      border: 'none',
      cursor: 'pointer',
      color: '#98A2B3',
      padding: width ? 3 : 2,
      display: 'flex',
    }}
  >
    <svg width={width ? 15 : 16} height={width ? 15 : 16} viewBox="0 0 24 24" fill="none">
      <path d="M6 6l12 12M18 6 6 18" stroke="currentColor" strokeWidth="2" strokeLinecap="round" />
    </svg>
  </button>
)

function ModTable({
  kind,
  mods,
  note,
  onField,
  onRemove,
  onAdd,
}: {
  kind: 'fixed' | 'variable'
  mods: ModEntry[]
  note: string
  onField: (kind: 'fixed' | 'variable', idx: number, field: keyof ModEntry, value: string) => void
  onRemove: (kind: 'fixed' | 'variable', idx: number) => void
  onAdd: (kind: 'fixed' | 'variable', preset: string) => void
}) {
  const cell = (extra: React.CSSProperties): React.CSSProperties => ({
    padding: '8px 10px',
    border: '1px solid #D7DBE0',
    borderRadius: 8,
    outline: 'none',
    boxSizing: 'border-box',
    ...extra,
  })

  return (
    <>
      <div
        style={{
          fontSize: 11,
          fontWeight: 700,
          letterSpacing: '0.05em',
          textTransform: 'uppercase',
          color: '#667085',
          marginBottom: 9,
        }}
      >
        {kind === 'fixed' ? 'Fixed' : 'Variable'}
      </div>
      {MOD_HEADERS}
      <div style={{ display: 'flex', flexDirection: 'column', gap: 8 }} data-key={`${kind}Mods`}>
        {mods.map((m, i) => (
          <div key={i} style={{ display: 'flex', gap: 8, alignItems: 'center' }}>
            <input
              className="pio-input"
              value={m.pattern}
              onChange={(e) => onField(kind, i, 'pattern', e.target.value)}
              placeholder={kind === 'fixed' ? 'C' : 'M'}
              style={cell({ width: 52, flex: 'none', padding: '8px 7px', font: "13px 'IBM Plex Mono'", textAlign: 'center' })}
            />
            <input
              className="pio-input"
              value={m.label}
              onChange={(e) => onField(kind, i, 'label', e.target.value)}
              placeholder={kind === 'fixed' ? 'Carbamidomethyl' : 'Oxidation'}
              style={cell({ flex: 1, minWidth: 0, font: "12.5px 'IBM Plex Sans'" })}
            />
            <input
              className="pio-input"
              value={m.name}
              onChange={(e) => onField(kind, i, 'name', e.target.value)}
              placeholder={kind === 'fixed' ? 'Unimod:4' : 'Unimod:35'}
              style={cell({ width: 120, flex: 'none', font: "12px 'IBM Plex Mono'" })}
            />
            <input
              className="pio-input"
              value={m.mass}
              onChange={(e) => onField(kind, i, 'mass', e.target.value)}
              placeholder={kind === 'fixed' ? '57.021464' : '15.99491'}
              style={cell({ width: 96, flex: 'none', padding: '8px 8px', font: "12px 'IBM Plex Mono'", textAlign: 'right' })}
            />
            {removeBtn(() => onRemove(kind, i), 'Remove', 21)}
          </div>
        ))}
      </div>
      {mods.length === 0 && (
        <div style={{ fontSize: 12, color: '#98A2B3', padding: '4px 0' }}>
          No {kind} modifications.
        </div>
      )}
      <select
        value=""
        onChange={(e) => onAdd(kind, e.target.value)}
        style={{
          marginTop: 9,
          width: '100%',
          padding: '9px 36px 9px 11px',
          border: '1px dashed #CBD2DA',
          borderRadius: 9,
          font: "600 12.5px 'IBM Plex Sans'",
          color: 'var(--pio-accent)',
          background:
            "#FAFBFC url(\"data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='16' height='16' viewBox='0 0 16 16' fill='none'%3E%3Cpath d='M4 6.5l4 4 4-4' stroke='%232E4D7E' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round'/%3E%3C/svg%3E\") no-repeat right 12px center",
          backgroundSize: '16px 16px',
          cursor: 'pointer',
          outline: 'none',
          appearance: 'none',
          WebkitAppearance: 'none',
        }}
      >
        <option value="">+ Add {kind} modification…</option>
        {Object.entries(MOD_PRESETS).map(([id, def]) => (
          <option key={id} value={id}>
            {def.menuLabel}
          </option>
        ))}
      </select>
      {note && (
        <div style={{ marginTop: 8, fontSize: 11.5, color: '#B45309' }}>{note}</div>
      )}
    </>
  )
}

interface Props {
  params: BuildParams
  fastaNotes: Note[]
  libNote: Note
  calibNote: Note
  modNote: { fixed: string; variable: string }
  onParam: (key: string, value: string) => void
  onToggle: (key: string) => void
  onOpenLoad: () => void
  onAddFasta: () => void
  onBrowseFasta: (idx: number) => void
  onFastaField: (idx: number, field: keyof FastaEntry, value: string) => void
  onFastaPreset: (idx: number, preset: HeaderPresetId) => void
  onFastaRegex: (idx: number, key: keyof FastaEntry['regex'], value: string) => void
  onToggleFastaRegex: (idx: number) => void
  onRemoveFasta: (idx: number) => void
  onBrowseLibPath: () => void
  onBrowseCalibration: () => void
  onModField: (kind: 'fixed' | 'variable', idx: number, field: keyof ModEntry, value: string) => void
  onRemoveMod: (kind: 'fixed' | 'variable', idx: number) => void
  onAddMod: (kind: 'fixed' | 'variable', preset: string) => void
}

export function BuildSpecLibForm({
  params,
  fastaNotes,
  libNote,
  calibNote,
  modNote,
  onParam,
  onToggle,
  onOpenLoad,
  onAddFasta,
  onBrowseFasta,
  onFastaField,
  onFastaPreset,
  onFastaRegex,
  onToggleFastaRegex,
  onRemoveFasta,
  onBrowseLibPath,
  onBrowseCalibration,
  onModField,
  onRemoveMod,
  onAddMod,
}: Props) {
  const selectedModel = predictionModelById(params.predictionModel)

  // Mirrors clamp_digest_length_to_model on the Julia side. Shown only when the
  // requested range actually exceeds the model's, so the note appears exactly
  // when the build would narrow it. Non-numeric input is left to the field's own
  // validation rather than second-guessed here.
  const lengthClamp = (() => {
    const lim = selectedModel.peptideLength
    if (!lim) return null
    const lo = Number(params.minLen)
    const hi = Number(params.maxLen)
    if (!Number.isFinite(lo) || !Number.isFinite(hi)) return null
    if (lo >= lim.min && hi <= lim.max) return null
    return { min: Math.max(lo, lim.min), max: Math.min(hi, lim.max) }
  })()

  const pill = (active: boolean): React.CSSProperties => ({
    flex: 1,
    padding: '7px 4px',
    border: 'none',
    borderRadius: 6,
    font: "600 11px 'IBM Plex Sans'",
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
          <h2 style={H2}>FASTA input</h2>
          <button
            type="button"
            className="pio-link"
            onClick={onOpenLoad}
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
        </div>

        {params.fastaFiles.map((f, i) => {
          const note = fastaNotes[i] ?? { level: '' as const, msg: '' }
          const presetLabel = HEADER_PRESETS[f.presetId].label
          return (
            <div
              key={i}
              style={{
                border: '1px solid #E7EAEE',
                borderRadius: 11,
                padding: '13px 15px',
                marginBottom: 10,
                background: '#fff',
              }}
            >
              <div
                style={{
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'space-between',
                  gap: 10,
                  marginBottom: 10,
                }}
              >
                <div
                  style={{
                    fontSize: 13,
                    fontWeight: 700,
                    color: '#1B2A4A',
                    overflow: 'hidden',
                    textOverflow: 'ellipsis',
                    whiteSpace: 'nowrap',
                  }}
                >
                  {(f.path && f.path.split(/[\\/]/).pop()) || 'New FASTA file'}
                </div>
                {removeBtn(() => onRemoveFasta(i), 'Remove')}
              </div>

              <div style={{ display: 'flex', gap: 8 }}>
                <input
                  className="pio-input"
                  value={f.path}
                  onChange={(e) => onFastaField(i, 'path', e.target.value)}
                  placeholder="/path/to/proteome.fasta"
                  style={{
                    flex: 1,
                    padding: '9px 11px',
                    borderRadius: 8,
                    font: "12.5px 'IBM Plex Mono'",
                    outline: 'none',
                    minWidth: 0,
                    border: `1px solid ${note.level === 'error' ? '#E5484D' : '#D7DBE0'}`,
                  }}
                />
                <button
                  type="button"
                  className="pio-browse"
                  onClick={() => onBrowseFasta(i)}
                  style={{
                    flex: 'none',
                    padding: '0 14px',
                    border: '1px solid #D7DBE0',
                    borderRadius: 8,
                    background: '#F8FAFB',
                    font: "600 12px 'IBM Plex Sans'",
                    color: '#344054',
                    cursor: 'pointer',
                  }}
                >
                  Browse
                </button>
              </div>
              {note.msg && (
                <div
                  style={{
                    marginTop: 7,
                    fontSize: 11.5,
                    color: note.level === 'error' ? '#C0392B' : '#B45309',
                  }}
                >
                  ⚠&nbsp; {note.msg}
                </div>
              )}

              <div
                style={{
                  display: 'grid',
                  gridTemplateColumns: '1fr 1.6fr',
                  gap: 12,
                  alignItems: 'end',
                  marginTop: 11,
                }}
              >
                <div>
                  <label
                    style={{
                      display: 'block',
                      fontSize: 11.5,
                      fontWeight: 600,
                      color: '#475467',
                      marginBottom: 5,
                    }}
                  >
                    Species name
                  </label>
                  <input
                    className="pio-input"
                    value={f.name}
                    onChange={(e) => onFastaField(i, 'name', e.target.value)}
                    placeholder="HUMAN"
                    style={SMALL_INPUT}
                  />
                </div>
                <div>
                  <label
                    style={{
                      display: 'block',
                      fontSize: 11.5,
                      fontWeight: 600,
                      color: '#475467',
                      marginBottom: 5,
                    }}
                  >
                    Header format
                  </label>
                  <div style={{ display: 'flex', gap: 3, padding: 3, background: '#EEF1F4', borderRadius: 8 }}>
                    {(['uniprot', 'ensembl', 'refseq', 'custom'] as const).map((id) => (
                      <button
                        key={id}
                        type="button"
                        onClick={() => onFastaPreset(i, id)}
                        style={pill(f.presetId === id)}
                      >
                        {HEADER_PRESETS[id].label}
                      </button>
                    ))}
                  </div>
                </div>
              </div>

              <div
                style={{
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'space-between',
                  marginTop: 10,
                }}
              >
                <span
                  style={{
                    fontSize: 11.5,
                    color: f.auto && f.presetId !== 'custom' ? 'var(--pio-accent)' : '#98A2B3',
                  }}
                >
                  {f.presetId === 'custom'
                    ? 'Custom header regex'
                    : `${f.auto ? 'Auto-detected: ' : 'Format: '}${presetLabel}`}
                </span>
                <button
                  type="button"
                  className="pio-link-underline"
                  onClick={() => onToggleFastaRegex(i)}
                  style={{
                    background: 'none',
                    border: 'none',
                    cursor: 'pointer',
                    font: "600 11.5px 'IBM Plex Sans'",
                    color: 'var(--pio-accent)',
                  }}
                >
                  {f.showRegex ? 'Hide header regex' : 'Header parsing regex'}
                </button>
              </div>

              {f.showRegex && (
                <div
                  style={{
                    marginTop: 12,
                    paddingTop: 12,
                    borderTop: '1px dashed #E2E6EA',
                    display: 'flex',
                    flexDirection: 'column',
                    gap: 9,
                  }}
                >
                  {(['accessions', 'genes', 'proteins', 'organisms'] as const).map((k) => (
                    <div key={k}>
                      <label
                        style={{
                          display: 'block',
                          fontSize: 11,
                          color: '#667085',
                          marginBottom: 3,
                          textTransform: 'capitalize',
                        }}
                      >
                        {k}
                      </label>
                      <input
                        className="pio-input"
                        value={f.regex[k]}
                        onChange={(e) => onFastaRegex(i, k, e.target.value)}
                        style={{
                          width: '100%',
                          padding: '7px 9px',
                          border: '1px solid #D7DBE0',
                          borderRadius: 7,
                          font: "11.5px 'IBM Plex Mono'",
                          outline: 'none',
                          boxSizing: 'border-box',
                        }}
                      />
                    </div>
                  ))}
                </div>
              )}
            </div>
          )
        })}

        {params.fastaFiles.length === 0 && (
          <div
            style={{
              padding: 16,
              border: '1px dashed #E2E6EA',
              borderRadius: 11,
              textAlign: 'center',
              fontSize: 12.5,
              color: '#98A2B3',
              marginBottom: 10,
            }}
          >
            No FASTA files yet — add one to begin.
          </div>
        )}

        <button
          type="button"
          className="pio-dashed-add"
          data-key="fastaAdd"
          onClick={onAddFasta}
          style={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            gap: 8,
            width: '100%',
            padding: '10px 14px',
            border: '1.5px dashed #CBD2DA',
            borderRadius: 10,
            background: 'none',
            cursor: 'pointer',
            font: "600 13px 'IBM Plex Sans'",
            color: 'var(--pio-accent)',
          }}
        >
          <svg width="15" height="15" viewBox="0 0 24 24" fill="none">
            <path d="M12 5v14M5 12h14" stroke="currentColor" strokeWidth="2" strokeLinecap="round" />
          </svg>
          Add FASTA file
        </button>
      </section>

      {/* Not in the design, but Pioneer's own simplified template carries this
          key as "optional but recommended" and warns on every build without
          it. */}
      <section style={CARD}>
        <div style={{ display: 'flex', alignItems: 'baseline', gap: 12, marginBottom: 4 }}>
          <h2 style={H2}>Calibration file</h2>
        </div>
        <p style={{ margin: '0 0 14px', fontSize: 12, color: '#98A2B3', lineHeight: 1.5 }}>
          One MS data file from the experiment this library is for. Pioneer reads it to detect
          fragment and precursor m/z bounds instead of assuming defaults.
        </p>
        <div style={{ display: 'flex', gap: 8 }}>
          <input
            className="pio-input"
            data-key="calibrationFile"
            value={params.calibrationFile}
            onChange={(e) => onParam('calibrationFile', e.target.value)}
            placeholder="/path/to/one_run.arrow  (optional)"
            style={{
              flex: 1,
              padding: '9px 12px',
              borderRadius: 9,
              font: "13px 'IBM Plex Mono'",
              color: '#1A2230',
              background: '#fff',
              outline: 'none',
              minWidth: 0,
              border: `1px solid ${calibNote.level === 'error' ? '#E5484D' : '#D7DBE0'}`,
            }}
          />
          <button
            type="button"
            className="pio-browse"
            onClick={onBrowseCalibration}
            style={{
              flex: 'none',
              padding: '0 14px',
              border: '1px solid #D7DBE0',
              borderRadius: 9,
              background: '#F8FAFB',
              font: "600 12.5px 'IBM Plex Sans'",
              color: '#344054',
              cursor: 'pointer',
            }}
          >
            Browse
          </button>
        </div>
        {calibNote.msg && (
          <div
            style={{
              marginTop: 9,
              fontSize: 12,
              lineHeight: 1.4,
              color: calibNote.level === 'error' ? '#C0392B' : '#B45309',
            }}
          >
            ⚠&nbsp; {calibNote.msg}
          </div>
        )}
      </section>

      <section style={CARD}>
        <div style={{ display: 'flex', alignItems: 'baseline', gap: 12, marginBottom: 14 }}>
          <h2 style={H2}>Library output</h2>
        </div>
        <div style={{ display: 'flex', gap: 8 }}>
          <input
            className="pio-input"
            data-key="libPath"
            value={params.libPath}
            onChange={(e) => onParam('libPath', e.target.value)}
            placeholder="/path/to/output/my_library"
            style={{
              flex: 1,
              padding: '9px 12px',
              borderRadius: 9,
              font: "13px 'IBM Plex Mono'",
              color: '#1A2230',
              background: '#fff',
              outline: 'none',
              minWidth: 0,
              border: `1px solid ${libNote.level === 'error' ? '#E5484D' : '#D7DBE0'}`,
            }}
          />
          <button
            type="button"
            className="pio-browse"
            onClick={onBrowseLibPath}
            style={{
              flex: 'none',
              padding: '0 14px',
              border: '1px solid #D7DBE0',
              borderRadius: 9,
              background: '#F8FAFB',
              font: "600 12.5px 'IBM Plex Sans'",
              color: '#344054',
              cursor: 'pointer',
            }}
          >
            Browse
          </button>
        </div>
        {libNote.msg ? (
          <div
            style={{
              marginTop: 9,
              fontSize: 12,
              lineHeight: 1.4,
              color: libNote.level === 'error' ? '#C0392B' : '#B45309',
            }}
          >
            ⚠&nbsp; {libNote.msg}
          </div>
        ) : (
          params.libPath.trim() && (
            <div style={{ marginTop: 9, fontSize: 11.5, color: '#98A2B3', fontFamily: "'IBM Plex Mono'" }}>
              Creates {params.libPath.trim().endsWith('.poin') ? params.libPath.trim() : `${params.libPath.trim()}.poin`}
            </div>
          )
        )}
      </section>

      <section style={CARD}>
        <h2 style={{ ...H2, marginBottom: 5 }}>Fragment prediction</h2>
        <p style={{ margin: '0 0 14px', fontSize: 12, color: '#98A2B3', lineHeight: 1.5 }}>
          Which model predicts fragment intensities. All are served by Koina, so a
          build needs network access.
        </p>
        <select
          value={params.predictionModel}
          onChange={(e) => onParam('predictionModel', e.target.value)}
          style={{
            width: '100%',
            padding: '9px 36px 9px 11px',
            border: '1px solid #CBD2DA',
            borderRadius: 9,
            font: "600 12.5px 'IBM Plex Sans'",
            color: '#1D2939',
            background:
              "#FFFFFF url(\"data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='16' height='16' viewBox='0 0 16 16' fill='none'%3E%3Cpath d='M4 6.5l4 4 4-4' stroke='%232E4D7E' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round'/%3E%3C/svg%3E\") no-repeat right 12px center",
            backgroundSize: '16px 16px',
            cursor: 'pointer',
            outline: 'none',
            appearance: 'none',
            WebkitAppearance: 'none',
          }}
        >
          {PREDICTION_MODELS.map((m) => (
            <option key={m.id} value={m.id}>
              {m.label}
            </option>
          ))}
        </select>
        <div style={{ marginTop: 9, fontSize: 11.5, color: '#98A2B3', lineHeight: 1.5 }}>
          {selectedModel.note}
          {selectedModel.peptideLength && (
            <>
              {' '}
              Accepts peptides {selectedModel.peptideLength.min}–
              {selectedModel.peptideLength.max} residues.
            </>
          )}
        </div>
      </section>

      <section style={CARD}>
        <div style={{ display: 'flex', alignItems: 'baseline', gap: 12, marginBottom: 14 }}>
          <h2 style={H2}>Digestion</h2>
        </div>
        <div
          style={{
            display: 'grid',
            gridTemplateColumns: '1fr 1fr 1fr',
            gap: 14,
            alignItems: 'start',
          }}
        >
          <NumField fieldKey="minLen" value={params.minLen} onChange={onParam} />
          <NumField fieldKey="maxLen" value={params.maxLen} onChange={onParam} />
          <NumField fieldKey="missedCleav" value={params.missedCleav} onChange={onParam} />
          <NumField fieldKey="minCharge" value={params.minCharge} onChange={onParam} />
          <NumField fieldKey="maxCharge" value={params.maxCharge} onChange={onParam} />
          <NumField fieldKey="maxVarMods" value={params.maxVarMods} onChange={onParam} />
        </div>
        {lengthClamp && (
          <div
            style={{
              marginTop: 14,
              padding: '10px 12px',
              borderRadius: 9,
              background: '#FFF7E6',
              border: '1px solid #FFE0A3',
              fontSize: 11.5,
              color: '#7A5A11',
              lineHeight: 1.5,
            }}
          >
            {selectedModel.label} accepts {selectedModel.peptideLength!.min}–
            {selectedModel.peptideLength!.max} residues. Pioneer will narrow this
            digest to {lengthClamp.min}–{lengthClamp.max} and log a warning;
            peptides outside that range are not built.
          </div>
        )}
      </section>

      <section style={CARD}>
        <h2 style={{ ...H2, marginBottom: 5 }}>Modifications</h2>
        <p style={{ margin: '0 0 16px', fontSize: 12, color: '#98A2B3', lineHeight: 1.5 }}>
          Fixed mods apply to every matching residue; variable mods are optional, up to{' '}
          <strong style={{ color: '#667085', fontWeight: 600 }}>Max var. mods</strong> per peptide.
        </p>
        <ModTable
          kind="fixed"
          mods={params.fixedMods}
          note={modNote.fixed}
          onField={onModField}
          onRemove={onRemoveMod}
          onAdd={onAddMod}
        />
        <div style={{ height: 18 }} />
        <ModTable
          kind="variable"
          mods={params.variableMods}
          note={modNote.variable}
          onField={onModField}
          onRemove={onRemoveMod}
          onAdd={onAddMod}
        />
      </section>

      <section style={CARD}>
        <div style={{ display: 'flex', alignItems: 'baseline', gap: 12, marginBottom: 14 }}>
          <h2 style={H2}>Options</h2>
        </div>
        <div style={{ display: 'flex', flexDirection: 'column', gap: 11 }}>
          {(
            [
              ['addDecoys', 'Add decoys', 'Decoy sequences for FDR control'],
              ['includeContaminants', 'Include contaminants', 'Append the common contaminants set'],
              ['predictFragments', 'Predict fragments', 'Deep model for fragment intensities'],
            ] as const
          ).map(([key, title, hint]) => (
            <div
              key={key}
              style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', gap: 14 }}
            >
              <div>
                <div style={{ fontSize: 13, fontWeight: 600, color: '#344054' }}>{title}</div>
                <div style={{ fontSize: 11.5, color: '#98A2B3' }}>{hint}</div>
              </div>
              <Toggle on={params[key]} fieldKey={key} onClick={() => onToggle(key)} />
            </div>
          ))}
        </div>
      </section>
    </div>
  )
}
