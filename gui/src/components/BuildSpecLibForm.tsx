/** The BuildSpecLib page: FASTA input, library output, digestion,
 *  modifications and options. Ported from the `isBuild` branch of the design.
 */
import { useState } from 'react'

import { InfoDot } from './InfoDot'
import { NumField } from './NumField'
import { Toggle } from './Toggle'
import { HEADER_PRESETS } from '../lib/fasta'
import {
  findMod,
  allowedSiteValues,
  initialSite,
  isRequiredFixedMod,
  modsForModel,
  occupiedResidues,
  siteAllowed,
  siteOptions,
  unimodId,
} from '../lib/koinaMods'
import {
  CUSTOM_ENZYME,
  ENZYMES,
  SAMPLE_SEQUENCE,
  enzymeByPattern,
  previewDigest,
} from '../lib/enzymes'
import { BROWSE, HINT, LABEL } from '../lib/styles'
import {
  PREDICTION_MODELS,
  isPrositModel,
  predictionModelById,
  unlocalizedMods,
} from '../lib/types'
import type { BuildParams, FastaEntry, HeaderPresetId, ModEntry } from '../lib/types'
import { conflictingResidues, modPatternResidues } from '../lib/validate'
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

/** The chevron the form's other selects draw, so the site picker matches them
 *  rather than falling back to the native control. */
const CHEVRON =
  "url(\"data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='16' height='16' viewBox='0 0 16 16' fill='none'%3E%3Cpath d='M4 6.5l4 4 4-4' stroke='%232E4D7E' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round'/%3E%3C/svg%3E\") no-repeat right 8px center"

/** Wide enough for the longest all-sites value the catalogue offers — Methyl's
 *  ten residues — so the column does not resize as modifications are added. */
const SITE_CELL: React.CSSProperties = {
  width: 112,
  flex: 'none',
  padding: '8px 10px',
  font: "12.5px 'IBM Plex Mono'",
  boxSizing: 'border-box',
}

const MOD_HEADERS = (
  <div style={{ display: 'flex', gap: 8, alignItems: 'center', marginBottom: 6 }}>
    <span style={{ width: 112, flex: 'none', fontSize: 10.5, fontWeight: 600, color: '#98A2B3' }}>
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
  modelId,
  modelLabel,
  note,
  onField,
  onRemove,
  occupied,
  conflicts,
  onAdd,
}: {
  kind: 'fixed' | 'variable'
  mods: ModEntry[]
  modelId: string
  modelLabel: string
  note: string
  onField: (kind: 'fixed' | 'variable', idx: number, field: keyof ModEntry, value: string) => void
  onRemove: (kind: 'fixed' | 'variable', idx: number) => void
  /** Residues the fixed modifications already hold. Empty for the fixed table,
   *  which is what defines them. */
  occupied: Set<string>
  /** Residues a fixed and a variable modification both claim. Both tables get
   *  the same set: the conflict is between two rows and either can resolve it,
   *  so marking only one would point at the wrong half as often as not. */
  conflicts: Set<string>
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

  // Label, accession and mass are all determined by the chosen modification, so
  // they are shown rather than edited. Only the site is a real choice.
  const readOnly = (extra: React.CSSProperties): React.CSSProperties =>
    cell({ background: '#F7F8FA', color: '#667085', ...extra })

  // Carbamidomethyl is pinned on C only. It stays offerable as a variable
  // modification wherever the model allows it on another residue -- K on the
  // PTM models -- and drops out only when C is the single site it has.
  // A modification whose every site is already held by a fixed one has no
  // variable form left to add, so it drops out of the menu rather than being
  // offered and then rejected.
  const available = modsForModel(modelId).filter(
    (d) => initialSite(modelId, kind, d, occupied) !== null,
  )
  const inList = new Set(mods.map((m) => unimodId(m.name)).filter((v) => v !== null))

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
        {mods.map((m, i) => {
          // A config loaded from disk can name a modification this model does
          // not accept. Keep the row and flag it rather than dropping data the
          // user did not ask us to touch.
          const def = findMod(modelId, m.name)
          const bad = def === null || !siteAllowed(def, m.pattern)
          // A residue cannot be fixed and variable at once: the fixed mod takes
          // every one of them, so the variable one lands on top. Pioneer rejects
          // the config, so this has to be visible while it is being made.
          const clash = [...modPatternResidues(m.pattern)].filter((r) => conflicts.has(r))
          const warn = bad
            ? `${modelLabel} does not accept ${m.label || m.name || 'this modification'}${
                def ? ` on ${m.pattern}` : ''
              }. Koina will reject the build.`
            : clash.length
              ? `${clash.sort().join(', ')} is claimed by both a fixed and a variable ` +
                `modification. A fixed one takes every matching residue, so the variable ` +
                `one would land on top of it. Remove one, or narrow a site.`
              : undefined
          // Red for the conflict, which stops the build; the existing amber
          // stays for a modification this model merely does not accept.
          const line = clash.length ? '#E5484D' : bad ? '#B45309' : null
          const ink = clash.length ? '#C0392B' : bad ? '#B45309' : null
          return (
            <div key={i} style={{ display: 'flex', gap: 8, alignItems: 'center' }} title={warn}>
              {def ? (
                <select
                  value={m.pattern}
                  onChange={(e) => onField(kind, i, 'pattern', e.target.value)}
                  style={{
                    ...SITE_CELL,
                    padding: '8px 24px 8px 10px',
                    border: `1px solid ${line ?? '#CBD2DA'}`,
                    borderRadius: 9,
                    color: ink ?? '#1D2939',
                    background: `#FFFFFF ${CHEVRON}`,
                    backgroundSize: '14px 14px',
                    cursor: 'pointer',
                    outline: 'none',
                    appearance: 'none',
                    WebkitAppearance: 'none',
                  }}
                >
                  {/* A stored pattern outside the model's sites still needs an
                      option, or the select would silently show a different one. */}
                  {!siteAllowed(def, m.pattern) && <option value={m.pattern}>{m.pattern || '—'}</option>}
                  {/* A fixed alkylation row may not drop C, and a variable one
                      may not take it, so those choices are not offered. */}
                  {siteOptions(def)
                    .filter(
                      (o) =>
                        allowedSiteValues(modelId, kind, m.name, [o.value], occupied).length > 0 ||
                        o.value === m.pattern,
                    )
                    .map((o) => (
                      <option key={o.value} value={o.value}>
                        {o.label}
                      </option>
                    ))}
                </select>
              ) : (
                <div
                  style={readOnly({
                    ...SITE_CELL,
                    borderColor: line ?? '#B45309',
                    color: ink ?? '#B45309',
                  })}
                >
                  {m.pattern || '—'}
                </div>
              )}
              <div
                style={readOnly({
                  flex: 1,
                  minWidth: 0,
                  font: "12.5px 'IBM Plex Sans'",
                  overflow: 'hidden',
                  textOverflow: 'ellipsis',
                  whiteSpace: 'nowrap',
                  ...(line ? { borderColor: line, color: ink ?? undefined } : null),
                })}
              >
                {def ? def.label : m.label || 'Not supported by this model'}
              </div>
              <div style={readOnly({ width: 120, flex: 'none', font: "12px 'IBM Plex Mono'" })}>
                {m.name || '—'}
              </div>
              <div
                style={readOnly({
                  width: 96,
                  flex: 'none',
                  padding: '8px 8px',
                  font: "12px 'IBM Plex Mono'",
                  textAlign: 'right',
                })}
              >
                {m.mass || '—'}
              </div>
              {kind === 'fixed' && isRequiredFixedMod(modelId, m.name, m.pattern) ? (
                // No control at all rather than a disabled one: there is no
                // state in which this becomes removable, and a greyed button
                // invites hunting for the condition that enables it.
                <div style={{ width: 21, flex: 'none' }} title="Required by this model" />
              ) : (
                removeBtn(() => onRemove(kind, i), 'Remove', 21)
              )}
            </div>
          )
        })}
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
        {available.map((d) => (
          <option key={d.unimod} value={String(d.unimod)} disabled={inList.has(d.unimod)}>
            {d.label} ({d.sites.join('')}) · {d.mass >= 0 ? '+' : '−'}
            {Math.abs(d.mass).toFixed(2)}
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
  onEnzyme: (id: string) => void
  cleavageNote: Note
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
  onEnzyme,
  cleavageNote,
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

  /** The peptides the current digestion settings yield from the sample
   *  sequence, or null when the cleavage rule will not compile.
   *
   *  A half-typed length reads as its default rather than blanking the
   *  preview. */
  const previewPeptides = (() => {
    const int = (raw: string, fallback: number) => {
      const v = Number(raw)
      return Number.isFinite(v) && v >= 0 ? Math.round(v) : fallback
    }
    return previewDigest(SAMPLE_SEQUENCE, {
      pattern: params.cleavageRegex,
      specificity: params.digestSpecificity,
      minLength: int(params.minLen, 7),
      maxLength: int(params.maxLen, 40),
      missedCleavages: int(params.missedCleav, 1),
    })
  })()

  // Four peptides is enough to see the shape of the digest; the rest are one
  // click away. Semi specificity runs to dozens even on a 32-residue sequence,
  // and reading them is a fair thing to want -- "+50 more" is a claim the user
  // should be able to check.
  const [peptidesExpanded, setPeptidesExpanded] = useState(false)
  const PREVIEW_SHOWN = 4
  const previewHidden = Math.max((previewPeptides?.length ?? 0) - PREVIEW_SHOWN, 0)

  // Recomputed on every keystroke, so the pair of rows turns red the moment the
  // clash is created rather than when Run is pressed.
  const modConflicts = conflictingResidues(params.fixedMods, params.variableMods)

  /** The fragment ceiling in words: what a precursor at each end of the range
   *  actually gets, and where the clamp starts.
   *
   *  Mirrors frag_bound_polynomials + FragBoundModel's clamp. A slope and an
   *  intercept do not tell you what window a precursor gets, and the clamp is
   *  the part people are surprised by -- 2.00 x 1250 + 10 is 2510 m/z, which no
   *  Orbitrap records. */
  const fragCeilingSummary = (() => {
    const n = (raw: string, fallback: number) => {
      const v = Number(raw)
      return Number.isFinite(v) ? v : fallback
    }
    const fragMax = n(params.fragMzMax, 2020)
    const precMin = n(params.precMzMin, 390)
    const precMax = n(params.precMzMax, 1010)
    const rule = params.fragBoundsRule
    if (rule === 'constant') {
      return `Every precursor gets fragments up to ${fragMax} m/z.`
    }
    const [slope, intercept] =
      rule === 'thermo_auto'
        ? [2.04, 24.1]
        : rule === 'thermo_auto_documented'
          ? [2.0, 10.0]
          : [n(params.fragCeilingSlope, 2), n(params.fragCeilingIntercept, 10)]
    const at = (p: number) => Math.min(slope * p + intercept, fragMax)
    const round = (v: number) => Math.round(v * 10) / 10
    const crossing = slope > 0 ? (fragMax - intercept) / slope : Infinity
    const clamped =
      crossing < precMax
        ? ` Held at ${fragMax} above precursor ${round(Math.max(crossing, precMin))}.`
        : ''
    return (
      `Ceiling ${round(at(precMin))} m/z at precursor ${round(precMin)}, ` +
      `${round(at(precMax))} at ${round(precMax)}.${clamped}`
    )
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
      <section style={CARD} data-drop="fastaAdd">
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
                  style={BROWSE}
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
        <div style={{ display: 'flex', alignItems: 'center', gap: 8, marginBottom: 4 }}>
          <h2 style={H2}>Reference MS file</h2>
          <InfoDot text="Any one run from the experiment this library is for. Pioneer reads its scan headers to detect the fragment and precursor m/z bounds, instead of assuming defaults that may not match your method. It is not used for calibration in the retention-time sense, and nothing from it ends up in the library." />
        </div>
        <p style={{ margin: '0 0 12px', fontSize: 12, color: '#98A2B3', lineHeight: 1.5 }}>
          Recommended. Without one, set the m/z bounds yourself below.
        </p>
        {/* The choice comes first: everything under it depends on it, and
            leaving the file field live while it is ignored invites filling it
            in and wondering why nothing changed. */}
        <div style={{ display: 'flex', gap: 8, marginBottom: 12 }}>
          <button
            type="button"
            onClick={() => !params.autoDetectFragBounds && onToggle('autoDetectFragBounds')}
            style={pill(params.autoDetectFragBounds)}
          >
            Detect from a file
          </button>
          <button
            type="button"
            onClick={() => params.autoDetectFragBounds && onToggle('autoDetectFragBounds')}
            style={pill(!params.autoDetectFragBounds)}
          >
            Set m/z bounds manually
          </button>
        </div>
        {params.autoDetectFragBounds && (
        <div data-drop="calibrationFile" style={{ display: 'flex', gap: 8 }}>
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
            style={BROWSE}
          >
            Browse
          </button>
        </div>
        )}
        {params.autoDetectFragBounds && calibNote.msg && (
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

        {!params.autoDetectFragBounds && (
          <>
            <div
              style={{
                display: 'grid',
                gridTemplateColumns: '1fr 1fr 1fr 1fr',
                gap: 14,
                alignItems: 'start',
              }}
            >
              <NumField fieldKey="fragMzMin" value={params.fragMzMin} onChange={onParam} />
              <NumField fieldKey="fragMzMax" value={params.fragMzMax} onChange={onParam} />
              <NumField fieldKey="precMzMin" value={params.precMzMin} onChange={onParam} />
              <NumField fieldKey="precMzMax" value={params.precMzMax} onChange={onParam} />
            </div>

            <div style={{ display: 'flex', alignItems: 'center', gap: 8, margin: '16px 0 5px' }}>
              <label
                htmlFor="frag-bounds-rule"
                style={{ fontSize: 11, color: '#667085' }}
              >
                Fragment ceiling
              </label>
              <InfoDot text="Thermo instruments in Scan Range Mode = Auto move the MS2 ceiling with the isolation window, so a flat ceiling keeps fragments the instrument never recorded. SCIEX and Bruker hold one fixed range, which is what Constant does. The gain is small — about 0.4% of intensity, and nothing at all at charge 2 — so this is a correctness fix rather than a way to find more peptides." />
            </div>
            <select
              id="frag-bounds-rule"
              className="pio-input"
              value={params.fragBoundsRule}
              onChange={(e) => onParam('fragBoundsRule', e.target.value)}
              style={{
                width: '100%',
                maxWidth: 340,
                padding: '8px 34px 8px 10px',
                border: '1px solid #D7DBE0',
                borderRadius: 8,
                font: "12.5px 'IBM Plex Sans'",
                color: '#1A2230',
                background: `#FFFFFF ${CHEVRON}`,
                backgroundSize: '16px 16px',
                cursor: 'pointer',
                outline: 'none',
                appearance: 'none',
                WebkitAppearance: 'none',
              }}
            >
              <option value="constant">Constant — the max above, at every precursor</option>
              <option value="thermo_auto_documented">Thermo Auto — 2.00 × precursor + 10</option>
              {/* Neither the measured variant nor explicit coefficients are
                  offered. Two Thermo rules differing by 2% invites a choice
                  nobody has grounds to make, and the config key still takes
                  both for a method that needs one. They are rendered only when
                  a loaded config already selected them, so opening and
                  re-saving it does not quietly rewrite the rule. */}
              {params.fragBoundsRule === 'thermo_auto' && (
                <option value="thermo_auto">
                  Thermo Auto measured — 2.04 × precursor + 24.1
                </option>
              )}
              {params.fragBoundsRule === 'custom' && (
                <option value="custom">Custom (from the loaded config)</option>
              )}
            </select>

            {params.fragBoundsRule === 'custom' && (
              <div
                style={{
                  display: 'grid',
                  gridTemplateColumns: '1fr 1fr 2fr',
                  gap: 14,
                  alignItems: 'start',
                  marginTop: 12,
                }}
              >
                <NumField
                  fieldKey="fragCeilingSlope"
                  value={params.fragCeilingSlope}
                  onChange={onParam}
                />
                <NumField
                  fieldKey="fragCeilingIntercept"
                  value={params.fragCeilingIntercept}
                  onChange={onParam}
                />
              </div>
            )}

            {/* The resolved rule, spelled out. A slope and intercept do not say
                what window a precursor actually gets, and the clamp is the part
                people are surprised by. */}
            <div style={{ marginTop: 10, fontSize: 11.5, color: '#98A2B3', lineHeight: 1.5 }}>
              {fragCeilingSummary}
            </div>
          </>
        )}
      </section>

      <section style={CARD}>
        <div style={{ display: 'flex', alignItems: 'baseline', gap: 12, marginBottom: 14 }}>
          <h2 style={H2}>Library output</h2>
        </div>
        <div data-drop="libPath" style={{ display: 'flex', gap: 8 }}>
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
            style={BROWSE}
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
            gridTemplateColumns: 'minmax(0, 2fr) minmax(0, 1fr)',
            gap: 14,
            alignItems: 'end',
          }}
        >
          <div>
        <label style={{ ...LABEL }}>Enzyme</label>
        <select
          data-key="enzyme"
          value={params.customEnzyme ? CUSTOM_ENZYME : (enzymeByPattern(params.cleavageRegex)?.id ?? CUSTOM_ENZYME)}
          onChange={(e) => onEnzyme(e.target.value)}
          style={{
            width: '100%',
            height: 36,
            padding: '0 9px',
            border: '1px solid #D6DAE1',
            borderRadius: 9,
            fontSize: 13,
            color: '#1B2A4A',
            background: '#fff',
            boxSizing: 'border-box',
          }}
        >
          {ENZYMES.map((e) => (
            <option key={e.id} value={e.id}>
              {e.label} — {e.rule}
            </option>
          ))}
          <option value={CUSTOM_ENZYME}>Custom cleavage rule…</option>
        </select>
          </div>
          <div>
            <label
              htmlFor="digest-specificity"
              style={{ display: 'block', fontSize: 11, color: '#667085', marginBottom: 5 }}
            >
              Specificity
            </label>
            <select
              id="digest-specificity"
              className="pio-input"
              value={params.digestSpecificity}
              onChange={(e) => onParam('digestSpecificity', e.target.value)}
              style={{
                width: '100%',
                padding: '8px 34px 8px 10px',
                border: '1px solid #D7DBE0',
                borderRadius: 8,
                font: "12.5px 'IBM Plex Sans'",
                color: '#1A2230',
                background: `#FFFFFF ${CHEVRON}`,
                backgroundSize: '16px 16px',
                cursor: 'pointer',
                outline: 'none',
                appearance: 'none',
                WebkitAppearance: 'none',
              }}
            >
              <option value="full">Full (both termini)</option>
              <option value="semi">Semi (either terminus)</option>
              <option value="semi-n">Semi-N (C terminus required)</option>
              <option value="semi-c">Semi-C (N terminus required)</option>
            </select>
          </div>
        </div>

        {(params.customEnzyme || !enzymeByPattern(params.cleavageRegex)) && (
          <div style={{ marginTop: 12 }}>
            <label style={LABEL}>Cleavage regex</label>
            <input
              data-key="cleavageRegex"
              value={params.cleavageRegex}
              onChange={(e) => onParam('cleavageRegex', e.target.value)}
              spellCheck={false}
              style={{
                width: '100%',
                height: 36,
                padding: '0 11px',
                border: `1px solid ${cleavageNote.level === 'error' ? '#F0B4A8' : '#D6DAE1'}`,
                borderRadius: 9,
                font: "12.5px 'IBM Plex Mono'",
                color: '#1B2A4A',
                background: '#fff',
                boxSizing: 'border-box',
              }}
            />
            <p style={{ ...HINT, marginTop: 6 }}>
              The first element is the residue cleaved after; anything following it is
              context. Use a lookahead to cut before a residue, as Asp-N does.
            </p>
          </div>
        )}

        {/* What the rule does, rather than a verdict on it: a pattern can be
            valid and still not be the digest you meant, and a compile check
            cannot tell you that. */}
        <div
          style={{
            marginTop: 10,
            padding: '9px 11px',
            borderRadius: 9,
            background: cleavageNote.level === 'error' ? '#FEF2F2' : '#F7F8FA',
            border: `1px solid ${cleavageNote.level === 'error' ? '#FECACA' : '#EDEFF3'}`,
          }}
        >
          <div style={{ fontSize: 11, color: '#98A2B3', marginBottom: 4 }}>
            {SAMPLE_SEQUENCE}
          </div>
          <div
            style={{
              font: "12.5px 'IBM Plex Mono'",
              color: cleavageNote.level === 'error' ? '#C0392B' : '#1B2A4A',
              wordBreak: 'break-word',
              // Expanded, a semi digest is dozens of peptides. Scroll them
              // rather than pushing the rest of the form down the page.
              ...(peptidesExpanded ? { maxHeight: 132, overflowY: 'auto' as const } : {}),
            }}
          >
            {cleavageNote.level === 'error' ? (
              cleavageNote.msg
            ) : !previewPeptides ? null : previewPeptides.length === 0 ? (
              'No peptides — nothing survives these length limits.'
            ) : (
              <>
                {previewPeptides.length} peptide{previewPeptides.length === 1 ? '' : 's'} ·{' '}
                {(peptidesExpanded
                  ? previewPeptides
                  : previewPeptides.slice(0, PREVIEW_SHOWN)
                ).join(' · ')}
                {previewHidden > 0 && (
                  <>
                    {' '}
                    <button
                      type="button"
                      className="pio-link-underline"
                      onClick={() => setPeptidesExpanded((v) => !v)}
                      style={{
                        background: 'none',
                        border: 'none',
                        padding: 0,
                        cursor: 'pointer',
                        font: "12.5px 'IBM Plex Mono'",
                        color: 'var(--pio-accent)',
                      }}
                    >
                      {peptidesExpanded ? 'show fewer' : `+${previewHidden} more`}
                    </button>
                  </>
                )}
              </>
            )}
          </div>
        </div>

        <div style={{ height: 14 }} />
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
          occupied={new Set<string>()}
          conflicts={modConflicts}
          mods={params.fixedMods}
          modelId={params.predictionModel}
          modelLabel={selectedModel.label}
          note={modNote.fixed}
          onField={onModField}
          onRemove={onRemoveMod}
          onAdd={onAddMod}
        />
        <div style={{ height: 18 }} />
        {isPrositModel(params.predictionModel) &&
          unlocalizedMods(params.variableMods).length > 0 && (
          <div
            style={{
              marginTop: 12,
              padding: '10px 12px',
              borderRadius: 9,
              background: '#FFFBEB',
              border: '1px solid #FDE68A',
              fontSize: 12,
              lineHeight: 1.5,
              color: '#92400E',
            }}
          >
            <strong style={{ fontWeight: 600 }}>Experimental.</strong> Searching a
            Prosit-predicted library that carries variable modifications is not yet a
            supported combination: Pioneer does not report site-localization confidence,
            so a modified residue is placed but the placement is not scored.
          </div>
        )}
        <ModTable
          kind="variable"
          occupied={occupiedResidues(params.fixedMods)}
          conflicts={modConflicts}
          mods={params.variableMods}
          modelId={params.predictionModel}
          modelLabel={selectedModel.label}
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
              ['debugLogging', 'Debug logging', 'Verbose console output'],
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
