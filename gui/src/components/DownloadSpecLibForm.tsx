/** The DownloadSpecLib page.
 *
 *  Complementary to BuildSpecLib: rather than predicting a library from a
 *  FASTA, fetch one that has already been built and published.
 *
 *  The catalog comes from `DownloadSpecLib --list --json`, so this component
 *  never talks to the network itself — the webview could not anyway, since the
 *  CSP allows no external origin. Everything shown here is derived by Pioneer
 *  from the repository, which keeps one description of a library rather than
 *  two that can disagree.
 */
import { useMemo, useState } from 'react'
import { BROWSE, HINT, LABEL } from '../lib/styles'
import { parseLibraryName } from '../lib/libraryName'
import { InfoDot } from './InfoDot'
import { Toggle } from './Toggle'
import type { DownloadParams, RemoteLibrary } from '../lib/types'
import { type Note } from '../lib/validate'

const CARD: React.CSSProperties = {
  background: '#fff',
  border: '1px solid #E7EAEE',
  borderRadius: 13,
  padding: '18px 20px',
}

const H2: React.CSSProperties = {
  margin: 0,
  fontSize: 14,
  fontWeight: 650,
  color: '#1B2A4A',
}

const noteStyle = (note: Note): React.CSSProperties => ({
  marginTop: 9,
  fontSize: 12,
  lineHeight: 1.4,
  color: note.level === 'error' ? '#C0392B' : note.level === 'warn' ? '#B45309' : '#98A2B3',
})

const INPUT: React.CSSProperties = {
  flex: 1,
  height: 36,
  padding: '0 11px',
  border: '1px solid #D6DAE1',
  borderRadius: 9,
  fontSize: 13,
  color: '#1B2A4A',
  background: '#fff',
  minWidth: 0,
}

interface Props {
  params: DownloadParams
  libraries: RemoteLibrary[]
  /** True while the catalog is being fetched. */
  loading: boolean
  /** Set when the catalog could not be read at all. */
  error: string
  destNote: Note
  onParam: (key: keyof DownloadParams, value: string) => void
  onToggleForce: () => void
  onSelect: (name: string) => void
  onBrowseDest: () => void
  onRefresh: () => void
}

/** A small turning ring, sized to sit inside a line of text. */
function Spinner() {
  return (
    <span
      aria-hidden="true"
      style={{
        display: 'inline-block',
        width: 12,
        height: 12,
        borderRadius: '50%',
        border: '2px solid #C7D2FE',
        borderTopColor: '#2563EB',
        animation: 'pio-spin 0.8s linear infinite',
        flex: 'none',
      }}
    />
  )
}

/** A bar sliding along a thin track, for work with no known length. */
function BusyBar() {
  return (
    <div
      role="progressbar"
      aria-label="Loading libraries"
      style={{
        position: 'relative',
        height: 3,
        marginTop: 12,
        borderRadius: 2,
        overflow: 'hidden',
        background: '#E7EAEE',
      }}
    >
      <div
        style={{
          position: 'absolute',
          inset: 0,
          width: '40%',
          borderRadius: 2,
          background: '#2563EB',
          animation: 'pio-slide 1.1s ease-in-out infinite',
        }}
      />
    </div>
  )
}

/** What a row shows. A curated manifest title wins over the parsed name; a
 *  name that follows no convention is shown whole. */
interface Row {
  library: RemoteLibrary
  organism: string
  latin: string
  proteomes: string[]
  model: string
  release: string
}

function toRow(library: RemoteLibrary): Row {
  const parsed = parseLibraryName(library.name)
  const curated = library.title && library.title !== library.name ? library.title : ''
  return {
    library,
    organism: curated || parsed?.common || library.name,
    latin: parsed?.latin ?? '',
    proteomes: parsed?.proteomes ?? [],
    model: library.model || parsed?.model || '',
    release: parsed ? parsed.release.replace('_', '-') : '',
  }
}

type SortKey = 'organism' | 'model' | 'release' | 'size'

const SORTERS: Record<SortKey, (a: Row, b: Row) => number> = {
  organism: (a, b) => a.organism.localeCompare(b.organism),
  model: (a, b) => a.model.localeCompare(b.model),
  release: (a, b) => a.release.localeCompare(b.release),
  size: (a, b) => a.library.total_bytes - b.library.total_bytes,
}

const TH: React.CSSProperties = {
  textAlign: 'left',
  padding: '7px 10px',
  fontSize: 11.5,
  fontWeight: 600,
  color: '#667085',
  borderBottom: '1px solid #E7EAEE',
  whiteSpace: 'nowrap',
  userSelect: 'none',
}

const TD: React.CSSProperties = {
  padding: '8px 10px',
  fontSize: 13,
  color: '#1B2A4A',
  borderBottom: '1px solid #F0F2F5',
  verticalAlign: 'middle',
}

const MONO: React.CSSProperties = { fontFamily: '"IBM Plex Mono", monospace', fontSize: 12 }

function LibraryTable({
  libraries,
  selected,
  onSelect,
}: {
  libraries: RemoteLibrary[]
  selected: string
  onSelect: (name: string) => void
}) {
  const [filter, setFilter] = useState('')
  const [sort, setSort] = useState<{ key: SortKey; desc: boolean }>({ key: 'organism', desc: false })

  const rows = useMemo(() => {
    const q = filter.trim().toLowerCase()
    const all = libraries.map(toRow).filter(
      (r) =>
        !q ||
        [r.organism, r.latin, r.model, r.library.name, ...r.proteomes].some((s) =>
          s.toLowerCase().includes(q),
        ),
    )
    all.sort(SORTERS[sort.key])
    if (sort.desc) all.reverse()
    return all
  }, [libraries, filter, sort])

  const header = (key: SortKey | null, label: string, align: 'left' | 'right' = 'left') => (
    <th
      style={{ ...TH, textAlign: align, cursor: key ? 'pointer' : 'default' }}
      onClick={
        key ? () => setSort((s) => ({ key, desc: s.key === key ? !s.desc : false })) : undefined
      }
    >
      {label}
      {key && sort.key === key ? (sort.desc ? ' ▾' : ' ▴') : ''}
    </th>
  )

  return (
    <>
      <input
        style={{ ...INPUT, width: '100%', marginTop: 12, boxSizing: 'border-box' }}
        value={filter}
        placeholder="Filter by organism, species or proteome ID"
        onChange={(e) => setFilter(e.target.value)}
      />
      <div style={{ marginTop: 10, border: '1px solid #E7EAEE', borderRadius: 10, overflow: 'hidden' }}>
        <table style={{ width: '100%', borderCollapse: 'collapse' }}>
          <thead style={{ background: '#F8F9FB' }}>
            <tr>
              {header('organism', 'Organism')}
              {header(null, 'Proteome')}
              {header('model', 'Model')}
              {header('release', 'UniProt')}
              {header('size', 'Size', 'right')}
            </tr>
          </thead>
          <tbody>
            {rows.map((r) => {
              const active = r.library.name === selected
              return (
                <tr
                  key={r.library.name}
                  onClick={() => onSelect(r.library.name)}
                  title={r.library.name}
                  style={{ cursor: 'pointer', background: active ? '#F5F8FF' : undefined }}
                >
                  <td style={TD}>
                    <label style={{ display: 'flex', alignItems: 'center', gap: 9, cursor: 'pointer' }}>
                      <input
                        type="radio"
                        name="library"
                        checked={active}
                        onChange={() => onSelect(r.library.name)}
                        style={{ accentColor: '#2563EB', margin: 0 }}
                      />
                      <span>
                        <span style={{ fontWeight: 600 }}>{r.organism}</span>
                        {r.latin && r.latin !== r.organism && (
                          <span style={{ display: 'block', fontSize: 11.5, color: '#667085', fontStyle: 'italic' }}>
                            {r.latin}
                          </span>
                        )}
                      </span>
                    </label>
                  </td>
                  <td style={{ ...TD, ...MONO, color: '#475467' }}>
                    {r.proteomes.map((p) => (
                      <div key={p}>{p}</div>
                    ))}
                  </td>
                  <td style={TD}>{r.model}</td>
                  <td style={{ ...TD, ...MONO }}>{r.release}</td>
                  <td style={{ ...TD, textAlign: 'right', whiteSpace: 'nowrap', color: '#475467' }}>
                    {r.library.size_human}
                  </td>
                </tr>
              )
            })}
            {rows.length === 0 && (
              <tr>
                <td colSpan={5} style={{ ...TD, ...HINT }}>
                  No library matches “{filter}”.
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>
    </>
  )
}

export function DownloadSpecLibForm({
  params,
  libraries,
  loading,
  error,
  destNote,
  onParam,
  onToggleForce,
  onSelect,
  onBrowseDest,
  onRefresh,
}: Props) {
  const selected = libraries.find((l) => l.name === params.selected)
  const selectedRow = selected && toRow(selected)

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 14 }}>
      <section style={CARD}>
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
          <h2 style={H2}>Available libraries</h2>
          <button
            type="button"
            style={{ ...BROWSE, display: 'inline-flex', alignItems: 'center', gap: 7 }}
            onClick={onRefresh}
            disabled={loading}
          >
            {loading && <Spinner />}
            {loading ? 'Loading…' : 'Refresh'}
          </button>
        </div>

        {loading && <BusyBar />}

        {error && <p style={noteStyle({ level: 'error', msg: error })}>{error}</p>}

        {!error && loading && libraries.length === 0 && (
          <p style={{ ...HINT, marginTop: 12 }}>Reading the library repository…</p>
        )}

        {!error && !loading && libraries.length === 0 && (
          <p style={{ ...HINT, marginTop: 12 }}>
            No libraries found in the repository.
          </p>
        )}

        {libraries.length > 0 && (
          <LibraryTable
            libraries={libraries}
            selected={params.selected}
            onSelect={onSelect}
          />
        )}
      </section>

      {selected && (
        <section style={CARD}>
          <h2 style={{ ...H2, overflowWrap: 'anywhere' }}>{selectedRow?.organism}</h2>
          <p style={{ ...HINT, marginTop: 6, overflowWrap: 'anywhere' }}>
            {selected.name} · {selected.size_human} in {selected.n_files} files
          </p>
          {selected.recommended_for && (
            <p style={{ fontSize: 12.5, color: '#344054', marginTop: 8 }}>
              {selected.recommended_for}
            </p>
          )}
          <dl
            style={{
              display: 'grid',
              gridTemplateColumns: 'max-content 1fr',
              gap: '5px 14px',
              margin: '12px 0 0',
              fontSize: 12.5,
            }}
          >
            {Object.entries(selected.details).map(([label, value]) => (
              <div key={label} style={{ display: 'contents' }}>
                <dt style={{ color: '#667085' }}>{label}</dt>
                <dd style={{ margin: 0, color: '#1B2A4A', minWidth: 0, overflowWrap: 'anywhere' }}>{value}</dd>
              </div>
            ))}
          </dl>
        </section>
      )}

      <section style={CARD}>
        <h2 style={H2}>Destination</h2>
        <label style={{ ...LABEL, marginTop: 12 }}>Download into</label>
        <div style={{ display: 'flex', gap: 8 }}>
          <input
            style={INPUT}
            value={params.dest}
            placeholder="Choose a folder"
            onChange={(e) => onParam('dest', e.target.value)}
          />
          <button type="button" style={BROWSE} onClick={onBrowseDest}>
            Browse
          </button>
        </div>
        {destNote.msg && <p style={noteStyle(destNote)}>{destNote.msg}</p>}

        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            gap: 16,
            marginTop: 14,
          }}
        >
          <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
            <div style={{ fontSize: 13, fontWeight: 600, color: '#344054' }}>
              Overwrite if already downloaded
            </div>
            <InfoDot text="Otherwise the download stops if this library is already in the folder." />
          </div>
          <Toggle on={params.force} fieldKey="force" onClick={onToggleForce} />
        </div>
      </section>
    </div>
  )
}
