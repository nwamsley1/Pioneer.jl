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
import { BROWSE, HINT, LABEL } from '../lib/styles'
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

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 14 }}>
      <section style={CARD}>
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
          <h2 style={H2}>Available libraries</h2>
          <button type="button" style={BROWSE} onClick={onRefresh} disabled={loading}>
            {loading ? 'Loading…' : 'Refresh'}
          </button>
        </div>

        {error && <p style={noteStyle({ level: 'error', msg: error })}>{error}</p>}

        {!error && loading && libraries.length === 0 && (
          <p style={{ ...HINT, marginTop: 12 }}>Reading the library repository…</p>
        )}

        {!error && !loading && libraries.length === 0 && (
          <p style={{ ...HINT, marginTop: 12 }}>
            No libraries found in the repository.
          </p>
        )}

        <div style={{ display: 'flex', flexDirection: 'column', gap: 8, marginTop: 12 }}>
          {libraries.map((library) => {
            const active = library.name === params.selected
            return (
              <button
                key={library.name}
                type="button"
                onClick={() => onSelect(library.name)}
                style={{
                  textAlign: 'left',
                  padding: '11px 13px',
                  borderRadius: 10,
                  cursor: 'pointer',
                  border: active ? '1.5px solid #2563EB' : '1px solid #E7EAEE',
                  background: active ? '#F5F8FF' : '#fff',
                }}
              >
                <div style={{ display: 'flex', justifyContent: 'space-between', gap: 12 }}>
                  <span style={{ fontSize: 13, fontWeight: 600, color: '#1B2A4A' }}>
                    {library.title || library.name}
                  </span>
                  <span style={{ fontSize: 12, color: '#667085', whiteSpace: 'nowrap' }}>
                    {library.size_human}
                  </span>
                </div>
                {library.description && (
                  <div style={{ ...HINT, marginTop: 4 }}>{library.description}</div>
                )}
                {library.model && (
                  <div style={{ ...HINT, marginTop: 2 }}>model: {library.model}</div>
                )}
              </button>
            )
          })}
        </div>
      </section>

      {selected && (
        <section style={CARD}>
          <h2 style={H2}>{selected.title || selected.name}</h2>
          <p style={{ ...HINT, marginTop: 6 }}>
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
                <dd style={{ margin: 0, color: '#1B2A4A' }}>{value}</dd>
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
        <p style={noteStyle(destNote)}>{destNote.msg}</p>

        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            gap: 16,
            marginTop: 14,
          }}
        >
          <div>
            <div style={{ fontSize: 13, fontWeight: 600, color: '#344054' }}>
              Replace an existing copy
            </div>
            <div style={{ fontSize: 11.5, color: '#98A2B3' }}>
              Without this, a library already at the destination stops the download
              rather than being overwritten
            </div>
          </div>
          <Toggle on={params.force} fieldKey="force" onClick={onToggleForce} />
        </div>
      </section>
    </div>
  )
}
