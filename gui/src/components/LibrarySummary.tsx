/** What the selected `.poin` says about itself, shown under the library field.
 *
 *  A library path is otherwise just a string, and the wrong one is not obvious
 *  until a search has run for an hour and found nothing. This turns it into
 *  something checkable: which model predicted it, what was digested, which
 *  modifications it covers.
 *
 *  Everything comes from the library's own config. Older libraries predate some
 *  of the keys — `prediction_model` in particular — so anything absent is
 *  omitted rather than guessed at.
 */
import type { LibraryInfo } from '../lib/backend'

const MODEL_LABELS: Record<string, string> = {
  altimeter: 'Altimeter',
  prosit_2020_hcd: 'Prosit 2020 HCD',
  prosit_2024_ptm: 'Prosit 2024 PTM',
  prosit_2025_40ptm: 'Prosit 2025 40-PTM',
}

function Row({ label, children }: { label: string; children: React.ReactNode }) {
  return (
    <div style={{ display: 'flex', gap: 8, alignItems: 'baseline' }}>
      <span style={{ flex: 'none', width: 96, fontSize: 11.5, color: '#98A2B3' }}>{label}</span>
      <span style={{ fontSize: 11.5, color: '#475467', minWidth: 0 }}>{children}</span>
    </div>
  )
}

export function LibrarySummary({ info }: { info: LibraryInfo | null }) {
  if (!info || !info.is_library) return null

  if (info.error) {
    return (
      <div style={{ marginTop: 9, fontSize: 11.5, color: '#B45309' }}>
        ⚠&nbsp; {info.error}
      </div>
    )
  }

  const digest = [
    info.length_range && `${info.length_range} residues`,
    info.charge_range && `charge ${info.charge_range}`,
    info.missed_cleavages && `${info.missed_cleavages} missed cleavage${info.missed_cleavages === '1' ? '' : 's'}`,
  ].filter(Boolean)

  const mods = [...info.fixed_mods.map((m) => `${m} (fixed)`), ...info.variable_mods]
  // Unrecorded counts as Altimeter: BuildSpecLib produced nothing else before
  // the key existed, and both libraries in hand -- the ecoli test one and the
  // published human one -- have no prediction_model while recording nce 26.0.
  // Requiring the key would leave precisely the libraries this row is for still
  // showing a number. The Model row above stays honest about not knowing, which
  // is a question about provenance rather than about what the number means.
  const isAltimeter = info.prediction_model === 'altimeter' || info.prediction_model === ''

  return (
    <div
      style={{
        marginTop: 10,
        padding: '10px 12px',
        borderRadius: 9,
        background: '#F7F8FA',
        border: '1px solid #E7EAEE',
        display: 'flex',
        flexDirection: 'column',
        gap: 5,
      }}
    >
      <Row label="Model">
        {info.prediction_model
          ? MODEL_LABELS[info.prediction_model] ?? info.prediction_model
          : // Not a failure: libraries built before the key existed simply do
            // not record it, and saying so beats implying a default.
            <span style={{ color: '#98A2B3' }}>not recorded in this library</span>}
      </Row>
      {digest.length > 0 && <Row label="Digest">{digest.join(' · ')}</Row>}
      {info.fastas.length > 0 && (
        <Row label="FASTA">
          {info.fastas.join(', ')}
          {info.include_contaminants && !info.fastas.includes('CONTAM') && ' · contaminants'}
        </Row>
      )}
      {mods.length > 0 && <Row label="Mods">{mods.join(', ')}</Row>}
      {/* Altimeter predicts spline coefficients across collision energy rather
          than at one setting, so the number a config happens to record is not
          the energy the library is good for. Saying "splines" states what the
          library actually is; printing 26.0 invites matching it to an
          instrument method that has nothing to do with it.

          Only when the model is recorded as Altimeter: an unrecorded model is
          left as the number the config gives, since the row above declines to
          guess at the model and this one should not either. */}
      {(isAltimeter || info.nce) && (
        <Row label="NCE">{isAltimeter ? 'splines' : info.nce}</Row>
      )}
      <Row label="Decoys">{info.has_decoys ? 'included' : 'none'}</Row>
    </div>
  )
}
