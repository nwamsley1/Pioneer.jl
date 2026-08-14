/** The optional name for the run about to be started.
 *
 *  Two call sites: inside SearchDIA's Essentials card, where it sits with the
 *  paths it will be remembered alongside, and in a card of its own above the
 *  other commands' forms, which have no equivalent place to put it. One
 *  component so the two cannot drift.
 */
import { LABEL } from '../lib/styles'

interface Props {
  value: string
  /** The name the run will actually get, once collisions with existing runs
   *  are resolved. Empty when nothing has been typed. */
  resolved: string
  onChange: (value: string) => void
}

export function JobNameField({ value, resolved, onChange }: Props) {
  const trimmed = value.trim()
  return (
    <div>
      <label htmlFor="pio-job-name" style={LABEL}>
        Job name
      </label>
      <input
        id="pio-job-name"
        data-key="jobName"
        value={value}
        placeholder="Optional"
        onChange={(e) => onChange(e.target.value)}
        // Geometry matches the path inputs it now sits under, so the Essentials
        // card reads as one set of fields. The face stays sans: those are paths,
        // this is a name.
        style={{
          width: '100%',
          padding: '9px 12px',
          border: '1px solid #D7DBE0',
          borderRadius: 9,
          fontSize: 13,
          color: '#1A2230',
          background: '#fff',
          outline: 'none',
          boxSizing: 'border-box',
        }}
      />
      {/* Only once there is something to say about the name. The blank case used
          to explain that one would be generated, which "Optional" implies. */}
      {trimmed && (
        <p style={{ fontSize: 11.5, color: '#98A2B3', margin: '8px 0 0' }}>
          {resolved === trimmed
            ? `This run will be called ${resolved}.`
            : `${trimmed} is already taken — this run will be called ${resolved}.`}
        </p>
      )}
    </div>
  )
}
