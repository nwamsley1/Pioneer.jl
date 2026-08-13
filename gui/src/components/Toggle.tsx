/** Port of the `toggle(on)` helper in the design's Component class. */
interface Props {
  on: boolean
  onClick: () => void
  fieldKey: string
}

export function Toggle({ on, onClick, fieldKey }: Props) {
  return (
    <button
      type="button"
      data-key={fieldKey}
      onClick={onClick}
      style={{ background: 'none', border: 'none', padding: 0, cursor: 'pointer', flex: 'none' }}
    >
      <span
        style={{
          // Trimmed from 38x22. The setting labels are 12px, so at the old size
          // the switch was the heaviest thing on the card and on/off state read
          // ahead of the values. Still a comfortable hit target.
          width: 34,
          height: 20,
          borderRadius: 10,
          position: 'relative',
          display: 'inline-block',
          flex: 'none',
          transition: 'background .15s',
          background: on ? 'var(--pio-accent)' : '#CDD2D9',
        }}
      >
        <span
          style={{
            position: 'absolute',
            top: 2,
            left: on ? 16 : 2,
            width: 16,
            height: 16,
            borderRadius: '50%',
            background: '#fff',
            boxShadow: '0 1px 2px rgba(0,0,0,.25)',
            transition: 'left .15s',
          }}
        />
      </span>
    </button>
  )
}
