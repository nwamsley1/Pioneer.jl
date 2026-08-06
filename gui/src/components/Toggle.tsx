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
          width: 38,
          height: 22,
          borderRadius: 11,
          position: 'relative',
          display: 'inline-block',
          flex: 'none',
          transition: 'background .15s',
          background: on ? '#2E4D7E' : '#CDD2D9',
        }}
      >
        <span
          style={{
            position: 'absolute',
            top: 2,
            left: on ? 18 : 2,
            width: 18,
            height: 18,
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
