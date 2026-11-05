/** The full license text, opened from Settings. The text is bundled into the
 *  app at build time from the repository's LICENSE, so it is always the terms
 *  this build was released under.
 */
import licenseText from '../../../LICENSE?raw'

export function LicenseModal({ onClose }: { onClose: () => void }) {
  return (
    <div
      onClick={onClose}
      style={{
        position: 'fixed',
        inset: 0,
        background: 'rgba(15,20,27,0.45)',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        zIndex: 60,
        padding: 40,
      }}
    >
      <div
        role="dialog"
        aria-label="License"
        onClick={(e) => e.stopPropagation()}
        style={{
          background: 'var(--pio-nav)',
          borderRadius: 14,
          width: 640,
          maxWidth: '100%',
          maxHeight: '100%',
          display: 'flex',
          flexDirection: 'column',
          overflow: 'hidden',
          boxShadow: '0 20px 60px rgba(0,0,0,0.4)',
        }}
      >
        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            padding: '16px 20px',
            borderBottom: '1px solid rgba(255,255,255,0.08)',
          }}
        >
          <span style={{ fontSize: 13, fontWeight: 600, color: '#E7EBF0' }}>License</span>
          <button
            type="button"
            onClick={onClose}
            aria-label="Close"
            style={{
              background: 'none',
              border: 'none',
              color: '#98A6BC',
              fontSize: 18,
              cursor: 'pointer',
              lineHeight: 1,
            }}
          >
            ×
          </button>
        </div>
        <pre
          style={{
            margin: 0,
            padding: '16px 20px 20px',
            overflow: 'auto',
            whiteSpace: 'pre-wrap',
            font: "12px/1.55 'IBM Plex Mono'",
            color: '#C8D0DA',
          }}
        >
          {licenseText}
        </pre>
      </div>
    </div>
  )
}
