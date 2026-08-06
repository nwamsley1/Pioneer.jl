/** "Load a previous configuration": point at a results folder or a config.json
 *  and fill the form from it. Ported from the LOAD CONFIG MODAL block.
 *
 *  Unlike the prototype, Browse opens a real dialog and the file is read from
 *  disk — the textarea shows what was actually read, and stays editable.
 */
import { useState } from 'react'
import { pickFile, pickFolder, readConfig } from '../lib/backend'

interface Props {
  onClose: () => void
  /** Returns an error message, or null on success. */
  onApply: (draft: string) => string | null
}

export function LoadConfigModal({ onClose, onApply }: Props) {
  const [path, setPath] = useState('')
  const [draft, setDraft] = useState('')
  const [error, setError] = useState('')

  const load = async (p: string) => {
    setPath(p)
    try {
      const text = await readConfig(p)
      setDraft(text)
      try {
        JSON.parse(text)
        setError('')
      } catch (e) {
        setError((e as Error).message)
      }
    } catch (e) {
      setDraft('')
      setError(String(e))
    }
  }

  const browseFile = async () => {
    const picked = await pickFile('Choose a config.json', 'Pioneer config', ['json'])
    if (picked) await load(picked)
  }

  const browseFolder = async () => {
    const picked = await pickFolder('Choose a results folder')
    if (picked) await load(picked)
  }

  const onDraft = (v: string) => {
    setDraft(v)
    try {
      JSON.parse(v)
      setError('')
    } catch (e) {
      setError((e as Error).message)
    }
  }

  const apply = () => {
    const err = onApply(draft)
    if (err) setError(err)
  }

  const browseBtn: React.CSSProperties = {
    flex: 'none',
    padding: '0 16px',
    border: '1px solid #D7DBE0',
    borderRadius: 9,
    background: '#F8FAFB',
    font: "600 12.5px 'IBM Plex Sans'",
    color: '#344054',
    cursor: 'pointer',
  }

  return (
    <div
      onClick={onClose}
      style={{
        position: 'absolute',
        inset: 0,
        background: 'rgba(15,20,27,0.45)',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        zIndex: 25,
        padding: 40,
      }}
    >
      <div
        onClick={(e) => e.stopPropagation()}
        style={{
          background: '#fff',
          borderRadius: 14,
          width: 600,
          maxWidth: '100%',
          maxHeight: '100%',
          display: 'flex',
          flexDirection: 'column',
          overflow: 'hidden',
          boxShadow: '0 20px 60px rgba(0,0,0,0.3)',
        }}
      >
        <div style={{ padding: '18px 22px', borderBottom: '1px solid #EEF1F4' }}>
          <div style={{ fontSize: 15, fontWeight: 700, color: '#1B2A4A' }}>
            Load a previous configuration
          </div>
          <div style={{ fontSize: 12.5, color: '#667085', marginTop: 3, lineHeight: 1.5 }}>
            Point to a results folder or a{' '}
            <code style={{ fontFamily: "'IBM Plex Mono'", fontSize: 12 }}>config.json</code> from an
            earlier run. Its parameters fill the form so you can rerun it or adjust a few settings.
          </div>
        </div>
        <div style={{ padding: '18px 22px', display: 'flex', flexDirection: 'column', gap: 16 }}>
          <div>
            <label
              style={{
                display: 'block',
                fontSize: 12,
                fontWeight: 600,
                color: '#344054',
                marginBottom: 7,
              }}
            >
              Folder or config.json
            </label>
            <div style={{ display: 'flex', gap: 8 }}>
              <input
                className="pio-input"
                value={path}
                onChange={(e) => setPath(e.target.value)}
                onBlur={() => path.trim() && load(path)}
                placeholder="/path/to/previous_run"
                style={{
                  flex: 1,
                  padding: '11px 13px',
                  border: '1px solid #D7DBE0',
                  borderRadius: 9,
                  font: "13px 'IBM Plex Mono'",
                  outline: 'none',
                  minWidth: 0,
                }}
              />
              <button type="button" className="pio-browse" onClick={browseFolder} style={browseBtn}>
                Folder
              </button>
              <button type="button" className="pio-browse" onClick={browseFile} style={browseBtn}>
                File
              </button>
            </div>
          </div>
          <div>
            <div
              style={{
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'space-between',
                marginBottom: 7,
              }}
            >
              <label style={{ fontSize: 12, fontWeight: 600, color: '#344054' }}>config.json</label>
              <span style={{ fontSize: 11, color: '#98A2B3' }}>read from file · editable</span>
            </div>
            <textarea
              className="pio-textarea"
              value={draft}
              onChange={(e) => onDraft(e.target.value)}
              spellCheck={false}
              style={{
                width: '100%',
                height: '32vh',
                resize: 'vertical',
                boxSizing: 'border-box',
                background: '#1B2A4A',
                color: '#C8D0DA',
                border: '1px solid #D7DBE0',
                borderRadius: 9,
                padding: 13,
                font: "12px/1.6 'IBM Plex Mono'",
                outline: 'none',
              }}
            />
            <div
              style={{
                minHeight: 18,
                padding: '6px 2px 0',
                font: "12px 'IBM Plex Mono'",
                color: error ? '#D9534F' : 'transparent',
              }}
            >
              {error ? `⚠  ${error}` : '.'}
            </div>
          </div>
        </div>
        <div
          style={{
            padding: '14px 22px',
            borderTop: '1px solid #EEF1F4',
            display: 'flex',
            justifyContent: 'flex-end',
            gap: 10,
          }}
        >
          <button
            type="button"
            className="pio-light-btn"
            onClick={onClose}
            style={{
              padding: '9px 18px',
              borderRadius: 9,
              border: '1px solid #D7DBE0',
              background: '#fff',
              color: '#344054',
              font: "600 13px 'IBM Plex Sans'",
              cursor: 'pointer',
            }}
          >
            Cancel
          </button>
          <button
            type="button"
            onClick={apply}
            disabled={!!error || !draft.trim()}
            style={{
              padding: '9px 18px',
              borderRadius: 9,
              border: 'none',
              font: "600 13px 'IBM Plex Sans'",
              ...(error || !draft.trim()
                ? { background: '#A7D8CF', color: '#fff', cursor: 'not-allowed' }
                : { background: '#2E4D7E', color: '#FFFFFF', cursor: 'pointer' }),
            }}
          >
            Load parameters
          </button>
        </div>
      </div>
    </div>
  )
}
