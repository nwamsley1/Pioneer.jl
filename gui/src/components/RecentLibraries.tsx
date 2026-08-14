/** Pick a spectral library you have searched with before.
 *
 *  The paths come from the run history; whether they are still there comes from
 *  the filesystem, checked when the menu opens rather than up front. Checking on
 *  load would stat every path at startup for a menu that is usually not opened,
 *  and would go stale the moment a library was moved.
 */
import { useEffect, useRef, useState } from 'react'

import { inspectPath } from '../lib/backend'
import { pathTail } from '../lib/recent'

interface Props {
  /** Candidate paths, most recent first. */
  paths: string[]
  onPick: (path: string) => void
}

export function RecentLibraries({ paths, onPick }: Props) {
  const [open, setOpen] = useState(false)
  const [live, setLive] = useState<string[] | null>(null)
  const box = useRef<HTMLDivElement>(null)

  useEffect(() => {
    if (!open) return
    let cancelled = false
    void (async () => {
      const checked = await Promise.all(
        paths.map(async (p) => {
          const info = await inspectPath(p).catch(() => null)
          // Still there, and still a library: a folder that has been emptied or
          // replaced is as useless as one that has been deleted.
          return info && info.exists && info.is_pion_library ? p : null
        }),
      )
      if (!cancelled) setLive(checked.filter((p): p is string => p !== null))
    })()
    return () => {
      cancelled = true
    }
  }, [open, paths])

  // Close on an outside click, as a menu should.
  useEffect(() => {
    if (!open) return
    const onDown = (e: MouseEvent) => {
      if (!box.current?.contains(e.target as Node)) setOpen(false)
    }
    document.addEventListener('mousedown', onDown)
    return () => document.removeEventListener('mousedown', onDown)
  }, [open])

  if (!paths.length) return null

  return (
    <div ref={box} style={{ position: 'relative', flex: 'none' }}>
      <button
        type="button"
        className="pio-browse"
        onClick={() => setOpen((v) => !v)}
        title="Libraries you have searched with before"
        style={{
          flex: 'none',
          padding: '0 12px',
          height: '100%',
          border: '1px solid #E3E8EC',
          borderRadius: 9,
          background: '#FBFCFD',
          font: "500 12.5px 'IBM Plex Sans'",
          color: '#667085',
          cursor: 'pointer',
          whiteSpace: 'nowrap',
        }}
      >
        Recent
      </button>
      {open && (
        <div
          style={{
            position: 'absolute',
            top: 'calc(100% + 6px)',
            right: 0,
            zIndex: 30,
            minWidth: 320,
            maxWidth: 520,
            maxHeight: 260,
            overflowY: 'auto',
            background: '#fff',
            border: '1px solid #E3E8EC',
            borderRadius: 10,
            boxShadow: '0 8px 24px rgba(15,20,27,0.14)',
            padding: 4,
          }}
        >
          {live === null ? (
            <div style={{ padding: '10px 12px', fontSize: 12, color: '#98A2B3' }}>
              Checking…
            </div>
          ) : live.length === 0 ? (
            // Distinguished from having no history at all: the button is not
            // rendered in that case, so reaching here means they have all moved.
            <div style={{ padding: '10px 12px', fontSize: 12, color: '#98A2B3', lineHeight: 1.5 }}>
              None of the libraries you have used are still at their old paths.
            </div>
          ) : (
            live.map((p) => {
              const { name, parent } = pathTail(p)
              return (
                <button
                  key={p}
                  type="button"
                  className="pio-row-hover"
                  onClick={() => {
                    onPick(p)
                    setOpen(false)
                  }}
                  style={{
                    display: 'block',
                    width: '100%',
                    textAlign: 'left',
                    padding: '7px 10px',
                    border: 'none',
                    borderRadius: 7,
                    background: 'none',
                    cursor: 'pointer',
                  }}
                >
                  <div style={{ font: "600 12.5px 'IBM Plex Sans'", color: '#1B2A4A' }}>
                    {name}
                  </div>
                  {parent && (
                    <div
                      style={{
                        font: "11px 'IBM Plex Mono'",
                        color: '#98A2B3',
                        overflow: 'hidden',
                        textOverflow: 'ellipsis',
                        whiteSpace: 'nowrap',
                        direction: 'rtl',
                        textAlign: 'left',
                      }}
                    >
                      {parent}
                    </div>
                  )}
                </button>
              )
            })
          )}
        </div>
      )}
    </div>
  )
}
