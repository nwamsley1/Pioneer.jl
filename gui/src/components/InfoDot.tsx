/** A hoverable ⓘ beside a label, for a field whose name cannot carry the whole
 *  explanation.
 *
 *  This used the native `title`, whose timing belongs to the operating system:
 *  around a second on macOS, which is a long wait for text you are deliberately
 *  reaching for, and no way to tune it. The tooltip is drawn here instead so
 *  the delay is ours — long enough that a pointer crossing the dot on its way
 *  somewhere else does not summon it, short enough that aiming at it feels
 *  answered.
 *
 *  `aria-label` stays, so the text is still reachable without a pointer.
 */
import { useEffect, useRef, useState } from 'react'

/** How long the pointer must rest before the description appears. */
const OPEN_DELAY_MS = 350

export function InfoDot({ text, tone = 'light' }: { text: string; tone?: 'light' | 'dark' }) {
  const [at, setAt] = useState<{ left: number; top: number } | null>(null)
  const ref = useRef<HTMLSpanElement | null>(null)
  const timer = useRef<number | null>(null)

  const cancel = () => {
    if (timer.current !== null) {
      window.clearTimeout(timer.current)
      timer.current = null
    }
  }

  const open = () => {
    cancel()
    timer.current = window.setTimeout(() => {
      const r = ref.current?.getBoundingClientRect()
      if (r) setAt({ left: r.left + r.width / 2, top: r.bottom + 8 })
    }, OPEN_DELAY_MS)
  }

  const close = () => {
    cancel()
    setAt(null)
  }

  // A dot can be unmounted while its timer is pending — switching tabs mid-hover
  // does exactly that — and the callback would then set state on a dead
  // component.
  useEffect(() => cancel, [])

  // Scrolling moves the dot out from under a tooltip positioned in viewport
  // coordinates, so the tooltip goes rather than drifting away from its anchor.
  useEffect(() => {
    if (!at) return
    const onScroll = () => close()
    window.addEventListener('scroll', onScroll, true)
    return () => window.removeEventListener('scroll', onScroll, true)
  }, [at])

  return (
    <>
      <span
        ref={ref}
        aria-label={text}
        tabIndex={0}
        onMouseEnter={open}
        onMouseLeave={close}
        onFocus={open}
        onBlur={close}
        onKeyDown={(e) => e.key === 'Escape' && close()}
        style={{
          display: 'inline-flex',
          alignItems: 'center',
          justifyContent: 'center',
          width: 14,
          height: 14,
          flex: 'none',
          borderRadius: '50%',
          // The sidebar is dark; the forms are not. Two fixed palettes rather
          // than currentColor, so neither has to be styled by its parent.
          border: `1px solid ${tone === 'dark' ? 'var(--pio-nav-hair-strong)' : '#B6BFC9'}`,
          color: tone === 'dark' ? 'var(--pio-nav-fg-faint)' : '#8A93A0',
          font: "italic 700 9.5px 'IBM Plex Serif', Georgia, serif",
          // The `font` shorthand does not cover these two, so without them the
          // dot inherits the sidebar heading's uppercase and tracking and renders
          // as a spaced capital I rather than the lowercase i the forms show.
          textTransform: 'none',
          letterSpacing: 'normal',
          cursor: 'help',
          userSelect: 'none',
        }}
      >
        i
      </span>
      {at && (
        <span
          role="tooltip"
          style={{
            // Fixed, so a card or the sidebar's own scroll box cannot clip it.
            position: 'fixed',
            left: at.left,
            top: at.top,
            transform: 'translateX(-50%)',
            // Clamped rather than centred blindly: a dot near the window edge
            // would otherwise place half the text off-screen.
            maxWidth: 'min(340px, calc(100vw - 24px))',
            zIndex: 60,
            padding: '8px 10px',
            borderRadius: 8,
            background: '#1B2A4A',
            color: '#F2F5F9',
            font: "12px/1.45 'IBM Plex Sans'",
            boxShadow: '0 6px 20px rgba(15,20,27,0.28)',
            pointerEvents: 'none',
            whiteSpace: 'normal',
          }}
        >
          {text}
        </span>
      )}
    </>
  )
}
