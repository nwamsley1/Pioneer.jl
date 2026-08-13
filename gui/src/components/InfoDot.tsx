/** A hoverable ⓘ beside a label, for a field whose name cannot carry the whole
 *  explanation.
 *
 *  Uses the native `title` rather than a custom popover: this is supplementary
 *  detail rather than something needed to use the field, and a title is
 *  keyboard- and screen-reader-reachable for free.
 */
export function InfoDot({ text, tone = 'light' }: { text: string; tone?: 'light' | 'dark' }) {
  return (
    <span
      title={text}
      aria-label={text}
      tabIndex={0}
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
  )
}
