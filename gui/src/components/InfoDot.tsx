/** A hoverable ⓘ beside a label, for a field whose name cannot carry the whole
 *  explanation.
 *
 *  Uses the native `title` rather than a custom popover: this is supplementary
 *  detail rather than something needed to use the field, and a title is
 *  keyboard- and screen-reader-reachable for free.
 */
export function InfoDot({ text }: { text: string }) {
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
        border: '1px solid #B6BFC9',
        color: '#8A93A0',
        font: "italic 700 9.5px 'IBM Plex Serif', Georgia, serif",
        cursor: 'help',
        userSelect: 'none',
      }}
    >
      i
    </span>
  )
}
