/** Styles shared across forms.
 *
 *  `LABEL` was defined identically in SearchDiaForm and ConvertRawForm, and
 *  NumField had drifted to its own lighter treatment (12px regular #475467).
 *  That put two different label styles on the same card — the q-value field
 *  read as less important than the toggles beside it, which inverts the real
 *  hierarchy. One definition, so they cannot diverge again.
 */
import type { CSSProperties } from 'react'

/** The name of a setting: a form field's label, or a toggle's title. */
export const LABEL: CSSProperties = {
  display: 'block',
  fontSize: 12,
  fontWeight: 600,
  color: '#344054',
  marginBottom: 6,
}

/** `LABEL` for a label whose container already provides the spacing. */
export const LABEL_TIGHT: CSSProperties = { ...LABEL, marginBottom: 0 }

/** The explanatory line under a setting's name. */
export const HINT: CSSProperties = {
  fontSize: 11.5,
  color: '#98A2B3',
}

/** The Browse button beside a path input.
 *
 *  Deliberately quieter than the input it sits next to. The input is the
 *  control you act on; Browse is the fallback for when typing a path is
 *  inconvenient. It previously carried a heavier border and semibold text than
 *  the field, so it won the row it was meant to support.
 *
 *  Was defined identically in SearchDiaForm and ConvertRawForm and inlined
 *  three more times in BuildSpecLibForm — five copies, already drifting on
 *  radius and font size.
 */
export const BROWSE: CSSProperties = {
  flex: 'none',
  padding: '0 14px',
  border: '1px solid #E3E8EC',
  borderRadius: 9,
  background: '#FBFCFD',
  font: "500 12.5px 'IBM Plex Sans'",
  color: '#667085',
  cursor: 'pointer',
}

/** One button of a segmented control -- the pill pair that picks between two
 *  shapes of the same field (single file vs folder, .raw vs .mzML, one folder
 *  vs chosen files). Was defined inside ConvertRawForm; the SearchDIA file-list
 *  toggle is the same control, so it lives here rather than being copied.
 *
 *  Wrap the buttons in a container with `padding: 3`, a `#EEF1F4` background
 *  and `borderRadius: 10` -- the raised look comes from the active button
 *  sitting on that inset track.
 */
export const seg = (active: boolean): CSSProperties => ({
  padding: '7px 14px',
  border: 'none',
  borderRadius: 8,
  font: "600 12.5px 'IBM Plex Sans'",
  cursor: 'pointer',
  transition: 'all .12s',
  ...(active
    ? { background: '#fff', color: '#1B2A4A', boxShadow: '0 1px 3px rgba(15,20,27,0.12)' }
    : { background: 'none', color: '#667085' }),
})

/** The track a `seg` pair sits in. */
export const SEG_TRACK: CSSProperties = {
  display: 'inline-flex',
  padding: 3,
  background: '#EEF1F4',
  borderRadius: 10,
}

/** Space reserved at the top of the window for the macOS traffic lights.
 *
 *  The window sets titleBarStyle "Overlay", so on macOS the webview starts at
 *  y=0 and the lights float over our own chrome. That setting is macOS-only —
 *  Windows and Linux keep an ordinary title bar above the webview, and
 *  reserving the strip there would just be dead space under it.
 */
export const TITLEBAR_H = navigator.userAgent.includes('Mac OS X') ? 28 : 0
