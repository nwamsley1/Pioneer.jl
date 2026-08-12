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

/** Space reserved at the top of the window for the macOS traffic lights.
 *
 *  The window sets titleBarStyle "Overlay", so on macOS the webview starts at
 *  y=0 and the lights float over our own chrome. That setting is macOS-only —
 *  Windows and Linux keep an ordinary title bar above the webview, and
 *  reserving the strip there would just be dead space under it.
 */
export const TITLEBAR_H = navigator.userAgent.includes('Mac OS X') ? 28 : 0
