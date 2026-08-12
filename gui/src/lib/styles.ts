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
