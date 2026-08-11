/** Chrome colour themes.
 *
 *  The palettes live in global.css as custom properties; this only decides which
 *  one is active. Navy is the default and carries no `data-theme` attribute, so
 *  a failed restore lands on the original look rather than an unstyled one.
 */

export type ThemeId = 'navy' | 'brown' | 'olive' | 'charcoal'

export interface ThemeDef {
  id: ThemeId
  label: string
  /** The swatch colour: the surface the theme is named after. */
  swatch: string
  /** Accent, so the swatch can hint at what the Run button will look like. */
  accent: string
}

export const THEMES: ThemeDef[] = [
  { id: 'navy', label: 'Navy', swatch: '#1B2A4A', accent: '#3E62A0' },
  { id: 'brown', label: 'Brown', swatch: '#2B2019', accent: '#A8703F' },
  { id: 'olive', label: 'Olive', swatch: '#262B1C', accent: '#74853F' },
  { id: 'charcoal', label: 'Charcoal', swatch: '#23262B', accent: '#52626F' },
]

const THEME_KEY = 'pioneerConsole.theme'

export function loadTheme(): ThemeId {
  try {
    const v = localStorage.getItem(THEME_KEY)
    if (v && THEMES.some((t) => t.id === v)) return v as ThemeId
  } catch {
    /* private mode — fall through to the default */
  }
  return 'navy'
}

export function applyTheme(id: ThemeId): void {
  const el = document.documentElement
  // Navy is the :root default, so it is expressed by the absence of the
  // attribute rather than by a rule of its own.
  if (id === 'navy') el.removeAttribute('data-theme')
  else el.setAttribute('data-theme', id)
  try {
    localStorage.setItem(THEME_KEY, id)
  } catch {
    /* the theme just will not persist */
  }
}
