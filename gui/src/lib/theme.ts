/** Chrome colour themes.
 *
 *  The palettes live in global.css as custom properties; this only decides which
 *  one is active. Charcoal is the default and carries no `data-theme` attribute,
 *  so a failed restore lands on it rather than on an unstyled shell.
 */

export type ThemeId = 'charcoal' | 'navy' | 'burgundy' | 'brown' | 'sand' | 'olive'

export interface ThemeDef {
  id: ThemeId
  label: string
  /** The swatch colour: the surface the theme is named after. */
  swatch: string
  /** The Run button's colour, so the swatch previews it. Not always the
   *  same as --pio-accent: the light palette needs a lighter button than the
   *  accent it uses for link text. */
  accent: string
}

export const THEMES: ThemeDef[] = [
  { id: 'charcoal', label: 'Charcoal', swatch: '#23262B', accent: '#52626F' },
  { id: 'navy', label: 'Navy', swatch: '#1B2A4A', accent: '#3E62A0' },
  { id: 'burgundy', label: 'Burgundy', swatch: '#2E1A20', accent: '#A8485E' },
  { id: 'brown', label: 'Brown', swatch: '#2B2019', accent: '#A8703F' },
  { id: 'sand', label: 'Sand', swatch: '#CDBE9C', accent: '#B08F4E' },
  { id: 'olive', label: 'Olive', swatch: '#262B1C', accent: '#74853F' },
]

const THEME_KEY = 'pioneerConsole.theme'

export function loadTheme(): ThemeId {
  try {
    const v = localStorage.getItem(THEME_KEY)
    if (v && THEMES.some((t) => t.id === v)) return v as ThemeId
  } catch {
    /* private mode — fall through to the default */
  }
  return 'charcoal'
}

export function applyTheme(id: ThemeId): void {
  const el = document.documentElement
  // Charcoal is the :root default, so it is expressed by the absence of the
  // attribute rather than by a rule of its own.
  if (id === 'charcoal') el.removeAttribute('data-theme')
  else el.setAttribute('data-theme', id)
  try {
    localStorage.setItem(THEME_KEY, id)
  } catch {
    /* the theme just will not persist */
  }
}
