/** Dropping a file or folder onto the field that wants it.
 *
 *  Not HTML5 drag-and-drop. Tauri's `dragDropEnabled` defaults to true, which
 *  means the OS drop is handled natively and `ondrop` never fires inside the
 *  webview — so this listens to Tauri's own event instead. The trade is that
 *  the event is window-level and carries a position rather than a target, so
 *  the field under the cursor has to be found by hit test.
 */
import { getCurrentWebview } from '@tauri-apps/api/webview'
import type { UnlistenFn } from '@tauri-apps/api/event'

/** What each droppable field accepts. Keyed by the `data-key` the field
 *  already carries for validation focus, so nothing new is threaded through
 *  the forms. */
export const DROP_FIELDS: Record<string, 'dir' | 'file' | 'either'> = {
  // A .poin library is a directory of tables, not a file.
  msData: 'dir',
  library: 'dir',
  results: 'dir',
  libPath: 'dir',
  convertOutput: 'dir',
  calibrationFile: 'file',
  // ConvertRAW takes a single .raw or a folder of them.
  convertInput: 'either',
  // The FASTA row's add button: dropping onto it appends rows.
  fastaAdd: 'file',
}

const HILITE = 'pio-droptarget'

/** The droppable element under a window position, or null.
 *
 *  The event reports physical pixels; `elementFromPoint` wants CSS pixels, so
 *  the ratio has to be divided out or the hit test lands somewhere else
 *  entirely on a Retina display.
 */
function targetAt(x: number, y: number): HTMLElement | null {
  const r = window.devicePixelRatio || 1
  const el = document.elementFromPoint(x / r, y / r)
  if (!el) return null
  // A row first — it covers the label, input, Browse button and padding. Then
  // data-key, for targets that are a single control (the FASTA add button).
  const row = el.closest('[data-drop]') as HTMLElement | null
  if (row) return row
  const hit = el.closest('[data-key]') as HTMLElement | null
  const key = hit?.getAttribute('data-key') ?? ''
  return key in DROP_FIELDS ? hit : null
}

function clearHighlight(): void {
  document.querySelectorAll(`.${HILITE}`).forEach((e) => e.classList.remove(HILITE))
}

export interface DropHandler {
  /** Called with the field's data-key and the dropped paths. */
  (key: string, kind: 'dir' | 'file' | 'either', paths: string[]): void
}

/** Start listening. Returns a promise for the unlisten function. */
export function listenForDrops(onDrop: DropHandler): Promise<UnlistenFn> {
  return getCurrentWebview().onDragDropEvent((event) => {
    const p = event.payload
    if (p.type === 'leave') {
      clearHighlight()
      return
    }
    if (p.type === 'drop') {
      const el = targetAt(p.position.x, p.position.y)
      clearHighlight()
      if (!el || !p.paths.length) return
      const key = (el.getAttribute('data-drop') ?? el.getAttribute('data-key')) as string
      onDrop(key, DROP_FIELDS[key], p.paths)
      return
    }
    // enter / over: show what would receive it.
    const el = targetAt(p.position.x, p.position.y)
    if (el?.classList.contains(HILITE)) return
    clearHighlight()
    el?.classList.add(HILITE)
  })
}
