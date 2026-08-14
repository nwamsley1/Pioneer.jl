/** Spectral libraries this install has searched with before.
 *
 *  No new storage: every submitted run already persists its whole parameter set
 *  in the history store's `snapshot` column, so a SearchDIA row carries the
 *  library it used. `history_load` returns every row, so the paths are already
 *  in memory by the time the sidebar has a queue — including from before this
 *  feature existed.
 */
import type { Job } from './types'

/**
 * Distinct library paths from SearchDIA runs, most recently used first.
 *
 * Every run counts, whatever became of it: a search that died in scoring still
 * used a real library, and one that was cancelled says nothing about the path.
 * A path that has since been deleted is not filtered here — that needs the
 * filesystem, and the caller does it when the list is shown, so it cannot go
 * stale in between.
 *
 * Ordered by run number rather than `finishedAt`, which is zero for runs
 * restored from before finish times were recorded.
 */
export function recentLibraries(jobs: Job[]): string[] {
  const byRecency = jobs
    .filter((j) => j.snapshot.cmd === 'searchdia')
    .slice()
    .sort((a, b) => b.runNo - a.runNo)

  const seen = new Set<string>()
  const out: string[] = []
  for (const j of byRecency) {
    if (j.snapshot.cmd !== 'searchdia') continue
    const path = j.snapshot.search.library.trim()
    if (!path || seen.has(path)) continue
    seen.add(path)
    out.push(path)
  }
  return out
}

/** The last path segment, for showing a library by name with its folder under it. */
export function pathTail(path: string): { name: string; parent: string } {
  const trimmed = path.replace(/[/\\]+$/, '')
  const cut = Math.max(trimmed.lastIndexOf('/'), trimmed.lastIndexOf('\\'))
  return cut === -1
    ? { name: trimmed, parent: '' }
    : { name: trimmed.slice(cut + 1), parent: trimmed.slice(0, cut) }
}
