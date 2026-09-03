/** Spectral libraries this install has searched with before.
 *
 *  Used, built or downloaded — all three are libraries you might reach for
 *  next, and a freshly built one would otherwise be missing precisely when it
 *  is most wanted.
 *
 *  No new storage: every submitted run already persists its whole parameter set
 *  in the history store's `snapshot` column, so each row carries the library it
 *  used or produced. `history_load` returns every row, so the paths are already
 *  in memory — including from before this feature existed.
 */
import type { Job } from './types'
import { downloadTargetPath, libraryTargetPath } from './validate'

/** The library a run either used or produced, if any. */
export function libraryOf(job: Job): string {
  switch (job.snapshot.cmd) {
    case 'searchdia':
      return job.snapshot.search.library.trim()
    // Built and downloaded libraries count too: having just made one, the next
    // thing anyone does is search with it, and it will not be in the list on
    // the strength of a search that has not happened yet.
    case 'buildspeclib':
      return libraryTargetPath(job.snapshot.build.libPath)
    case 'downloadspeclib':
      return downloadTargetPath(job.snapshot.download.dest, job.snapshot.download.selected)
    default:
      return ''
  }
}

/**
 * Distinct library paths, most recently used or produced first.
 *
 * Every run counts, whatever became of it: a search that died in scoring still
 * used a real library, and one that was cancelled says nothing about the path.
 * A build that failed is the one exception worth noting — it may never have
 * written the library — but the existence check the caller runs covers that
 * without needing to reason about it here.
 *
 * A path that has since been deleted is not filtered here: that needs the
 * filesystem, and the caller does it when the list is shown, so it cannot go
 * stale in between.
 *
 * Ordered by run number rather than `finishedAt`, which is zero for runs
 * restored from before finish times were recorded.
 */
export function recentLibraries(jobs: Job[]): string[] {
  const byRecency = jobs.slice().sort((a, b) => b.runNo - a.runNo)
  const seen = new Set<string>()
  const out: string[] = []
  for (const j of byRecency) {
    const path = libraryOf(j)
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
