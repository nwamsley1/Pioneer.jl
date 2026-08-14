/** Reordering the run queue.
 *
 *  The scheduler takes the first job in `jobs` whose status is `queued`, so a
 *  job's position in this array is its place in line. Reordering is therefore a
 *  rewrite of the array, not a separate ordering stored beside it.
 */
import type { Job } from './types'

/**
 * Move the queued job `dragId` into the position held by `dropId`.
 *
 * Only queued jobs move. A running job has already started and a finished one
 * has no position left to change, so their slots in the array are held fixed
 * and the queued entries are rewritten into the remaining ones — a running job
 * keeps the row it is displayed at rather than being shuffled by a drag that
 * had nothing to do with it.
 *
 * Run numbers travel with the position, not the job: the numbers the queued
 * jobs hold are redistributed so they still ascend down the queue. Dragging the
 * last of 42, 43, 44 to the front leaves the queue reading 42, 43, 44 rather
 * than 44, 42, 43. A number is therefore a place in line only while a run is
 * waiting; from the moment it starts it is fixed for life, which is what the
 * history depends on.
 *
 * Returns the original array unchanged when either id is not a queued job, or
 * when the move is a no-op, so React can skip the re-render.
 */
export function moveQueuedJob(jobs: Job[], dragId: string, dropId: string): Job[] {
  const slots: number[] = []
  jobs.forEach((j, i) => {
    if (j.status === 'queued') slots.push(i)
  })
  const queued = slots.map((i) => jobs[i])
  const from = queued.findIndex((j) => j.id === dragId)
  const to = queued.findIndex((j) => j.id === dropId)
  if (from === -1 || to === -1 || from === to) return jobs

  // Captured before the move: the same numbers come back out, in the same
  // ascending order, just held by different jobs.
  const numbers = queued.map((j) => j.runNo).sort((a, b) => a - b)

  const [moved] = queued.splice(from, 1)
  queued.splice(to, 0, moved)

  const next = [...jobs]
  queued.forEach((j, k) => {
    if (j.runNo !== numbers[k]) queued[k] = { ...j, runNo: numbers[k] }
  })
  slots.forEach((slot, k) => {
    next[slot] = queued[k]
  })
  return next
}
