/** Human-readable names for runs.
 *
 *  Every job used to be titled after its command, so a queue of four searches
 *  read "SearchDIA / SearchDIA / SearchDIA / SearchDIA" and the only way to tell
 *  them apart was the output path. An adjective-noun pair is short enough for a
 *  narrow sidebar, easy to say out loud ("did brisk-otter finish?"), and needs no
 *  model to generate.
 *
 *  40 x 40 gives 1,600 pairs, so collisions are unlikely but not impossible;
 *  `generateRunName` takes the names already in use and retries.
 */

const ADJECTIVES = [
  'amber', 'arctic', 'bold', 'brisk', 'bronze', 'calm', 'clever', 'cobalt',
  'copper', 'crisp', 'dapper', 'deft', 'eager', 'fleet', 'gentle', 'gilded',
  'glacial', 'golden', 'harbor', 'ivory', 'jade', 'keen', 'lively', 'lucid',
  'mellow', 'nimble', 'noble', 'polar', 'quiet', 'rapid', 'rustic', 'sage',
  'scarlet', 'silent', 'sleek', 'solar', 'stellar', 'sturdy', 'sunny', 'swift',
]

const NOUNS = [
  'acorn', 'anchor', 'badger', 'beacon', 'bison', 'canyon', 'cedar', 'comet',
  'coral', 'crane', 'delta', 'ember', 'falcon', 'fjord', 'garnet', 'harbor',
  'heron', 'ibex', 'juniper', 'kestrel', 'lantern', 'lynx', 'maple', 'meadow',
  'nectar', 'onyx', 'osprey', 'otter', 'pebble', 'quartz', 'raven', 'ridge',
  'sable', 'summit', 'tamarind', 'thistle', 'vireo', 'walnut', 'willow', 'zephyr',
]

const pick = <T,>(xs: readonly T[]): T => xs[Math.floor(Math.random() * xs.length)]

/**
 * A name not already in `taken`.
 *
 * Gives up after a bounded number of attempts and appends a counter rather than
 * looping forever — with 1,600 pairs that only matters for absurd queue lengths,
 * but an unbounded retry on a full namespace would hang the UI.
 */
export function generateRunName(taken: Iterable<string> = []): string {
  const used = new Set(taken)
  for (let i = 0; i < 40; i++) {
    const name = `${pick(ADJECTIVES)}-${pick(NOUNS)}`
    if (!used.has(name)) return name
  }
  let n = 2
  const base = `${pick(ADJECTIVES)}-${pick(NOUNS)}`
  while (used.has(`${base}-${n}`)) n++
  return `${base}-${n}`
}

/** The name a run should carry.
 *
 *  Empty input keeps the generated adjective-noun name. A supplied name is used
 *  as typed unless it is already taken, in which case it gains the next free
 *  suffix -- the bare name counts as the first, so a second `analysis` becomes
 *  `analysis-2`. The suffix style is deliberately the one `generateRunName`
 *  already falls back to, rather than a second convention.
 *
 *  `taken` spans persisted history, not just this session, because the run list
 *  does.
 */
export function resolveRunName(desired: string, taken: Iterable<string> = []): string {
  const name = desired.trim()
  if (!name) return generateRunName(taken)
  const used = new Set(taken)
  if (!used.has(name)) return name
  let n = 2
  while (used.has(`${name}-${n}`)) n++
  return `${name}-${n}`
}
