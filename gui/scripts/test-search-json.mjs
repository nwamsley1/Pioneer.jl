import assert from 'node:assert/strict'
import { mkdtempSync, rmSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { fileURLToPath } from 'node:url'
import { createRequire } from 'node:module'
import { spawnSync } from 'node:child_process'

// SearchDIA config written by the form, and read back from a saved config. Compiled with the existing
// TypeScript dependency, as test-convert-args.mjs does.
const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..')
const out = mkdtempSync(path.join(tmpdir(), 'pioneer-search-json-'))
try {
  const compile = spawnSync(process.execPath, [
    path.join(root, 'node_modules/typescript/bin/tsc'),
    '--target', 'ES2020', '--module', 'commonjs', '--moduleResolution', 'node',
    '--skipLibCheck', '--esModuleInterop', '--resolveJsonModule',
    '--outDir', out, 'src/lib/config.ts', 'src/lib/types.ts',
  ], { cwd: root, stdio: 'inherit' })
  assert.equal(compile.status, 0, 'search config builder must compile')
  const require = createRequire(path.join(out, 'test.cjs'))
  const { buildSearchJson, searchConfigToState } = require('./config.js')
  const { SEARCH_DEFAULTS } = require('./types.js')

  // MS1 quantification (maxLFQ.ms1_quantification): off by default, always written, read back
  assert.equal(SEARCH_DEFAULTS.ms1Quant, false)
  assert.equal(buildSearchJson(SEARCH_DEFAULTS, null).maxLFQ.ms1_quantification, false)
  const on = buildSearchJson({ ...SEARCH_DEFAULTS, ms1Quant: true }, null)
  assert.equal(on.maxLFQ.ms1_quantification, true)
  assert.equal(on.maxLFQ.run_to_run_normalization, SEARCH_DEFAULTS.runToRunNorm)
  assert.equal(searchConfigToState(on).ms1Quant, true)
  // a config written before the option existed leaves the toggle as it stands
  const old = buildSearchJson(SEARCH_DEFAULTS, null)
  delete old.maxLFQ.ms1_quantification
  assert.equal(searchConfigToState(old).ms1Quant, undefined)

  console.log('SearchDIA config checks passed')
} finally {
  rmSync(out, { recursive: true, force: true })
}
