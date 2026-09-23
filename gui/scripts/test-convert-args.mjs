import assert from 'node:assert/strict'
import { mkdtempSync, rmSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { fileURLToPath } from 'node:url'
import { createRequire } from 'node:module'
import { spawnSync } from 'node:child_process'

// Compile the production builder with the existing TypeScript dependency.
// No test runner or additional package is required.
const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..')
const out = mkdtempSync(path.join(tmpdir(), 'pioneer-convert-args-'))
try {
  const compile = spawnSync(process.execPath, [
    path.join(root, 'node_modules/typescript/bin/tsc'),
    '--target', 'ES2020', '--module', 'commonjs', '--moduleResolution', 'node',
    '--skipLibCheck', '--esModuleInterop', '--resolveJsonModule',
    '--outDir', out, 'src/lib/config.ts', 'src/lib/types.ts',
  ], { cwd: root, stdio: 'inherit' })
  assert.equal(compile.status, 0, 'argument builder must compile')
  const require = createRequire(path.join(out, 'test.cjs'))
  const { buildConvertArgs, convertCommandLine } = require('./config.js')
  const { CONVERT_DEFAULTS } = require('./types.js')
  assert.equal(CONVERT_DEFAULTS.batchSize, '1000')
  const raw = { ...CONVERT_DEFAULTS, input: ' /data/input.raw ',
    outputDir: ' /ssd/output ', skipExisting: true, threadsPerFile: '4',
    batchSize: '50', scanChunkSize: '17', concurrentFiles: '8' }
  assert.deepEqual(buildConvertArgs(raw), [
    '/data/input.raw', '--output-dir', '/ssd/output', '--skip-existing',
    '--threads-per-file', '4', '--batch-size', '50', '--scan-chunk-size', '17',
  ])
  for (const inputMode of ['files', 'folder']) {
    const args = buildConvertArgs({ ...raw, inputMode })
    assert(!args.includes('--concurrent-files'))
    assert(!args.includes('-n'))
  }
  assert(!convertCommandLine(raw).includes('--concurrent-files'))
  assert.deepEqual(buildConvertArgs({ ...raw, format: 'mzml' }), [
    '/data/input.raw', '--output-dir', '/ssd/output', '--skip-existing',
    '--concurrent-files', '8', '--skip-header',
  ])
  console.log('RAW and mzML converter argument regression checks passed')
} finally {
  rmSync(out, { recursive: true, force: true })
}
