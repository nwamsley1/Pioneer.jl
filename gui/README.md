# Pioneer Console

A desktop front end for the packaged Pioneer CLI, implementing the
**"Pioneer GUI - 1 Console"** design from the Claude Design project
`cd3be228-6186-470a-a0b9-7cff100d78e4`.

Tauri v2 (Rust) + React 18 + TypeScript.

## How this ships

The GUI lives in `gui/` inside Pioneer.jl so that one tag builds both halves and
the front end can never drift from the binaries it drives.

It does **not** bundle a copy of the Pioneer distribution. Installer builds
embed a stable version key in the GUI, so each GUI resolves its matching CLI
even when multiple versions are installed:

| platform | CLI installed to | GUI installed to |
|---|---|---|
| macOS `.pkg` | `/usr/local/lib/pioneer/<version>` | `/Applications/Pioneer <version>.app` |
| Linux `.deb` | `/opt/pioneer/<version>` | version-labelled desktop entry |
| Windows setup | `%ProgramFiles%\Pioneer\<version>` | version-labelled Start Menu and optional desktop shortcuts |

Bundling would have duplicated ~600 MB inside a payload that already contains it.

`resolve_home` in `src-tauri/src/pioneer.rs` tries, in order: `PIONEER_HOME`, the
versioned platform install location compiled into installer builds, then
`<resources>/pioneer` (only present if someone deliberately builds a
self-contained bundle). If none match, the error names every path it tried.
An unversioned macOS development build follows
`/usr/local/lib/pioneer/current`, the same active installation selected by the
unqualified `/usr/local/bin/pioneer` command.

### Coverage

No codecov change is needed. `tests.yml` runs `julia-actions/julia-processcoverage`
and uploads only the resulting `lcov.info`, which is built from Julia `.cov`
files under `src/` — TypeScript and Rust produce nothing it reads, so `gui/` is
invisible to coverage automatically.

## Running it

```sh
npm install
npm run tauri:dev
```

For a release bundle:

```sh
npm run tauri:build
```

## Finding the Pioneer binaries

`src-tauri/src/pioneer.rs` resolves the distribution root in this order:

1. **`PIONEER_HOME`** — an unpacked distribution. This is the dev loop, and the
   escape hatch when the bundled copy is wrong.
2. **The installed distribution** — a release GUI uses its embedded version
   key; an unversioned macOS development GUI uses the active `current` link.
3. **`<app resources>/pioneer`** — used only by a deliberately self-contained
   app bundle.

```sh
export PIONEER_HOME=/path/to/Pioneer
npm run tauri:dev
```

The resolved root and its source are shown under the form, so a wrong pickup is
visible rather than mysterious.

### Three traps this app has to handle

All three are properties of the packaged distribution, not of this app:

- **`convert-raw` and `convert-mzml` are different programs.** Only
  `PioneerConverter` reads Thermo `.raw`.
- **`bin/SearchDIA` does not set `JULIA_NUM_THREADS`.** Calling the executable
  directly leaves it single-threaded (~3.7x slower). The `pioneer` wrapper
  defaults it to `auto`, but only when unset (`if [ -z "$JULIA_NUM_THREADS" ]`) —
  it does *not* override an inherited value. This app invokes `bin/SearchDIA`
  directly and sets both variables itself:
  - `JULIA_NUM_THREADS` — the sidebar thread count.
  - `JULIA_NUM_GC_THREADS` — `(threads + 1) / 2` then `,1`, matching the
    wrapper's own formula. Note the name: `JULIA_GC_THREADS` is **not** a Julia
    variable and is silently ignored, leaving GC threads at the full worker
    count. Verified with `Threads.ngcthreads()`.
- **ConvertRAW is not a Julia program and has no params file.** It is a separate
  .NET 8 binary taking a positional RAW path plus flags
  (`-o/--output-dir`, `--skip-existing`, `-n/--concurrent-files`,
  `-t/--threads-per-file`, `-b/--batch-size`, `--scan-chunk-size`). The design
  invented a `{input:{mode,path},output:{path}}` JSON for it; no such file
  exists. `runner::Invocation` models the two shapes, and the Julia thread
  settings are skipped for it.
- **A spectral library is a directory, not a file — and the extension is
  `.poin`.** It holds `precursors_table.arrow`, `detailed_fragments.jls`, the
  fragment index tables and so on. Browse therefore opens a *folder* picker, and
  validation checks for those marker files: an empty folder named `lib.poin`
  passes an extension test and then fails deep inside Pioneer with
  `Fragment file not found`. The design's `.pion` placeholder was a typo —
  there are zero hits for `.pion` in the Pioneer source and `BuildSpecLib`
  appends `.poin`. Libraries named `.pion` do exist in the wild, so the marker
  files decide and an unexpected extension is only a warning.

## FASTA header presets are verified against real files

The design shipped four header-format presets. Only **UniProt** was correct —
its patterns are byte-identical to the ones Pioneer's own `GetBuildLibParams`
emits. The Ensembl and RefSeq patterns were written against invented headers and
matched nothing (or worse, matched wrongly). Measured over 100 real headers each:

| preset | field | design | corrected |
|---|---|---|---|
| Ensembl | accessions | 100/100 | 100/100 |
| Ensembl | genes | 100/100 | 100/100 |
| Ensembl | proteins | **0/100** — no ` protein:` field exists | 100/100 |
| Ensembl | organisms | **0/100** — no `taxon:`/`species:` field | *(empty — absent from the format)* |
| RefSeq | accessions | 100/100 | 100/100 |
| RefSeq | genes | **0/100** — no `GN=` field | *(empty — absent from the format)* |
| RefSeq | proteins | **wrong**: `' (.*?) '` captured only `"thr"` of `"thr operon leader peptide"` | 100/100 |
| RefSeq | organisms | **0/100** — no `OS=` field | 100/100 |

The RefSeq `proteins` case was the dangerous one: it *matched*, so nothing looked
broken, and silently wrote a one-word protein name into the library.

An empty pattern is deliberate, not a placeholder — `chronologer_prep.jl` maps
empty strings to `nothing` ("this format has no such field"), which is better
than a wrong capture. Verified by building real libraries with empty patterns.

Verified against Ensembl release-current peptide FASTAs (yeast R64-1-1, human
GRCh38) and an NCBI RefSeq protein FASTA (E. coli K-12 ASM584v2), then confirmed
end-to-end by reading `proteins_table.arrow` out of the built libraries.

One accepted limitation: for Ensembl, `genes` captures the stable gene ID
(`ENSG…`) rather than the friendlier `gene_symbol:`. `gene_symbol:` is absent in
some species (yeast) and appears *after* `gene:` in the header, so a single
pattern cannot prefer it.

## One field added beyond the design

**Calibration file** (`calibration_raw_file`). The design has no field for it,
but Pioneer's own *simplified* template carries it as "optional but recommended"
and every build without one emits two warnings. It is not cosmetic — building the
same yeast library with and without it:

| | fragment m/z bounds | precursor m/z range |
|---|---|---|
| without | fixed default 150–2020 | fixed default 390–1010 |
| with | fitted `FragBoundModel` polynomial | measured 495.5–601.5 |

So the form has the field, and leaving it blank shows a standing amber warning
naming the defaults you will get instead. A path that is set but wrong blocks the
run; absent does not.

## PioneerConverter is grafted in, not built

It never passes through the Julia build. `.github/workflows/build_app_*.yml`
downloads the latest release of the separate `nwamsley1/PioneerConverter` repo,
picks the asset matching the platform (`osx-x64`, `osx-arm64`, …), and copies
`bin/*` and `lib/*` into the distribution. To reproduce locally:

```sh
unzip PioneerConverter-osx-x64-<ver>.zip -d converter
cp -r converter/PioneerConverter-osx-x64/bin/* "$PIONEER_HOME/bin/"
cp -r converter/PioneerConverter-osx-x64/lib/* "$PIONEER_HOME/lib/"
```

It parallelises on two nested levels, so the knobs multiply
(`Parallel.ForEach` over files, each file split across `ScanReaderWorker`s from
Thermo's `RawFileReaderFactory.CreateThreadManager`):

| flag | default | effect |
|---|---|---|
| `--concurrent-files` | 2 | files converted at once |
| `--threads-per-file` | 3 | scan-reader threads within each file (1 = sequential, no thread manager) |

Defaults give ~6 threads total. The form shows the product live, and the Julia
thread picker is hidden on this page because it does nothing here.

Because the two knobs *multiply*, individually reasonable values can still
exceed the machine — so unlike the Julia picker (clamped at its own control),
this is enforced at run time: `convertThreadsNote` blocks the run when
`concurrent-files × threads-per-file` exceeds the logical core count. The form
and the validator share `convertTotalThreads`, so the readout and the block
cannot disagree.

## The log drawer emulates a terminal

Julia's progress bars are terminal-driven, so the log pane has to understand a
little terminal control or a single bar becomes dozens of lines. Measured from
a real SearchDIA run:

- **stderr** carries progress bars *only* — 54 writes, 54 `\n`, 54 lone `\r`,
  and 41 `ESC[1A` (cursor up). ProgressBars.jl commits each redraw with a
  newline, then moves the cursor back over it.
- **stdout** carries `[ Info: …` and `┌ Warning: …` — no escapes, no carriage
  returns.

`LineSplitter` in `src-tauri/src/runner.rs` interprets `\r`, `CRLF` and
`ESC[<n>A`, discards other CSI codes, and tells the frontend how many committed
lines a segment replaces. `applyLine` in `src/App.tsx` applies that, restricted
to the stream the control came from so interleaved output is never eaten. On the
run above this collapses 54 events to 13 lines — one finished bar per stage.

`cargo test --lib` covers the splitter.

## What is implemented

| Area | State |
|---|---|
| SearchDIA form (Essentials / Confidence & output / Advanced) | done |
| Real path validation against the filesystem | done |
| Native folder & file pickers | done |
| Params JSON view + hand-edit, preserving unknown keys | done |
| Load a previous `config.json` | done |
| Job queue, one run at a time, with live streamed output | done |
| Cancel a run (kills the process group) | done |
| BuildSpecLib (FASTA rows, digestion, mods, options) | done |
| ConvertRAW (file/folder, output, parallelism, advanced) | done |
| Viewer tabs ("Analysis" section) | not yet — separate design file |

## Notes on the port

The design file is a `.dc.html`: a `<x-dc>` template compiled to React by
`support.js`, with the logic in a `class Component extends DCLogic`. Porting it
meant translating `{{ expr }}` → `{expr}`, `<sc-for>` → `.map()`, `<sc-if>` →
`&&`, inline style strings → style objects, and `style-hover=` attributes → real
CSS in `src/global.css`.

Two things in the design deliberately did **not** carry over:

- **The stepper.** The Component class holds `currentStep`, `buildSteps()`,
  `next`/`back`/`goStep` and `stepValid()` — inherited from the Guided variant —
  but the Console template never renders any of it. SearchDIA is one page.
- **The simulated backend.** The prototype advanced a job every 620 ms through a
  hardcoded array of fake log lines, filled paths from `data-sample` attributes,
  compared results folders against a hardcoded list of three directories, and
  failed any run whose paths matched `/fail/i`. All of that is replaced by real
  process execution, real dialogs and real `stat` calls.

Fonts are bundled via `@fontsource` rather than fetched from Google Fonts, so the
app renders identically offline and under the app's CSP.
