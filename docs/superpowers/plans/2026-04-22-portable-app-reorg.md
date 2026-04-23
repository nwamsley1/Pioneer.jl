# Portable App Reorganization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reorganize Pioneer so portable PackageCompiler apps can keep all safe code in the sysimage, while risky or optional code moves into separate runtime-loadable packages that can persist native caches in a writable user depot.

**Architecture:** Convert the current monolithic `Pioneer` package into a workspace-style monorepo with a thin umbrella compatibility package plus command-family app packages. Build small, purpose-specific app sysimages from those app packages. Keep SIMD and plotting out of the search sysimage by moving them behind explicit package boundaries, then warm their package caches into a stable runtime depot so they land in category 2 instead of recompiling on every launch.

**Tech Stack:** Julia 1.12, PackageCompiler, Pkg/package caches/pkgimages, app-specific `Project.toml` environments, optional runtime packages, GitHub Actions installer workflows.

---

## Current Diagnosis

The current build is structurally hostile to portable app compilation:

- `create_app(".", ...)` builds one sysimage from the entire `Pioneer` package for all executables in [`.github/workflows/build_app_linux.yml`](/Users/dennisgoldfarb/Programming/Pioneer.jl/.github/workflows/build_app_linux.yml:143), [`.github/workflows/build_app_windows.yml`](/Users/dennisgoldfarb/Programming/Pioneer.jl/.github/workflows/build_app_windows.yml:130), and [`.github/workflows/build_app_macos.yml`](/Users/dennisgoldfarb/Programming/Pioneer.jl/.github/workflows/build_app_macos.yml:201).
- [`src/Pioneer.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Pioneer.jl:19) eagerly imports `LoopVectorization`, `Plots`, `StatsPlots`, `GR`, `LaTeXStrings`, and `Measures`.
- [`src/importScripts.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/importScripts.jl:17) recursively includes most of `src/`, so `using Pioneer` drags almost the entire project graph into the app build.
- SIMD code is mixed into common structs and search logic in files such as [`src/structs/SparseArray.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/structs/SparseArray.jl:172), [`src/Routines/SearchDIA/PSMs/spectralDistanceMetrics.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/SearchDIA/PSMs/spectralDistanceMetrics.jl:231), and [`src/utils/ML/probitRegression.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/utils/ML/probitRegression.jl:51).
- Plotting code is mixed into search result types and workflow code in files such as [`src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/types.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/types.jl:179), [`src/Routines/SearchDIA/WriteOutputs/qcPlots.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/SearchDIA/WriteOutputs/qcPlots.jl:72), and [`src/utils/pdfUtils.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/utils/pdfUtils.jl:1).

This means “safe” executables like `GetSearchParams` and `convertMzML` still inherit the portable-build risk of `SearchDIA`, plotting, and host-sensitive SIMD code.

## Recommended End State

Adopt a three-layer structure:

1. **Portable workspace packages that are safe to sysimage**
   - `PioneerCommon`
   - `PioneerPredict`
   - `PioneerConvert`
   - `PioneerSearch`

2. **Runtime-only optional packages that should not be in the portable search sysimage**
   - `PioneerPlots`
   - `PioneerSIMD`

3. **Thin app packages used only for PackageCompiler**
   - `apps/PioneerParamsApp`
   - `apps/PioneerPredictApp`
   - `apps/PioneerConvertApp`
   - `apps/PioneerSearchApp`

Keep a top-level `Pioneer` umbrella package for developer ergonomics, docs, and backward-compatible imports. `Pioneer` should reexport stable APIs from the workspace packages, but the installer workflows should stop compiling the root package directly.

## Why This Architecture

This is the best fit for your stated priorities:

- **Portability first:** `PioneerSIMD` becomes optional runtime code instead of portable sysimage code.
- **Maximize sysimage coverage:** each app package gets its own smaller sysimage that still contains all safe code for that command family.
- **Maximize category 2 coverage:** `PioneerPlots` and `PioneerSIMD` remain bundled as app dependencies, but are loaded at runtime from a stable writable depot so Julia can persist pkgimages and reuse them across launches.
- **Minimize user-facing breakage:** `using Pioneer` can continue to work while the app build path uses purpose-specific packages.

## Packages And Responsibilities

### `packages/PioneerCommon`

**Responsibility:** Safe shared code used by almost every command.

**Move here:**
- Logging and version helpers from [`src/Pioneer.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Pioneer.jl:19)
- Common structs from [`src/structs`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/structs)
- Non-plotting/non-SIMD utilities from [`src/utils`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/utils)
- Parameter parsing primitives and shared config types
- File operations and Arrow/CSV helpers

**Do not depend on:** `LoopVectorization`, `HostCPUFeatures`, `Plots`, `StatsPlots`, `GR`, `Measures`, `LaTeXStrings`

### `packages/PioneerPredict`

**Responsibility:** Spectral-library generation and related parameter generation.

**Move here:**
- [`src/Routines/BuildSpecLib.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/BuildSpecLib.jl:1)
- [`src/Routines/GenerateParams.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/GenerateParams.jl:1) build/parse-spec-lib pieces
- [`src/Routines/BuildSpecLib/**/*`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/BuildSpecLib.jl:1)
- Koina/Chronologer code

**Notes:** `estimate_collision_ev.jl` currently plots. Split the plotting bits into `PioneerPlots` and keep the numeric logic in `PioneerPredict`.

### `packages/PioneerConvert`

**Responsibility:** mzML-to-Arrow conversion and its CLI.

**Move here:**
- [`src/Routines/mzmlConverter/convertMzML.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/mzmlConverter/convertMzML.jl:1)
- Any converter-specific parsing helpers

**Target outcome:** This package should be fully sysimage-safe and independent of search, plotting, and SIMD.

### `packages/PioneerSearch`

**Responsibility:** SearchDIA orchestration, scoring pipeline, tuning logic, and outputs that do not require plotting or host-specific SIMD.

**Move here:**
- [`src/Routines/SearchDIA.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/SearchDIA.jl:1)
- Portable parts of [`src/Routines/SearchDIA/**/*`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/SearchDIA)
- ML and scoring abstractions from [`src/utils/ML`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/utils/ML)

**Key rule:** `PioneerSearch` owns the public search APIs and interfaces, but not plotting implementations or `@turbo` kernels.

### `packages/PioneerPlots`

**Responsibility:** Search/build QC plots and PDF composition.

**Move here:**
- [`src/Routines/SearchDIA/WriteOutputs/qcPlots.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/SearchDIA/WriteOutputs/qcPlots.jl:1)
- [`src/Routines/SearchDIA/WriteOutputs/plotRTAlignment.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/SearchDIA/WriteOutputs/plotRTAlignment.jl:1)
- Plotting helpers embedded in tuning-search `utils.jl` files
- [`src/Routines/SearchDIA/SearchMethods/ProteinScoringSearch/qc_plots.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/SearchDIA/SearchMethods/ProteinScoringSearch/qc_plots.jl:1)
- [`src/utils/pdfUtils.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/utils/pdfUtils.jl:1)

**Design rule:** `PioneerSearch` should store portable diagnostics data, not `Plots.Plot` objects. `PioneerPlots` should translate those diagnostics into plot objects on demand.

### `packages/PioneerSIMD`

**Responsibility:** Optional accelerated kernels and their dispatch layer.

**Move here:**
- `@turbo` kernels from:
  - [`src/structs/ArrayDict.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/structs/ArrayDict.jl:34)
  - [`src/structs/Counter.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/structs/Counter.jl:77)
  - [`src/structs/SparseArray.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/structs/SparseArray.jl:172)
  - [`src/Routines/SearchDIA/PSMs/spectralDistanceMetrics.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/SearchDIA/PSMs/spectralDistanceMetrics.jl:231)
  - [`src/utils/ML/probitRegression.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/utils/ML/probitRegression.jl:51)

**Design rule:** `PioneerSearch` must define portable fallback implementations and call named kernel functions. `PioneerSIMD` adds accelerated methods by overloading those kernel entrypoints.

## App Packaging Strategy

Stop building the root package directly with `create_app(".")`.

Instead, create command-family app packages:

- `apps/PioneerParamsApp` for `GetSearchParams`, `GetBuildLibParams`, `GetParseSpecLibParams`
- `apps/PioneerConvertApp` for `convertMzML`
- `apps/PioneerPredictApp` for `BuildSpecLib`
- `apps/PioneerSearchApp` for `SearchDIA`

Each app package gets its own `Project.toml` and `src/*.jl` with only the executables it owns. This lets PackageCompiler build different sysimages for different command families while the installer still ships a single product with multiple binaries.

### Search App Rules

`PioneerSearchApp` should:

- include `PioneerSearch` in the sysimage
- keep `PioneerPlots` and `PioneerSIMD` in the app project as dependencies
- not import `PioneerPlots` or `PioneerSIMD` during `create_app`
- load those packages at runtime only when needed

This is the part that turns plotting/SIMD into category 2 instead of “compile at every launch.”

## Runtime Cache Strategy

PackageCompiler apps already use an app-local project and depot under `share/julia`, but durable category 2 behavior needs a stable writable depot that survives reinstalls and upgrades.

### Required behavior

- On startup, launcher scripts or wrappers should prepend a user-specific writable depot before the bundled `share/julia` depot.
- Suggested platform paths:
  - Linux: `~/.cache/pioneer/julia`
  - macOS: `~/Library/Caches/Pioneer/julia`
  - Windows: `%LOCALAPPDATA%\\Pioneer\\julia`

### Required warm-up flow

Add a runtime warm-up script, for example:

- `src/build/warm_runtime_cache.jl`

This script should:

- activate the app project
- import `PioneerPlots` and `PioneerSIMD`
- run a small representative workload for each
- call `Pkg.precompile()` or a purpose-built warm-up entrypoint

Use that script in one or both of these ways:

- installer post-install validation
- first-run lazy warm-up when cache files are absent

## Migration Phases

### Phase 1: Freeze Behavior And Add Guardrails

**Files:**
- Create: `test/packaging/test_app_command_smokes.jl`
- Create: `scripts/report_app_build_metrics.jl`
- Modify: [`.github/workflows/build_app_linux.yml`](/Users/dennisgoldfarb/Programming/Pioneer.jl/.github/workflows/build_app_linux.yml:1)
- Modify: [`.github/workflows/build_app_windows.yml`](/Users/dennisgoldfarb/Programming/Pioneer.jl/.github/workflows/build_app_windows.yml:1)
- Modify: [`.github/workflows/build_app_macos.yml`](/Users/dennisgoldfarb/Programming/Pioneer.jl/.github/workflows/build_app_macos.yml:1)

**Outcome:** Measure current app size, build time, and portable-build failures before restructuring.

### Phase 2: Replace `importScripts.jl` With Explicit Package Boundaries

**Files:**
- Modify: [`src/Pioneer.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Pioneer.jl:19)
- Delete eventually: [`src/importScripts.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/importScripts.jl:17)
- Create: `packages/PioneerCommon/Project.toml`
- Create: `packages/PioneerCommon/src/PioneerCommon.jl`
- Create: `packages/PioneerPredict/Project.toml`
- Create: `packages/PioneerPredict/src/PioneerPredict.jl`
- Create: `packages/PioneerConvert/Project.toml`
- Create: `packages/PioneerConvert/src/PioneerConvert.jl`
- Create: `packages/PioneerSearch/Project.toml`
- Create: `packages/PioneerSearch/src/PioneerSearch.jl`

**Outcome:** The build graph becomes explicit and reviewable.

### Phase 3: Extract Plotting Out Of Search Data Structures

**Files:**
- Create: `packages/PioneerPlots/Project.toml`
- Create: `packages/PioneerPlots/src/PioneerPlots.jl`
- Modify: [`src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/types.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/types.jl:166)
- Modify: [`src/Routines/SearchDIA/SearchMethods/NceTuningSearch/NceTuningSearch.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/SearchDIA/SearchMethods/NceTuningSearch/NceTuningSearch.jl:49)
- Modify: [`src/Routines/SearchDIA/SearchMethods/QuadTuningSearch/QuadTuningSearch.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/SearchDIA/SearchMethods/QuadTuningSearch/QuadTuningSearch.jl:54)
- Move plotting writers/helpers listed above into `packages/PioneerPlots`

**Outcome:** `PioneerSearch` results become portable data containers, not plotting objects.

### Phase 4: Extract SIMD Behind Explicit Kernel Interfaces

**Files:**
- Create: `packages/PioneerSIMD/Project.toml`
- Create: `packages/PioneerSIMD/src/PioneerSIMD.jl`
- Modify: [`src/structs/ArrayDict.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/structs/ArrayDict.jl:1)
- Modify: [`src/structs/Counter.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/structs/Counter.jl:1)
- Modify: [`src/structs/SparseArray.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/structs/SparseArray.jl:1)
- Modify: [`src/Routines/SearchDIA/PSMs/spectralDistanceMetrics.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Routines/SearchDIA/PSMs/spectralDistanceMetrics.jl:1)
- Modify: [`src/utils/ML/probitRegression.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/utils/ML/probitRegression.jl:1)

**Outcome:** Search can run portably with fallback kernels, while `PioneerSIMD` becomes an acceleration layer instead of a build requirement.

### Phase 5: Introduce App Packages And Stop Sysimaging The Monolith

**Files:**
- Create: `apps/PioneerParamsApp/Project.toml`
- Create: `apps/PioneerParamsApp/src/PioneerParamsApp.jl`
- Create: `apps/PioneerPredictApp/Project.toml`
- Create: `apps/PioneerPredictApp/src/PioneerPredictApp.jl`
- Create: `apps/PioneerConvertApp/Project.toml`
- Create: `apps/PioneerConvertApp/src/PioneerConvertApp.jl`
- Create: `apps/PioneerSearchApp/Project.toml`
- Create: `apps/PioneerSearchApp/src/PioneerSearchApp.jl`
- Modify installer workflows under `.github/workflows/`

**Outcome:** Each binary gets a right-sized sysimage instead of one shared monolith.

### Phase 6: Add Durable Runtime Cache Warming

**Files:**
- Create: `src/build/warm_runtime_cache.jl`
- Modify launcher/wrapper scripts in installer assets
- Modify search-app runtime startup path
- Add validation coverage in installer smoke tests

**Outcome:** `PioneerPlots` and `PioneerSIMD` become category 2 across launches instead of paying JIT every time.

### Phase 7: Preserve Developer Ergonomics

**Files:**
- Modify: [`src/Pioneer.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/src/Pioneer.jl:19)
- Modify: [`docs/make.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/docs/make.jl:18)
- Modify: [`test/runtests.jl`](/Users/dennisgoldfarb/Programming/Pioneer.jl/test/runtests.jl:19)
- Modify docs under [`docs/src`](/Users/dennisgoldfarb/Programming/Pioneer.jl/docs/src)

**Outcome:** `using Pioneer` still works for docs and most tests while the internal app build uses split packages.

## Verification Criteria

The reorganization is successful only if all of these become true:

- Linux, Windows, and macOS app builds succeed on Julia 1.12 with portable CPU targets.
- `PioneerSearchApp` sysimage does not import `LoopVectorization`, `HostCPUFeatures`, `Plots`, `StatsPlots`, `GR`, `Measures`, or `LaTeXStrings`.
- A first launch after install populates cache files for `PioneerPlots` and `PioneerSIMD` in the user depot.
- A second launch reuses those caches and shows materially lower latency for the same workload.
- Docs and the default developer test path still work with `using Pioneer`.

## Recommended Rollout Order

1. Guardrails/metrics
2. `PioneerCommon`, `PioneerPredict`, `PioneerConvert`, `PioneerSearch`
3. `PioneerPlots`
4. `PioneerSIMD`
5. app packages + workflow switch
6. runtime cache warming
7. docs/tests cleanup

Do **not** start with `PioneerSIMD`. Without explicit package/app boundaries first, that extraction will be noisy and hard to verify.

## Main Risks

- The current dynamic include system hides many implicit dependencies. Expect the first package split to surface missing imports and ordering bugs.
- Tests currently assume a very broad `using Pioneer` namespace. The umbrella package must remain until package-level tests are in place.
- Plot objects are currently stored inside search result structs. That must be converted to portable diagnostics data before plotting can move cleanly.
- Some tuning/search files mix numeric logic and plotting in the same file. Those files will need surgical splitting, not blanket moves.

## Recommendation

Implement the package split and app-package split together as the architectural spine of the project. Treat plotting and SIMD extraction as two deliberate follow-on slices. This gives you the only credible path to all three goals at once:

- portable installers
- maximum safe sysimage coverage
- persistent non-sysimage caches across launches
