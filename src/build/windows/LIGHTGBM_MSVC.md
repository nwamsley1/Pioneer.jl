# LightGBM on Windows: use the MSVC build, not the shipped MinGW one

**Status:** implemented on branch `feature/windows-lightgbm-msvc`.
**Applies to:** Windows only. Linux/macOS are unaffected and unchanged.

## TL;DR

The `LightGBM_jll` binary Julia ships on Windows is a **mingw-w64** build whose OpenMP
runtime (GNU `libgomp` on winpthreads) does **not** parallelize the LightGBM tree-learner
`FIT`: on an 18-core box it runs at **~1 core**. Microsoft's **MSVC** build of the *same*
LightGBM version (3.3.5) saturates all cores and produces **bit-identical** results
(Pioneer trains LightGBM with `seed=1776, deterministic=true, force_row_wise=true`).

This branch applies that swap in two places:

| Channel | Who | Mechanism |
|---|---|---|
| **Packaged Windows app** | end users | build-time DLL swap + a WiX bootstrapper (`PioneerSetup.exe`) that installs the VC++ redistributable, then Pioneer. One UAC prompt at install; nothing at first launch. |
| **Source / `Pkg.add` installs** | developers | `Pioneer.setup_windows_lightgbm()` writes a `LightGBM_jll` artifact override. Run once, restart Julia. |

Both are pinned to a SHA-256-verified MSVC `lib_lightgbm.dll` v3.3.5.

---

## 1. Symptom and root cause

The same Pioneer search on identical hardware spends far more wall time in the LightGBM
stages on Windows than on Linux. The cause is **which OpenMP runtime the binary links**:

- `LightGBM.jl` v2.2.2 → `LightGBM_jll` v3.3.5+1. The Windows artifact is built by
  Yggdrasil with mingw-w64 and takes OpenMP from GNU `libgomp` (on winpthreads). There is
  no MSVC build in the JLL.
- OpenMP **is** present (the DLL imports `GOMP_parallel`, `omp_set_num_threads`, …). This is
  **not** a "no threading" problem — it is a *threading-efficiency* problem. LightGBM's tree
  learner opens a very large number of **short-lived** parallel regions (per split, per
  bin-construction pass). `libgomp`'s fork/join and barriers on winpthreads cost far more per
  region than a Linux futex or MSVC's `vcomp140`; when per-region work is small relative to
  barrier cost, worker threads sit in the barrier instead of computing.
- Diagnostic signature: **low CPU utilization during FIT**, not busy cores.

LightGBM upstream documents exactly this (Installation Guide / FAQ #8: Visual Studio "may be
10x faster than MinGW" for many-core systems).

## 2. Evidence (measured on an i9-10980XE, 18 physical / 36 logical, Julia 1.12.6, `-t 18`)

**Isolation benchmark** (LightGBM C library only, 1,000,000 × 50; cores = process
CPU-seconds / wall):

| build | FIT wall | FIT cores busy | PREDICT wall |
|---|---|---|---|
| MinGW (shipped) | 4.56 s | **1.1 / 18** | 0.54 s |
| MSVC (override) | 2.27 s | **14.9 / 18** | 0.57 s |

Free env-var tuning did **not** help FIT (all stayed ~1 core): `OMP_PROC_BIND=true` +
`OMP_PLACES=cores` → 4.61 s; 9 threads → 4.26 s; 4 threads → 6.08 s.

**End-to-end, 3-file warm-start `SearchDIA`** (YEAST_KO, MBR on, q=0.01; warmup search then a
timed search in the same session so JIT is excluded; only the DLL differs):

| | MinGW | MSVC |
|---|---|---|
| Precursor Scoring stage (heavy LightGBM) | 60.5 s | **26.0 s** (−57%) |
| Main Search stage (mostly non-LightGBM) | 40.9 s | 38.4 s |
| Protein Scoring stage | 5.7 s | 5.7 s |
| **warm total** | 2.59 min | **1.98 min** (−24%) |
| unique precursors | 83,316 | **83,316** |
| protein-group rows | 13,827 | **13,827** |

IDs are **bit-identical** — the correctness gate. The faster binary is numerically the same.

**Why Main Search barely moves (it uses LightGBM too — the fix IS applied there).** The
override is process-wide: exactly one `lib_lightgbm.dll` is dlopened, so every LightGBM call
uses MSVC. Main Search's model is just deliberately *lightweight*, and the stage is dominated
by non-LightGBM spectral deconvolution (Julia threads, not OpenMP). Microbenchmark of the two
real model shapes:

| Pioneer model | shape | MinGW FIT | MSVC FIT |
|---|---|---|---|
| MainSearch per-file (`MAINSEARCH_LGBM_HP`) | 250k × 30, depth 4, 15 leaves, 50 iters | 0.79 s | 0.45 s |
| Precursor Scoring (`SCORING_LGBM_HP`) | 1M × 40, depth 8, 64 leaves, 200 iters | 10.94 s | 5.00 s |

So the MainSearch model *does* get faster (1.76×), but at ~0.5 s/file it is a small slice of a
~14 s/file stage; the LightGBM time — and the win — concentrates in Precursor Scoring.

## 3. The fix

The verified binary is Microsoft's MSVC `lib_lightgbm.dll` v3.3.5, distributed in the official
PyPI wheel (a zip). It imports `VCOMP140.DLL` / `VCRUNTIME140.dll` / `MSVCP140.dll` and is
2.4 MB (vs the 21 MB mingw DLL).

- Wheel: `https://files.pythonhosted.org/packages/6a/31/13f80e5abac627c53375c1dc49404d622d929272a231c05b2f4f7bb98b9e/lightgbm-3.3.5-py3-none-win_amd64.whl`
- Wheel SHA-256: `02a40745c1972cf4a2cde764c7739228f45178c2237af2df40fde7063a58ac6a`
- DLL SHA-256: `a52a9494c894f9987f9be24c97e0dec7c37fd23b1184aa0ac39d3200ffa51699`

### 3a. Packaged app (end users) — build-time swap

`PackageCompiler.create_app` bundles a private copy of the LightGBM artifact into the app tree.
Two workflow steps in `.github/workflows/build_app_windows.yml` (added by this branch) do the work.

**"Swap in MSVC LightGBM DLL"** runs right after `create_app` and:

1. downloads the pinned wheel and verifies its SHA-256;
2. extracts `lib_lightgbm.dll` and verifies its SHA-256;
3. overwrites every bundled `lib_lightgbm.dll` under the app tree (and **fails the build** if
   none is found, so a future `create_app` layout change can never silently skip the fix).

**"Download VC++ redistributable"** + **"Build installer bundle"** then wrap the Pioneer MSI in a
WiX Burn bootstrapper (`src/build/windows/Bundle.wxs`) that chains Microsoft's official
`vc_redist.x64.exe` ahead of the MSI. The shipped artifact is a single `PioneerSetup-*.exe`; the
bare `.msi` is still emitted as an internal component (see §4). The MSVC DLL needs **no**
`Overrides.toml` — the bundled artifact simply *is* the MSVC build — and its runtime is provided
by the redistributable the bootstrapper installs.

### 3b. Source / `Pkg.add` installs (developers) — `setup_windows_lightgbm()`

For users who run Pioneer from a normal Julia session (not the packaged app), call once:

```julia
using Pioneer
setup_windows_lightgbm()   # no-op off Windows; idempotent
# then RESTART Julia
```

It downloads + SHA-256-verifies the MSVC DLL, lays it out with a `lightgbm.exe` copied from the
current artifact (the JLL declares two products, so both must exist), and merges a
`LightGBM_jll` override into `~/.julia/artifacts/Overrides.toml` (UUID
`0e4427ef-1ff7-5cd7-8faa-8ff0877bb2ec`). It has **no new Julia dependencies** — download,
SHA-256, and unzip go through PowerShell; the TOML merge uses `Base.TOML`. Implementation:
`src/utils/install/windows_lightgbm.jl`.

## 4. Visual C++ Redistributable

The MSVC `lib_lightgbm.dll` depends on the VC++ 2015-2022 runtime. Its imports are **`VCOMP140.DLL`
(the MSVC OpenMP runtime — the whole reason the swap works), `VCRUNTIME140.dll`,
`VCRUNTIME140_1.dll`, and `MSVCP140.dll`** (and `msvcp140.dll` itself also imports
`vcruntime140_1.dll`). Note `vcruntime140_1.dll` — a *separate* file from `vcruntime140.dll`,
introduced with VS 2019 for the C++ exception-handling ABI. It is easy to miss when hand-listing
DLLs, and omitting it makes LightGBM fail to load on a clean machine.

**Packaged app — the installer installs the official redistributable (implemented).** The Burn
bootstrapper chains `vc_redist.x64.exe` (fetched at build time from Microsoft's official
`https://aka.ms/vs/17/release/vc_redist.x64.exe`) ahead of the Pioneer MSI, with a registry
`DetectCondition` so it is skipped when the runtime is already present and a `3010` exit-code
handler for the reboot-required case. This runs **at install time** (one UAC prompt, before the
app is ever launched — not lazily on first use, which would be a chicken-and-egg trap since the
missing runtime is exactly what the app needs to start). It is:

- **complete by construction** — installs the entire runtime, so `vcruntime140_1` and any future
  dependency are covered automatically, with no per-DLL bookkeeping;
- **correct on provenance** — `vc_redist.x64.exe` is Microsoft's redistributable, explicitly
  licensed for distribution, unlike copying the same-named DLLs out of `System32` (covered by the
  Windows licence, not the "Distributable Code" grant);
- **self-updating** — tracks the latest security-patched runtime;
- **independent of the build machine** — fetched from Microsoft, not scraped from the runner's
  `System32` or a VS `VC\Redist` folder (that folder is an *optional* VS component and is not
  guaranteed present even when Visual Studio / Build Tools is installed).

It requires the installing user to have admin rights (the MSI already installs per-machine into
`Program Files`, so this is no new constraint).

**Rejected: app-local runtime DLLs.** Copying `vcomp140`/`vcruntime140`/`vcruntime140_1`/`msvcp140`
next to the app avoids elevation and keeps a portable copy, but we would own keeping that DLL set
complete (the first cut missed `vcruntime140_1`), and copying from `System32` takes the file from
the source that does not grant redistribution. If a no-admin/portable install is ever needed, the
provenance-correct app-local route is to **extract** the DLLs from the downloaded
`vc_redist.x64.exe` rather than from `System32`.

**Source installs:** `setup_windows_lightgbm()` **detects** the runtime (looks for `vcomp140.dll`
etc. in `System32`) and, if missing, prints the official download link. A Julia function cannot
install it silently (that needs elevation), and dev machines almost always already have it.

## 5. Maintenance: version coupling (important)

The MSVC DLL version **must** match `LightGBM_jll`'s C API, which is pinned by
`LightGBM = "2.2"` in `Project.toml` (→ LightGBM 3.3.5). If that compat is ever bumped:

1. find the matching Microsoft MSVC wheel for the new LightGBM version,
2. update the URL + both SHA-256 hashes in **both**
   `.github/workflows/build_app_windows.yml` **and** `src/utils/install/windows_lightgbm.jl`.

The CI step fails the build if the pinned hash does not match, so a stale/incompatible pin is
caught rather than shipped. Consider adding an assertion in `validate_installers.yml` that the
packaged app's `lib_lightgbm.dll` is the MSVC build and loads.

## 6. Verify / roll back

**Verify (either channel):**

```julia
using LightGBM, Libdl
filter(h -> occursin("lightgbm", lowercase(h)), Libdl.dllist())
# source install: path should be under <depot>/pioneer_lgbm_msvc/bin
# packaged app: the bundled DLL should be ~2.4 MB (MSVC) not ~21 MB (mingw)
```

A quick FIT on ~1M×50 should now saturate the cores (Task Manager) instead of sitting near
one core.

**CI coverage and its limit.** `validate_installers.yml` (Windows job) installs the
`PioneerSetup-*.exe` bootstrapper and runs a smoke search, so the bundle is exercised end to
end on every validation run. But that runner already has the VC++ runtime, so the redist's
`DetectCondition` finds it present and **skips** installation — CI therefore proves the bundle
builds/installs/runs, **not** that the redist installs correctly on a bare machine. That case
must be checked **manually once on a runtime-free Windows VM** before a release: on a clean image
with no VC++ redistributable, run `PioneerSetup-*.exe`, confirm the redist actually installs, and
confirm a search runs (LightGBM loads). Hosted runners cannot cover this because they all ship the
runtime.

**Roll back a source install:** delete the `0e4427ef-1ff7-5cd7-8faa-8ff0877bb2ec` entry from
`~/.julia/artifacts/Overrides.toml` (and optionally remove `<depot>/pioneer_lgbm_msvc`), then
restart Julia. **Roll back the app:** revert the CI workflow steps (the swap step and the bundle
build/upload).

## 7. Known dead end (do not retry)

Raising `num_threads` for LightGBM calls made *from Julia worker threads* (e.g. inside the
per-file `parallel_foreach!` in MainSearch / the scoring passes) is ~8× faster but silently
**corrupts results (−28% IDs)** — a multi-threaded OpenMP region launched from a Julia worker
thread is unsafe here. Those calls stay single-threaded; the DLL swap is orthogonal and
preserved exact IDs.
