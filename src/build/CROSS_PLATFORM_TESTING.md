# Cross-platform build testing

Instructions for validating the PackageCompiler binary on hardware other than the machine it was
diagnosed on (x86_64 macOS, Julia 1.12.5). Written to be handed to a fresh session with no context.

## Background: what was broken and what was fixed

Four independent breakages, all introduced by `6bab0dadb` ("Consolidate SearchDIA pipeline..."),
which changed the package without updating the build:

| # | bug | effect |
|---|---|---|
| 1 | workflows declared `"ParseSpecLib"=>"main_ParseSpecLib"`, removed from the package | `create_app` aborts before compiling — **this is what broke CI** |
| 2 | workflows installed Julia 1.11; `Project.toml` requires `1.12.5` | Pioneer cannot resolve |
| 3 | `dev/Project.toml` pinned `PackageCompiler = "=2.2.3"` | builds fine on Julia ≥1.12, then **dies at startup** on `Base.PROGRAM_FILE` (JuliaLang/julia#56933) |
| 4 | Zenodo `yeast.poin` predates the `.poin` format change | `SearchDIA` precompile target threw; `try/catch` swallowed it, so a **green build shipped with SearchDIA uncompiled** |

Bug 4 is the important one for testing: it was invisible from outside. `snoop.jl` now collects
failures and rethrows at the end, so a broken target fails the build instead of degrading the
artifact. Do **not** set `PIONEER_SNOOP_ALLOW_FAILURES=1` when validating.

All four are fixed on branch `fix/restore-binary-build`. What is *not* yet verified is whether any of
this holds on ARM macOS or Windows.

## Tier 1 — contract check (~20 s, run this first)

Catches workflow/package drift without invoking PackageCompiler.

```
julia --project=. scratchpad/check_build_contract.jl     # path may differ; see the branch
```

Expect `CONTRACT OK`. If an entrypoint reports `MISSING`, the workflows and the package have drifted
again — fix that before spending an hour on a build.

## Tier 2 — fast build (~15 min, no snoop)

```
julia --threads auto --project=dev -e '
  using PackageCompiler
  create_app(".", "build_local/Pioneer";
    incremental=true, force=true,
    executables=[
      "GetSearchParams"=>"main_GetSearchParams",
      "GetBuildLibParams"=>"main_GetBuildLibParams",
      "GetParseSpecLibParams"=>"main_GetParseSpecLibParams",
      "BuildSpecLib"=>"main_BuildSpecLib",
      "SearchDIA"=>"main_SearchDIA",
      "convertMzML"=>"main_convertMzML"])'
```

**Then actually run the binary.** A completed build proves nothing — the 2.2.3 build also completed
and only failed at startup:

```
build_local/Pioneer/bin/GetSearchParams \
  data/ecoli_test/altimeter_ecoli.poin data/ecoli_test/raw /tmp/out \
  --params-path /tmp/out/params.json
echo "exit=$?"        # must be 0, and /tmp/out/params.json must be valid JSON
```

## Tier 3 — full build (30–60 min, what CI does)

Same as tier 2 but `incremental=false` and `precompile_execution_file="src/build/snoop.jl"`.
Requires the precompile data; `snoop.jl` builds its own library via the `BuildSpecLib` target, so no
Zenodo download is needed for `SearchDIA` any more. `BuildSpecLib` needs network access (Koina).

## Platform-specific things to check

### ARM macOS (Apple Silicon)

**The `cpu_target` default is much weaker on ARM.** From `PackageCompiler/src/PackageCompiler.jl`:

```julia
Sys.ARCH === :x86_64  ? "generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)" :
Sys.ARCH === :aarch64 ? "generic"   #= is this really the best here? =#
```

x86_64 gets three multiversioned variants; **aarch64 gets a single `generic` build with no
multiversioning at all** — and that comment is upstream's, not ours. Worth measuring whether an
Apple-Silicon-specific target helps.

**Hand-written SIMD is the thing to watch.** Three `Core.Intrinsics.llvmcall` sites emit explicit
LLVM vector IR:

- `src/structs/SpectralLibrary/PartitionedFragmentIndex/search.jl:31` (`_vcmpge_mask`)
- `src/Routines/SearchDIA/CommonSearchUtils/fusedScan.jl:89,96` (`_fused_vcmpge_mask`, `_fused_vcmple_mask`)

All three are `fcmp oge <8 x float>` → `bitcast <8 x i1> to i8`, i.e. a fixed **256-bit** vector.
This should be portable — LLVM legalises `<8 x float>` to whatever the target has, so on NEON
(128-bit) it becomes two 128-bit operations: correct, just more instructions. **Verify rather than
assume**, because two things could bite:

1. **Alignment.** `unsafe_load(Ptr{F32x8}(pointer(arr, i)))` loads 32 bytes from an arbitrary element
   offset. x86 tolerates unaligned vector loads; confirm AArch64 does here too. A misaligned fault
   would show as a crash or bus error inside fragment-index search, not as a wrong answer.
2. **The `i8` mask bitcast.** `<8 x i1>` → `i8` is target-independent in IR but check the lowering is
   sane rather than a scalarised mess.

Quick check without a full build:

```
julia --project=. -e '
  using Pioneer
  a = rand(Float32, 1000); sort!(a)
  # exercise both SIMD scanners against a scalar reference
  for t in (0.1f0, 0.5f0, 0.9f0, 2.0f0)
    simd = Pioneer.simd_find_first_ge(a, 1, length(a), t)
    ref  = something(findfirst(>=(t), a), length(a)+1)
    @assert simd == ref "simd_find_first_ge disagrees at $t: $simd vs $ref"
  end
  println("SIMD scanners agree with scalar reference")'
```

### Windows (Intel)

- The workflow builds `.exe` wrappers and an installer; check `src/build/windows/`.
- Path separators: `snoop.jl` and the precompile configs use `./data/...` relative paths — confirm
  they resolve when the working directory differs.
- `build_app_windows.yml` had the same three workflow bugs and the same fixes; verify the Julia
  version there is `1.12.6`.

## What to report back

1. Which tier reached, and the **exit code and output of an actual binary invocation** — not just
   whether the build completed.
2. Full text of the first error, if any. `create_app` failures are usually near the top of a long log.
3. `versioninfo()` and `Sys.CPU_NAME` from the build machine.
4. For ARM: whether the SIMD assertion above passes, and whether `SearchDIA` completes on
   `data/precompile/search_ecoli_altimeter.json`.
5. Build wall time (x86_64 reference: 882.8 s for tier 2, incremental, no snoop).
