# macOS Packaging — Signing/Notarization Streamlining Plan

## Goal

Cut macOS Build & Package time by only signing what’s necessary, avoiding deep re‑sign passes, and pruning non‑runtime files. Keep the output fully notarizable and identical in behavior at install/run time.

## Current State (Summary)

- Workflow: `.github/workflows/build_app_macos.yml` runs on push/PR and on release via `workflow_call`.
- Packaging steps:
  1) Build app layout under `pkgroot/usr/local/Pioneer`.
  2) Codesign all frameworks, then run a broad Mach‑O scan and sign every file that looks like a binary, then deep‑sign the entire tree again.
  3) `pkgbuild` → `productsign` → `notarytool --wait` → `stapler` → upload artifacts.
- Pain points:
  - The “sign everything” loop scans huge artifact trees (Qt, GR, Julia artifacts) and re‑signs many non‑runtime files (e.g., `*.o`) and duplicates signing work due to a final deep sign.
  - Six‑hour job timeout indicates signing is a major bottleneck even before/alongside notarization.

## Principles

- Only sign code objects that ship and execute:
  - App bundle main binaries (`*/Contents/MacOS/*`)
  - Framework bundles (`*.framework`)
  - CLI executables under `bin/` (executable bit)
  - Dynamic libraries/plugins (`*.dylib`, `*.so`, Qt/GR plugins)
- Do not sign build debris and non‑runtime files:
  - Object files `*.o`, static archives `*.a`, intermediate `objects-Release` trees
- Avoid deep re‑sign at the end; keep targeted signing passes and a single verification step.
- Reduce per‑file timestamp usage to lower network overhead; notarization of the final `.pkg` is the authoritative ticket.

## Implementation Steps

1) Restrict packaging triggers (optional but recommended)
   - Make heavy packaging/signing/notarization release‑only:
     - Remove `push`/`pull_request` triggers from `build_app_macos.yml` and keep `workflow_call` + `workflow_dispatch`.
     - Or keep push/PR but gate the four heavy steps (codesign, notarize, staple, upload signed pkg) with `if: startsWith(github.ref, 'refs/tags/v')`.

2) Pre‑prune non‑runtime files
   - Before signing, delete object/static files and intermediate build dirs within the payload:
     ```bash
     # Remove object files and archives
     find "$PKGROOT/usr/local/$APP" -type f \( -name '*.o' -o -name '*.a' \) -delete
     # Remove known intermediate dirs
     find "$PKGROOT/usr/local/$APP" -type d -name 'objects-Release' -prune -exec rm -rf {} +
     ```

3) Targeted signing passes (replace broad scan + deep sign)
   - Frameworks (first):
     ```bash
     find "$PKGROOT/usr/local/$APP" -type d -name '*.framework' -print0 | \
       while IFS= read -r -d '' fw; do
         echo "Signing framework: $fw"
         codesign --force --options runtime \
           --entitlements src/build/osx/entitlements.plist \
           --sign "$CODESIGN_IDENTITY" "$fw"
       done
     ```
   - App bundle binaries (GR apps):
     ```bash
     find "$PKGROOT/usr/local/$APP" -type f -path '*/Contents/MacOS/*' -print0 | \
       while IFS= read -r -d '' f; do
         echo "Signing app binary: $f"
         codesign --force --options runtime \
           --entitlements src/build/osx/entitlements.plist \
           --sign "$CODESIGN_IDENTITY" "$f"
       done
     ```
   - CLI executables in bin/:
     ```bash
     if [ -d "$PKGROOT/usr/local/$APP/bin" ]; then
       find "$PKGROOT/usr/local/$APP/bin" -type f -perm -111 -print0 | \
         while IFS= read -r -d '' f; do
           echo "Signing CLI: $f"
           codesign --force --options runtime \
             --entitlements src/build/osx/entitlements.plist \
             --sign "$CODESIGN_IDENTITY" "$f"
         done
     fi
     ```
   - Dynamic libraries and plugins (avoid frameworks):
     ```bash
     find "$PKGROOT/usr/local/$APP" \( -name '*.dylib' -o -name '*.so' \) \
       -not -path '*/.framework/*' -print0 | \
       while IFS= read -r -d '' f; do
         echo "Signing dylib/plugin: $f"
         codesign --force --options runtime \
           --entitlements src/build/osx/entitlements.plist \
           --sign "$CODESIGN_IDENTITY" "$f"
       done
     ```
   - Notes:
     - Deliberately omit `--timestamp` on per‑file signing to reduce latency (final notarized pkg is authoritative).
     - No final `codesign --deep` over the whole tree.

4) Verify signatures once
   ```bash
   codesign --verify --deep --strict --verbose=2 "$PKGROOT/usr/local/$APP"
   ```

5) Build, sign, and notarize package (unchanged)
   - `pkgbuild` → `productsign` → `notarytool --wait` → `stapler`.
   - Consider adding a readable timeout to the Notarize step (e.g., 120 minutes) for earlier feedback, while recognizing hosted runner hard limit is 6h.

6) Logging improvements
   - Print counts per signing pass (e.g., frameworks signed, app binaries signed, dylibs signed) to make hotspots visible in run logs.

## Optional Cleanup (further reductions)

- Clear extended attributes (if present) before signing to avoid quarantine flag propagation:
  ```bash
  xattr -rd com.apple.quarantine "$PKGROOT/usr/local/$APP" || true
  ```
- Remove unused Qt build artifacts (if safe): e.g., `plugins/*/objects-*` directories — already covered by the generic prune above.
- If verification fails due to a missed plugin path in a specific release, add a targeted mini‑pass for that subpath rather than re‑introducing `--deep` globally.

## Risk Assessment and Mitigations

- Risk: A required code object is missed by the targeted passes.
  - Mitigation: Keep a strict `codesign --verify --deep --strict` check. In failure cases, inspect the path and add a specific signing rule for that subtree.
- Risk: Removing `--timestamp` on individual files could reduce trust chain detail.
  - Mitigation: Final notarization of the `.pkg` is the requirement for Gatekeeper; per‑file timestamp is not required.
- Risk: Accidental deletion of needed files during prune.
  - Mitigation: Restrict prune to `*.o`, `*.a`, and `objects-Release` only; do not remove `.dylib`/`.so` or app resources.

## Validation Plan

1) Dry run on `workflow_dispatch` with the new signing passes; capture per‑step durations.
2) Confirm `codesign --verify --deep --strict` passes.
3) Confirm notarization status Accepted and stapling succeeds.
4) Smoke test: install `.pkg` on a clean macOS VM for both x64 and arm64; run CLI and launch GR apps.
5) Compare job total time vs. baseline; target >50% reduction.

## Rollout and Backout

- Rollout: Merge workflow change, test on a release tag candidate.
- Backout: Revert to previous workflow or temporarily add a final `--deep` pass for a subset path if a corner case appears. Keep the old logic in VCS for quick rollback.

## Follow‑ups (Optional)

- Make packaging fully release‑only to keep CI lightweight on pushes/PRs.
- Add summaries to job logs with itemized counts and durations per signing phase.
