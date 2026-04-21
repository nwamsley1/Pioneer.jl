# Release Runbook

This file is the operator-facing checklist for Pioneer releases.

The internal workflow reference at `docs/ci_cd_workflow.md` explains the full GitHub Actions map. This runbook is intentionally shorter and answers one question: what do we do when we are about to ship a release?

## Release Ritual

The release process is:

1. Run `prepare_release.yml` and merge the release-preparation PR.
2. Let `pre-release.yml` run automatically on the merge commit.
3. Review the installer validation jobs for:
   - Windows
   - Linux
   - macOS Intel
   - macOS Apple Silicon
4. If the pre-release checks look good, create the `vX.Y.Z` tag on the same commit.
5. Let `release.yml` publish the release artifacts.

`pre-release.yml` remains manually dispatchable for reruns, but the normal path is the automatic run after the release PR merges. `nightly.yml` runs the same installed-artifact checks on `develop` as an early warning lane, but it does not replace the merged-PR `pre-release.yml` gate before tagging.

Those installed-artifact checks now include a real `convert-raw` smoke test using the stable `ci-fixtures-v1` RAW fixture release asset, plus checksum verification before conversion.

## Artifact Policy

The release should publish only native installer artifacts:

- Windows: `.msi`
- Linux: `.deb`
- macOS Intel: `.pkg`
- macOS Apple Silicon: `.pkg`

Raw zipped binaries should not be treated as supported release downloads.

`release.yml` publishes the GitHub release installer assets only. Installer-validation logs and smoke outputs live in the `pre-release.yml` run artifacts for 3 days.

## Update Rule

If the release ritual changes, update this file and `docs/ci_cd_workflow.md` in the same pull request as the workflow change.
