# CI/CD Workflow Overview

This is internal developer documentation for the GitHub Actions workflows in Pioneer.jl.

This file is intentionally kept outside `docs/src` so it is not included in the compiled Documenter site.

For the operator-facing release checklist, see `docs/release.md`.

## Workflow Inventory

The repository currently contains these GitHub Actions workflows:

- `tests.yml` - runs the source test suite on Ubuntu.
- `docs.yml` - builds docs for pull requests, branch pushes, tags, and manual runs; deploys docs from `develop` and `v*` tags.
- `build_app_linux.yml` - reusable/manual Linux installer build workflow that produces `.deb` packages.
- `build_app_macos.yml` - reusable/manual macOS installer build workflow for Intel and Apple Silicon that produces notarized `.pkg` packages.
- `build_app_windows.yml` - reusable/manual Windows installer build workflow that produces `.msi` packages.
- `validate_installers.yml` - reusable installer-validation workflow shared by manual pre-release runs and nightly validation.
- `pre-release.yml` - automatic post-merge and manual rerun workflow for installer validation before tagging.
- `nightly.yml` - nightly installer-validation workflow that runs on `develop` only when `develop` has advanced since the last successful nightly run.
- `prepare_release.yml` - updates `Project.toml` and opens or refreshes the release PR.
- `release.yml` - validates the tagged version, rebuilds the installer artifacts for the tag, publishes the GitHub release, and pushes Docker images.
- `regression_slurm.yml` - runs the large regression suite on the cluster and publishes regression reports.
- `cancel_regression_on_merge.yml` - cancels in-flight regression runs when a pull request merges.
- `CompatHelper.yml` - opens dependency compatibility update PRs on a daily schedule.
- `registrator.yml` - manually invokes Julia Registrator.
- `TagBot.yml` - handles Julia package tagging automation.

## Current Trigger Summary

| Workflow | Primary triggers | Notes |
| --- | --- | --- |
| `tests.yml` | `push` to `main`/`develop`, pull requests, manual dispatch | Source tests only. |
| `docs.yml` | `push` to `main`/`develop`, `v*` tags, pull requests, manual dispatch | Builds on all triggers; deploys only from `develop` and version tags. |
| `build_app_*` | manual dispatch, `workflow_call` | Reusable installer builders used by `pre-release.yml` and `release.yml`. |
| `validate_installers.yml` | `workflow_call` | Shared installer build-and-validate body. |
| `pre-release.yml` | merged release PRs, manual dispatch | Runs automatically on merged release PRs and stays available for manual reruns. |
| `nightly.yml` | daily schedule on `develop`, manual dispatch | Runs heavy installer validation only when `develop` has moved since the last successful nightly. |
| `prepare_release.yml` | manual dispatch | Bumps `Project.toml` and updates the release PR. |
| `release.yml` | `push` of `v*.*.*` tags | Publication workflow. |
| `regression_slurm.yml` | `develop` pushes, release events, PRs, manual dispatch | Long-running cluster regressions. |
| `cancel_regression_on_merge.yml` | merged pull requests | Cleanup helper. |
| `CompatHelper.yml` | daily schedule, manual dispatch | Dependency maintenance. |
| `registrator.yml` | manual dispatch | Julia package registration. |
| `TagBot.yml` | tag pushes, issue comments, manual dispatch | Registry/tag automation. |

## Current Behavior By Stage

### Source Validation

- `tests.yml` is the main source-validation workflow.
- It currently runs on Ubuntu only.
- It intentionally stays focused on package tests, not installer validation.

### Documentation

- `docs.yml` always builds docs when triggered.
- It only deploys documentation for `develop` and tagged releases.
- Versioned docs are published through Documenter using the `stable`, `v#.#`, and `dev` channels defined in `docs/make.jl`.

### Packaging

- The three `build_app_*` workflows are reusable building blocks for release automation.
- They are no longer part of the default branch/PR validation flow.
- They publish installer artifacts only.
- Raw zipped binaries are no longer treated as supported release deliverables.

### Installed-Artifact Validation

- `validate_installers.yml` is the single reusable workflow for native-installer validation.
- `pre-release.yml` uses it as the release gate before tagging.
- Merging the PR opened by `prepare_release.yml` automatically triggers `pre-release.yml` on the merge commit.
- `pre-release.yml` also remains manually dispatchable for reruns against a specific ref or SHA.
- `nightly.yml` uses it on `develop` as an early warning lane.
- The nightly run skips automatically when the current `develop` HEAD already has a successful nightly result.

### Release Orchestration

The current release-related automation is split across several workflows:

1. `prepare_release.yml` updates `Project.toml` and opens or refreshes the release PR.
2. Merging that release PR automatically triggers `pre-release.yml` on the merge commit.
3. `nightly.yml` gives `develop` a recurring installed-artifact check between releases.
4. `registrator.yml` is dispatched manually to register the Julia package version.
5. `TagBot.yml` handles registry-driven tagging automation.
6. `release.yml` responds to the pushed release tag, rebuilds the platform installer artifacts, publishes the GitHub release, and pushes Docker images.
7. `regression_slurm.yml` also reacts to release publication and generates the cluster-backed regression report.

## Pre-release Validation Scope

`pre-release.yml` validates the installed CLI on clean OS runners by:

- installing the native package for each supported OS target
- checking `pioneer --help`
- downloading the stable RAW fixture from the `ci-fixtures-v1` release asset, verifying its checksum, and converting it with `pioneer convert-raw`
- generating search and library-build parameter templates
- converting the checked-in mzML smoke input
- building a small predicted library
- running a search against the checked-in Arrow and `.poin` test data

Smoke outputs and installer artifacts from the pre-release run are retained for 3 days.

`nightly.yml` runs the same validation scope on `develop`. Scheduled nightly runs only proceed when `develop` has advanced since the last successful nightly run. Manual dispatch of `nightly.yml` always runs the full validation path.

## Release Flow

The release sequence is:

1. Merge the release-preparation PR onto the target release branch.
2. Let `pre-release.yml` run automatically on the merge commit.
3. Review the installer validation jobs:
   - `validate windows`
   - `validate linux`
   - `validate osx intel`
   - `validate osx apple silicon`
4. If the pre-release workflow passes, create the `vX.Y.Z` tag on that same commit.
5. Let `release.yml` rebuild the installers for publication and publish only the supported installer artifacts.

This intentionally separates "prove the installers work" from "publish the release".

## Maintenance Rule

If workflow triggers, artifact policy, or the release ritual change, update this file and `docs/release.md` in the same pull request.
