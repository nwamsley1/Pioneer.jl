# CI/CD Workflow Overview

This document outlines the GitHub Actions workflows that power continuous integration and delivery for Pioneer.jl.

## Workflows

The repository uses several workflows:

- `tests.yml` – run the test suite on Ubuntu
- `docs.yml` – build docs for pull requests, tags, and manual runs; deploy docs from `main` and `v*` tags
- `build_app_linux.yml`, `build_app_macos.yml`, `build_app_windows.yml` – build and package applications
- `CompatHelper.yml` – update package compatibility constraints (scheduled daily)
- `TagBot.yml` – tag releases and interact with the Julia package registry
- `registrator.yml` – submit releases to the Julia General registry

## Event Matrix

| Event | Condition | Tests | Build Docs  | Deploy Docs | Compile | Release | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| push (tag) | `v*.*.*` | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | Release new version |
| push | `main`/`develop` | :white_check_mark: | :white_check_mark: | Only if `develop` | :white_check_mark: | :x: | Merge or hotfix |
| pull request | `main`/`develop` | :white_check_mark: | :white_check_mark: | :x: | :white_check_mark: | :x: | Completed feature |
| pull request | not `main`/`develop` | :white_check_mark: | :white_check_mark: | :x: | :x: | :x: | WIP feature |
| push | not `main`/`develop` | :white_check_mark: | :white_check_mark: | :x: | :x: | :x: | WIP feature |

## Manual Workflow Dispatches

Manual dispatches allow running workflows on demand.

| Dispatch | Condition | Tests | Build Docs | Deploy Docs | Compile | Release | Purpose |
|----------|-----------|-------|------|--------|---------|---------|---------|
| `tests.yml`, `docs.yml`, `build_app_*` | `v*.*.*` | :white_check_mark: | :white_check_mark: | :white_check_mark:  | :white_check_mark: | :white_check_mark: | Release new version |
| `tests.yml`, `docs.yml`, `build_app_*` | no tag | :white_check_mark: | :white_check_mark: | `main`/`develop` | :white_check_mark: | :x: | Test/build dev version |
| `registrator.yml` | | :x: | :x: | :x: | :x: | :white_check_mark: | Register release with Julia registry |

## Additional Notes

- The docs workflow builds documentation for all pull requests, tags, and manual runs, but only deploys to `gh-pages` when invoked on `main`, `develop`, or a `v*` tag.
- Tag-based releases (`v*.*.*`) trigger tests, docs, builds, publishing of GitHub releases, and versioned documentation.
- `CompatHelper` runs nightly to propose dependency updates; `TagBot` publishes tags to the Julia registry.
- Documentation is published for `stable`, semantic version tags (`v#.#.#`), and `dev`.