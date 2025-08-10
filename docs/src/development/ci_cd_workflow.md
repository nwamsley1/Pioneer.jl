# CI/CD Workflow Overview

This document outlines the GitHub Actions workflows that power continuous integration and delivery for Pioneer.jl.

## Workflows

The repository uses several workflows:

- `tests.yml` – run the test suite on Ubuntu
- `docs.yml` – build docs for pull requests, tags, and manual runs; deploy docs from `main` and `v*` tags
- `build_app_linux.yml`, `build_app_macos.yml`, `build_app_windows.yml` – build and package applications; reusable via `workflow_call`
- `release.yml` – orchestrate cross-platform builds and publish GitHub releases for tags
- `CompatHelper.yml` – update package compatibility constraints (scheduled daily)
- `TagBot.yml` – tag releases and interact with the Julia package registry
- `registrator.yml` – run Registrator to open a registration PR in the Julia General registry

## Event Matrix

| Event | Condition | Tests | Build Docs  | Deploy Docs | Compile | Release | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| push (tag) | `v*.*.*` | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | Release new version |
| push | `main`/`develop` | :white_check_mark: | :white_check_mark: | if `develop` | :white_check_mark: | :x: | Merge or hotfix |
| pull request | `main`/`develop` | :white_check_mark: | :white_check_mark: | :x: | :white_check_mark: | :x: | Completed feature |
| pull request | not `main`/`develop` | :white_check_mark: | :white_check_mark: | :x: | :x: | :x: | WIP feature |
| push | not `main`/`develop` | :white_check_mark: | :white_check_mark: | :x: | :x: | :x: | WIP feature |

## Manual Workflow Dispatches

Manual dispatches allow running workflows on demand.

| Workflow | Condition | Tests | Build Docs | Deploy Docs | Compile | Release | Purpose |
|----------|-----------|-------|-----------|-------------|---------|---------|---------|
| `registrator.yml` | Julia registration success | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | Release new version |
| `release.yml` | `v*.*.*` | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | Build & publish release |
| `tests.yml`, `docs.yml`, `build_app_*` | `v*.*.*` | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x: | Manually rerun or troubleshoot |
| `tests.yml`, `docs.yml`, `build_app_*` | no tag | :white_check_mark: | :white_check_mark: | `develop` | :white_check_mark: | :x: | Test/build dev version |
  
## Release Process

To cut a new release:

1. Update the version number in `Project.toml`.
2. Merge that change into the `main` branch.
3. Manually trigger `registrator.yml`, which runs `JuliaRegistries/Registrator@v1` to open or update a registration pull request in the Julia General registry.
4. After the registry PR is merged, `TagBot` creates a tag and GitHub release.
5. The tag triggers the docs workflow to build and deploy versioned documentation and the `release.yml` workflow. The release workflow calls the platform build workflows (`build_app_linux.yml`, `build_app_macos.yml`, `build_app_windows.yml`) to compile the application and attach the artifacts to the GitHub release.

## Pre-releases

- Tags containing a hyphen (e.g., `v1.2.0-rc1`) are treated as pre-releases.
- TagBot marks such tags as GitHub prereleases.
- Build workflows still run, and packaging steps strip pre-release/build metadata when generating Windows installers.

## Additional Notes

- The docs workflow builds documentation for all pull requests, tags, and manual runs, but only deploys to `gh-pages` when invoked on `main`, `develop`, or a `v*` tag.
- Tag-based releases (`v*.*.*`) trigger tests, docs, builds, publishing of GitHub releases, and versioned documentation.
- `CompatHelper` runs nightly to propose dependency updates; `TagBot` publishes tags to the Julia General registry.
- Documentation is published for `stable`, semantic version tags (`v#.#.#`), and `dev`.