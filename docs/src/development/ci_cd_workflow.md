# CI/CD Workflow Overview

This document outlines the GitHub Actions workflows that power continuous integration and delivery for Pioneer.jl.

## Workflows

The repository uses several workflows:

- `tests.yml` – run the test suite on Ubuntu with Julia 1.11
- `docs.yml` – build docs for pull requests, tags, and manual runs; deploy docs from `main`
- `build_app_linux.yml`, `build_app_macos.yml`, `build_app_windows.yml` – build and package applications
- `CompatHelper.yml` – update package compatibility constraints (scheduled daily)
- `TagBot.yml` – tag releases and interact with the Julia package registry

## Event Matrix

| Event | Condition | Tests | Build Docs  | Deploy Docs | Compile | Release | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| workflow dispatch | tag provided (`v*.*.*`) | :white_check_mark: | :white_check_mark: | Only if `main` | :white_check_mark: | :white_check_mark: | Release new version |
| workflow dispatch | no tag provided | :white_check_mark: | :white_check_mark: | Only if `main` | :white_check_mark: | :x: | Test/build dev version |
| push (tag) | `v*.*.*` | :white_check_mark: | :white_check_mark: | :x: | :white_check_mark: | :white_check_mark: | ? |
| push | `main`/`develop` | :white_check_mark: | :white_check_mark: | :x:| :white_check_mark: | :x: | Hotfix |
| pull request | `main`/`develop` | :white_check_mark: | :white_check_mark: | :x: | :white_check_mark: | :x: | Completed feature |
| pull request | not `main`/`develop` | :white_check_mark: | :white_check_mark: | :x: | :x: | :x: | WIP new feature |
| push | not `main`/`develop` | :white_check_mark: | :white_check_mark: | :x: | :x: | :x: | WIP new feature |

## Additional Notes

- The docs workflow builds documentation for all pull requests, tags, and manual runs, but only deploys to `gh-pages` when invoked on `main`.
- Build workflows use concurrency controls to cancel in-progress runs when new commits arrive on the same reference.
- Tag-based releases (`v*.*.*`) trigger tests, docs, builds, and publishing of GitHub releases.
- `CompatHelper` runs nightly to propose dependency updates; `TagBot` publishes tags to the Julia registry.