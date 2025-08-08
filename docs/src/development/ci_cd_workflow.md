# CI/CD Workflow Overview

This document outlines the GitHub Actions workflows that power continuous integration and delivery for Pioneer.jl.

## Workflows

The repository uses several workflows:

- `tests.yml` – run the test suite on Ubuntu with Julia 1.11
- `docs.yml` – build docs for pull requests, tags, and manual runs; deploy docs from `main`
- `build_app_linux.yml`, `build_app_macos.yml`, `build_app_windows.yml` – build and package applications for Linux, macOS, and Windows
- `CompatHelper.yml` – update package compatibility constraints (scheduled daily)
- `TagBot.yml` – tag releases and interact with the Julia package registry

## Event Matrix

| Event | Ref / Condition | Tests | Docs – Build | Docs – Deploy | Build & Package (Linux/macOS/Windows) | Publish Release (Linux/macOS/Windows) |
| --- | --- | --- | --- | --- | --- | --- |
| push | branch == `main` | Yes | Yes | Yes | Yes | No |
| push | branch == `develop` | Yes | No | No | Yes | No |
| push | branch ∉ {`main`, `develop`} | Yes | No | No | No | No |
| push (tag) | tag matches `v*.*.*` | Yes | Yes | No | Yes | Yes |
| pull_request | base in {`main`, `develop`} | Yes | Yes | No | Yes | No |
| pull_request | base ∉ {`main`, `develop`} | Yes | Yes | No | No | No |
| workflow_dispatch | tag provided | Yes | Yes | Only if ref == `main` | Yes | Yes |
| workflow_dispatch | no tag provided | Yes | Yes | Only if ref == `main` | Yes | No |

## Additional Notes

- The docs workflow builds documentation for all pull requests, tags, and manual runs, but only deploys to `gh-pages` when invoked on `main`.
- Build workflows use concurrency controls to cancel in-progress runs when new commits arrive on the same reference.
- Tag-based releases (`v*.*.*`) trigger tests, docs, builds, and publishing of GitHub releases.
- `CompatHelper` runs nightly to propose dependency updates; `TagBot` publishes tags to the Julia registry.