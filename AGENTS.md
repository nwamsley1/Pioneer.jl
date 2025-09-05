# Repository Guidelines

## Project Structure & Module Organization
- `src/`: Package source. Entry point `src/Pioneer.jl` defines the `Pioneer` module; subpackages live in `Routines/`, `structs/`, and `utils/`.
- `test/`: Test suite (`runtests.jl`, topic files like `test_multi_score_tuning.jl`, and subfolders). Place new tests near related code paths.
- `docs/`: Documenter.jl setup (`make.jl`, `Project.toml`, `src/`).
- `assets/`, `figures/`, `data/`: Static assets, plots/images, and sample inputs (do not commit large datasets).
- `scripts/`: Helper scripts (e.g., `branch_cleanup.sh`).

## Build, Test, and Development Commands
- Setup deps (first run): `julia --project -e 'using Pkg; Pkg.instantiate()'`
- Run tests: `julia --project -e 'using Pkg; Pkg.test()'`
- Run specific tests: `julia --project test/runtests.jl`
- REPL dev: `julia --project` then `using Pioneer`
- Build docs: `julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'`

## Coding Style & Naming Conventions
- Julia style with 4-space indentation; no trailing whitespace.
- Names: modules/types `CamelCase`; functions/variables `lowercase_with_underscores`; constants `UPPER_CASE`.
- Docstrings: use triple-quoted docstrings with examples and argument descriptions.
- Keep functions focused; prefer explicit imports and clear return types when helpful.

## Testing Guidelines
- Framework: `Test` stdlib; use `@testset` with descriptive names.
- Place tests under `test/` mirroring directories (e.g., `test/Routines/...`).
- Strive to keep or increase coverage (Codecov badge tracks CI). Add regression tests for bug fixes.

## Commit & Pull Request Guidelines
- Branching: work on `feature/*` (enhancements) or `hotfix/*` (critical fixes). Merge to `develop` (or `main` for hotfixes as needed).
- Conventional Commits: `feat:`, `fix:`, `docs:`, `refactor:`, `test:`, `chore:`.
- PRs: clear description, link issues (`Fixes #123`), include tests, update docs where applicable, and ensure CI is green.

## Security & Configuration Tips
- Do not commit secrets or large raw datasets; prefer paths/config files.
- Parameters are managed via JSON for CLI workflows; keep example configs small and anonymized.
- Use `.gitignore` for generated outputs; keep `data/` minimal and reproducible.

