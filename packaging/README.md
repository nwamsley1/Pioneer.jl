# Pioneer application packaging

Development, tests, and the packaged runtime use the repository-root
`Project.toml` and committed `Manifest.toml`. The `compiler/` environment
independently locks PackageCompiler and build-only dependencies so those tools
do not become part of the shipped runtime graph. The optional `dev/` environment
remains available for interactive tooling, but it does not drive app builds.

Build from the repository root with:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=packaging/compiler -e 'using Pkg; Pkg.instantiate()'
julia --project=. packaging/check_app_contract.jl
julia --threads auto --project=packaging/compiler \
    packaging/build_app.jl build/Pioneer/Applications/Pioneer
```

`build_app.jl` is the single source of truth for executables, precompile inputs,
and reviewed CPU targets. `PIONEER_APP_CPU_TARGET` is available only for
diagnostic overrides.

When dependencies change, resolve and instantiate the root and compiler
environments, review the manifest diffs, run the contract check, and validate
installed applications on Windows x64, Linux x64, macOS x64, and macOS arm64.
Do not edit either manifest by hand.
