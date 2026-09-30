# Decisions

Scope and design decisions.

## 2026-09-24

- **Scope: stepped SWATH only.** ZT Scan DIA (scanning Q1) is out of scope.
  Only `.wiff` + `.wiff.scan` file sets are supported. `.wiff2` is an encrypted
  database and is not supported.
- **Mac first.** Development and the msConvert oracle run on the Apple Silicon
  Mac. Move to a Windows PC only if the oracle cannot be made to work here.
- **New package**, `SciexWiff.jl`, not a Pioneer submodule.
- **CFB/OLE2 parser.** No pure-Julia CFB reader found in the General registry
  (only LibXLS, which wraps a C library), so we write a minimal read-only one.
