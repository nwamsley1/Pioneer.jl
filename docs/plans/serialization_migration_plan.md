# Implementation Plan: JLD2 → Base Julia Serialization

**Branch**: `feature/serialization-for-fragments`
**Created**: 2026-02-02
**Status**: In Progress

---

## Overview

Convert spectral library fragment storage from JLD2 format to Julia's base `Serialization` module with Zlib compression (using existing CodecZlib dependency).

### Files Being Converted

| Current File | New File | Contents |
|-------------|----------|----------|
| `detailed_fragments.jld2` | `detailed_fragments.jls` | `Vector{DetailedFrag}` or `Vector{SplineDetailedFrag}` |
| `precursor_to_fragment_indices.jld2` | `precursor_to_fragment_indices.jls` | `Vector{UInt64}` |
| `spline_knots.jld2` | `spline_knots.jls` | `NTuple{M, T}` |

---

## Motivation

### JLD2 (Current)

**Pros:**
- HDF5-compatible format readable by other languages (Python h5py, MATLAB)
- Named keys for stored objects (dictionary-like access)
- Partial/lazy loading capability
- Compression support built-in

**Cons:**
- External dependency (`JLD2.jl`)
- Slower than native serialization for simple data structures
- HDF5 complexity overhead for simple use cases
- Occasional breaking changes between JLD2 versions
- Type reconstruction can be fragile across package versions

### Base Julia Serialization (Proposed)

**Pros:**
- **Zero additional dependencies** - Serialization is in Base, CodecZlib already a dependency
- **Faster read/write** for native Julia types
- Simpler API: `serialize(io, obj)` / `deserialize(io)`
- More stable across Julia versions for simple struct types
- Better suited for single-object-per-file patterns

**Cons:**
- **Julia-only** - cannot read files from Python/R/MATLAB
- No partial/lazy loading
- No named keys (one object per serialize call)

### Why This is a Good Fit

1. The data is only consumed by Pioneer.jl (Julia-only)
2. We're storing a single `Vector{DetailedFrag}` per file (no need for named keys)
3. The structs are simple, fixed-layout types (ideal for Serialization)
4. Performance improvement for large libraries

---

## Compression (CodecZlib - already a dependency)

All three formats use the same DEFLATE algorithm internally:

| Format | Stream Types | Header | Checksum | Notes |
|--------|--------------|--------|----------|-------|
| **Gzip** | `GzipCompressor/DecompressorStream` | 10 bytes | CRC32 | Most portable, `.gz` compatible |
| **Zlib** | `ZlibCompressor/DecompressorStream` | 2 bytes | Adler32 | Slightly smaller than Gzip |
| **Deflate** | `DeflateCompressor/DecompressorStream` | None | None | Smallest, but no integrity check |

Compression levels 1-9 available (1=fastest, 6=balanced, 9=best compression).

**Decision**: Use **Zlib compression (level 6)** as default - good balance, smaller header than Gzip, includes checksum.

---

## Files to Modify

| File | Changes |
|------|---------|
| `src/utils/serialization.jl` | **NEW** - Utility functions for serialize/deserialize with compression |
| `src/structs/LibraryIon.jl` | Update `save_detailed_frags` and `load_detailed_frags` |
| `src/Routines/BuildSpecLib/build/build_poin_lib.jl` | Update save calls, change `.jld2` → `.jls` |
| `src/Routines/BuildSpecLib.jl` | Update `spline_knots` save |
| `src/Routines/SearchDIA/ParseInputs/loadSpectralLibrary.jl` | Update load calls with backwards compatibility |

---

## Implementation

### Phase 1: Create Serialization Utilities

Create `src/utils/serialization.jl`:
```julia
using Serialization
using CodecZlib

function serialize_to_jls(filepath::String, data; level::Int=6)
    open(filepath, "w") do file_io
        stream = ZlibCompressorStream(file_io; level=level)
        try
            serialize(stream, data)
        finally
            close(stream)
        end
    end
end

function deserialize_from_jls(filepath::String)
    open(filepath, "r") do file_io
        stream = ZlibDecompressorStream(file_io)
        try
            deserialize(stream)
        finally
            close(stream)
        end
    end
end
```

### Phase 2: Update Save Functions

- `LibraryIon.jl`: Update `save_detailed_frags` to use `serialize_to_jls`
- `build_poin_lib.jl`: Change file extensions and use new serialization
- `BuildSpecLib.jl`: Update `spline_knots` save

### Phase 3: Update Load Functions (with backwards compatibility)

- Try `.jls` first
- Fall back to `.jld2` with deprecation warning
- Eventually remove JLD2 fallback

---

## Backwards Compatibility

```
Load Flow:
1. Check for .jls file → deserialize_from_jls()
2. If not found, check for .jld2 → JLD2.load() + warning
3. If neither found → error
```

---

## Benchmarking

After implementation, run `benchmarks/compression_benchmark.jl` to compare:
- JLD2 (baseline)
- Serialization + No compression
- Serialization + Gzip (levels 1, 6, 9)
- Serialization + Zlib (levels 1, 6, 9)
- Serialization + Deflate (levels 1, 6, 9)
