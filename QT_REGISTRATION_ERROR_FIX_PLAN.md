# Qt Registration Error Fix Plan for Pioneer.jl Linux Binaries

## Problem Description

When running Pioneer.jl compiled binaries on Linux (downloaded as .zip), the application crashes with:
```
Cannot add multiple registrations for QtQml.Models
[3702804] signal 6 (-6): Aborted
```

This error **only** occurs with:
- ✅ Linux binary .zip downloads
- ❌ Linux .deb installer (works fine)
- ❌ Windows installer/binaries (works fine)  
- ❌ macOS installer/binaries (works fine)

## Root Cause Analysis

### Stack Trace Analysis
The error originates from:
1. `Qt6Declarative_jll` initialization during Julia startup
2. Attempting to register QML modules that are already registered
3. Qt6Core's `qAbort()` being called due to duplicate registration

### Why This Happens
1. **GR.jl Dependency Chain**: 
   - Pioneer.jl → Plots/StatsPlots → GR.jl → Qt6 libraries (Qt6Base_jll, Qt6Declarative_jll)
   
2. **PackageCompiler Bundling**:
   - PackageCompiler bundles all JLL artifacts into the compiled application
   - Qt libraries get loaded multiple times during initialization
   - QML module registration happens multiple times, causing the fatal error

3. **Platform Differences**:
   - The .deb package properly isolates libraries with correct paths
   - The raw .zip binary doesn't handle Qt library loading correctly
   - Windows/macOS have different dynamic library loading mechanisms that avoid this issue

## Proposed Solutions

### Solution 1: Lazy GR Initialization (RECOMMENDED)
**Approach**: Delay GR/Plots initialization to prevent Qt conflicts at startup

**Implementation**:
```julia
# src/Pioneer.jl modifications
module Pioneer
# Remove from top-level imports:
# using Plots, StatsPlots

# Add lazy loading function:
function load_plotting_backend()
    @eval using Plots
    @eval using StatsPlots
    ENV["GKSwstype"] = "100"  # Use non-interactive backend
    ENV["PLOTS_DEFAULT_BACKEND"] = "GR"
end

# Modify __init__():
function __init__()
    # Set environment to prevent Qt auto-initialization
    ENV["QT_QPA_PLATFORM"] = "offscreen"
    ENV["GKSwstype"] = "100"
    # Don't load GR immediately
end
```

**Pros**: 
- Minimal code changes
- Maintains all functionality
- Works across all platforms

**Cons**:
- Slight delay when first plot is created
- Need to ensure plotting functions call `load_plotting_backend()`

### Solution 2: Environment Variable Workaround
**Approach**: Set environment variables to control Qt loading behavior

**Implementation**:
1. Update Linux wrapper script (`src/build/CLI/pioneer`):
```bash
#!/bin/bash
# Set Qt environment to prevent duplicate registration
export QT_PLUGIN_PATH=""
export QT_QPA_PLATFORM="offscreen"
export GKSwstype="100"
export LD_LIBRARY_PATH="$SCRIPT_DIR/lib:$LD_LIBRARY_PATH"

# Rest of script...
```

2. Add to GitHub Actions (`build_app_linux.yml`):
```yaml
- name: Add wrapper scripts with environment setup
  shell: bash
  run: |
    cat << 'EOF' > build/Pioneer_${{ matrix.identifier }}/Applications/Pioneer/pioneer
    #!/bin/bash
    export QT_PLUGIN_PATH=""
    export QT_QPA_PLATFORM="offscreen"
    export GKSwstype="100"
    SCRIPT_DIR="$(cd -P "$(dirname "$0")" && pwd)"
    export LD_LIBRARY_PATH="$SCRIPT_DIR/lib:$LD_LIBRARY_PATH"
    exec "$SCRIPT_DIR/bin/$SUBCOMMAND" "${SUBCOMMAND_ARGS[@]}"
    EOF
    chmod +x build/Pioneer_${{ matrix.identifier }}/Applications/Pioneer/pioneer
```

**Pros**: 
- No Julia code changes needed
- Quick to implement

**Cons**:
- Linux-specific fix
- May affect plotting output quality

### Solution 3: Switch Plotting Backend
**Approach**: Replace GR with a backend that doesn't use Qt

**Options**:
1. **PlotlyJS**: Web-based, no Qt dependencies
2. **UnicodePlots**: Terminal-based, minimal dependencies
3. **PGFPlotsX**: LaTeX-based, no Qt

**Implementation**:
```julia
# Project.toml - Replace GR with:
PlotlyJS = "f0f68f2c-4968-5e81-91da-67840de0976a"

# src/Pioneer.jl
ENV["PLOTS_DEFAULT_BACKEND"] = "PlotlyJS"
```

**Pros**: 
- Completely eliminates Qt issues
- May reduce binary size

**Cons**:
- Different plot appearance
- May require code changes for plot customization

### Solution 4: Improve Compilation Strategy
**Approach**: Modify PackageCompiler settings to handle Qt properly

**Implementation in `build_app_linux.yml`**:
```julia
create_app(".", "build/Pioneer_${{ matrix.identifier }}/Applications/Pioneer/",
    incremental=false, 
    force=true,
    filter_stdlibs=true,  # Add this
    include_lazy_artifacts=false,  # Add this
    cpu_target="generic;sandybridge,-xsaveopt,clone_all",
    executables=[...],
    precompile_execution_file="src/build/snoop.jl"
)
```

**Additional step**: Deduplicate Qt libraries
```bash
- name: Deduplicate Qt libraries
  run: |
    QT_DIR="build/Pioneer_${{ matrix.identifier }}/Applications/Pioneer/share/julia/artifacts"
    # Find and remove duplicate Qt libraries
    find "$QT_DIR" -name "libQt6*.so*" -type f | while read -r file; do
      basename=$(basename "$file")
      # Keep only first occurrence
      find "$QT_DIR" -name "$basename" -type f | tail -n +2 | xargs rm -f
    done
```

**Pros**: 
- Addresses root cause
- Cleaner package structure

**Cons**:
- Requires testing to ensure nothing breaks
- May affect other JLL packages

### Solution 5: Package Structure Fix
**Approach**: Restructure how Qt libraries are packaged

**Implementation**:
1. Create separate Qt directory:
```bash
- name: Reorganize Qt libraries
  run: |
    BUILD_DIR="build/Pioneer_${{ matrix.identifier }}/Applications/Pioneer"
    mkdir -p "$BUILD_DIR/qt6"
    
    # Move Qt libraries to dedicated directory
    find "$BUILD_DIR" -name "libQt6*.so*" -exec mv {} "$BUILD_DIR/qt6/" \;
    
    # Fix RPATH with patchelf
    find "$BUILD_DIR/bin" -type f -executable | while read -r exe; do
      patchelf --set-rpath '$ORIGIN/../qt6:$ORIGIN/../lib' "$exe"
    done
```

2. Update wrapper to set Qt paths:
```bash
export LD_LIBRARY_PATH="$SCRIPT_DIR/qt6:$LD_LIBRARY_PATH"
```

**Pros**: 
- Clean separation of Qt libraries
- Prevents conflicts with system Qt

**Cons**:
- Requires `patchelf` tool
- More complex packaging

## Testing Strategy

1. **Local Testing**:
   ```bash
   # Test with problematic environment
   unset QT_PLUGIN_PATH
   ./bin/SearchDIA test_params.json
   ```

2. **CI Testing**:
   - Add test step in GitHub Actions after packaging
   - Download and test .zip on Ubuntu container

3. **Verification Matrix**:
   - [ ] Linux .zip binary works
   - [ ] Linux .deb installer still works
   - [ ] Plotting functionality intact
   - [ ] No performance regression

## Recommended Implementation Order

1. **Immediate Fix**: Implement Solution 2 (Environment Variables)
   - Quick workaround for users
   - Can be deployed immediately

2. **Short Term**: Implement Solution 1 (Lazy Loading)
   - Proper fix that maintains functionality
   - Test thoroughly before release

3. **Long Term**: Consider Solution 4 or 5
   - Clean up packaging process
   - Optimize binary size and structure

## Alternative Approaches

If the above solutions don't work:

1. **Conditional Compilation**: Build separate binaries with/without plotting
2. **Dynamic Loading**: Load GR.jl only when needed via `@require`
3. **Container Solution**: Provide Docker/Singularity images instead of raw binaries
4. **Static Linking**: Investigate static linking of Qt libraries

## Resources and References

- [PackageCompiler.jl Issue #381](https://github.com/JuliaLang/PackageCompiler.jl/issues/381) - Duplicate Artifacts
- [GR.jl Issue #314](https://github.com/jheinen/GR.jl/issues/314) - BinaryBuilder support
- [Qt Documentation on QML Registration](https://www.qt.io/blog/qml-type-registration-in-qt-5.15)
- [Julia Discourse: GR_jll errors](https://discourse.julialang.org/t/standalone-packagecompiler-app-with-plots-fails-with-gr-jll-error/81507)

## Notes for Implementation

- The `__init__()` function already has a comment about preventing Qt registration conflicts
- GR is only used for plotting/visualization, not core functionality
- Consider making plotting optional or plugin-based in future versions
- The .deb package works because it uses system paths and proper library isolation

## Contact

For questions or issues with implementation, contact:
- Maintainer: n.t.wamsley@wustl.edu
- GitHub Issues: [Pioneer.jl repository]