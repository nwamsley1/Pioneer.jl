# Solution 1: Lazy GR Initialization - Detailed Implementation Plan

## Overview
Implement lazy loading of GR/Plots to prevent Qt6Declarative_jll registration conflicts at startup. This solution delays the loading of plotting libraries until they are actually needed, avoiding the Qt initialization conflict that occurs during binary startup on Linux.

## Branch Strategy
- Branch name: `fix/qt-registration-lazy-gr-loading`
- Base branch: Current working directory (appears to be a release tag)
- Commits will be atomic and well-documented

## Implementation Steps

### Step 1: Create and Push Branch
```bash
git checkout -b fix/qt-registration-lazy-gr-loading
git push -u origin fix/qt-registration-lazy-gr-loading
```

### Step 2: Modify src/Pioneer.jl

#### 2.1 Remove Direct Plotting Imports
**Current code (lines to modify):**
```julia
# Line ~36
using Plots, Polynomials, ProgressBars, Printf
# Line ~38
using StatsPlots, SentinelArrays
```

**Change to:**
```julia
# Line ~36
using Polynomials, ProgressBars, Printf
# Line ~38
using SentinelArrays
# Note: Plots and StatsPlots will be loaded lazily
```

#### 2.2 Add Lazy Loading Infrastructure
**Add after the export statement (before module end):**
```julia
# Lazy loading state
const PLOTTING_LOADED = Ref(false)

"""
    ensure_plotting_loaded()

Lazily loads Plots and StatsPlots packages with GR backend.
This prevents Qt registration conflicts at startup on Linux.
"""
function ensure_plotting_loaded()
    if !PLOTTING_LOADED[]
        @info "Loading plotting libraries..."
        
        # Set environment to prevent Qt conflicts
        ENV["GKSwstype"] = "100"  # Use non-interactive backend
        ENV["QT_QPA_PLATFORM"] = "offscreen"  # Prevent Qt GUI initialization
        
        # Load plotting packages
        @eval Main using Plots
        @eval Main using StatsPlots
        
        # Import into Pioneer module
        @eval using Plots
        @eval using StatsPlots
        
        # Ensure GR backend is set
        ENV["PLOTS_DEFAULT_BACKEND"] = "GR"
        
        PLOTTING_LOADED[] = true
        @info "Plotting libraries loaded successfully"
    end
end

# Re-export plotting functions with lazy loading
for func in [:plot, :scatter, :histogram, :heatmap, :savefig, :plot!, :scatter!, 
             :histogram!, :violin, :violin!, :boxplot, :boxplot!, :gr]
    @eval begin
        function $func(args...; kwargs...)
            ensure_plotting_loaded()
            Plots.$func(args...; kwargs...)
        end
        export $func
    end
end

# Handle StatsPlots specific functions
for func in [:groupedbar, :groupedbar!, :groupedhist, :groupedhist!]
    @eval begin
        function $func(args...; kwargs...)
            ensure_plotting_loaded()
            StatsPlots.$func(args...; kwargs...)
        end
        export $func
    end
end
```

#### 2.3 Modify __init__ Function
**Current code:**
```julia
function __init__()
    # Don't initialize gr() immediately - let it be initialized when first used
    # This prevents Qt registration conflicts
    ENV["PLOTS_DEFAULT_BACKEND"] = "GR"
end
```

**Change to:**
```julia
function __init__()
    # Prevent Qt registration conflicts by setting environment early
    ENV["GKSwstype"] = "100"  # Non-interactive GR backend
    ENV["QT_QPA_PLATFORM"] = "offscreen"  # No Qt GUI
    # Note: PLOTS_DEFAULT_BACKEND will be set when plotting is loaded
    
    # Ensure we're not in a conflicting Qt environment
    if haskey(ENV, "QT_PLUGIN_PATH")
        @warn "QT_PLUGIN_PATH is set, this may cause conflicts" QT_PLUGIN_PATH=ENV["QT_PLUGIN_PATH"]
    end
end
```

### Step 3: Find and Update Plotting Usage

#### 3.1 Search for Plotting Functions
Files to check for plotting usage:
- `src/Routines/SearchDIA/WriteOutputs/qcPlots.jl`
- `src/Routines/SearchDIA/WriteOutputs/plotRTAlignment.jl`
- Any file using `plot`, `scatter`, `histogram`, etc.

#### 3.2 Update Each Plotting File
**Add at the beginning of each file that uses plotting:**
```julia
# Ensure plotting is loaded before using
Pioneer.ensure_plotting_loaded()
```

OR, if the functions are called through Pioneer's re-exported functions, no changes needed.

### Step 4: Update Tests

#### 4.1 Modify Test Files
For any test files that use plotting directly:
```julia
using Pioneer
Pioneer.ensure_plotting_loaded()  # Before any plotting tests
```

### Step 5: Create Compatibility Layer

Create new file: `src/utils/PlottingCompat.jl`
```julia
"""
    PlottingCompat.jl
    
Compatibility layer for lazy-loaded plotting functionality.
Ensures backward compatibility while preventing Qt conflicts.
"""
module PlottingCompat

using ..Pioneer: ensure_plotting_loaded

# Create a macro for conditional plotting
macro with_plots(expr)
    quote
        ensure_plotting_loaded()
        $(esc(expr))
    end
end

export @with_plots

end # module
```

Include this in `src/importScripts.jl`:
```julia
include(joinpath(@__DIR__, "utils", "PlottingCompat.jl"))
using .PlottingCompat
```

### Step 6: Handle Edge Cases

#### 6.1 Precompilation Script Updates
Modify `src/build/snoop.jl` to handle lazy loading:
```julia
# Add after using Pioneer
if get(ENV, "PIONEER_PRECOMPILE_PLOTS", "true") == "true"
    Pioneer.ensure_plotting_loaded()
end
```

#### 6.2 Environment Detection
Add platform-specific handling:
```julia
function should_use_offscreen_backend()
    # Always use offscreen on Linux binaries
    if Sys.islinux() && isfile(joinpath(@__DIR__, "..", "..", "bin", "SearchDIA"))
        return true
    end
    # Check if we're in a compiled app
    if ccall(:jl_generating_output, Cint, ()) == 1
        return true
    end
    return false
end
```

### Step 7: Testing Plan

#### 7.1 Local Testing
```bash
# Test 1: Basic functionality without plotting
julia --project=. -e 'using Pioneer; println("No plots: OK")'

# Test 2: Lazy loading
julia --project=. -e 'using Pioneer; Pioneer.ensure_plotting_loaded(); println("Plots loaded: OK")'

# Test 3: Full workflow
julia --project=. test/runtests.jl
```

#### 7.2 Compilation Testing
```julia
using PackageCompiler
# Create test app with new changes
create_app(".", "test_build/",
    executables=["SearchDIA" => "main_SearchDIA"],
    precompile_execution_file="src/build/snoop.jl",
    force=true
)
```

#### 7.3 Linux Binary Testing
```bash
# In Linux environment or Docker
cd test_build/bin
./SearchDIA test_params.json
# Should not crash with Qt error
```

### Step 8: Documentation Updates

#### 8.1 Update CLAUDE.md
Add section about lazy loading:
```markdown
### Plotting and Visualization
Pioneer uses lazy loading for plotting libraries to prevent Qt conflicts on Linux:
- Plots are not loaded at startup
- First plot call triggers loading
- Use `Pioneer.ensure_plotting_loaded()` to manually load if needed
```

#### 8.2 Update README if needed
Add note about plotting being lazily loaded.

## Commit Plan

### Commit 1: Remove direct plotting imports
```bash
git add src/Pioneer.jl
git commit -m "fix: Remove direct Plots/StatsPlots imports to prevent Qt conflicts

- Removed Plots and StatsPlots from top-level using statements
- Preparing for lazy loading implementation
- Part of fix for Linux binary Qt registration error"
```

### Commit 2: Implement lazy loading infrastructure
```bash
git add src/Pioneer.jl
git commit -m "feat: Add lazy loading for plotting libraries

- Added ensure_plotting_loaded() function
- Re-exported plotting functions with lazy loading
- Set appropriate environment variables for Qt conflict prevention
- Fixes Qt6Declarative_jll registration error on Linux binaries"
```

### Commit 3: Update plotting usage files
```bash
git add src/Routines/SearchDIA/WriteOutputs/*.jl
git commit -m "refactor: Update plotting files for lazy loading

- Ensure plotting is loaded before use in QC plots
- Maintain backward compatibility
- All plotting now goes through lazy loading mechanism"
```

### Commit 4: Add compatibility layer
```bash
git add src/utils/PlottingCompat.jl src/importScripts.jl
git commit -m "feat: Add plotting compatibility layer

- Created PlottingCompat module for easier migration
- Added @with_plots macro for conditional plotting
- Improved code organization for plotting functionality"
```

### Commit 5: Update precompilation and tests
```bash
git add src/build/snoop.jl test/
git commit -m "test: Update tests and precompilation for lazy plotting

- Modified snoop.jl to handle lazy loading during precompilation
- Updated test files to ensure plotting is loaded when needed
- Added environment variable control for precompilation"
```

### Commit 6: Documentation
```bash
git add CLAUDE.md README.md SOLUTION_1_LAZY_GR_IMPLEMENTATION_PLAN.md
git commit -m "docs: Document lazy loading solution for Qt conflicts

- Added documentation about lazy plotting loading
- Explained the Qt conflict resolution
- Updated technical documentation"
```

## Rollback Plan

If issues arise:
1. `git checkout main` (or previous branch)
2. `git branch -D fix/qt-registration-lazy-gr-loading`
3. Implement Solution 2 (environment variables) as a quicker fix

## Success Criteria

1. ✅ Linux binary runs without Qt registration error
2. ✅ Plotting still works when needed
3. ✅ No performance regression in non-plotting workflows
4. ✅ All existing tests pass
5. ✅ Compiled binaries work on Ubuntu, Debian, RHEL
6. ✅ .deb installer continues to work as before

## Potential Issues and Mitigations

### Issue 1: Plotting functions not found
**Mitigation**: Ensure all plotting functions are properly re-exported

### Issue 2: Module loading conflicts
**Mitigation**: Use `@eval Main using` then `@eval using` pattern

### Issue 3: Precompilation failures
**Mitigation**: Add conditional precompilation based on environment

### Issue 4: Threading issues with lazy loading
**Mitigation**: Add thread safety with locks if needed:
```julia
const PLOTTING_LOCK = ReentrantLock()
function ensure_plotting_loaded()
    lock(PLOTTING_LOCK) do
        # ... loading code
    end
end
```

## Alternative Approaches Within Solution 1

If the main approach has issues:

1. **Conditional Compilation**: Use `@static if` to conditionally include plotting
2. **Requires.jl Pattern**: Use `@require` for conditional loading
3. **Separate Plotting Module**: Move all plotting to a submodule that's loaded on demand

## Files to Monitor

Critical files that may need updates:
- `src/Pioneer.jl` - Main module file
- `src/Routines/SearchDIA/WriteOutputs/qcPlots.jl` - QC plotting
- `src/Routines/SearchDIA/WriteOutputs/plotRTAlignment.jl` - RT alignment plots
- `src/build/snoop.jl` - Precompilation
- `src/importScripts.jl` - File imports

## Testing Checklist

- [ ] Local Julia REPL works without plotting
- [ ] Local Julia REPL works with plotting
- [ ] Tests pass locally
- [ ] Compilation succeeds
- [ ] Linux binary runs without Qt error
- [ ] Linux binary can create plots
- [ ] Windows binary still works
- [ ] macOS binary still works
- [ ] .deb package still works
- [ ] Performance is not degraded

## Next Steps After Implementation

1. Test on multiple Linux distributions
2. Get user confirmation that the issue is resolved
3. Consider implementing Solution 2 as additional safety
4. Plan for v0.2.3 release with this fix
5. Update GitHub Actions to test for this specific issue