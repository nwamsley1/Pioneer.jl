# CLAUDE.md - Pioneer LoggingSystem Technical Documentation

This file provides comprehensive technical guidance for Claude Code when working with the Pioneer logging system.

## System Architecture

The Pioneer logging system implements a three-tier parallel logging architecture:

```
User Code → Logging Macros → WarningCapturingLogger → TeeLogger → [Console, Simple File, Full File]
```

### Core Components

1. **LoggingSystem.jl** - Main orchestration
   - Manages logger lifecycle
   - Stores global logger state in `LOGGER_STATE`
   - Handles initialization and finalization

2. **Loggers.jl** - Logger factories
   - `create_console_logger()` - Interactive console output with verbosity control
   - `create_simplified_logger()` - Clean file output for user consumption
   - `create_full_logger()` - Comprehensive debug output with all metadata

3. **Macros.jl** - Custom logging macros
   - User-facing: `@user_info`, `@user_warn`, `@user_error`
   - Debug levels: `@debug_l1`, `@debug_l2`, `@debug_l3`, `@trace`
   - Progress tracking: `@progress_start`, `@progress_update`, `@progress_end`

4. **WarningTracker.jl** & **WarningCapturingLogger.jl**
   - Captures and categorizes warnings
   - Tracks warning counts and patterns
   - Generates warning reports at session end

5. **Configuration.jl** - Config management
   - `LoggingConfig` struct with all settings
   - JSON-based configuration (NOT environment variables)
   - Runtime adjustable settings

## Critical Implementation Details

### Log Levels (Integer Values)

```julia
# Standard Julia levels
Error = 2000
Warn = 1000  
Info = 0
Debug = -1000

# Custom Pioneer levels
USER_INFO_LEVEL = -100    # Between Info and Debug
DEBUG_L1_LEVEL = -1000    # Standard debug
DEBUG_L2_LEVEL = -1100    # More detailed
DEBUG_L3_LEVEL = -1200    # Very detailed
TRACE_LEVEL = -2000       # Everything

# Console verbosity mappings
:silent => 2001    # Only critical errors
:minimal => 1000   # Warnings and above
:normal => -100    # User info and above
:verbose => -500   # More details
:debug => -1000    # All debug messages
```

### Load Order Issues ⚠️

**CRITICAL**: The logging system is loaded AFTER many other modules due to importScripts.jl ordering:

1. ParseInputs (line 182) loads BEFORE LoggingSystem (line 275)
2. This means files in ParseInputs CANNOT use logging macros
3. loadSpectralLibrary.jl cannot use `@user_info` - must use println or nothing

**Current load order in importScripts.jl:**
```
1. Core utilities
2. ParseInputs (includes loadSpectralLibrary.jl)
3. ... many other modules ...
4. LoggingSystem (at the end)
```

### The TeeLogger Pattern

The system uses LoggingExtras.TeeLogger to broadcast messages to multiple loggers:

```julia
TeeLogger(
    console_logger,      # User-facing console output
    simplified_logger,   # Clean file log
    full_logger         # Debug file log
)
```

Each logger can have different:
- Minimum levels (MinLevelLogger wrapper)
- Formatting (TransformerLogger wrapper)
- Filtering (EarlyFilteredLogger wrapper)

### TransformerLogger Usage

TransformerLoggers modify log records in-flight:

```julia
TransformerLogger(base_logger) do log
    # Modify log record
    clean_kwargs = filter(k -> k ∉ [:internal_metadata], log.kwargs)
    return merge(log, (kwargs = clean_kwargs,))
end
```

## Known Issues and Solutions

### 1. The Hanging Issue (Resolved)

**Problem**: Complex logging setup caused system to hang at "Loading Spectral Library..."

**Cause**: Unknown interaction between multiple TransformerLogger layers or file I/O

**Solution**: Temporarily simplified to just console logger for testing

### 2. Console Verbosity Not Working

**Problem**: All messages printed regardless of verbosity setting

**Cause**: Direct println statements bypassing logging system

**Files with println to check**:
- loadSpectralLibrary.jl (can't be fixed due to load order)
- QuadTuningSearch/utils.jl (debug statements - removed)

### 3. Metadata Leaking to Console

**Problem**: Console showed "is_user_message = true" and other metadata

**Solution**: Clean kwargs in TransformerLogger before display

## Integration with SearchDIA

SearchDIA.jl initializes logging at the start:

```julia
# Get config from params.logging
config = LoggingConfig(
    console_level = Symbol(params.logging.console_level),
    ...
)

# Currently using simplified setup
console_logger = create_console_logger(config)
global_logger(console_logger)

# Legacy compatibility functions
function dual_println(args...)
    msg = join(string.(args))
    @user_info msg
end
```

## Testing the Logging System

### Check Console Filtering
```julia
# Should show at :normal level
@user_info "This should appear"
@user_warn "This warning too"

# Should NOT show at :normal level  
@debug_l1 "This should be hidden"
```

### Check File Output
- Simple log: `results_dir/pioneer_search_log.txt`
- Full log: `results_dir/pioneer_debug.log`
- Warnings: `results_dir/warnings.txt`

### Verify Warning Tracking
```julia
@user_warn "Test warning" category=:test
# Should appear in warnings.txt after finalize_logging()
```

## Best Practices for Modifications

1. **Always test with actual SearchDIA runs** - Unit tests don't catch integration issues

2. **Be careful with TransformerLogger chains** - Too many layers can cause issues

3. **Remember load order** - Not all files can use logging macros

4. **Keep console logger simple** - Complex formatting can interfere with progress bars

5. **Test all verbosity levels** - Ensure filtering works correctly

## Configuration Through JSON

The system reads from `params.logging`:

```json
{
    "logging": {
        "console_level": "normal",    // :silent, :minimal, :normal, :verbose, :debug
        "enable_progress": true,       // Progress bar integration
        "enable_warnings": true,       // Warning tracking
        "rotation_size_mb": 100.0,     // Log rotation threshold
        "max_warnings": 10000          // Max warnings to track
    }
}
```

## Memory and Performance Considerations

1. **File I/O buffering** - FormatLogger handles buffering automatically
2. **Warning accumulation** - Limited by max_warnings setting
3. **Log rotation** - Automatic rotation when size exceeds threshold
4. **Thread safety** - Logging is thread-safe but may cause contention

## Future Improvements to Consider

1. **Fix load order** - Move LoggingSystem earlier in importScripts.jl
2. **Async file writing** - Reduce I/O blocking
3. **Structured logging** - JSON output for machine parsing
4. **Log aggregation** - Combine logs from parallel processes
5. **Performance metrics** - Track logging overhead

## Debugging Tips

1. **If logging causes hangs**: Simplify to ConsoleLogger only
2. **If messages don't appear**: Check log level calculations
3. **If metadata leaks**: Review TransformerLogger chains
4. **If files aren't created**: Check directory permissions and mkpath calls
5. **If load errors occur**: Remember macro availability based on load order

## Related Files

- `/src/importScripts.jl` - Controls load order (LINE 275-286 for LoggingSystem)
- `/src/Routines/SearchDIA.jl` - Main integration point (LINE 128-148)
- `/src/Routines/SearchDIA/ParseInputs/parseParams.jl` - Parses logging config
- Example config: `/data/ecoli_test/ecoli_test_params.json`