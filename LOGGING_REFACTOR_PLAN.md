
# Pioneer.jl Logging System Refactor Plan

## Executive Summary
Refactor the debug logging system to:
1. Remove line numbers from `@debug_l1` output (both console and file)
2. Keep line numbers for `@debug_l2` and `@debug_l3`
3. Change file writing behavior: only write debug messages to file when their level ≤ DEBUG_CONSOLE_LEVEL
4. Ensure consistency between console and file output based on debug level

## Current Behavior Analysis

### Current Implementation
- **Location**: `/src/Pioneer.jl` lines 250-350 (functions), 385-420 (macros)
- **Debug Level Setting**: Set from `params.logging.debug_console_level` in SearchDIA.jl (line 137)
- **Default Level**: 0 (no debug output)

### Current Behavior
1. **Console Output**: 
   - Shows debug messages when `DEBUG_CONSOLE_LEVEL[] >= message_level`
   - All debug levels show line numbers on console
   
2. **File Output**:
   - ALL debug messages are ALWAYS written to debug file regardless of DEBUG_CONSOLE_LEVEL
   - All debug levels include line numbers in file

3. **Line Number Display**:
   - Console: Shows as `└ @ Module file:line` for all debug levels
   - File: Shows as ` @ Module file:line` appended to message for all debug levels

## Proposed Changes

### 1. Update Debug Functions

#### A. `debug_l1` Function (lines 252-275)
**Changes Required**:
- Remove line number display from both console and file output
- Only write to file when `DEBUG_CONSOLE_LEVEL[] >= 1`

```julia
function debug_l1(msg::String, file::String="", line::String="", mod::String="")
    # Only process if debug level allows
    if DEBUG_CONSOLE_LEVEL[] >= 1
        # Console output WITHOUT line numbers
        printstyled("┌ ", color=:blue)
        printstyled("Debug:", bold=true, color=:blue)
        println(" ", msg)
        # NO LINE NUMBER OUTPUT
        
        # File output WITHOUT line numbers
        if DEBUG_FILE[] !== nothing
            debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
            println(DEBUG_FILE[], "[$debug_timestamp] [DEBUG1] $msg")
            flush(DEBUG_FILE[])
        end
    end
    # If debug level < 1, no output at all
end
```

#### B. `debug_l2` Function (lines 277-300)
**Changes Required**:
- Keep line number display for both console and file
- Only write to file when `DEBUG_CONSOLE_LEVEL[] >= 2`

```julia
function debug_l2(msg::String, file::String="", line::String="", mod::String="")
    # Only process if debug level allows
    if DEBUG_CONSOLE_LEVEL[] >= 2
        # Console output WITH line numbers
        printstyled("┌ ", color=:blue)
        printstyled("Debug:", bold=true, color=:blue)
        println(" ", msg)
        
        if !isempty(file) && !isempty(line)
            printstyled("└ ", color=:blue)
            println("@ $mod $file:$line")
        end
        
        # File output WITH line numbers
        if DEBUG_FILE[] !== nothing
            debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
            source_loc = ""
            if !isempty(file) && !isempty(line)
                source_loc = " @ $mod $file:$line"
            end
            println(DEBUG_FILE[], "[$debug_timestamp] [DEBUG2] $msg$source_loc")
            flush(DEBUG_FILE[])
        end
    end
    # If debug level < 2, no output at all
end
```

#### C. `debug_l3` Function (lines 302-325)
**Changes Required**:
- Keep line number display for both console and file
- Only write to file when `DEBUG_CONSOLE_LEVEL[] >= 3`

```julia
function debug_l3(msg::String, file::String="", line::String="", mod::String="")
    # Only process if debug level allows
    if DEBUG_CONSOLE_LEVEL[] >= 3
        # Console output WITH line numbers
        printstyled("┌ ", color=:blue)
        printstyled("Debug:", bold=true, color=:blue)
        println(" ", msg)
        
        if !isempty(file) && !isempty(line)
            printstyled("└ ", color=:blue)
            println("@ $mod $file:$line")
        end
        
        # File output WITH line numbers
        if DEBUG_FILE[] !== nothing
            debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
            source_loc = ""
            if !isempty(file) && !isempty(line)
                source_loc = " @ $mod $file:$line"
            end
            println(DEBUG_FILE[], "[$debug_timestamp] [DEBUG3] $msg$source_loc")
            flush(DEBUG_FILE[])
        end
    end
    # If debug level < 3, no output at all
end
```

#### D. `trace_msg` Function (lines 327-350)
**Changes Required**:
- Only write to file when `DEBUG_CONSOLE_LEVEL[] >= 4`
- Keep line numbers (as trace is for detailed debugging)

```julia
function trace_msg(msg::String, file::String="", line::String="", mod::String="")
    # Only process if debug level allows (4+)
    if DEBUG_CONSOLE_LEVEL[] >= 4
        # Console output WITH line numbers
        printstyled("┌ ", color=:blue)
        printstyled("Debug:", bold=true, color=:blue)
        println(" ", msg)
        
        if !isempty(file) && !isempty(line)
            printstyled("└ ", color=:blue)
            println("@ $mod $file:$line")
        end
        
        # File output WITH line numbers
        if DEBUG_FILE[] !== nothing
            debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
            source_loc = ""
            if !isempty(file) && !isempty(line)
                source_loc = " @ $mod $file:$line"
            end
            println(DEBUG_FILE[], "[$debug_timestamp] [TRACE] $msg$source_loc")
            flush(DEBUG_FILE[])
        end
    end
    # If debug level < 4, no output at all
end
```

### 2. No Changes Required for Macros
The macros (@debug_l1, @debug_l2, @debug_l3, @trace) don't need changes as they just pass through the source location information. The functions will decide whether to use it.

### 3. Update Documentation

#### Update CLAUDE.md Logging Section
```markdown
#### Configuration
JSON parameters control debug verbosity:
```json
"logging": {
    "debug_console_level": 0  // 0-4, controls which debug levels show
}
```

- `0` (default): No debug messages (console or file)
- `1`: Show DEBUG1 messages (console and file)
- `2`: Show DEBUG1 and DEBUG2 messages (console and file)
- `3`: Show DEBUG1, DEBUG2, and DEBUG3 messages (console and file)
- `4+`: Show all debug messages including TRACE (console and file)

**Important**: Debug messages are only written to files when their level ≤ debug_console_level

#### Line Number Behavior
- `@debug_l1`: No line numbers (simplified output)
- `@debug_l2`: Includes line numbers (for detailed debugging)
- `@debug_l3`: Includes line numbers (for detailed debugging)
- `@trace`: Includes line numbers (for trace-level debugging)
```

## Implementation Steps

1. **Backup Current Code**
   - Save current Pioneer.jl logging functions

2. **Update Functions** (in order)
   - Modify `debug_l1` function (remove line numbers, add level check for file)
   - Modify `debug_l2` function (add level check for file)
   - Modify `debug_l3` function (add level check for file)
   - Modify `trace_msg` function (add level check for file)

3. **Test Each Level**
   - Test with `debug_console_level: 0` - Verify NO debug output
   - Test with `debug_console_level: 1` - Verify only DEBUG1 (no line numbers)
   - Test with `debug_console_level: 2` - Verify DEBUG1 (no lines) + DEBUG2 (with lines)
   - Test with `debug_console_level: 3` - Verify DEBUG1-3 with appropriate line numbers
   - Test with `debug_console_level: 4` - Verify all debug including TRACE

4. **Update Documentation**
   - Update CLAUDE.md with new behavior
   - Add examples of each debug level output

## Testing Plan

### Test Cases

1. **Level 0 (Default)**
   - No debug messages in console
   - No debug messages in debug file
   - Only @user_info, @user_warn, @user_error appear

2. **Level 1**
   - @debug_l1 messages appear in console (no line numbers)
   - @debug_l1 messages appear in file (no line numbers)
   - @debug_l2, @debug_l3, @trace do not appear anywhere

3. **Level 2**
   - @debug_l1 messages appear (no line numbers)
   - @debug_l2 messages appear (with line numbers)
   - @debug_l3, @trace do not appear anywhere

4. **Level 3**
   - @debug_l1 messages appear (no line numbers)
   - @debug_l2, @debug_l3 messages appear (with line numbers)
   - @trace does not appear

5. **Level 4+**
   - All debug messages appear with appropriate line number behavior

### Test Script
```julia
# Test file: test_logging.jl
using Pioneer

function test_logging(level::Int)
    # Set debug level
    Pioneer.DEBUG_CONSOLE_LEVEL[] = level
    
    println("\n=== Testing with debug_console_level = $level ===")
    
    @user_info "This is an info message"
    @debug_l1 "This is a DEBUG1 message (should have no line numbers)"
    @debug_l2 "This is a DEBUG2 message (should have line numbers)"
    @debug_l3 "This is a DEBUG3 message (should have line numbers)"
    @trace "This is a TRACE message (should have line numbers)"
end

# Test all levels
for level in 0:4
    test_logging(level)
end
```

## Benefits of This Change

1. **Cleaner Output**: DEBUG1 messages will be cleaner without line numbers
2. **Reduced Log Size**: Debug file won't contain messages above the configured level
3. **Consistency**: Console and file outputs will be synchronized
4. **Performance**: Less I/O when debug level is low (no unnecessary file writes)
5. **Flexibility**: Users can control verbosity for both console AND file

## Risks and Mitigation

### Risk 1: Breaking Existing Debug Workflows
- **Mitigation**: Document changes clearly, provide migration guide

### Risk 2: Loss of Debug Information
- **Mitigation**: Users can increase debug level if they need more detail

### Risk 3: Confusion About Line Numbers
- **Mitigation**: Clear documentation about which levels include line numbers

## Timeline

- **Phase 1** (30 min): Update debug functions
- **Phase 2** (20 min): Test all debug levels
- **Phase 3** (10 min): Update documentation
- **Total**: ~1 hour

## Success Criteria

1. ✅ @debug_l1 shows no line numbers (console and file)
2. ✅ @debug_l2 and @debug_l3 show line numbers
3. ✅ Debug messages only written to file when level ≤ DEBUG_CONSOLE_LEVEL
4. ✅ All tests pass with different debug levels
5. ✅ Documentation updated to reflect new behavior