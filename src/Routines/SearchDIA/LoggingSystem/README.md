# Pioneer Logging System - User Guide

## Overview

The Pioneer logging system provides parallel output to console and files with configurable verbosity levels, warning tracking, and progress bar integration.

## For Users: Configuring Logging

### Configuration in JSON Parameters

Add a `logging` section to your parameters JSON file:

```json
{
    "logging": {
        "console_level": "normal",    
        "enable_progress": true,       
        "enable_warnings": true,       
        "rotation_size_mb": 100.0,     
        "max_warnings": 10000          
    }
}
```

### Console Verbosity Levels

- **`"silent"`** - Only critical errors (almost nothing)
- **`"minimal"`** - Warnings and errors only
- **`"normal"`** - Standard operation messages (recommended)
- **`"verbose"`** - Detailed operation information
- **`"debug"`** - Everything including debug traces

### Output Files

**Note**: File logging is currently disabled. All output goes to console only.

When file logging is re-enabled, logs will be in your results directory:

- **`pioneer_search_log.txt`** - Clean, timestamped log of the run
- **`pioneer_debug.log`** - Comprehensive debug information (if issues occur)
- **`warnings.txt`** - Summary of all warnings encountered

### Example Configurations

**For normal operation:**
```json
{
    "logging": {
        "console_level": "normal",
        "enable_progress": true,
        "enable_warnings": true
    }
}
```

**For debugging issues:**
```json
{
    "logging": {
        "console_level": "debug",
        "enable_progress": false,
        "enable_warnings": true
    }
}
```

**For automated pipelines:**
```json
{
    "logging": {
        "console_level": "minimal",
        "enable_progress": false,
        "enable_warnings": true
    }
}
```

## For Developers: Using Logging Macros

### When to Use Each Macro

#### User-Facing Messages

**`@user_info`** - Important operational messages users should see
```julia
@user_info "Loading Spectral Library..."
@user_info "Processing $(n_files) MS files"
@user_info "Search completed in $(elapsed_time) seconds"
```

**`@user_print`** - Decorative output without any prefix (for tables, separators)
```julia
@user_print "=" ^ 80  # Separator line
@user_print "DIA Search Performance Report"  # Header without prefix
@user_print rpad("Step", 30) * " " * rpad("Time", 10)  # Table formatting
```

**`@user_warn`** - Issues users should be aware of but aren't fatal
```julia
@user_warn "Only one MS file detected, disabling match-between-runs"
@user_warn "Low peptide count: $(count)" category=:identification
@user_warn "Missing isotope data for precursor $(prec_id)" category=:quantification
```

**`@user_error`** - Critical issues that prevent operation
```julia
@user_error "MS data directory not found: $(path)"
@user_error "Invalid parameter value: $(param_name) = $(value)"
@user_error "Library loading failed: $(error_message)"
```

#### Debug Messages (Hidden from Users by Default)

**`@debug_l1`** - High-level debugging for developers
```julia
@debug_l1 "Starting parameter tuning with tolerance range: $(range)"
@debug_l1 "Thread $(thread_id) processing batch $(batch_id)"
@debug_l1 "Model training completed with $(n_features) features"
```

**`@debug_l2`** - Detailed debugging information
```julia
@debug_l2 "PSM scoring details" psm_id=id score=score features=features
@debug_l2 "RT calibration point" observed=obs_rt expected=exp_rt
@debug_l2 "Fragment match" mz=mz intensity=intensity error_ppm=error
```

**`@debug_l3`** - Very detailed/verbose debugging
```julia
@debug_l3 "Matrix multiplication step" matrix_size=size result=result
@debug_l3 "Memory allocation" bytes=bytes location=location
@debug_l3 "Lock acquired" thread=thread_id resource=resource
```

**`@trace`** - Extreme detail (rarely used)
```julia
@trace "Entering function" function_name=name args=args
@trace "Loop iteration" i=i value=value
```

### Best Practices

#### DO Use Logging Macros For:
- Status updates about long-running operations
- Warnings about suboptimal conditions
- Errors that affect results
- Debug information for troubleshooting

#### DON'T Use Logging For:
- Temporary debug statements (use `@show` instead)
- High-frequency operations in loops (performance impact)
- Binary data or very large strings
- Sensitive information (passwords, paths with user info)

### Adding Context to Messages

Include relevant variables for debugging:

```julia
# Good - provides context
@user_warn "Peptide coverage below threshold" coverage=cov threshold=min_cov peptide_count=n

# Less helpful
@user_warn "Low peptide coverage"
```

### Categories for Warnings

Use categories to group related warnings:

```julia
@user_warn "Issue with RT alignment" category=:rt_alignment
@user_warn "Issue with mass calibration" category=:mass_calibration
@user_warn "Missing MS1 data" category=:quantification
```

Common categories:
- `:identification` - PSM/peptide identification issues
- `:quantification` - Integration/quantification problems  
- `:rt_alignment` - Retention time calibration
- `:mass_calibration` - Mass accuracy issues
- `:library` - Spectral library problems
- `:input_data` - Issues with input files

### Progress Tracking

For long operations with known iterations:

```julia
@progress_start "Processing files" total_files

for (i, file) in enumerate(files)
    @progress_update "Processing files" i total_files
    # ... process file ...
end

@progress_end "Processing files"
```

### Performance Considerations

1. **Minimize string interpolation in hot loops** - Build message once outside loop
2. **Use appropriate levels** - Don't use @user_info for debug information
3. **Batch related messages** - Combine multiple related warnings into one
4. **Check log level before expensive operations**:

```julia
if get_logger_level() <= DEBUG_L2_LEVEL
    expensive_debug_info = compute_debug_details()
    @debug_l2 "Detailed info" info=expensive_debug_info
end
```

## Important Limitations

### Cannot Use Logging Macros In These Files

Due to module load order, these files load BEFORE the logging system:

- `ParseInputs/loadSpectralLibrary.jl`
- `ParseInputs/paramsChecks.jl` 
- Any file loaded before line 275 in `importScripts.jl`

In these files, either:
- Use `println` sparingly for critical messages
- Avoid output entirely
- Return status/errors for the caller to log

## Examples

### Typical Search Method Implementation

```julia
function performSearch!(method::MySearchMethod)
    @user_info "Starting MySearch phase..."
    
    # Get parameters
    params = getParameters(method)
    @debug_l1 "Parameters loaded" param_count=length(params)
    
    # Process data
    results = []
    @progress_start "Processing scans" n_scans
    
    for (i, scan) in enumerate(scans)
        @progress_update "Processing scans" i n_scans
        
        try
            result = process_scan(scan)
            push!(results, result)
            @debug_l2 "Scan processed" scan_id=scan.id result_score=result.score
            
        catch e
            @user_warn "Failed to process scan $(scan.id)" error=e category=:processing
            continue
        end
    end
    
    @progress_end "Processing scans"
    
    # Validate results
    if length(results) < minimum_results
        @user_warn "Low result count" count=length(results) minimum=minimum_results
    end
    
    @user_info "MySearch completed: $(length(results)) results"
    @debug_l1 "Memory usage" bytes=Base.summarysize(results)
end
```

### Error Handling Pattern

```julia
function load_data(path::String)
    @debug_l1 "Attempting to load data" path=path
    
    if !isfile(path)
        @user_error "Data file not found: $path"
        throw(ArgumentError("File not found"))
    end
    
    try
        data = Arrow.Table(path)
        @debug_l1 "Data loaded successfully" rows=length(data) 
        return data
        
    catch e
        @user_error "Failed to load data file" path=path error=e
        throw(e)
    end
end
```

## Troubleshooting

### Messages Not Appearing

1. Check your console_level setting in JSON
2. Verify you're using the right macro for the level
3. Ensure the file can use logging macros (see Limitations)

### Too Much Output

- Change console_level to "minimal" or "silent"
- Disable progress bars: `"enable_progress": false`
- Check for debug macros that should be removed

### Log Files Not Created

- Verify results directory exists and is writable
- Check that logging system is initialized (SearchDIA does this)
- Look for errors during logger initialization

### Performance Impact

If logging is slowing down execution:
1. Reduce console_level to "minimal"  
2. Remove debug logging from inner loops
3. Disable progress tracking if not needed
4. Check for expensive string interpolation