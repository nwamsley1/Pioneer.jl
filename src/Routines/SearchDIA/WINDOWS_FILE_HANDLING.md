# Windows File Handling in SearchDIA

## Overview

Windows has more restrictive file locking compared to Unix-based systems (Linux/macOS). This document describes the specific challenges and solutions implemented in Pioneer.jl for cross-platform compatibility.

## Known Issues

### 1. File Locking
- Windows prevents deletion of files that have open handles
- Arrow files can maintain file handles even after operations complete
- Julia's garbage collector may not immediately release file handles

### 2. Permission Errors
Common error: `IOError: unlink("file.arrow"): permission denied (EACCES)`

This occurs when:
- Another process has the file open
- Antivirus software is scanning the file
- The file handle hasn't been released by Julia's GC
- The file is memory-mapped

### 3. Antivirus Interference
- Real-time scanning can temporarily lock files
- Large Arrow files may trigger extended scans
- File operations may fail intermittently

## Solutions Implemented

### 1. Safe File Removal (`safeRm`)
Located in `src/utils/safeFileOps.jl`

Features:
- Multiple retry attempts with exponential backoff
- Forces garbage collection before each attempt
- Fallback to Windows `del` command
- Final fallback: rename file instead of delete
- Cross-platform: simple `rm()` on Unix

Usage:
```julia
# Replace this:
rm(filepath)
rm(filepath, force=true)

# With this:
safeRm(filepath)
safeRm(filepath, force=true)
```

### 2. Safe Arrow Writing (`writeArrow`)
Located in `src/utils/writeArrow.jl`

Features:
- Writes to temporary file first
- Handles existing file removal with retries
- Uses move operation to minimize lock time
- Fallback to copy+delete if move fails

Usage:
```julia
# Replace this:
Arrow.write(filepath, df)

# With this:
writeArrow(filepath, df)
```

## Best Practices

### 1. Always Use Safe Functions
- Use `safeRm()` instead of `rm()` for all file deletions
- Use `writeArrow()` instead of `Arrow.write()`
- These functions are no-ops on Unix but essential on Windows

### 2. File Handle Management
```julia
# Bad: File handle may persist
table = Arrow.Table(filepath)
# ... use table ...
rm(filepath)  # May fail on Windows

# Good: Explicit handle management
table = Arrow.Table(filepath)
# ... use table ...
table = nothing  # Release reference
GC.gc()         # Force garbage collection
safeRm(filepath)  # Now safe to remove
```

### 3. Antivirus Configuration
For production use, consider:
- Adding output directories to antivirus exclusions
- Disabling real-time scanning for Arrow file extensions
- Using a dedicated output drive with minimal scanning

### 4. Temporary Files
```julia
# Use unique temporary names to avoid conflicts
temp_path = tempname() * ".arrow"
# ... write to temp_path ...
mv(temp_path, final_path, force=true)
```

## Error Handling

### Permission Denied Errors
If you encounter permission errors:

1. Check if antivirus is scanning the file
2. Ensure no other Julia processes have the file open
3. Try manual GC: `GC.gc()`
4. Use Process Explorer (Windows) to find file handles

### Debugging File Locks
```julia
# Enable verbose mode to see retry attempts
ENV["PIONEER_DEBUG_FILES"] = "1"

# Check if file is locked
try
    open(filepath, "r+") do io
        # File is not locked
    end
catch e
    @warn "File is locked" filepath error=e
end
```

## Performance Considerations

### 1. File Operations are Slower on Windows
- Expect 10-30% slower I/O operations
- Plan for longer runtimes with many files
- Consider batching operations

### 2. Memory Usage
- Windows file caching is less efficient
- May need to increase Julia heap size
- Monitor memory usage with Task Manager

### 3. Concurrent Access
- Avoid multiple processes writing to same directory
- Use file-based locking if needed
- Consider process synchronization

## Platform-Specific Code

Always check platform when implementing file operations:

```julia
if Sys.iswindows()
    # Windows-specific handling
    safeRm(filepath)
else
    # Unix handling
    rm(filepath, force=true)
end
```

## Troubleshooting Checklist

1. ✓ Using `safeRm()` instead of `rm()`?
2. ✓ Using `writeArrow()` instead of `Arrow.write()`?
3. ✓ Releasing file handles before deletion?
4. ✓ Antivirus exclusions configured?
5. ✓ Sufficient disk space?
6. ✓ No other processes accessing files?
7. ✓ Temporary directory accessible?

## Future Improvements

1. Implement file-based locking mechanism
2. Add automatic retry for Arrow.append operations
3. Create Windows-specific performance benchmarks
4. Investigate memory-mapped file alternatives