# Windows GC Crash Fix Plan

## Problem Summary

**Symptom**: Intermittent `EXCEPTION_ACCESS_VIOLATION` in Julia's garbage collector on Windows during parallel file operations.

**Root Cause**: Multiple threads calling `GC.gc()` concurrently during parallel file sorting in MaxLFQ search, creating race conditions in Julia's memory allocator.

**Location**: `src/utils/writeArrow.jl:49` (called from `src/utils/FileOperations/io/ArrowOperations.jl:85` during parallel processing)

**Call Stack**:
```
MaxLFQSearch.jl:183 → sort_file_by_keys!(parallel=true)
  ↓
ArrowOperations.jl:124 → Threads.@threads for ref in refs
  ↓ (multiple threads)
ArrowOperations.jl:85 → writeArrow(file_path(ref), df)
  ↓ (each thread)
writeArrow.jl:49 → GC.gc()  ⚠️ CONCURRENT CALLS CRASH
  ↓
Julia GC → EXCEPTION_ACCESS_VIOLATION
```

---

## Solution Options

### Option 1: Thread-Safe GC with Lock ⭐ (Recommended)

**Approach**: Add a global lock to ensure only one thread calls GC at a time.

**Implementation**:
```julia
# At module level (top of writeArrow.jl)
const GC_LOCK = ReentrantLock()

function writeArrow(fpath::String, df::AbstractDataFrame)
    fpath = normpath(fpath)
    if Sys.iswindows()
        tpath = tempname() * ".arrow"
        Arrow.write(tpath, df)

        if isfile(fpath)
            max_retries = 5
            for i in 1:max_retries
                try
                    # Thread-safe GC call
                    lock(GC_LOCK) do
                        GC.gc()
                    end

                    rm(fpath, force=true)
                    break
                catch e
                    if i == max_retries
                        # ... existing fallback logic
                    else
                        sleep(0.1 * i)
                    end
                end
            end
        end

        mv(tpath, fpath, force=true)
    else
        Arrow.write(fpath, df)
    end
    return nothing
end
```

**Pros**:
- Minimal code change
- Guarantees thread safety
- No change to algorithm logic
- Low overhead (lock only held during GC)

**Cons**:
- Serializes GC calls (could create bottleneck with many threads)
- Doesn't address whether explicit GC is actually needed
- GC is still global (affects all threads even outside lock)

**Performance Impact**:
- Minimal if GC calls are infrequent
- Could become bottleneck if many files need retry logic
- Threads will wait sequentially for GC access

---

### Option 2: Use finalize() on Arrow Objects

**Approach**: Instead of calling GC explicitly, use `finalize()` on specific Arrow objects to release file handles.

**Implementation**:
```julia
function writeArrow(fpath::String, df::AbstractDataFrame)
    fpath = normpath(fpath)
    if Sys.iswindows()
        tpath = tempname() * ".arrow"

        # Write to temp, capturing the Arrow writer object
        Arrow.write(tpath, df)

        if isfile(fpath)
            max_retries = 5
            for i in 1:max_retries
                try
                    # Option 2a: Try to finalize any existing Arrow.Table at fpath
                    # (Would need to track Arrow objects globally - not practical)

                    # Option 2b: Just rely on Julia's finalizers without explicit GC
                    # Add a small delay to let finalizers run
                    sleep(0.05)

                    # Option 2c: Force finalization of the file system
                    # (Not directly possible - finalizers are object-based)

                    rm(fpath, force=true)
                    break
                catch e
                    if i == max_retries
                        # ... existing fallback logic
                    else
                        sleep(0.1 * i)
                    end
                end
            end
        end

        mv(tpath, fpath, force=true)
    else
        Arrow.write(fpath, df)
    end
    return nothing
end
```

**Challenge with finalize() Approach**:

The problem is that we don't have a reference to the Arrow object that's holding the file handle. The file handle issue occurs when:
1. Some *other* code has read the file (e.g., `Arrow.Table(fpath)`)
2. That Arrow.Table object is still in scope somewhere
3. The file handle is still open
4. We try to delete/overwrite the file

**Why finalize() is difficult here**:
- We'd need to track which Arrow.Table objects refer to which files (global state)
- Other parts of codebase might have references we don't know about
- Arrow.jl internally manages file handles - we can't easily finalize them
- `finalize()` only works if you have the object reference

**Possible Enhanced Implementation**:
```julia
# Create a global registry of Arrow objects by file path
const ARROW_FILE_REGISTRY = Dict{String, Vector{WeakRef}}()
const REGISTRY_LOCK = ReentrantLock()

function register_arrow_object(path::String, obj)
    lock(REGISTRY_LOCK) do
        if !haskey(ARROW_FILE_REGISTRY, path)
            ARROW_FILE_REGISTRY[path] = WeakRef[]
        end
        push!(ARROW_FILE_REGISTRY[path], WeakRef(obj))
    end
end

function finalize_arrow_objects_for_path(path::String)
    lock(REGISTRY_LOCK) do
        if haskey(ARROW_FILE_REGISTRY, path)
            for weak_ref in ARROW_FILE_REGISTRY[path]
                obj = weak_ref.value
                if obj !== nothing
                    finalize(obj)
                end
            end
            delete!(ARROW_FILE_REGISTRY, path)
        end
    end
end

function writeArrow(fpath::String, df::AbstractDataFrame)
    fpath = normpath(fpath)
    if Sys.iswindows()
        tpath = tempname() * ".arrow"
        Arrow.write(tpath, df)

        if isfile(fpath)
            max_retries = 5
            for i in 1:max_retries
                try
                    # Try to finalize any registered Arrow objects for this path
                    finalize_arrow_objects_for_path(fpath)
                    sleep(0.05)  # Let finalizers complete

                    rm(fpath, force=true)
                    break
                catch e
                    if i == max_retries
                        # ... existing fallback logic
                    else
                        sleep(0.1 * i)
                    end
                end
            end
        end

        mv(tpath, fpath, force=true)
    else
        Arrow.write(fpath, df)
    end
    return nothing
end
```

**Pros**:
- More targeted than full GC
- Thread-safe (no concurrent GC issues)
- Only affects specific objects

**Cons**:
- Requires global object tracking infrastructure
- Need to modify all Arrow.Table() calls to register objects
- Complex to implement correctly
- May not catch all file handle holders
- Relies on external code cooperating with registry
- WeakRefs add overhead and complexity

**Verdict**: This is **much more complex** than the lock approach and requires invasive changes throughout the codebase.

---

### Option 3: Remove GC Call Entirely

**Approach**: Remove `GC.gc()` and rely on longer retry delays + more aggressive fallbacks.

**Implementation**:
```julia
function writeArrow(fpath::String, df::AbstractDataFrame)
    fpath = normpath(fpath)
    if Sys.iswindows()
        tpath = tempname() * ".arrow"
        Arrow.write(tpath, df)

        if isfile(fpath)
            max_retries = 10  # Increased retries
            for i in 1:max_retries
                try
                    # No GC call - just retry with exponential backoff
                    if i > 1
                        sleep(0.05 * (2^(i-1)))  # Exponential backoff
                    end

                    rm(fpath, force=true)
                    break
                catch e
                    if i == max_retries
                        # Final attempt: rename instead of delete
                        backup_path = fpath * ".backup_" * string(time_ns())
                        try
                            mv(fpath, backup_path, force=true)
                        catch
                            error("Unable to remove or rename existing file: $fpath")
                        end
                    end
                end
            end
        end

        mv(tpath, fpath, force=true)
    else
        Arrow.write(fpath, df)
    end
    return nothing
end
```

**Pros**:
- Completely eliminates the race condition
- Simple implementation
- Natural finalizers may work with longer delays

**Cons**:
- May increase failure rate on Windows (file handle not released)
- Longer total retry time
- Less deterministic (relies on Julia's automatic GC timing)
- Could accumulate backup files if retries fail

**Testing needed**:
- Does Windows file handle issue actually occur without GC?
- What's the minimum sleep time needed?
- How often do retries fail?

---

### Option 4: Conditional Parallelism on Windows

**Approach**: Disable parallel sorting on Windows to avoid concurrent GC calls.

**Implementation in ArrowOperations.jl**:
```julia
function sort_file_by_keys!(refs::Vector{<:FileReference}, keys::Symbol...;
                           reverse::Union{Bool, Vector{Bool}}=false,
                           parallel::Bool=true,
                           show_progress::Bool=true)
    # Validate reverse vector length
    if reverse isa Vector{Bool} && length(reverse) != length(keys)
        error("Length of reverse vector must match number of sort keys")
    end

    # Disable parallelism on Windows to avoid GC race conditions
    use_parallel = parallel && !Sys.iswindows() && length(refs) > 1

    if use_parallel
        if show_progress
            Threads.@threads for ref in refs
                if exists(ref)
                    sort_file_by_keys!(ref, keys...; reverse=reverse, show_progress=false)
                end
            end
        else
            # ... parallel code
        end
    else
        # Sequential processing
        if show_progress
            p = Progress(length(refs), desc="Sorting files...")
        end
        for ref in refs
            if exists(ref)
                sort_file_by_keys!(ref, keys...; reverse=reverse, show_progress=false)
                show_progress && next!(p)
            end
        end
    end

    return refs
end
```

**Pros**:
- Completely avoids the race condition
- No changes to writeArrow needed
- Simple and safe

**Cons**:
- Significant performance penalty on Windows (sequential vs parallel)
- Doesn't fix the underlying issue
- Windows users get degraded performance

**Performance Impact**: With 312 files and 8 threads, could be 4-8x slower for sort operations on Windows.

---

### Option 5: Hybrid Lock + Reduced Frequency

**Approach**: Add lock AND reduce how often GC is called.

**Implementation**:
```julia
const GC_LOCK = ReentrantLock()

function writeArrow(fpath::String, df::AbstractDataFrame)
    fpath = normpath(fpath)
    if Sys.iswindows()
        tpath = tempname() * ".arrow"
        Arrow.write(tpath, df)

        if isfile(fpath)
            max_retries = 5
            for i in 1:max_retries
                try
                    # Only call GC on later retries (not first attempt)
                    if i >= 3
                        lock(GC_LOCK) do
                            GC.gc()
                        end
                    elseif i > 1
                        # Just sleep on early retries
                        sleep(0.05 * i)
                    end

                    rm(fpath, force=true)
                    break
                catch e
                    if i == max_retries
                        # ... existing fallback logic
                    else
                        sleep(0.1 * i)
                    end
                end
            end
        end

        mv(tpath, fpath, force=true)
    else
        Arrow.write(fpath, df)
    end
    return nothing
end
```

**Pros**:
- Combines safety of lock with reduced contention
- Only uses GC when really needed (after initial retries fail)
- Best of both worlds

**Cons**:
- More complex logic
- Still has lock contention (but less frequent)

---

## Recommendation

**Primary Fix: Option 1 (Thread-Safe GC with Lock)**

Reasoning:
1. **Minimal code change** - Single lock addition
2. **Guaranteed safety** - No race conditions possible
3. **No algorithm changes** - Maintains current logic
4. **Low risk** - Lock overhead is minimal for infrequent operations
5. **Easy to test** - Clear before/after behavior

**Secondary Enhancement: Option 5 (Hybrid)**

If lock contention becomes an issue, upgrade to Option 5:
- Only call GC on 3rd+ retry
- Maintains thread safety
- Reduces lock contention

**Reject Options**:
- **Option 2 (finalize)**: Too complex, requires invasive changes
- **Option 3 (remove GC)**: Risk of more Windows failures
- **Option 4 (disable parallel)**: Unacceptable performance penalty

---

## Implementation Steps

### Step 1: Add Lock (Option 1)
```julia
# In src/utils/writeArrow.jl at top of file
const GC_LOCK = ReentrantLock()

# In writeArrow function, line 49:
lock(GC_LOCK) do
    GC.gc()
end
```

### Step 2: Test on Windows
- Run `BuildSpecLib` test multiple times
- Monitor for crashes
- Check performance impact

### Step 3: If Performance Issues, Upgrade to Option 5
- Only call locked GC on 3rd+ retry
- Profile lock contention

### Step 4: Monitor Other GC Calls
Search for other parallel contexts with GC.gc():
```bash
grep -r "GC.gc()" src/ | grep -i "thread\|parallel"
```

Current known locations:
- `src/Routines/BuildSpecLib.jl:220` (inside threaded loop)
- `src/Routines/BuildSpecLib.jl:341` (inside threaded loop)
- `src/utils/FileOperations/streaming/MergeOperations.jl:325` (inside threaded loop)

These may need similar fixes.

---

## Testing Strategy

### 1. Reproduction Test (Windows)
```julia
# Create test that triggers parallel sorting
using Pioneer
include("test/UnitTests/BuildSpecLib/test_build_spec_lib.jl")

# Run 10 times to catch intermittent crash
for i in 1:10
    @info "Iteration $i"
    test_build_spec_lib()
end
```

### 2. Performance Benchmark
```julia
# Compare with/without lock
@time sort_file_by_keys!(refs, :score; parallel=true)
```

### 3. Stress Test
```julia
# Many threads, many files
Threads.@threads for i in 1:1000
    writeArrow("test_$i.arrow", df)
end
```

### 4. Validate Other Platforms
- Ensure lock doesn't break macOS/Linux
- Verify no performance regression on non-Windows

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Lock contention bottleneck | Low | Medium | Upgrade to Option 5 (conditional GC) |
| Lock overhead | Very Low | Low | Lock is only held during GC (~ms) |
| Other GC race conditions | Medium | High | Audit all GC.gc() calls in parallel code |
| Windows-specific regression | Low | High | Extensive Windows testing |

---

## Alternative Investigations

If Option 1 doesn't fully resolve the issue, investigate:

1. **Arrow.jl file handle management**
   - Check if Arrow.jl has built-in cleanup methods
   - File issue with Arrow.jl about Windows file handles

2. **Julia GC settings**
   - Try different GC configurations (`--gcthreads` flag)
   - Check if newer Julia versions handle concurrent GC better

3. **Windows-specific file APIs**
   - Use Win32 APIs directly via `ccall`
   - Implement more robust file deletion with proper handle closing

4. **Process-level GC control**
   - Investigate if GC can be disabled temporarily during critical sections
   - Use `GC.enable(false)` / `GC.enable(true)` (but this affects whole process)

---

## Future Considerations

1. **Remove Windows special case entirely**
   - Investigate if Arrow.jl 3.0+ handles Windows file handles better
   - Test if modern Julia versions (1.11+) have better GC thread safety

2. **Use immutable file operations**
   - Write to unique temp files, never overwrite
   - Clean up old files separately in non-critical path

3. **Lazy file cleanup**
   - Queue files for deletion
   - Single-threaded cleanup worker
   - Avoids concurrent GC issues entirely
