# Log Message Truncation — Implementation Plan

## Goal

Prevent any single log write from ballooning the log files (e.g., when an error message accidentally embeds a giant array). Add a configurable per‑message byte cap so long messages are truncated safely before being written to console or files.

## Scope and Non‑Goals

- In scope: Truncate individual log messages to a maximum size (bytes) across all logging sinks (console, essential, debug, warnings).
- Not in scope (but considered later): Log rotation, global file size caps, duplicate‑message throttling.

## Configuration

- Add a new parameter: `logging.max_message_bytes`.
  - Type: integer (bytes)
  - Default: `4096` (4 KB)
  - Min allowed: `1024` (1 KB) — guards against unrealistic settings.
  - Max allowed: e.g. `1048576` (1 MB) — sanity cap.
  - Environment override: if `ENV["PIONEER_MAX_LOG_MSG_BYTES"]` is set, it takes precedence (validated and clamped to [min, max]).

### Where to put it

- Update the default JSON configuration (example configs):
  - `assets/example_config/defaultSearchParams.json` (or the currently recommended default JSON file) under the existing `logging` block:
    ```json
    "logging": {
      "debug_console_level": 0,
      "max_message_bytes": 4096
    }
    ```
- Update the user docs at `docs/src/user_guide/parameters.md` to document the new field (default, meaning, and ENV override).

## Design

### Helper function: truncate_for_log

Add a small helper that ensures UTF‑8 safe truncation by bytes, with a clear suffix annotation:

```julia
const MAX_LOG_MSG_BYTES = Ref{Int}(4096)  # set at startup from params/env

function truncate_for_log(msg::String; max_bytes::Int=MAX_LOG_MSG_BYTES[])
    # Fast path
    nb = ncodeunits(msg)
    if nb <= max_bytes
        return msg
    end

    # Walk codepoints and accumulate bytes until adding the next char would exceed max_bytes
    n = 0
    lastidx = 0
    @inbounds for (i, ch) in enumerate(eachindex(msg))
        cu = ncodeunits(msg[i])  # bytes for this char
        if n + cu > max_bytes
            break
        end
        n += cu
        lastidx = i
    end

    # Fallback if something odd happens
    lastidx == 0 && return msg[1:prevind(msg, end)] * " … [truncated]"

    prefix = @views String(msg[1:lastidx])
    omitted = nb - n
    return string(prefix, " … [truncated ", omitted, " bytes]")
end
```

Notes:
- Uses byte accounting (`ncodeunits`) for accuracy.
- Ensures we never split a multi‑byte codepoint by advancing at character boundaries (`eachindex`).
- Adds an explicit suffix with the number of omitted bytes.

### Integration points

- Centralize all message writes through `truncate_for_log`:
  - `user_info`, `user_warn`, `user_error`, `user_print` in `src/Pioneer.jl`.
  - Any place that builds a message via `sprint(showerror, e, bt)` must pass the result through `truncate_for_log` before printing.

#### Example wiring (conceptual, not exact code):

```julia
function user_error(msg::String, file::String="", line::String="", mod::String="")
    msg_trunc = truncate_for_log(msg)
    # console
    println(msg_trunc)
    # files (ESSENTIAL/CONSOLE/DEBUG/WARNINGS): also use msg_trunc
end

# For exception printing
err_str = sprint(showerror, e, bt)
println(truncate_for_log(err_str))
```

### Parameter loading

- Extend the parameter parsing to read `params.logging.max_message_bytes` into the runtime config.
- On SearchDIA start, set `MAX_LOG_MSG_BYTES[]` to:
  - ENV override if present and valid; else
  - params value; else default 4096.
- Clamp to `[1024, 1048576]`.

## Backward Compatibility

- Default behavior (4096 bytes) ensures most messages are unaffected. Only unusually long messages will be truncated.
- Existing configs without the new field continue to work (default applied).
- ENV override offers a runtime escape hatch for ops (e.g., set lower during CI).

## Testing Plan

- Unit‑style tests (lightweight):
  - Truncate a short ASCII message → unchanged.
  - Truncate a long ASCII message → length limited; suffix present.
  - Truncate a long UTF‑8 string with multibyte characters → no invalid UTF‑8; suffix present.
  - Truncate a gigantic error string created via `sprint(showerror, …)` → limited, suffix present.

- Integration sanity:
  - Simulate a failure that previously dumped huge arrays (e.g., `BoundsError` with big vector) and verify the debug/warnings logs contain a short, truncated message instead of megabytes of data.

## Rollout Steps

1. Add `logging.max_message_bytes` to example/default JSON(s).
2. Update parameter parsing (where logging settings are read) to store the value.
3. Add `MAX_LOG_MSG_BYTES` ref and `truncate_for_log` helper in `src/Pioneer.jl`.
4. Apply truncation in `user_info`, `user_warn`, `user_error`, `user_print`, and any direct `sprint(showerror, …)` prints.
5. Clamp and initialize from params/env at SearchDIA startup.
6. Update documentation (`docs/src/user_guide/parameters.md`).
7. (Optional) Add a minimal test that asserts truncation occurs.

## Risks and Mitigations

- Risk: Truncation could hide critical context.
  - Mitigation: Preserve head of message; append explicit “truncated N bytes”. Operators can raise the cap via config/env if needed.
- Risk: Performance overhead of character iteration.
  - Mitigation: Only executed when message exceeds threshold; normal messages take the fast path.

## Future Extensions (Optional)

- Duplicate warning throttling: per‑window max repeats to prevent millions of identical lines.
- Global log size guard: after each file/phase, abort if logs exceed a threshold (e.g., 200 MB).
- Simple rotation: when a log exceeds limit, rename to `*.YYYYMMDD-HHMMSS.log` and start fresh (keep at most N).

