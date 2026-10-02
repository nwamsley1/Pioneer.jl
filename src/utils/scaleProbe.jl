# Opt-in profiling probes for large searches.
#
# PIONEER_SCALE_PROBE selects the mode, read once at the first probe:
#   unset or "0"  off: each probe is one branch on a cached Symbol
#   "1"           record wall time, current RSS, Julia live heap and cumulative allocation
#   "gc"          as "1", but run a full GC first, so `live_gb` is memory still referenced
# Records append to PIONEER_SCALE_PROBE_FILE (default `scale_probe.tsv` in the working directory),
# one tab-separated line per probe:
#   unix_time  stage  event  file_idx  seconds  rss_gb  live_gb  allocd_gb  value
# `seconds` is the duration the caller measured (NaN when none); `value` is an extra number
# (for scale_probe_size, the summarysize in GB).

const _SCALE_PROBE_MODE = Ref{Symbol}(:unread)
const _SCALE_PROBE_LOCK = ReentrantLock()

function scale_probe_mode()
    if _SCALE_PROBE_MODE[] === :unread
        v = lowercase(strip(get(ENV, "PIONEER_SCALE_PROBE", "0")))
        _SCALE_PROBE_MODE[] = v == "gc" ? :gc : v in ("1", "true", "yes") ? :on : :off
    end
    return _SCALE_PROBE_MODE[]
end
scale_probe_on() = scale_probe_mode() !== :off

"Current resident set size in GB (not the peak): /proc on Linux, `ps` elsewhere."
function current_rss_gb()
    if Sys.islinux()
        pages = parse(Int, split(read("/proc/self/statm", String))[2])
        return pages * 4096 / 2^30
    end
    kb = tryparse(Int, strip(read(`ps -o rss= -p $(getpid())`, String)))
    return kb === nothing ? NaN : kb / 2^20
end

function scale_probe(stage::AbstractString, event::AbstractString;
                     file_idx::Integer = 0, seconds::Real = NaN, value::Real = NaN, measure::Bool = true)
    mode = scale_probe_mode()
    mode === :off && return nothing
    # measure = false records only the caller's timing (no GC, no memory readings)
    measure && mode === :gc && GC.gc(true)
    live = measure ? Base.gc_live_bytes() / 2^30 : NaN
    allocd = Base.gc_total_bytes(Base.gc_num()) / 2^30
    rss = measure ? current_rss_gb() : NaN
    line = join((round(time(); digits = 3), stage, event, file_idx, round(seconds; digits = 3),
                 round(rss; digits = 3), round(live; digits = 3), round(allocd; digits = 3),
                 round(value; digits = 6)), '\t')
    path = get(ENV, "PIONEER_SCALE_PROBE_FILE", "scale_probe.tsv")
    lock(_SCALE_PROBE_LOCK) do
        open(io -> println(io, line), path, "a")
    end
    return nothing
end

"Record `Base.summarysize(obj)` in GB as `value` (only when probing is on: summarysize walks the object)."
scale_probe_size(stage::AbstractString, event::AbstractString; file_idx::Integer = 0, obj) =
    scale_probe_size(stage, event, obj; file_idx)
function scale_probe_size(stage::AbstractString, event::AbstractString, obj; file_idx::Integer = 0)
    scale_probe_on() || return nothing
    scale_probe(stage, event; file_idx, value = Base.summarysize(obj) / 2^30)
end
