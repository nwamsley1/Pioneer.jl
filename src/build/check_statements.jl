# Guard against `src/build/precompile_statements.jl` silently rotting.
#
# `create_app` wraps every precompile statement in try/catch and skips anything that fails to parse,
# resolve, or eval (PackageCompiler.jl, create_sysimg_object_file). A statements file can therefore
# go 90% dead -- because a dependency bumped a version, or because Pioneer code moved -- and the
# build stays green while shipping a binary that JIT-compiles at startup again. `Manifest.toml` is
# gitignored, so CI resolves fresh dependency versions on every run and this drift is expected over
# time, not hypothetical.
#
# This replays the same parse -> eval pipeline and fails loudly when the survival rate drops. It is
# the analogue of what b18ac0176 did for snoop.jl, which had the identical failure mode: a green
# build shipping an uncompiled SearchDIA.
#
# Run:  julia --project=. src/build/check_statements.jl [min_survival_fraction]
# Costs about a minute; run it BEFORE create_app so a failure is cheap.

using Base.Meta
using Pioneer

const FILE = joinpath(@__DIR__, "precompile_statements.jl")

# 99.5%, not 95%. Dropping the file-count parameter from ArrowTableReference invalidated exactly 120
# statements -- the ones covering the whole search pipeline -- and at a 95% floor that scored 97.95%
# and passed. A percentage is a crude proxy: 120 pipeline statements matter far more than 120
# arbitrary Base ones. 0.5% of ~5,800 leaves ~29 statements of slack for ordinary dependency drift,
# which is enough to avoid false alarms while still catching a structural break.
const MIN_SURVIVAL = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 0.995

isfile(FILE) || error("missing $FILE -- create_app references it, so the build would ship without " *
                      "the real-workload statements")

const Staging = Module()

# Mirrors PackageCompiler's UndefVarError recovery: pull top-level modules into the staging area by
# name and retry, which is how statements naming a package Pioneer does not re-export resolve.
function try_eval!(ps)
    while true
        try
            Core.eval(Staging, ps)
            return nothing
        catch e
            if e isa UndefVarError
                dep = string(e.var)
                mods = filter(p -> p.first.name == dep, Base.loaded_modules)
                length(mods) == 1 || return "unresolved module $dep"
                Base.eval(Staging, :($(Symbol(dep)) = $(only(mods).second)))
            else
                return first(split(sprint(showerror, e), '\n'))
            end
        end
    end
end

function survey(file)
    total = 0
    kept = 0
    reasons = Dict{String,Int}()
    for raw in eachline(file)
        stmt = strip(raw)
        (isempty(stmt) || startswith(stmt, "#")) && continue
        total += 1
        ps = try
            Meta.parse(stmt)
        catch
            reasons["parse failed"] = get(reasons, "parse failed", 0) + 1
            continue
        end
        if !isexpr(ps, :call)
            reasons["not a call"] = get(reasons, "not a call", 0) + 1
            continue
        end
        popfirst!(ps.args)
        ps.head = :tuple
        why = try_eval!(ps)
        if why === nothing
            kept += 1
        else
            reasons[why] = get(reasons, why, 0) + 1
        end
    end
    return total, kept, reasons
end

total, kept, reasons = survey(FILE)
survival = total == 0 ? 0.0 : kept / total
println("precompile_statements.jl: $kept / $total survive ($(round(100survival; digits=2))%)")
if !isempty(reasons)
    println("top failure reasons:")
    for (r, n) in first(sort(collect(reasons); by = last, rev = true), 10)
        println("  ", rpad(n, 6), r)
    end
end

if survival < MIN_SURVIVAL
    error("""
          Only $(round(100survival; digits=2))% of precompile statements still resolve \
          (floor is $(round(100MIN_SURVIVAL; digits=2))%).

          The shipped binary would lose the startup coverage this file provides. Regenerate it -- \
          see the header of src/build/precompile_statements.jl for the trace command -- and commit \
          the result.
          """)
end
println("OK")
