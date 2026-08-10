# Refresh a generated precompile corpus after dependencies or method shapes
# change without throwing away statements harvested from expensive real-data
# searches. Existing statements that no longer resolve are removed, then valid
# statements from one or more new trace files are appended and deduplicated.
#
# Run:
#   julia --project=. src/build/refresh_precompile_statements.jl \
#       src/build/precompile_statements.jl /tmp/precompile_statements.jl \
#       /tmp/new_trace.jl [...]

using Base.Meta
using Pioneer

length(ARGS) >= 3 || error(
    "usage: refresh_precompile_statements.jl INPUT OUTPUT TRACE [TRACE ...]",
)

const INPUT_FILE = abspath(ARGS[1])
const OUTPUT_FILE = abspath(ARGS[2])
const TRACE_FILES = abspath.(ARGS[3:end])
const Staging = Module()

function resolves(statement::AbstractString)
    expression = try
        Meta.parse(statement)
    catch
        return false
    end
    isexpr(expression, :call) || return false
    popfirst!(expression.args)
    expression.head = :tuple

    while true
        try
            Core.eval(Staging, expression)
            return true
        catch error_value
            error_value isa UndefVarError || return false
            dependency = string(error_value.var)
            modules = filter(pair -> pair.first.name == dependency, Base.loaded_modules)
            length(modules) == 1 || return false
            Base.eval(Staging, :($(Symbol(dependency)) = $(only(modules).second)))
        end
    end
end

input_lines = readlines(INPUT_FILE)
first_statement = findfirst(line -> startswith(strip(line), "precompile("), input_lines)
isnothing(first_statement) && error("no precompile statements found in $INPUT_FILE")
header = input_lines[1:(first_statement - 1)]

candidates = String[]
append!(candidates, input_lines[first_statement:end])
for trace_file in TRACE_FILES
    isfile(trace_file) || error("missing trace file $trace_file")
    append!(candidates, readlines(trace_file))
end

seen = Set{String}()
survivors = String[]
dropped = Ref(0)
for raw_line in candidates
    statement = strip(raw_line)
    startswith(statement, "precompile(") || continue
    if statement in seen
        continue
    elseif resolves(statement)
        push!(seen, statement)
        push!(survivors, statement)
    else
        dropped[] += 1
    end
end

open(OUTPUT_FILE, "w") do io
    for line in header
        println(io, line)
    end
    for statement in survivors
        println(io, statement)
    end
end

println(
    "precompile refresh: kept $(length(survivors)) unique statements, " *
    "dropped $(dropped[]) unresolved statements",
)
