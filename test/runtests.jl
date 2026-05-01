# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# Pioneer test suite — driver.
#
# Runs each test "part" in its own Julia subprocess so each starts with
# a clean compiler / heap state. We tried adding hoists individually for
# the order-dependent crashes (signal-11 in inference passes, GC heap
# corruption SIGABRT) but it became a treadmill — moving one victim
# surfaces the next. The actual root cause is upstream in Julia's
# compiler/GC; subprocess isolation sidesteps it cleanly.
#
# Compatibility:
#   - `Pkg.test`: invokes `runtests.jl` automatically; this driver
#     spawns subprocesses inside that test process. Standard.
#   - Code coverage (julia-actions/julia-processcoverage): when CI
#     passes `--code-coverage=user|all`, we propagate that flag into
#     each subprocess so each one writes its own .cov files. The action
#     merges them across PIDs.
#   - Threads: propagate `--threads=$(Threads.nthreads())`.

const TEST_PARTS = String[
    "runtests_part1_heavy.jl",
    "runtests_part2_units.jl",
    "runtests_part3_units.jl",
]

const PROJECT_ROOT = abspath(joinpath(@__DIR__, ".."))

"""
Translate `Base.JLOptions().code_coverage` (an integer enum) back into
the corresponding command-line flag(s) for a child julia process.

  0 = none
  1 = user            → `--code-coverage=user`
  2 = all             → `--code-coverage=all`
  3 = @<tracked_path> → `--code-coverage=@<JLOptions.tracked_path>`

`Pkg.test(coverage=true)` on Julia 1.12 uses form 3 — it pins coverage
collection to the package's source tree. Without handling 3 we
silently drop the flag and the subprocesses generate zero `.cov`
files (= the codecov "0% coverage" failure this branch hit).
"""
function _coverage_flags()
    opts = Base.JLOptions()
    cov = opts.code_coverage
    cov == 1 && return String["--code-coverage=user"]
    cov == 2 && return String["--code-coverage=all"]
    if cov == 3
        path_ptr = opts.tracked_path
        path_ptr == C_NULL && return String[]
        return String["--code-coverage=@" * unsafe_string(path_ptr)]
    end
    return String[]
end

function _run_part(rel_path::AbstractString)
    script = joinpath(@__DIR__, rel_path)
    isfile(script) || error("Test part not found: $script")

    cov_flags = _coverage_flags()
    nthreads_flag = "--threads=$(Threads.nthreads())"
    cmd = `$(Base.julia_cmd()) --project=$PROJECT_ROOT $nthreads_flag $cov_flags $script`

    println("\n", "="^78)
    println("[runtests.jl] Running part: ", rel_path)
    isempty(cov_flags) || println("[runtests.jl]   coverage: $(join(cov_flags, ' '))")
    println("[runtests.jl]   threads:  $(Threads.nthreads())")
    println("="^78)

    proc = run(pipeline(cmd; stdout = stdout, stderr = stderr); wait = true)
    return proc.exitcode
end

let
    results = Tuple{String, Int32}[]
    for part in TEST_PARTS
        push!(results, (part, _run_part(part)))
    end

    println("\n", "="^78)
    println("[runtests.jl] Test part summary")
    println("="^78)
    for (part, code) in results
        status = code == 0 ? "PASS" : "FAIL (exit=$code)"
        println("  ", rpad(part, 36), status)
    end

    failed = [p for (p, code) in results if code != 0]
    if !isempty(failed)
        error("Test parts failed: ", join(failed, ", "))
    end
end
