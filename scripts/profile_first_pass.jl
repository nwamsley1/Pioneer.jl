#!/usr/bin/env julia
#
# Profile FirstPassSearch deconvolution on real data.
#
# The Profile.@profile is placed directly on perform_second_pass_search()
# inside FirstPassSearch.jl, so only that call is captured.
#
# Usage: julia --project=. --threads=10 scripts/profile_first_pass.jl <params.json>
#

using Profile
using PProf
using Pioneer

if length(ARGS) < 1
    error("Usage: julia --threads=10 scripts/profile_first_pass.jl <params.json>")
end
params_path = ARGS[1]

println("=" ^ 60)
println("  FirstPassSearch Profiler (targeted)")
println("  Threads: $(Threads.nthreads()) (FirstPass deconv forced to 1)")
println("=" ^ 60)
println()

# Warmup run — JIT compile everything
println("Warmup run...")
t0 = time()
try SearchDIA(params_path) catch end
println("  Done in $(round(time() - t0, digits=1))s")
println()

# Clear any warmup samples, use default Profile settings
Profile.clear()

println("Profiled run...")
t1 = time()
try SearchDIA(params_path) catch end
println("  Done in $(round(time() - t1, digits=1))s")
println()

# Print flat profile to terminal
println("=" ^ 70)
println("  Profile: top functions by sample count")
println("=" ^ 70)
Profile.print(; mincount=50, noisefloor=2, maxdepth=30)
println()

# Save pprof and start viewer
out_path = "first_pass_profile.pb.gz"
pprof(; out=out_path, web=false)
println("\nProfile written to $out_path")
PProf.refresh(file=out_path)
println("PProf at http://localhost:57599")
sleep(3600)
