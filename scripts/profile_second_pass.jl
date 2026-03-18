# Profile SecondPassSearch with PProf
#
# Usage:
#   julia --threads=auto scripts/profile_second_pass.jl [config_path]
#
# Defaults to ecoli test config for warmup, then profiles the real config.
# Produces a flamegraph at second_pass_profile.pb.gz viewable with:
#   using PProf; PProf.refresh(file="second_pass_profile.pb.gz")

using Pioneer
using Profile
using PProf

# Parse config path from command line, or use default
config_path = length(ARGS) >= 1 ? ARGS[1] : "./data/ecoli_test/ecoli_test_params.json"
warmup_path = "./data/ecoli_test/ecoli_test_params.json"

println("=== SecondPassSearch Profiling ===")
println("Config: $config_path")
println("Threads: $(Threads.nthreads())")

# Step 1: JIT warmup on small test data
println("\n--- Step 1: JIT warmup on ecoli test ---")
SearchDIA(warmup_path)
println("Warmup complete.")

# Step 2: Profile the real run
println("\n--- Step 2: Profiling $config_path ---")
Profile.clear()
Profile.init(n = 10_000_000, delay = 0.0001)  # 100μs sampling, large buffer

@profile SearchDIA(config_path)

# Step 3: Save flamegraph
outfile = "second_pass_profile.pb.gz"
println("\n--- Step 3: Saving flamegraph to $outfile ---")
pprof(out=outfile, web=false)
println("Done. View with:")
println("  using PProf; PProf.refresh(file=\"$outfile\")")
