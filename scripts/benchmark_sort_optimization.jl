#!/usr/bin/env julia
# Benchmark script: measure FirstPassSearch timing and verify correctness
# Usage: julia --threads 10 --gcthreads 5,1 --project=. scripts/benchmark_sort_optimization.jl

using Pioneer
using Arrow, DataFrames

const PARAMS_PATH = "./data/OlsenAstralThreeProteome200ng_globalq05pct.json"
const RESULTS_DIR = "/Users/nathanwamsley/Data/For_Figures/OlsenAstralThreeProteome200ng/OlsenAstralThreeProteome200ng_sort_benchmark"

# Modify results path to avoid overwriting existing results
params_json = read(PARAMS_PATH, String)
params_json = replace(params_json,
    r"\"results\":.*" => "\"results\": \"$RESULTS_DIR\"")

# Write temp params
const TEMP_PARAMS = "/tmp/sort_benchmark_params.json"
write(TEMP_PARAMS, params_json)

println("=" ^ 60)
println("SORT OPTIMIZATION BENCHMARK")
println("=" ^ 60)
println("Threads: $(Threads.nthreads())")
println("GC threads: $(ccall(:jl_n_gcthreads, Int32, ()))")
println("Commit: ", strip(read(`git rev-parse --short HEAD`, String)))
println()

# Run full pipeline
t_total = @elapsed begin
    SearchDIA(TEMP_PARAMS)
end

println()
println("=" ^ 60)
println("TOTAL PIPELINE TIME: $(round(t_total, digits=1))s")
println("=" ^ 60)

# Count unique targets from final results
psms_dir = joinpath(RESULTS_DIR, "passing_psms")
if isdir(psms_dir)
    psm_files = filter(f -> endswith(f, ".arrow"), readdir(psms_dir, join=true))
    all_prec_ids = Set{UInt32}()
    n_targets = 0
    for f in psm_files
        df = DataFrame(Arrow.Table(f))
        if "is_decoy" in names(df) && "precursor_idx" in names(df)
            targets = filter(r -> !r.is_decoy, df)
            for pid in targets.precursor_idx
                push!(all_prec_ids, UInt32(pid))
            end
        end
    end
    println("UNIQUE TARGET PRECURSORS: $(length(all_prec_ids))")
else
    # Try temp_data passing_psms
    temp_psms = joinpath(RESULTS_DIR, "temp_data", "scoring_search")
    if isdir(temp_psms)
        psm_files = filter(f -> endswith(f, ".arrow"), readdir(temp_psms, join=true))
        all_prec_ids = Set{UInt32}()
        for f in psm_files
            df = DataFrame(Arrow.Table(f))
            if "is_decoy" in names(df) && "precursor_idx" in names(df)
                targets = filter(r -> !r.is_decoy, df)
                for pid in targets.precursor_idx
                    push!(all_prec_ids, UInt32(pid))
                end
            end
        end
        println("UNIQUE TARGET PRECURSORS: $(length(all_prec_ids))")
    end
end
println("=" ^ 60)
