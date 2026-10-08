using Pioneer, LightGBM, Arrow, Tables, JSON, SHA, Dates

const ROOT = @__DIR__
const INPUT_DIR = "/scratch2/fs1/d.goldfarb/n.t.wamsley/runs/lgbm_threads_20261008T174117Z/inputs"
const STATS = mkpath(joinpath(ROOT, "stats"))
Pioneer.DEBUG_CONSOLE_LEVEL[] = 0

# Process-local prototype: current Pass-1 explicitly uses Threads.nthreads()
# for prediction even when the fitted classifier has a different budget.
# Keep its conversions/clamps identical and use that classifier's budget.
# No source file or package cache is modified by this override.
@eval Pioneer function _predict_pass1_block(cls, X::Matrix{Float32})
    n = size(X, 1)
    n == 0 && return Float32[]
    if cls isa NamedTuple && cls.kind == :constant
        return fill(clamp(Float32(cls.value), 1f-6, 1f0 - 1f-4), n)
    end
    raw = LightGBM.predict(cls, X; num_threads=cls.num_threads)
    scores = ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
    return clamp.(Float32.(scores), 1f-6, 1f0 - 1f-4)
end

function prepare_paths(name)
    directory = mkpath(joinpath(ROOT, name))
    paths = String[]
    for file in ("SWATH_r01.arrow", "SWATH_r02.arrow")
        path = joinpath(directory, file)
        cp(joinpath(INPUT_DIR, file), path)
        push!(paths, path)
    end
    return paths
end

function signature(paths)
    return [
        [(name=String(column), eltype=string(eltype(Tables.getcolumn(table, column))),
          length=length(Tables.getcolumn(table, column)),
          sha256=bytes2hex(sha256(reinterpret(UInt8, collect(Tables.getcolumn(table, column))))))
         for column in Tables.columnnames(table)]
        for table in (Arrow.Table(path * Pioneer.PASS1_SIDECAR_SUFFIX) for path in paths)
    ]
end

function one_pass(io, threads, repeat)
    paths = prepare_paths("threads_$(threads)_repeat_$(repeat)")
    GC.gc()
    cpu_start = Float64(ccall(:clock, Clong, ())) / 1_000_000
    measured = @timed Pioneer.train_and_predict_pass1_oom!(paths;
        features=Pioneer.model_features(Pioneer.ADVANCED_FEATURE_SET, false),
        compute_infold=true,
        lgbm_hp=merge(Pioneer.SCORING_LGBM_HP, (num_threads=threads,)),
        semisupervised=true)
    cpu_seconds = Float64(ccall(:clock, Clong, ())) / 1_000_000 - cpu_start
    state = measured.value
    metrics = (n_total=state.n_total, n_fold0=state.n_fold0, n_fold1=state.n_fold1,
               n_pool_fold0=state.n_pool_fold0, n_pool_fold1=state.n_pool_fold1,
               selected_iteration=state.semisupervised_iter,
               pool_targets_q01=state.target_q01, pool_decoys_q01=state.decoy_q01)
    values = signature(paths)
    row = (threads=threads, repeat=repeat, seconds=measured.time,
           cpu_seconds=cpu_seconds, allocated_bytes=measured.bytes,
           gc_seconds=measured.gctime, process_peak_rss_bytes=Sys.maxrss(),
           metrics..., sidecar_values=values)
    println(io, JSON.json(row)); flush(io)
    println(JSON.json(row)); flush(stdout)
    for cls in (state.cls_trained_on.fold0, state.cls_trained_on.fold1)
        cls isa LightGBM.LGBMClassification && Base.finalize(cls.booster)
    end
    return (metrics=metrics, values=values)
end

write(joinpath(STATS, "environment.json"), JSON.json(Dict(
    "utc" => string(Dates.now(Dates.UTC)), "job" => ENV["SLURM_JOB_ID"],
    "Julia_threads" => Threads.nthreads(), "caller_thread" => Threads.threadid(),
    "source_input" => INPUT_DIR, "source_rows" => 589057,
    "prototype" => "Only process-local Pass-1 prediction override: cls.num_threads instead of Threads.nthreads(). Training thread counts supplied through existing lgbm_hp.",
    "scope" => "Actual two-file semi-supervised Pass-1 pool, both folds, full OOF/in-fold Float32 sidecars. No replicated training rows. No whole-search/Windows validation.",
    "sha256" => bytes2hex(sha256(read(@__FILE__))),
), 2))

open(joinpath(STATS, "pass1_samples.jsonl"), "w") do io
    reference = one_pass(io, 24, 0) # JIT/cache warmup, also original-count reference
    for repeat in 1:3
        for threads in (repeat == 1 ? (24, 48, 64) : repeat == 2 ? (64, 48, 24) : (48, 24, 64))
            observed = one_pass(io, threads, repeat)
            @assert observed.metrics == reference.metrics "Pass-1 selected iteration/counts changed"
            @assert observed.values == reference.values "Pass-1 stored sidecar values changed"
        end
    end
end
write(joinpath(STATS, "COMPLETED"), "All nine Pass-1 trials match warmup metrics and complete sidecar values.\n")
println("PASS-1 THREAD EXPERIMENT COMPLETED")
