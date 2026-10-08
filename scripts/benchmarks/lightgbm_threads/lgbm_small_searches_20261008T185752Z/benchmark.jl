using Pioneer, LightGBM, Arrow, Tables, JSON, SHA, Dates, Random, Statistics

const ROOT = @__DIR__
const STATS = mkpath(joinpath(ROOT, "stats"))
const FIT_TIMINGS = NamedTuple[]
const PREDICT_TIMINGS = NamedTuple[]
const REPEATS = 5

process_cpu_seconds() = Float64(ccall(:clock, Clong, ())) / 1_000_000

# Process-local instrumentation preserves production fit behavior. No package
# file, cache or existing output is edited. Prediction alone needs a thread
# argument override: production Pass-1 explicitly resets it to Julia's count.
@eval Pioneer function _fit_pass1_booster(
    X::Matrix{Float32}, y::Vector{Bool}, lgbm_hp::NamedTuple,
)
    if isempty(y) || length(unique(y)) == 1
        constant_score = isempty(y) || !y[1] ? 0.0f0 : 1.0f0
        return (kind = :constant, value = constant_score, model = nothing)
    end
    cls = build_lightgbm_classifier(; _pass1_small_pool_hp(lgbm_hp, length(y))...)
    @assert Threads.threadid() == 1
    cpu_start = Main.process_cpu_seconds()
    started = time_ns()
    LightGBM.fit!(cls, X, _prepare_labels(y); verbosity = -1)
    _detach_lightgbm_training_data!(cls)
    push!(Main.FIT_TIMINGS, (threads=cls.num_threads, rows=size(X, 1),
        columns=size(X, 2), seconds=(time_ns() - started) / 1e9,
        cpu_seconds=Main.process_cpu_seconds() - cpu_start))
    return cls
end

@eval Pioneer function _predict_pass1_block(cls, X::Matrix{Float32})
    n = size(X, 1)
    n == 0 && return Float32[]
    if cls isa NamedTuple && cls.kind == :constant
        return fill(clamp(Float32(cls.value), 1f-6, 1f0 - 1f-4), n)
    end
    @assert Threads.threadid() == 1
    cpu_start = Main.process_cpu_seconds()
    started = time_ns()
    raw = LightGBM.predict(cls, X; num_threads=cls.num_threads)
    scores = ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
    result = clamp.(Float32.(scores), 1f-6, 1f0 - 1f-4)
    push!(Main.PREDICT_TIMINGS, (threads=cls.num_threads, rows=n,
        columns=size(X, 2), seconds=(time_ns() - started) / 1e9,
        cpu_seconds=Main.process_cpu_seconds() - cpu_start))
    return result
end

function emit(io, row)
    println(io, JSON.json(row)); flush(io)
    # The complete signatures are retained in JSONL; keep SLURM logs readable.
    brief = Dict(String(k) => v for (k, v) in pairs(row)
        if !(k in (:sidecar_values, :fit_calls, :predict_calls)))
    println(JSON.json(brief)); flush(stdout)
end

function sidecar_signatures(paths)
    return [
        [(name=String(column), eltype=string(eltype(Tables.getcolumn(table, column))),
          length=length(Tables.getcolumn(table, column)),
          sha256=bytes2hex(sha256(reinterpret(UInt8, collect(Tables.getcolumn(table, column))))))
         for column in Tables.columnnames(table)]
        for table in (Arrow.Table(path * Pioneer.PASS1_SIDECAR_SUFFIX) for path in paths)
    ]
end

function private_paths(source_paths, name)
    directory = mkpath(joinpath(ROOT, "trials", name))
    paths = String[]
    for source in source_paths
        path = joinpath(directory, basename(source))
        cp(source, path)
        push!(paths, path)
    end
    return paths
end

function one_pass(io, source_paths, case, threads, repeat; reference=nothing)
    paths = private_paths(source_paths, "$(case)_t$(threads)_r$(repeat)")
    GC.gc()
    empty!(FIT_TIMINGS); empty!(PREDICT_TIMINGS)
    cpu_start = process_cpu_seconds()
    measured = @timed Pioneer.train_and_predict_pass1_oom!(paths;
        features=Pioneer.model_features(Pioneer.ADVANCED_FEATURE_SET, false),
        compute_infold=true,
        lgbm_hp=merge(Pioneer.SCORING_LGBM_HP, (num_threads=threads,)),
        semisupervised=true)
    cpu_seconds = process_cpu_seconds() - cpu_start
    state = measured.value
    metrics = (n_total=state.n_total, n_fold0=state.n_fold0, n_fold1=state.n_fold1,
        n_pool_fold0=state.n_pool_fold0, n_pool_fold1=state.n_pool_fold1,
        selected_iteration=state.semisupervised_iter,
        pool_targets_q01=state.target_q01, pool_decoys_q01=state.decoy_q01,
        features=String.(state.available_features))
    values = sidecar_signatures(paths)
    row = (kind="pass1", case=case, files=length(paths), threads=threads, repeat=repeat,
        seconds=measured.time, cpu_seconds=cpu_seconds, allocated_bytes=measured.bytes,
        gc_seconds=measured.gctime, process_peak_rss_bytes=Sys.maxrss(), metrics...,
        fit_seconds=sum(x.seconds for x in FIT_TIMINGS; init=0.0),
        prediction_seconds=sum(x.seconds for x in PREDICT_TIMINGS; init=0.0),
        fit_calls=copy(FIT_TIMINGS), predict_calls=copy(PREDICT_TIMINGS),
        metrics_identical=reference === nothing || metrics == reference.metrics,
        sidecars_identical=reference === nothing || values == reference.values,
        sidecar_values=values)
    emit(io, row)
    for cls in (state.cls_trained_on.fold0, state.cls_trained_on.fold1)
        cls isa LightGBM.LGBMClassification && Base.finalize(cls.booster)
    end
    return (metrics=metrics, values=values)
end

function seed_search(io, dataset)
    config_path = joinpath(ROOT, "configs", dataset * ".json")
    config = JSON.parsefile(config_path)
    println("GENERATING CURRENT-DEVELOP PSMs: ", dataset); flush(stdout)
    empty!(FIT_TIMINGS); empty!(PREDICT_TIMINGS)
    GC.gc()
    measured = @timed Pioneer.SearchDIA(config_path)
    results = config["paths"]["results"]
    psmdir = joinpath(results, "temp_data", "main_search_psms")
    paths = sort(filter(p -> endswith(p, ".arrow") && !occursin("sidecar", basename(p)),
        readdir(psmdir; join=true)))
    @assert length(paths) == 3 "Expected three original PSM files"
    input_files = sort(filter(p -> endswith(p, ".arrow"),
        readdir(config["paths"]["ms_data"]; join=true)))
    emit(io, (kind="seed_search", case=dataset, threads=24, files=length(paths),
        seconds=measured.time, allocated_bytes=measured.bytes, gc_seconds=measured.gctime,
        process_peak_rss_bytes=Sys.maxrss(), config=config,
        raw_input_files=input_files,
        raw_input_sha256=[bytes2hex(sha256(read(p))) for p in input_files],
        psm_files=paths, psm_sha256=[bytes2hex(sha256(read(p))) for p in paths]))
    # Retain only small reports/logs when results are downloaded.
    log_dir = mkpath(joinpath(STATS, dataset))
    for file in readdir(results)
        if startswith(file, "pioneer_") && (endswith(file, ".log") || endswith(file, ".txt"))
            cp(joinpath(results, file), joinpath(log_dir, file))
        end
    end
    return paths
end

function main()
    @assert Threads.nthreads() == 24
    @assert Threads.threadid() == 1
    write(joinpath(STATS, "environment.json"), JSON.json(Dict(
        "utc" => string(Dates.now(Dates.UTC)), "job" => ENV["SLURM_JOB_ID"],
        "Julia_threads" => Threads.nthreads(), "caller_thread" => Threads.threadid(),
        "native_threads" => [24, 48], "measured_repeats" => REPEATS,
        "script_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "source_commit" => "cfb45759dccbbe209c9b83f46f4005d96117763b",
        "scope" => "Current-develop raw searches generate actual PSMs; timed paired Pass-1 trials on identical PSMs for 1-file and 3-file cohorts. No replicated observations or reduced training cap; no full-search 48-thread timing or Windows validation.",
        "prototype" => "Instrument production _fit_pass1_booster; Pass-1 prediction uses classifier.num_threads instead of Julia count; all native calls serial startup thread.",
        "memory" => "Process RSS high-water across cases and variants; Julia allocated bytes are cumulative, not native RSS.",
    ), 2))
    open(joinpath(STATS, "samples.jsonl"), "w") do io
        for dataset in ("SCP_Astral_250pg_3ms", "Olsen_Exploris_3P")
            paths = seed_search(io, dataset)
            Pioneer.DEBUG_CONSOLE_LEVEL[] = 0
            for count in (1, 3)
                case = dataset * "_$(count)file"
                selected = paths[1:count]
                reference = one_pass(io, selected, case, 24, 0)
                # Also warm the high-thread path; exclude both warmups.
                one_pass(io, selected, case, 48, 0; reference=reference)
                for repeat in 1:REPEATS
                    order = isodd(repeat) ? (24, 48) : (48, 24)
                    for threads in order
                        one_pass(io, selected, case, threads, repeat; reference=reference)
                    end
                end
            end
        end
    end
    write(joinpath(STATS, "COMPLETED"), "All seed searches and 40 measured Pass-1 trials completed. Check recorded parity booleans.\n")
    println("SMALL-SEARCH THREAD BENCHMARK COMPLETED")
end

main()
