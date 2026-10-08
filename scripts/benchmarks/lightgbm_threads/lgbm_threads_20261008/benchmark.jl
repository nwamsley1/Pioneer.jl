using Pioneer, LightGBM, Arrow, DataFrames, Tables, Random, Statistics, JSON, SHA, Dates, Sockets

const ROOT = @__DIR__
const STATS = mkpath(joinpath(ROOT, "stats"))
const THREAD_COUNTS = [16, 24, 32, 48, 64]
const REPEATS = 3

process_cpu_seconds() = Float64(ccall(:clock, Clong, ())) / 1_000_000

function timed_call(f)
    GC.gc()
    cpu_start = process_cpu_seconds()
    measured = @timed f()
    return (value=measured.value, seconds=measured.time, bytes=measured.bytes,
            gc_seconds=measured.gctime, cpu_seconds=process_cpu_seconds() - cpu_start)
end

function emit(io, row)
    println(io, JSON.json(row))
    flush(io)
    println(JSON.json(row))
    flush(stdout)
end

function differences(predictions, reference)
    @assert length(predictions) == length(reference)
    return (different=count(i -> predictions[i] != reference[i], eachindex(predictions)),
            max_abs=maximum(abs.(predictions .- reference); init=0.0))
end

function fit_one(X, labels, hp, threads, mode)
    model = Pioneer.build_lightgbm_classifier(; hp..., num_threads=threads)
    if mode == "col"
        model.force_row_wise = false
        model.force_col_wise = true
    end
    parameters = LightGBM.stringifyparams(model)
    dataset = timed_call() do
        ds = LightGBM.dataset_constructor(X, parameters, false)
        LightGBM.LGBM_DatasetSetField(ds, "label", labels)
        ds
    end
    training = timed_call() do
        LightGBM.fit!(model, dataset.value; verbosity=-1)
    end
    detach = timed_call() do
        Pioneer._detach_lightgbm_training_data!(model)
    end
    return model, (
        dataset_seconds=dataset.seconds, training_seconds=training.seconds,
        detach_seconds=detach.seconds,
        total_seconds=dataset.seconds + training.seconds + detach.seconds,
        cpu_seconds=dataset.cpu_seconds + training.cpu_seconds + detach.cpu_seconds,
        allocated_bytes=dataset.bytes + training.bytes + detach.bytes,
    )
end

function load_inputs()
    paths = [joinpath(ROOT, "inputs", name) for name in ("SWATH_r01.arrow", "SWATH_r02.arrow")]
    tables = Arrow.Table.(paths)
    available = intersect((Set(Symbol.(Tables.columnnames(table))) for table in tables)...)
    main = filter(f -> f in available, Pioneer.model_features(Pioneer.PRESCORE_FEATURES, false))
    advanced = filter(f -> f in available, Pioneer.model_features(Pioneer.ADVANCED_FEATURE_SET, false))
    wanted = unique(vcat(main, advanced, [:target, :cv_fold, :precursor_idx]))
    frame = vcat((DataFrame([f => collect(Tables.getcolumn(table, f)) for f in wanted]) for table in tables)...)
    for features in (main, advanced)
        if :num_enzymatic_termini in features && all(==(first(frame.num_enzymatic_termini)), frame.num_enzymatic_termini)
            filter!(!=(:num_enzymatic_termini), features)
        end
    end
    # Preserve precursor fold membership across files and stress-test repeats.
    membership = Dict{UInt32, UInt8}()
    for (precursor, fold) in zip(frame.precursor_idx, frame.cv_fold)
        @assert fold in (0, 1)
        prior = get!(membership, precursor, fold)
        @assert prior == fold "Precursor fold differs between source files"
    end
    rng = MersenneTwister(1776)
    train_pool = shuffle(rng, findall(==(UInt8(1)), frame.cv_fold))
    test_pool = shuffle(rng, findall(==(UInt8(0)), frame.cv_fold))
    @assert !isempty(train_pool) && !isempty(test_pool)
    information = Dict(
        "source_paths" => paths, "source_rows" => nrow(frame),
        "source_hashes" => [bytes2hex(open(sha256, path)) for path in paths],
        "train_fold" => 1, "test_fold" => 0,
        "unique_source_train_rows" => length(train_pool),
        "unique_source_test_rows" => length(test_pool),
        "main_features" => string.(main), "advanced_features" => string.(advanced),
        "missing_main_features" => string.(setdiff(Pioneer.model_features(Pioneer.PRESCORE_FEATURES, false), main)),
        "missing_advanced_features" => string.(setdiff(Pioneer.model_features(Pioneer.ADVANCED_FEATURE_SET, false), advanced)),
        "note" => "Larger matrices repeat saved observations within the same precursor CV fold; these are throughput stress cases, not independent cohorts or identification validation.",
    )
    write(joinpath(STATS, "inputs.json"), JSON.json(information, 2))
    return frame, main, advanced, train_pool, test_pool
end

function training_sweep(io, name, frame, features, train_pool, test_pool, n, hp)
    println("START TRAINING SCENARIO ", name)
    flush(stdout)
    columns = AbstractVector[frame[!, f] for f in features]
    rows = [train_pool[mod1(i, length(train_pool))] for i in 1:n]
    X = Pioneer.rows_matrix(columns, rows)
    labels = Pioneer._prepare_labels(frame.target[rows])
    check_rows = first(test_pool, min(100_000, length(test_pool)))
    Xcheck = Pioneer.rows_matrix(columns, check_rows)
    @assert length(unique(labels)) == 2
    baseline_predictions = nothing
    baseline_model = nothing
    mode_baselines = Dict{String, Vector{Float64}}()
    # Baseline fitted with actual production HP; native threads are explicit.
    baseline_model, baseline_metrics = fit_one(X, labels, hp, 24, "row")
    baseline_predictions = vec(LightGBM.predict(baseline_model, Xcheck; num_threads=24, verbosity=-1))
    mode_baselines["row"] = baseline_predictions
    emit(io, (kind="training_warmup", scenario=name, mode="row", threads=24,
              rows=n, features=length(features), repeat_rows=n > length(train_pool),
              baseline_metrics...))

    for repeat in 1:REPEATS
        # Rotate then reverse the sweep to limit ordering bias on one node.
        order = circshift(THREAD_COUNTS, repeat - 1)
        iseven(repeat) && reverse!(order)
        for threads in order
            model, metrics = fit_one(X, labels, hp, threads, "row")
            predictions = vec(LightGBM.predict(model, Xcheck; num_threads=24, verbosity=-1))
            parity = differences(predictions, baseline_predictions)
            emit(io, (kind="training", scenario=name, mode="row", threads=threads,
                      repeat=repeat, rows=n, features=length(features),
                      repeat_rows=n > length(train_pool), metrics...,
                      prediction_differences=parity.different, prediction_max_abs=parity.max_abs,
                      process_peak_rss_bytes=Sys.maxrss()))
            Base.finalize(model.booster)
        end
    end

    # Exploratory histogram layout test is separate from increasing threads.
    if n == 2_500_000
        for repeat in 1:REPEATS
            for threads in (isodd(repeat) ? (24, 64) : (64, 24))
                model, metrics = fit_one(X, labels, hp, threads, "col")
                predictions = vec(LightGBM.predict(model, Xcheck; num_threads=24, verbosity=-1))
                if !haskey(mode_baselines, "col")
                    mode_baselines["col"] = predictions
                end
                parity = differences(predictions, mode_baselines["col"])
                row_parity = differences(predictions, baseline_predictions)
                emit(io, (kind="training", scenario=name, mode="col", threads=threads,
                          repeat=repeat, rows=n, features=length(features), metrics...,
                          prediction_differences=parity.different, prediction_max_abs=parity.max_abs,
                          row_mode_prediction_differences=row_parity.different,
                          row_mode_prediction_max_abs=row_parity.max_abs,
                          process_peak_rss_bytes=Sys.maxrss()))
                Base.finalize(model.booster)
            end
        end
    end
    return baseline_model, features
end

function prediction_sweep(io, scenario, model, frame, features, test_pool)
    columns = AbstractVector[frame[!, f] for f in features]
    # Native matrix prediction and the complete Pioneer fill/predict/copy path.
    for n in (50_000, 500_000, 1_000_000)
        rows = [test_pool[mod1(i, length(test_pool))] for i in 1:n]
        matrix = Pioneer.rows_matrix(columns, rows)
        reference = vec(LightGBM.predict(model, matrix; num_threads=24, verbosity=-1))
        buffer = Float32[]
        output = Vector{Float64}(undef, n)
        for threads in THREAD_COUNTS
            LightGBM.predict(model, matrix; num_threads=threads, verbosity=-1)
            predictor = X -> LightGBM.predict(model, X; num_threads=threads, verbosity=-1)
            Pioneer.predict_rows!(output, predictor, buffer, columns, rows; batch_rows=n)
        end
        for repeat in 1:REPEATS
            order = circshift(THREAD_COUNTS, repeat - 1)
            iseven(repeat) && reverse!(order)
            for threads in order
                native = timed_call() do
                    vec(LightGBM.predict(model, matrix; num_threads=threads, verbosity=-1))
                end
                native_parity = differences(native.value, reference)
                emit(io, (kind="prediction_native", scenario=scenario, threads=threads,
                          repeat=repeat, rows=n, features=length(features),
                          seconds=native.seconds, cpu_seconds=native.cpu_seconds,
                          allocated_bytes=native.bytes,
                          prediction_differences=native_parity.different,
                          prediction_max_abs=native_parity.max_abs))
                predictor = X -> LightGBM.predict(model, X; num_threads=threads, verbosity=-1)
                complete = timed_call() do
                    Pioneer.predict_rows!(output, predictor, buffer, columns, rows; batch_rows=n)
                end
                complete_parity = differences(output, reference)
                emit(io, (kind="prediction_complete", scenario=scenario, threads=threads,
                          repeat=repeat, rows=n, features=length(features),
                          seconds=complete.seconds, cpu_seconds=complete.cpu_seconds,
                          allocated_bytes=complete.bytes,
                          prediction_differences=complete_parity.different,
                          prediction_max_abs=complete_parity.max_abs))
                native = nothing
                complete = nothing
            end
        end
    end
    Base.finalize(model.booster)
end

function main()
    environment = Dict(
        "utc" => string(Dates.now(Dates.UTC)), "host" => gethostname(),
        "Julia" => string(VERSION), "Julia_default_threads" => Threads.nthreads(),
        "caller_julia_thread_id" => Threads.threadid(),
        "LightGBM_jl" => string(Base.pkgversion(LightGBM)),
        "job" => get(ENV, "SLURM_JOB_ID", "unknown"),
        "cpu_allowed_list" => match(r"Cpus_allowed_list:\s*([^\n]+)", read("/proc/self/status", String)).captures[1],
        "parameters" => "Pioneer production hyperparameters; seed 1776, deterministic=true; row-wise primary, col-wise exploratory",
        "num_threads" => THREAD_COUNTS,
        "execution" => "All native calls sequential on startup Julia thread; Julia remains at 24; SLURM allocates 64 CPUs exclusively; no simultaneous model fits.",
        "cpu_seconds" => "Linux process CPU time; includes all native threads; GC outside timed regions",
        "memory" => "Sys.maxrss is process high-water across configurations, not variant-isolated RSS",
        "script_sha256" => bytes2hex(sha256(read(@__FILE__))),
    )
    write(joinpath(STATS, "environment.json"), JSON.json(environment, 2))
    frame, main_features, advanced, train_pool, test_pool = load_inputs()
    open(joinpath(STATS, "samples.jsonl"), "w") do io
        model, features = training_sweep(io, "main_250k", frame, main_features,
            train_pool, test_pool, 250_000, Pioneer.MAINSEARCH_LGBM_HP)
        prediction_sweep(io, "main_250k", model, frame, features, test_pool)
        for n in (1_000_000, 2_500_000)
            model, features = training_sweep(io, "scoring_$(n)", frame, advanced,
                train_pool, test_pool, n, Pioneer.SCORING_LGBM_HP)
            prediction_sweep(io, "scoring_$(n)", model, frame, features, test_pool)
        end
    end
    write(joinpath(STATS, "COMPLETED"), string(Dates.now(Dates.UTC)) * "\n")
    println("THREAD BENCHMARK COMPLETED")
end

main()
