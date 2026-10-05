# Utility helpers for working with LightGBM directly through its native API.

struct LightGBMModel
    booster::Union{LightGBM.LGBMClassification, LightGBM.LGBMRanking, Nothing}
    features::Vector{Symbol}
    constant_prediction::Union{Float32, Nothing}
end


# Serializes every entry into LightGBM's C API. LightGBM's bundled OpenMP
# (libomp) aborts in __kmpc_end_serialized_parallel when predict is entered
# concurrently from multiple Julia threads (pass1_oom.jl's parallel_foreach!).
const LGBM_C_LOCK = ReentrantLock()

"""
    feature_matrix(df, features) -> Matrix{Float32}

Construct a dense matrix with the columns in `features` converted to `Float32`.
Missing values are imputed with `0.0f0`. Columns are filled in parallel via a
typed inner helper to avoid Float64 intermediates and per-column broadcast
temporaries.
"""
function _fill_column!(M::Matrix{Float32}, j::Int, col::AbstractVector{T}) where {T<:Real}
    @inbounds @simd for i in eachindex(col)
        M[i, j] = Float32(col[i])
    end
end

function _fill_column!(M::Matrix{Float32}, j::Int, col::AbstractVector{<:Union{Missing, T}}) where {T<:Real}
    @inbounds for i in eachindex(col)
        v = col[i]
        M[i, j] = v === missing ? 0.0f0 : Float32(v)
    end
end

_fill_column!(::Matrix{Float32}, ::Int, col::AbstractVector) =
    throw(ArgumentError("Unsupported feature type $(eltype(col)) for LightGBM"))

function feature_matrix(df::AbstractDataFrame, features::Vector{Symbol})
    n = nrow(df)
    m = length(features)
    matrix = Matrix{Float32}(undef, n, m)
    cols = AbstractVector[df[!, f] for f in features]
    Threads.@threads for j in 1:m
        _fill_column!(matrix, j, cols[j])
    end
    return matrix
end

#############################################################################
# Reusable feature-matrix buffers
#############################################################################

"""
    LGBMMatrixBuffers

Reusable backing stores for the LightGBM feature matrices built per MS file.

MainSearch's file loop is sequential (`execute_search`), so one buffer set is
reused across every file instead of allocating a fresh `Matrix{Float32}` each
time. No matrix ever holds a whole file: rows are filled straight from the PSM
table's columns (`rows_matrix!`), and prediction runs in fixed-size batches
(`predict_rows!`), so the working set stays bounded however many PSMs a file has.

One buffer per distinct matrix, so no two simultaneously-live matrices ever
share memory:
- `train` — sub-sampled training rows (≤ max_train × nfeat)
- `batch` — one prediction batch (≤ `LGBM_PREDICT_BATCH_ROWS` × nfeat)

Buffers grow monotonically (`resize!` only when a bigger matrix is needed). Size
a buffer *before* wrapping any matrix from it: `resize!` may move the data,
which would leave a live wrap dangling.
"""
struct LGBMMatrixBuffers
    train::Vector{Float32}
    batch::Vector{Float32}
end

LGBMMatrixBuffers() = LGBMMatrixBuffers(Float32[], Float32[])

"Rows per LightGBM prediction batch (`predict_rows!`): ~100 MB at 49 features."
const LGBM_PREDICT_BATCH_ROWS = 500_000

"""
    _size_matrix_buffer!(buf, n, m)

Grow `buf` so it can back an `n × m` `Float32` matrix. Must run before any
matrix is wrapped from `buf` (see `LGBMMatrixBuffers`).
"""
@inline function _size_matrix_buffer!(buf::Vector{Float32}, n::Int, m::Int)
    need = n * m
    length(buf) < need && resize!(buf, need)
    return nothing
end

"""
    _wrap_matrix_buffer(buf, n, m) -> Matrix{Float32}

Wrap the first `n*m` elements of `buf` as an `n × m` column-major matrix without
allocating. `buf` must already be sized by `_size_matrix_buffer!` and must be
kept alive — rooted or `GC.@preserve`d — for as long as the matrix is read.
"""
@inline function _wrap_matrix_buffer(buf::Vector{Float32}, n::Int, m::Int)
    n == 0 && return Matrix{Float32}(undef, 0, m)
    length(buf) >= n * m ||
        throw(ArgumentError("matrix buffer holds $(length(buf)) elements, need $(n * m)"))
    return unsafe_wrap(Matrix{Float32}, pointer(buf), (n, m))
end

"""
    feature_matrix!(buf, df, features) -> Matrix{Float32}

`feature_matrix` into a caller-owned buffer: identical column-parallel fill, no
matrix allocation. Keep `buf` alive while the returned matrix is in use.
"""
function feature_matrix!(buf::Vector{Float32}, df::AbstractDataFrame, features::Vector{Symbol})
    n = nrow(df)
    m = length(features)
    _size_matrix_buffer!(buf, n, m)
    matrix = _wrap_matrix_buffer(buf, n, m)
    cols = AbstractVector[df[!, f] for f in features]
    Threads.@threads for j in 1:m
        _fill_column!(matrix, j, cols[j])
    end
    return matrix
end

#############################################################################
# Row-indexed fills and batched prediction
#############################################################################

# Same per-type conversion as `_fill_column!` (Float32(v); missing -> 0.0f0), with the rows chosen by an index
# vector, so a row matrix holds exactly the values `feature_matrix(df)[rows, :]` would.
function _fill_rows!(M::Matrix{Float32}, j::Int, col::AbstractVector{T},
                     rows::AbstractVector{<:Integer}) where {T<:Real}
    @inbounds for k in eachindex(rows)
        M[k, j] = Float32(col[rows[k]])
    end
end

function _fill_rows!(M::Matrix{Float32}, j::Int, col::AbstractVector{<:Union{Missing, T}},
                     rows::AbstractVector{<:Integer}) where {T<:Real}
    @inbounds for k in eachindex(rows)
        v = col[rows[k]]
        M[k, j] = v === missing ? 0.0f0 : Float32(v)
    end
end

_fill_rows!(::Matrix{Float32}, ::Int, col::AbstractVector, ::AbstractVector{<:Integer}) =
    throw(ArgumentError("Unsupported feature type $(eltype(col)) for LightGBM"))

"""
    fill_rows!(M, cols, rows) -> M

Write row `k` of `M` from row `rows[k]` of each feature column (`cols[j]` -> column `j`).
`M` must be `length(rows) × length(cols)`. Columns are filled in parallel, as in `feature_matrix!`.
"""
function fill_rows!(M::Matrix{Float32}, cols::Vector{AbstractVector}, rows::AbstractVector{<:Integer})
    size(M) == (length(rows), length(cols)) ||
        throw(DimensionMismatch("matrix is $(size(M)), rows × features is $((length(rows), length(cols)))"))
    # Validate the row indices once so the fill loops can be @inbounds.
    if !isempty(rows) && !isempty(cols)
        lo, hi = extrema(rows)
        n_src = length(cols[1])
        (1 <= lo && hi <= n_src) || throw(BoundsError(cols[1], lo < 1 ? lo : hi))
        all(c -> length(c) == n_src, cols) || throw(DimensionMismatch("feature columns differ in length"))
    end
    Threads.@threads for j in eachindex(cols)
        _fill_rows!(M, j, cols[j], rows)
    end
    return M
end

"""
    rows_matrix!(buf, cols, rows) -> Matrix{Float32}

The feature rows `rows` as a `length(rows) × length(cols)` matrix wrapped from `buf` (grown if needed).
Keep `buf` alive (`GC.@preserve`) while the returned matrix is in use.
"""
function rows_matrix!(buf::Vector{Float32}, cols::Vector{AbstractVector}, rows::AbstractVector{<:Integer})
    n, m = length(rows), length(cols)
    _size_matrix_buffer!(buf, n, m)
    return fill_rows!(_wrap_matrix_buffer(buf, n, m), cols, rows)
end

"`rows_matrix!` into a freshly allocated matrix."
rows_matrix(cols::Vector{AbstractVector}, rows::AbstractVector{<:Integer}) =
    fill_rows!(Matrix{Float32}(undef, length(rows), length(cols)), cols, rows)

"""
    predict_rows!(out, predict, buf, cols, rows; batch_rows = LGBM_PREDICT_BATCH_ROWS) -> out

Score the PSMs `rows` without materializing them all at once: fill at most `batch_rows` of them into one matrix
wrapped from `buf`, call `predict(matrix)`, and write prediction `k` of a batch to `out[i]` where `rows[i]` is
that batch's row `k`. So `out[i]` is the score of PSM `rows[i]` (`out` aligned with `rows`), whatever the batch
size. `predict` must return one value per matrix row, in row order (a vector or an n × 1 matrix).

Each batch is wrapped with its own exact shape, `n × m` over the front of `buf`: a view of the first `n` rows of
a larger matrix would not be contiguous, and LightGBM's predict needs a dense `Matrix`.
"""
function predict_rows!(out::AbstractVector{Float64}, predict, buf::Vector{Float32},
                       cols::Vector{AbstractVector}, rows::AbstractVector{<:Integer};
                       batch_rows::Int = LGBM_PREDICT_BATCH_ROWS)
    length(out) == length(rows) || throw(DimensionMismatch("out has $(length(out)) entries for $(length(rows)) rows"))
    batch_rows > 0 || throw(ArgumentError("batch_rows must be positive"))
    isempty(rows) && return out
    m = length(cols)
    _size_matrix_buffer!(buf, min(batch_rows, length(rows)), m)    # once, before any wrap
    GC.@preserve buf begin
        for lo in 1:batch_rows:length(rows)
            hi = min(lo + batch_rows - 1, length(rows))
            batch = view(rows, lo:hi)
            X = fill_rows!(_wrap_matrix_buffer(buf, length(batch), m), cols, batch)
            pred = predict(X)
            length(pred) == length(batch) ||
                throw(DimensionMismatch("predict returned $(length(pred)) values for $(length(batch)) rows"))
            @inbounds for k in eachindex(batch)
                out[lo + k - 1] = pred[k]
            end
        end
    end
    return out
end

function build_lightgbm_classifier(; num_iterations::Integer = 100,
                                    max_depth::Integer = -1,
                                    num_leaves::Integer = 31,
                                    learning_rate::Real = 0.05,
                                    feature_fraction::Real = 0.5,
                                    bagging_fraction::Real = 1.0,
                                    bagging_freq::Integer = 0,
                                    min_data_in_leaf::Integer = 500,
                                    min_gain_to_split::Real = 1.0,
                                    lambda_l1::Real = 0.0,
                                    lambda_l2::Real = 0.0,
                                    max_bin::Integer = 255,
                                    monotone_constraints::AbstractVector{<:Integer} = Int[],
                                    monotone_constraints_method::AbstractString = "basic",
                                    num_threads::Integer = Threads.nthreads(),
                                    metric = ["binary_logloss"],
                                    objective::AbstractString = "binary",
                                    is_unbalance = false,
                                    verbosity::Integer = -1)
    return LightGBM.LGBMClassification(
        objective = objective,
        metric = metric,
        learning_rate = float(learning_rate),
        num_iterations = Int(num_iterations),
        num_leaves = Int(num_leaves),
        max_depth = Int(max_depth),
        feature_fraction = float(feature_fraction),
        bagging_fraction = float(bagging_fraction),
        bagging_freq = Int(bagging_freq),
        min_data_in_leaf = Int(min_data_in_leaf),
        min_gain_to_split = float(min_gain_to_split),
        lambda_l1 = float(lambda_l1),
        lambda_l2 = float(lambda_l2),
        max_bin = Int(max_bin),
        monotone_constraints = Int.(monotone_constraints),
        monotone_constraints_method = String(monotone_constraints_method),
        num_threads = Int(num_threads),
        num_class = 1,
        verbosity = Int(verbosity),
        is_unbalance = is_unbalance,
        seed = 1776, # potentialy needed for stable results
        deterministic = true, # potentialy needed for stable results
        force_row_wise = true # potentialy needed for stable results
    )
end

function _prepare_labels(labels)
    label_vec = collect(labels)
    if isempty(label_vec)
        throw(ArgumentError("LightGBM requires at least one training example"))
    end

    if eltype(label_vec) <: Bool
        label_vec = Int.(label_vec)
    elseif eltype(label_vec) <: Integer
        label_vec = Int.(label_vec)
    elseif eltype(label_vec) <: AbstractFloat
        label_vec = Int.(round.(label_vec))
    else
        throw(ArgumentError("Unsupported label type $(eltype(label_vec)) for LightGBM"))
    end

    if any(x -> x ∉ (0, 1), label_vec)
        throw(ArgumentError("LightGBM requires binary labels encoded as 0 or 1."))
    end

    return label_vec
end

# Only for completed matrix fits whose training datasets are owned by this booster.
function _detach_lightgbm_training_data!(model::LightGBM.LGBMClassification)
    trained = model.booster
    isempty(trained.datasets) && return model
    # LightGBM.jl reloads only the trees when copying a booster.
    model.booster = deepcopy(trained)
    Base.finalize(trained)
    foreach(Base.finalize, trained.datasets)
    empty!(trained.datasets)
    return model
end

function fit_lightgbm_model(model::LightGBM.LGBMClassification,
                            feature_data::AbstractDataFrame,
                            labels::AbstractVector;
                            positive_label = true)
    features = Symbol.(names(feature_data))
    X = feature_matrix(feature_data, features)
    y_int = _prepare_labels(labels)

    unique_labels = unique(y_int)
    if length(unique_labels) == 1
        constant_prob = unique_labels[1] == 0 ? 0.0f0 : 1.0f0
        return LightGBMModel(nothing, features, constant_prob)
    end

    LightGBM.fit!(model, X, y_int; verbosity = -1)
    _detach_lightgbm_training_data!(model)
    return LightGBMModel(model, features, nothing)
end

function lightgbm_predict(model::LightGBMModel,
                          feature_data::AbstractDataFrame;
                          output_type = Float64)
    if model.booster === nothing
        n = nrow(feature_data)
        prob = model.constant_prediction === nothing ? 0.0f0 : model.constant_prediction
        return fill(convert(output_type, prob), n)
    end

    X = feature_matrix(feature_data, model.features)
    raw = LightGBM.predict(model.booster, X)
    ŷ = ndims(raw) == 2 ? dropdims(raw; dims = 2) : raw
    return convert.(output_type, ŷ)
end

function lightgbm_feature_importances(model::LightGBMModel)
    if model.booster === nothing
        return nothing
    end

    try
        return LightGBM.gain_importance(model.booster)
    catch
        return nothing
    end
end

predict(model::LightGBMModel, df::AbstractDataFrame) =
    lightgbm_predict(model, df; output_type = Float32)

function importance(model::LightGBMModel)
    gains = lightgbm_feature_importances(model)
    return gains === nothing ? nothing : collect(zip(model.features, gains))
end
