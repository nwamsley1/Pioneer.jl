const PIONEER_SIMD_PKGID = Base.PkgId(Base.UUID("05f169b2-1f58-4f3e-a5d2-7bfdd2b4f8e3"), "PioneerSIMD")
const PIONEER_SIMD_LOADED = Ref(false)

function _default_arraydict_reset!(keys, vals, size)
    @inbounds for i in 1:size
        vals[keys[i]] = zero(eltype(vals))
        keys[i] = zero(eltype(keys))
    end
    return nothing
end

function _default_counter_reset!(ids, counts, last_index)
    @inbounds for i in 1:last_index
        counts[ids[i]] = zero(eltype(counts))
        ids[i] = zero(eltype(ids))
    end
    return nothing
end

function _default_sparsearray_reset!(colval, rowval, x, nzval, matched, colptr, isotope, n_vals)
    @inbounds for i in 1:n_vals
        colval[i] = zero(eltype(colval))
        rowval[i] = zero(eltype(rowval))
        x[i] = zero(eltype(x))
        nzval[i] = zero(eltype(nzval))
        matched[i] = true
        colptr[i] = zero(eltype(colptr))
        isotope[i] = zero(eltype(isotope))
    end
    return nothing
end

function _default_fill_zero_chunk!(dest, chunk)
    @inbounds for i in chunk
        dest[i] = zero(eltype(dest))
    end
    return nothing
end

function _default_scaled_add_chunk!(dest, src, alpha, chunk)
    @inbounds @fastmath for i in chunk
        dest[i] += src[i] * alpha
    end
    return nothing
end

function _default_weighted_triple_sum(x, w, y)
    total = zero(promote_type(eltype(w), eltype(y)))
    @inbounds @fastmath for i in eachindex(w)
        total += x[i] * w[i] * y[i]
    end
    return total
end

function _default_sum_chunk(v, chunk)
    total = zero(eltype(v))
    @inbounds @fastmath for i in chunk
        total += v[i]
    end
    return total
end

function _default_sparse_axpy_rows!(dest, rowval, values, alpha, chunk)
    @inbounds @fastmath for n in chunk
        dest[rowval[n]] += alpha * values[n]
    end
    return nothing
end

const ARRAYDICT_RESET_IMPL = Ref{Function}(_default_arraydict_reset!)
const COUNTER_RESET_IMPL = Ref{Function}(_default_counter_reset!)
const SPARSEARRAY_RESET_IMPL = Ref{Function}(_default_sparsearray_reset!)
const FILL_ZERO_CHUNK_IMPL = Ref{Function}(_default_fill_zero_chunk!)
const SCALED_ADD_CHUNK_IMPL = Ref{Function}(_default_scaled_add_chunk!)
const WEIGHTED_TRIPLE_SUM_IMPL = Ref{Function}(_default_weighted_triple_sum)
const SUM_CHUNK_IMPL = Ref{Function}(_default_sum_chunk)
const SPARSE_AXPY_ROWS_IMPL = Ref{Function}(_default_sparse_axpy_rows!)

arraydict_reset_kernel!(keys, vals, size) = ARRAYDICT_RESET_IMPL[](keys, vals, size)
counter_reset_kernel!(ids, counts, last_index) = COUNTER_RESET_IMPL[](ids, counts, last_index)
sparsearray_reset_kernel!(colval, rowval, x, nzval, matched, colptr, isotope, n_vals) =
    SPARSEARRAY_RESET_IMPL[](colval, rowval, x, nzval, matched, colptr, isotope, n_vals)
fill_zero_chunk_kernel!(dest, chunk) = FILL_ZERO_CHUNK_IMPL[](dest, chunk)
scaled_add_chunk_kernel!(dest, src, alpha, chunk) = SCALED_ADD_CHUNK_IMPL[](dest, src, alpha, chunk)
weighted_triple_sum_kernel(x, w, y) = WEIGHTED_TRIPLE_SUM_IMPL[](x, w, y)
sum_chunk_kernel(v, chunk) = SUM_CHUNK_IMPL[](v, chunk)
sparse_axpy_rows_kernel!(dest, rowval, values, alpha, chunk) =
    SPARSE_AXPY_ROWS_IMPL[](dest, rowval, values, alpha, chunk)

function install_pioneer_simd_kernels!(;
    arraydict_reset_impl=nothing,
    counter_reset_impl=nothing,
    sparsearray_reset_impl=nothing,
    fill_zero_chunk_impl=nothing,
    scaled_add_chunk_impl=nothing,
    weighted_triple_sum_impl=nothing,
    sum_chunk_impl=nothing,
    sparse_axpy_rows_impl=nothing,
)
    !isnothing(arraydict_reset_impl) && (ARRAYDICT_RESET_IMPL[] = arraydict_reset_impl)
    !isnothing(counter_reset_impl) && (COUNTER_RESET_IMPL[] = counter_reset_impl)
    !isnothing(sparsearray_reset_impl) && (SPARSEARRAY_RESET_IMPL[] = sparsearray_reset_impl)
    !isnothing(fill_zero_chunk_impl) && (FILL_ZERO_CHUNK_IMPL[] = fill_zero_chunk_impl)
    !isnothing(scaled_add_chunk_impl) && (SCALED_ADD_CHUNK_IMPL[] = scaled_add_chunk_impl)
    !isnothing(weighted_triple_sum_impl) && (WEIGHTED_TRIPLE_SUM_IMPL[] = weighted_triple_sum_impl)
    !isnothing(sum_chunk_impl) && (SUM_CHUNK_IMPL[] = sum_chunk_impl)
    !isnothing(sparse_axpy_rows_impl) && (SPARSE_AXPY_ROWS_IMPL[] = sparse_axpy_rows_impl)
    return nothing
end

function load_pioneer_simd!()
    PIONEER_SIMD_LOADED[] && return true
    ccall(:jl_generating_output, Cint, ()) == 1 && return false
    try
        Base.require(PIONEER_SIMD_PKGID)
        PIONEER_SIMD_LOADED[] = true
        return true
    catch err
        if err isa ArgumentError
            return false
        end
        @warn "Failed to load PioneerSIMD; continuing with generic kernels" exception=(err, catch_backtrace())
        return false
    end
end
