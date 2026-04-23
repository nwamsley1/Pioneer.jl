module PioneerSIMD

using LoopVectorization
using PioneerCommon

const PIONEER_PKGID = Base.PkgId(Base.UUID("585db54e-1982-41c9-ae07-9e9eb56c7d61"), "Pioneer")

function _simd_arraydict_reset!(keys, vals, size)
    @turbo for i in 1:size
        vals[keys[i]] = zero(eltype(vals))
        keys[i] = zero(eltype(keys))
    end
    return nothing
end

function _simd_counter_reset!(ids, counts, last_index)
    @turbo for i in 1:last_index
        counts[ids[i]] = zero(eltype(counts))
        ids[i] = zero(eltype(ids))
    end
    return nothing
end

function _simd_sparsearray_reset!(colval, rowval, x, nzval, matched, colptr, isotope, n_vals)
    @turbo for i in 1:n_vals
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

function _simd_fill_zero_chunk!(dest, chunk)
    @turbo for i in chunk
        dest[i] = zero(eltype(dest))
    end
    return nothing
end

function _simd_scaled_add_chunk!(dest, src, alpha, chunk)
    @turbo for i in chunk
        dest[i] += src[i] * alpha
    end
    return nothing
end

function _simd_weighted_triple_sum(x, w, y)
    total = zero(promote_type(eltype(w), eltype(y)))
    @turbo for i in eachindex(w)
        total += x[i] * w[i] * y[i]
    end
    return total
end

function _simd_sum_chunk(v, chunk)
    total = zero(eltype(v))
    @turbo for i in chunk
        total += v[i]
    end
    return total
end

function _simd_sparse_axpy_rows!(dest, rowval, values, alpha, chunk)
    @turbo for n in chunk
        dest[rowval[n]] += alpha * values[n]
    end
    return nothing
end

function _install_into!(mod::Module)
    getfield(mod, :install_pioneer_simd_kernels!)(
        arraydict_reset_impl=_simd_arraydict_reset!,
        counter_reset_impl=_simd_counter_reset!,
        sparsearray_reset_impl=_simd_sparsearray_reset!,
        fill_zero_chunk_impl=_simd_fill_zero_chunk!,
        scaled_add_chunk_impl=_simd_scaled_add_chunk!,
        weighted_triple_sum_impl=_simd_weighted_triple_sum,
        sum_chunk_impl=_simd_sum_chunk,
        sparse_axpy_rows_impl=_simd_sparse_axpy_rows!,
    )
    return nothing
end

function __init__()
    _install_into!(PioneerCommon)
    haskey(Base.loaded_modules, PIONEER_PKGID) && _install_into!(Base.root_module(PIONEER_PKGID))
    return nothing
end

end
