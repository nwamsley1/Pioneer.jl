module PioneerAOTCanary

export simd_find_first_ge, weighted_sum

# Mirror the fixed-width LLVM operation used by Pioneer's fragment-index and
# fused-scan kernels. Its 256-bit vector must lower correctly in every CPU tier.
const F32x8 = NTuple{8, Core.VecElement{Float32}}

@inline _vbroadcast8(x::Float32)::F32x8 =
    ntuple(_ -> Core.VecElement(x), Val(8))

@inline _vload8(values::Vector{Float32}, i::Int)::F32x8 =
    unsafe_load(Ptr{F32x8}(pointer(values, i)))

@inline _vcmpge_mask(a::F32x8, b::F32x8)::UInt8 =
    Core.Intrinsics.llvmcall("""
        %cmp = fcmp oge <8 x float> %0, %1
        %mask = bitcast <8 x i1> %cmp to i8
        ret i8 %mask
    """, UInt8, Tuple{F32x8, F32x8}, a, b)

@inline function simd_find_first_ge(
    values::Vector{Float32},
    first_index::Int,
    last_index::Int,
    target::Float32,
)
    target_vector = _vbroadcast8(target)
    i = first_index
    @inbounds while i + 7 <= last_index
        mask = _vcmpge_mask(_vload8(values, i), target_vector)
        mask != 0x00 && return i + trailing_zeros(mask)
        i += 8
    end
    @inbounds while i <= last_index
        values[i] >= target && return i
        i += 1
    end
    return last_index + 1
end

# Representative of the Base @simd reductions that replaced Pioneer's
# LoopVectorization kernels. LLVM can specialize this method per sysimage tier.
function weighted_sum(x::Vector{Float64}, y::Vector{Float64})
    total = 0.0
    @inbounds @simd for i in eachindex(x, y)
        total += x[i] * y[i]
    end
    return total
end

end
