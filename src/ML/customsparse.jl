#test = StructArray{FRAG}(([1, 2, 3], [1, 2, 3], Float32[0, 0, 0], [1, 2, 3]))
#function testLOOP(test::StructArray{FRAG, 1, NamedTuple{(:row, :col, :nz_val, :x), Tuple{Array{Int64, 1}, Array{Int64, 1}, Array{Float32, 1}, Array{Int64, 1}}}, Int64})

mutable struct SparseArray{Ti<:Integer,T<:AbstractFloat}
    n_vals::Int64
    m::Int64
    n::Int64
    rowval::Vector{Ti}
    colval::Vector{Ti}
    nzval::Vector{T}
    mask::Vector{Bool}
    matched::Vector{Bool}
    x::Vector{T}
    colptr::Vector{Ti}
end

function setMask!(sa::SparseArray{Ti, T}, i::Int64, v::Bool) where {Ti<:Integer,T<:AbstractFloat}
    sa.mask[i] = v
end

@inline function selectpivot!(v::SparseArray{Ti, T}, lo::Integer, hi::Integer, o::Ordering) where {Ti<:Integer, T<:AbstractFloat}
    @inbounds begin
        mi = Base.midpoint(lo, hi)

        # sort v[mi] <= v[lo] <= v[hi] such that the pivot is immediately in place
        if lt(o, v.colval[lo], v.colval[mi])
            v.colval[mi], v.colval[lo] = v.colval[lo], v.colval[mi]
            v.rowval[mi], v.rowval[lo] = v.rowval[lo], v.rowval[mi]
            v.nzval[mi], v.nzval[lo] = v.nzval[lo], v.nzval[mi]
            v.x[mi], v.x[lo] = v.x[lo], v.x[mi]
            v.matched[mi], v.matched[lo] = v.matched[lo], v.matched[mi]
        end

        if lt(o, v.colval[hi], v.colval[lo])
            if lt(o, v.colval[hi], v.colval[mi])
                #v[hi], v[lo], v[mi] = v[lo], v[mi], v[hi]

                v.colval[hi], v.colval[lo], v.colval[mi] = v.colval[lo], v.colval[mi], v.colval[hi]
                v.rowval[hi], v.rowval[lo], v.rowval[mi] = v.rowval[lo], v.rowval[mi], v.rowval[hi]
                v.nzval[hi], v.nzval[lo], v.nzval[mi] = v.nzval[lo], v.nzval[mi], v.nzval[hi]
                v.x[hi], v.x[lo], v.x[mi] = v.x[lo], v.x[mi], v.x[hi]
                v.matched[hi], v.matched[lo], v.matched[mi] = v.matched[lo], v.matched[mi], v.matched[hi]

            else
                #v[hi], v[lo] = v[lo], v[hi]

                v.colval[hi], v.colval[lo] = v.colval[lo], v.colval[hi]
                v.rowval[hi], v.rowval[lo] = v.rowval[lo], v.rowval[hi]
                v.nzval[hi], v.nzval[lo] = v.nzval[lo], v.nzval[hi] 
                v.x[hi], v.x[lo] = v.x[hi], v.x[lo] 
                v.matched[hi], v.matched[lo] = v.matched[hi], v.matched[lo] 
            end
        end

        # return the pivot
        return v.colval[lo], v.rowval[lo], v.nzval[lo], v.x[lo], v.matched[lo]
    end
end

function partition!(v::SparseArray{Ti, T}, lo::Integer, hi::Integer, o::Ordering) where {Ti<:Integer, T<:AbstractFloat}
    pivot = selectpivot!(v, lo, hi, o)
    # pivot == v[lo], v[hi] > pivot
    i, j = lo, hi
    @inbounds while true
        i += 1; j -= 1
        while lt(o, v.colval[i], pivot[1]); i += 1; end;
        while lt(o, pivot[1], v.colval[j]); j -= 1; end;
        i >= j && break
        v.colval[i], v.colval[j] = v.colval[j], v.colval[i]
        v.rowval[i], v.rowval[j] = v.rowval[j], v.rowval[i]
        v.nzval[i], v.nzval[j] = v.nzval[j], v.nzval[i]
        v.x[i], v.x[j] = v.x[j], v.x[i]
        v.matched[i], v.matched[j] = v.matched[j], v.matched[i]
    end
    v.colval[j], v.colval[lo] = pivot[1], v.colval[j]
    v.rowval[j], v.rowval[lo] = pivot[2], v.rowval[j]
    v.nzval[j], v.nzval[lo] = pivot[3], v.nzval[j]
    v.x[j], v.x[lo] = pivot[4], v.x[j]
    v.matched[j], v.matched[lo] = pivot[5], v.matched[j]

    # v[j] == pivot
    # v[k] >= pivot for k > j
    # v[i] <= pivot for i < j
    return j
end

function specialsort!(v::SparseArray{Ti, T}, lo::Integer, hi::Integer, o::Ordering) where {Ti<:Integer, T<:AbstractFloat}
    @inbounds while lo < hi
        hi-lo <= 20 && return smallsort!(v, lo, hi, o)
        j = partition!(v, lo, hi, o)
        if j-lo < hi-j
            # recurse on the smaller chunk
            # this is necessary to preserve O(log(n))
            # stack space in the worst case (rather than O(n))
            lo < (j-1) && specialsort!(v, lo, j-1, o)
            lo = j+1
        else
            j+1 < hi && specialsort!(v, j+1, hi, o)
            hi = j-1
        end
    end
    return v
end

function smallsort!(v::SparseArray{Ti, T}, lo::Int64, hi::Int64, o::Ordering) where {Ti<:Integer, T<:AbstractFloat}
    #getkw lo hi scratch
    lo_plus_1 = (lo + 1)::Int64
    @inbounds for i = lo_plus_1:hi
        j = i
        col_x = v.colval[i]
        row_x = v.rowval[i]
        nzval_x = v.nzval[i]
        x_x = v.x[i]
        matched_x = v.matched[i]
        while j > lo
            #y = v[j-1]
            col_y = v.colval[j - 1]
            row_y = v.rowval[j - 1]
            nzval_y = v.nzval[j - 1]
            x_y = v.x[j - 1]
            matched_y = v.matched[j - 1]
            if !(lt(o, col_x, col_y)::Bool)
                break
            end
            v.colval[j] = col_y
            v.rowval[j] = row_y
            v.nzval[j] = nzval_y
            v.x[j] = x_y
            v.matched[j] = matched_y
            j -= 1
        end
        v.colval[j] = col_x
        v.rowval[j] = row_x
        v.nzval[j] = nzval_x
        v.x[j] = x_x
        v.matched[j] = matched_x
    end
    #scratch
end

function getRowVal(sa::SparseArray{Ti, T}, i::Int64) where {Ti<:Integer,T<:AbstractFloat}
    return sa.rowval[i]
end

function getColVal(sa::SparseArray{Ti, T}, i::Int64) where {Ti<:Integer,T<:AbstractFloat}
    return sa.colval[i]
end

function getNzVal(sa::SparseArray{Ti, T}, i::Int64) where {Ti<:Integer,T<:AbstractFloat}
    return sa.nzval[i]
end

function getRowColNZ(sa::SparseArray{Ti, T}, i::Int64) where {Ti<:Integer,T<:AbstractFloat}
    return sa.row_col_nzval_x[i]
end

SparseArray(N::Int) = SparseArray(
                    0,
                    0,
                    0,
                    Vector{Int64}(undef, N), #rowval 
                    Vector{Int64}(undef, N), #colval 
                    Vector{Float32}(undef, N), #nzval
                    ones(Bool, N), #mask
                    ones(Bool, N), #matched,
                    Vector{Float32}(undef, N), #x
                    Vector{Int64}(undef, N), #colptr
                    
)

function reset!(sa::SparseArray{Ti,T}) where {Ti<:Integer,T<:AbstractFloat}
    @turbo for i in range(1, sa.n_vals)
        sa.colval[i] = 0
        sa.rowval[i] = 0
        sa.x[i] = zero(T)
        sa.nzval[i] = zero(T)
        sa.mask[i] = true
        sa.matched[i] = true
        sa.colptr[i] = 0
    end
    sa.n_vals = 0
    sa.m = 0
    sa.n = 0
    return 
end

function sortSparse!(sa::SparseArray{Ti,T}) where {Ti<:Integer,T<:AbstractFloat}
    #Get sorted indices by column
    #partialsort!(sa.row_col_nzval[1:sa.n_vals], 
    #                    1:sa.n_vals, 
    #                    by = x -> x[2])
    specialsort!(sa, 1, sa.n_vals, Base.Order.Forward);
    max_col = 1
    max_row  = 1
    sa.colptr[1] = 1
    for i in range(1, sa.n_vals - 1)
        #If row greater than max row
        if sa.rowval[i + 1] > max_row
            max_row = sa.rowval[i + 1][1] 
        end 

        if sa.colval[i + 1] == sa.colval[i]
            continue
        else
            max_col += 1
            sa.colptr[max_col] = i + 1
        end
    end
    sa.colptr[max_col + 1] = sa.n_vals
    #Number of columns
    sa.n = max_col
    sa.m = max_row
end

function initResiduals!( r::Vector{T}, sa::SparseArray{Ti,T}, w::Vector{T}) where {Ti<:Integer,T<:AbstractFloat}
    #resize if necessary
    if length(r) < sa.n_vals
        append!(r, zeros(T, sa.n_vals - length(r)))
    end


    for col in range(1, sa.n)
        start = sa.colptr[col]
        stop = sa.colptr[col+1] - 1
        @turbo for n in start:stop
            r[sa.rowval[n]] = -sa.x[n]
        end
    end

    for col in range(1, sa.n)
        start = sa.colptr[col]
        stop = sa.colptr[col+1] - 1
        @turbo for n in start:stop
            r[sa.rowval[n]] += w[col]*sa.nzval[n]
        end
    end
end

function multiply!(sa::SparseArray{Ti,T}, w::Vector{T}, r::Vector{T}) where {Ti<:Integer,T<:AbstractFloat}
    #resize if necessary
    if length(r) < length(sa.n_vals)
        append!(r, zeros(T, length(sa.n_vals) - length(r)))
    end

    for col in range(1, sa.n)
        start = sa.colptr[col]
        stop = sa.colptr[col+1] - 1
        @turbo for n in start:stop
            r[sa.rowval[n]] += w[col]*sa.nzval[n]
        end
    end
end

struct SpectralScores{T<:AbstractFloat}
    scribe::T
    scribe_corrected::T
    scribe_fitted::T
    city_block::T
    city_block_fitted::T
    spectral_contrast::T
    spectral_contrast_corrected::T
    matched_ratio::T
    entropy_score::T
end

#=
Vector{SpectralScores{Float16}}(undef, 10)
sa_test = SparseArray(1000)
X_test, precID_to_col, H_ncol = buildDesignMatrix(ionMatches, ionMisses, nmatches, nmisses, sa_test)

include("src/ML/sparseNNLS.jl")
wt = ones(Float32, sa_test.n);
r = zeros(Float32, 10000)
for i in range(1, sa_test.m)
    r[i] = r[i] - X[i]
end
multiply!(sa_test, wt, r)
@time solveHuber!(sa_test, r, wt, Float32(10000), max_iter_outer = 100, max_iter_inner = 20, tol = sa_test.n);

test_scores = Vector{SpectralScores{Float16}}(undef, 100000)
getDistanceMetrics(wt, sa_test, test_scores)
test_unscored_psms = [LXTandem(Float32) for x in 1:10000]
ScoreFragmentMatches!(test_unscored_psms, precID_to_col, ionMatches, nmatches, Laplace(-1, 3.9))
include("src/PSM_TYPES/LibraryIntensity.jl")
test_scored_psms = Vector{LibPSM{Float32, Float16}}(undef, 1)
Score!(test_scored_psms,
       test_unscored_psms, 
       test_scores,
       wt,
       Float64(0.1),
       sa_test.n,
       Float32(10000))
test_scored_psms[1:10]

unscored_PSMs = UnorderedDictionary{UInt32, XTandem{Float32}}()

test_tandem = [LXTandem(Float32) for x in 1:10000]

test_IDtoCOL = UnorderedDictionary{UInt32, Int64}()
for (key, val) in pairs(precID_to_col)
    insert!(test_IDtoCOL, key, first(val))
end


for (key, val) in pairs(precID_to_col)
    delete!(precID_to_col, key)
end

ScoreFragmentMatches!(test_tandem, test_IDtoCOL, ionMatches, nmatches, err_dist)


[XTandem(Float32) for x in 1:10]



weights = ones(Float32, sa_test.n);
r = Hs*weights .- X;
@time solveHuber!(Hs, r, weights, Float32(10000), max_iter_outer = 100, max_iter_inner = 20, tol = Hs.n);

plot!(log10.(wt), log10(weights))

@time for i in (1, 10000)
    weights = ones(Float32, sa_test.n);
    solveHuber!(Hs, Hs*weights .- X, weights, huber_δ, max_iter_outer = 100, max_iter_inner = 20, tol = Hs.n);
    a += 1
end


a = 0
@time for i in (1, 10000)
    wt = ones(Float32, sa_test.n)
    residuals = zeros(Float32, 10000)
    for i in range(1, sa_test.m)
        residuals[i] = residuals[i] - X[i]
    end
    multiply!(sa_test, wt, residuals)
    old_wt = copy(wt)
    @time solveHuber!(sa_test, residuals, wt, huber_δ, max_iter_outer = 100, max_iter_inner = 10, tol = sa_test.n);
    plot(log2.(old_wt), log2.(wt), seriestype = :scatter)
    a += 1
end

weights = ones(Float32, sa_test.n);
residuals = Hs*weights .- X;

old_wt = copy(weights)
solveHuber!(Hs, Hs*weights .- X, weights, huber_δ, max_iter_outer = 1, max_iter_inner = 20, tol = Hs.n);
max_diff = 0.0
for i in range(1, length(old_wt))
    if old_wt != false
        if (old_wt[i] - weights[i])/old_wt[i] > max_diff
            max_diff = (old_wt[i] - weights[i])/old_wt[i]
        end
    end
end
println(max_diff)
argmax((old_wt .- weights)./max.(old_wt, 100))

max_diff = 0.0
for i in range(1, length(old_wt))
    if iszero
argmax((old_wt .- weights)./max.(old_wt, 100))
plot(log2.(old_wt), log2.(wt), seriestype = :scatter)
Hs[:,23]
X[Hs[:,23].!=0.0]

struct LibPSM{H,L<:Real} <: PSM
    #H is "high precision"
    #L is "low precision"

    #Ion Count Statistics
    best_rank::UInt8 #Highest ranking predicted framgent that was observed
    topn::UInt8 #How many of the topN predicted fragments were observed. 
    longest_y::UInt8
    longest_b::UInt8
    b_count::UInt8
    y_count::UInt8

    #Basic Metrics 
    poisson::L
    hyperscore::L
    log2_intensity_explained::L
    error::H

    #Spectral Simmilarity
    scribe::L
    scribe_corrected::L
    city_block::L
    spectral_contrast::L
    spectral_contrast_corrected::L
    matched_ratio::L
    entropy_score::L
    weight::H

    #Non-scores/Labels
    precursor_idx::UInt32
    ms_file_idx::UInt32
    scan_idx::UInt32
end
#=

jldsave("/Users/n.t.wamsley/Desktop/test.jld2";  X, Hs, IDtoCOL, last_matched_col, ionMatches, ionMisses, nmatches, nmisses)
@load "/Users/n.t.wamsley/Desktop/test.jld2" X
@load "/Users/n.t.wamsley/Desktop/test.jld2" Hs
@load "/Users/n.t.wamsley/Desktop/test.jld2" IDtoCOL
@load "/Users/n.t.wamsley/Desktop/test.jld2" last_matched_col
@load "/Users/n.t.wamsley/Desktop/test.jld2" ionMatches
@load "/Users/n.t.wamsley/Desktop/test.jld2" ionMisses
@load "/Users/n.t.wamsley/Desktop/test.jld2" nmatches
@load "/Users/n.t.wamsley/Desktop/test.jld2" nmisses


    sa_test = SparseArray(1000)
    X_test, precID_to_col, H_ncol = buildDesignMatrix(ionMatches, ionMisses, nmatches, nmisses, sa_test)
    sortSparse!(sa_test)
    wt = ones(Float32, sa_test.n)
    residuals = zeros(Float32, 10000)

    multiply!(sa_test, wt, residuals)


    multiply!(sa_test, )
    mult(x)
solveHuber!(Hs, Hs*weights .- X, weights, Float32(10000), max_iter_outer = 100, max_iter_inner = 20, tol = Hs.n);

errs = []
for i in range(1, sa_test.n_vals)
    row, col, nzval,x  = sa_test.row_col_nzval_x[i]
    push!(errs, abs(Hs[row, col] - nzval))
end

=#
err = []
for (key, value) in pairs(precID_to_row)
    if first(value) == 17
        println(key)
    end
    push!(err, first(value))
    #println(Int64(first(value)))
end

4942222
1954710

precID_to_row[4942222]

precID_to_row[1954710]



struct LibPSM{H,L<:Real} <: PSM
    #H is "high precision"
    #L is "low precision"

    #Ion Count Statistics
    best_rank::UInt8 #Highest ranking predicted framgent that was observed
    topn::UInt8 #How many of the topN predicted fragments were observed. 
    longest_y::UInt8
    longest_b::UInt8
    b_count::UInt8
    y_count::UInt8

    #Basic Metrics 
    poisson::L
    hyperscore::L
    log2_intensity_explained::L
    error::H

    #Spectral Simmilarity
    scribe::L
    scribe_corrected::L
    city_block::L
    spectral_contrast::L
    spectral_contrast_corrected::L
    matched_ratio::L
    entropy_score::L
    weight::H

    #Non-scores/Labels
    precursor_idx::UInt32
    ms_file_idx::UInt32
    scan_idx::UInt32
end

=#