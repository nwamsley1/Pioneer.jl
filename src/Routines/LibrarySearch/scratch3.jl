abstract type Node end
getVal(nd::Node) = nd.val
getBinID(nd::Node) = nd.bin_id
getTopID(nd::Node) = nd.top_id
getLastID(nd::Node) = nd.last_id

struct TTreeNode{I<:Unsigned,T<:Real} <: Node
    val::T
    bin_id::I
    top_id::I
    last_id::I
end

mutable struct TournamentTree{I<:Unsigned,T<:Real}
    nodes::Vector{TTreeNode{I,T}}
    n_bins::I
end
getNode(ttree::TournamentTree, node_id::UInt32) = ttree.nodes[node_id]
getNodes(ttree::TournamentTree) = ttree.nodes
setNode!(ttree::TournamentTree, node::TTreeNode, node_id::UInt32) = ttree.nodes[node_id] = node


function fillNode!(ttree::TournamentTree{I,T}, node_id::UInt32, val::T, bin_id::I, top_id::I, bin_len::I) where {I<:Unsigned, T<:Real}
    setNode!(ttree, 
                TTreeNode(val, bin_id, top_id, bin_len),
                node_id
    )
end

#Grow tree to accomadate N bins
function growTTree!(ttree::TournamentTree{I,T}, N::I) where {I<:Unsigned, T<:Real}
    ttree.nodes = Vector{TTreeNode{I, T}}(undef, N + N - 1)
    ttree.n_bins = N
end

#Find the smallest power of 2 greater than or equal
#to the query
function getNearestPowerOf2(n::UInt32)
    n -= one(UInt32)
    n |= n >> 1
    n |= n >> 2
    n |= n >> 4
    n |= n >> 8
    n |= n >> 16
    n += one(UInt32)
    return n
end


function buildTTree!(
    ttree::TournamentTree{I, T}, #Tournament Tree
    sub_bin_ranges::Vector{UnitRange{I}}, #Ranges of `values` that are each sorted in ascending order
    n_sub_bins::UInt32, #Number of sub bin ranges to consider
    values::Vector{T} #Values 
    ) where {I<:Unsigned,T<:Real}

    #Number of leaves needs to be the next highest power of 2. 
    N = getNearestPowerOf2(n_sub_bins)
    if N+N-one(UInt32) > length(getNodes(ttree))
        growTTree!(ttree, N)
    end
    #Fill leaves (first-layer nodes) 
    leaf_idx = one(UInt32)
    for sub_bin_range in sub_bin_ranges
        top_id = first(sub_bin_range)
        last_id = last(sub_bin_range)
        fillNode!(ttree, 
                    leaf_idx, #node_id
                    values[top_id], #val
                    leaf_idx, #bin_id
                    top_id,
                    last_id
                    )
        leaf_idx += one(UInt32)
    end
    for leaf_idx in range(n_sub_bins+one(UInt32), N) #Fill placeholder leaves (if last_sub_bin_idx is not a power of 2)
        fillNode!(ttree, 
                leaf_idx,
                typemax(T),
                leaf_idx, #leaf_id
                one(I),
                one(I)
        )
    end
    #Fill remaining nodes 
    lower_node_id = 1
    upper_node_id = N + 1
    n = N >> 1
    while n > 1#node_id <= (N - 1) #While there are nodes left to fill
        start = upper_node_id
        n = n >> 1
        while upper_node_id <= start + n + 1 #Fill current layer of nodes 

            node_a, node_b = ttree.nodes[lower_node_id], ttree.nodes[lower_node_id+1]

            if getVal(node_a) <= getVal(node_b)
                ttree.nodes[upper_node_id] = node_a
            else
                ttree.nodes[upper_node_id] = node_b
            end
            lower_node_id += 2
            upper_node_id += 1
        end
    end

end

function removeSmallestElement!(
    ttree::TournamentTree{I, T},
    values_in::Vector{T},
    values_out::Vector{T},
    out_idx::Int64
) where {I<:Unsigned,T<:Real}
    #Number of bins in the merge 
    n = ttree.n_bins
    #Get node containing minimum value 
    leading_node = getNode(ttree, n + n - one(UInt32))
    #Write minimum value to output array and increase index 
    values_out[out_idx] = getVal(leading_node)
    out_idx += 1
    #Repair tree 
    top_id,last_id,bin_id = getTopID(leading_node),getLastID(leading_node),getBinID(leading_node)
    if top_id === last_id #No more elements left in the bin
        setNode!(ttree, 
                TTreeNode(
                typemax(T), #No more elements left, so make sure this bin always loses 
                bin_id,
                top_id,
                last_id),
                bin_id)
    else
        top_id += one(I)
        setNode!(ttree, 
                TTreeNode(
                values_in[top_id],
                bin_id,
                top_id,
                last_id),
                bin_id)
    end
    #If the bin_id is even, subtract one
    lower_node_id = bin_id - (bin_id%I(2) === zero(I))
    N = ttree.n_bins
    n = (lower_node_id >> 1)
    while N > 1
        node_a, node_b = getNode(ttree, lower_node_id), getNode(ttree, lower_node_id + one(I))
        lower_node_id += N - n
        if getVal(node_a) <= getVal(node_b)
            setNode!(ttree, node_a,  lower_node_id)
        else
            setNode!(ttree, node_b,  lower_node_id)
        end
        lower_node_id = lower_node_id - (lower_node_id%I(2) === zero(I))
        N = N >> 1
        n = n >> 1
    end
end
using PProf, Profile
Profile.clear()
MS_TABLE = Arrow.Table(MS_TABLE_PATH)
@time begin
@profile RESULT = quantitationSearch(MS_TABLE, 
prosit_lib["precursors"],
prosit_lib["f_det"],
RT_INDICES[file_id_to_parsed_name[ms_file_idx]],
UInt32(ms_file_idx), 
RT_iRT[file_id_to_parsed_name[ms_file_idx]],
frag_err_dist_dict[ms_file_idx],
irt_err,#irt_errs[ms_file_idx]/3,
params_,  
ionMatches,
ionMisses,
IDtoCOL,
ionTemplates,
iso_splines,
chromatograms,
complex_scored_PSMs,
complex_unscored_PSMs,
complex_spectral_scores,
precursor_weights,
);

end
pprof(;webport=58603)


using PProf, Profile
Profile.clear()
test_spline = BSplineApprox(test_shapes, 
test_intensities, 3, 4, :ArcLen, :Uniform, extrapolate = true)
@profile @btime test_spline(10.0f0) 
pprof(;webport=58603)


N = 100
a = rand(100)
p = zeros(UInt32, N)
sortperm!(@view(p[1:20]),@view(a[1:20]),alg = PartialQuickSort(1:4))
issorted([a[i] for i in p[1:20]])

PSMs[!,:prec_mz] = [precursors[:mz][x] for x in PSMs[!,:precursor_idx]]
println("max observed mz ", maximum(PSMs[!,:prec_mz]))
println("min observed mz ", minimum(PSMs[!,:prec_mz]))
MS_TABLE[:centerMass][PSMs[!,:scan_idx][1]] -  MS_TABLE[:isolationWidth][PSMs[!,:scan_idx][1]]/2
MS_TABLE[:centerMass][PSMs[!,:scan_idx][1]] +  MS_TABLE[:isolationWidth][PSMs[!,:scan_idx][1]]/2

ms2_scans = MS_TABLE[:msOrder].==2
MS_TABLE[:lowMass][ms2_scans]

H, Hs = RESULT[1]
n = H.n_vals
H.idx[1:30]
Hs.idx[1:30]

MS2_CHROMS = groupby(PSMS, [:precursor_idx]);


MS2_CHROMS[N][!,[:RT,:weight,:b_count,:y_count,:scribe_fitted,:city_block_fitted,:spectral_contrast,:matched_ratio,:iso_rank]]
plot(MS2_CHROMS[N][!,:RT],
MS2_CHROMS[N][!,:weight], seriestype=:scatter)
N += 1
b = iso_splines.splines[1][1].knots[1:end - 1]
A = hcat(ones(length(b)), b)
x = A\[x for x in range(0, length(b)-1)]
ceil(Int, (1075 - 75)/1000)

y = 0.0
ceil(Int, first(x)+ y*last(x))


y = 80.0
ceil(Int, first(x)+ y*last(x))

y = 1080.0
ceil(Int, first(x)+ y*last(x))

plot(collect(LinRange(0, 10000.0, 10000)), [iso_splines.splines[1][1](x) for x in LinRange(0, 15000.0, 10000)])
plot!(collect(LinRange(0, 10000.0, 10000)), [iso_splines.splines[1][2](x) for x in LinRange(0, 15000.0, 10000)])
plot!(collect(LinRange(0, 10000.0, 10000)), [iso_splines.splines[1][3](x) for x in LinRange(0, 15000.0, 10000)])

plot(collect(LinRange(0, 3000.0, 10000)), [iso_splines.splines[1][1](x) for x in LinRange(0, 3000.0, 10000)])
plot!(collect(LinRange(0, 3000.0, 10000)), [iso_splines.splines[1][2](x) for x in LinRange(0, 3000.0, 10000)])
plot!(collect(LinRange(0, 3000.0, 10000)), [iso_splines.splines[1][3](x) for x in LinRange(0, 3000.0, 10000)])
vline!([75.0])

plot(collect(LinRange(0, 3000.0, 10000)), [iso_splines.splines[3][1](x) for x in LinRange(0, 3000.0, 10000)])
plot!(collect(LinRange(0, 3000.0, 10000)), [iso_splines.splines[3][2](x) for x in LinRange(0, 3000.0, 10000)])
plot!(collect(LinRange(0, 3000.0, 10000)), [iso_splines.splines[3][3](x) for x in LinRange(0, 3000.0, 10000)])
vline!([75.0])

[iso_splines.splines[3][2](1973) ]
iso_splines.splines[1][1].polynomials


plot(collect(LinRange(0, 5000.0, 10000)), [iso_splines.splines[1][1].polynomials[1](x) for x in LinRange(0, 5000.0, 10000)])
plot!(collect(LinRange(0, 5000.0, 10000)), [iso_splines.splines[1][1].polynomials[2](x) for x in LinRange(0, 5000.0, 10000)])
plot!(collect(LinRange(0, 5000.0, 10000)), [iso_splines.splines[1][1].polynomials[3](x) for x in LinRange(0, 5000.0, 10000)])
plot!(collect(LinRange(0, 5000.0, 10000)), [iso_splines.splines[1][1].polynomials[4](x) for x in LinRange(0, 5000.0, 10000)])
plot!(collect(LinRange(0, 5000.0, 10000)), [iso_splines.splines[1][1].polynomials[5](x) for x in LinRange(0, 5000.0, 10000)])

bins = LinRange(0, 5, 100)
histogram(PSMS[PSMS[!,:decoy],:entropy_score], alpha = 0.5, bins = bins)
histogram!(PSMS[PSMS[!,:target],:entropy_score], alpha = 0.5, bins = bins)


bins = LinRange(-0.5, 2, 100)
histogram(PSMS[PSMS[!,:decoy],:city_block_fitted], alpha = 0.5, bins = bins)
histogram!(PSMS[PSMS[!,:target],:city_block_fitted], alpha = 0.5, bins = bins)

test = quantitationSearch(MS_TABLE, 
prosit_lib["precursors"],
prosit_lib["f_det"],
RT_INDICES[file_id_to_parsed_name[ms_file_idx]],
UInt32(ms_file_idx), 
frag_err_dist_dict[ms_file_idx],
irt_errs[ms_file_idx],
params_,  
ionMatches,
ionMisses,
IDtoCOL,
ionTemplates,
iso_splines,
complex_scored_PSMs,
complex_unscored_PSMs,
complex_spectral_scores,
precursor_weights,
)

testHX = test[10]

hcat([testHX.nzval[range(testHX.colptr[104], testHX.colptr[105]-1)],
    testHX.matched[range(testHX.colptr[104], testHX.colptr[105]-1)],
    testHX.mask[range(testHX.colptr[104], testHX.colptr[105]-1)],
    testHX.x[range(testHX.colptr[104], testHX.colptr[105]-1)],
    testHX.rowval[range(testHX.colptr[104], testHX.colptr[105]-1)]])


    bins = LinRange(0, 2, 100)
    histogram(best_psms[best_psms[!,:target].&(best_psms[!,:q_value].<=0.01), :entropy_score], alpha = 0.5, bins = bins, normalize = :pdf)
    histogram!(best_psms[best_psms[!,:target].==false, :entropy_score], alpha = 0.5, bins = bins, normalize=:pdf)
    
PSM_FIRST = copy(PSMs_Dict[""]) # sum(PSM_FIRST[!,:q_value].<=0.01)
#=
julia> sum(PSM_FIRST[!,:q_value].<=0.01)
95119
=#
best_psms_a = copy(best_psms)
#=
julia> value_counts(best_psms_a[(best_psms_a[:,:q_value].<=0.01) .& (best_psms_a[:,:decoy].==false),:], [:file_name])
1×2 DataFrame
Row │ file_name  nrow  
    │ String     Int64 
─────┼──────────────────
1 │            97628
=#
post_quant_set = Set(best_psms_a[(best_psms_a[:,:q_value].<=0.01) .& (best_psms_a[:,:decoy].==false),:precursor_idx])
pre_quant_set = Set(PSM_FIRST[PSM_FIRST[!,:q_value].<=0.01,:precursor_idx])

setdiff(post_quant_set, pre_quant_set) #13974

missing_from_qaunt = collect(setdiff(pre_quant_set, post_quant_set)) #11465 Why are these missing
PSMS[!,:prec_mz] = [MS_TABLE[:centerMass][i] for i in PSMS[!,:scan_idx]]
N = 1000

MS2_CHROMS[(precursor_idx = missing_from_qaunt[N],iso_rank = 1)][!,
[:b_count,:y_count,:isotope_count,:scribe,:spectral_contrast,:matched_ratio,:city_block_fitted,:entropy_score,:max_entropy,:scan_idx,:prec_mz,:RT,:weight,:peak_area,:target]]
plot(MS2_CHROMS[(precursor_idx = missing_from_qaunt[N],iso_rank = 1)][!,:RT], 
MS2_CHROMS[(precursor_idx = missing_from_qaunt[N],iso_rank = 1)][!,:weight], seriestype=:scatter)

dtype = Float32;
gx, gw = gausslegendre(100);
state = GD_state(
    HuberParams(zero(dtype), zero(dtype),zero(dtype),zero(dtype)), #Initial params
    zeros(dtype, N), #t
    zeros(dtype, N), #y
    zeros(dtype, N), #data
    falses(N), #mask
    0, #number of iterations
    N #max index
    );
integratePrecursorMS2(MS2_CHROMS[(precursor_idx = missing_from_qaunt[N],iso_rank = 1)],
state,
gx::Vector{Float64},
gw::Vector{Float64},
intensity_filter_fraction =  Float32(params_[:integration_params]["intensity_filter_threshold"]),
α = 0.001f0,
half_width_at_α = 0.15f0,
isplot = true
);
best_psms_a[best_psms_a[!,:precursor_idx] .== missing_from_qaunt[N],[:precursor_idx,:scribe,:entropy_score,:weight,:q_value]]
N += 1

intensities, shapes = ModelMassErrs(
           frag_ppm_intensities,
           frag_ppm_errs,
           Float64(max_ppm),#params_[:presearch_params]["frag_tol_ppm"],
           n_intensity_bins = length(frag_ppm_errs)÷1500,#Int64(params_[:presearch_params]["samples_per_mass_err_bin"]),
           frag_err_quantile = 0.975,#params_[:frag_tol_params]["frag_tol_quantile"],
           out_fdir = mass_err_estimation_folder,
           out_fname = out_fname
       )


m6 = rlm(Float64.(intensities), Float64.(shapes), TauEstimator{TukeyLoss}(); initial_scale=:mad)


m, b = RobustModels.coef(m6)


#m, b = intensities[10:20,:]\(shapes[10:20])

plot(2 .^intensities[:,1], 2 .^shapes, seriestype=:scatter)
bins = LinRange(minimum(intensities[:,1]), maximum(intensities[:,1]), 1000)
plot!(2 .^bins, 2 .^[m*x + b for x in bins])


plot(intensities[:,1], 2 .^shapes, seriestype=:scatter)
bins = LinRange(minimum(intensities[:,1]), maximum(intensities[:,1]), 1000)
plot!(bins, 2 .^[m*x + b for x in bins])
hline!([2 .^shapes[end]])
hline!([2 .^shapes[1]])



S_interp = approximate(f, B, ApproxByInterpolation(B))  # or simply approximate(f, B)

using BSplineKit
x_interval = 0
ξs = range(0, 1; length = 15)
B = BSplineBasis(BSplineOrder(4), ξs)


test_interp = LinearInterpolation(ξs, 
                                    sin.(ξs),
                                    #savitzky_golay(ys, w, 3).y, 
                                    extrapolation_bc = Line())


S_minL2 = approximate(test_interp, B, MinimiseL2Error())

function testSplineFunc(a::BSplineKit.SplineApproximations.SplineApproximation,b::Float32)
    return a(b)
end

quantile(d::Laplace, p::Real) = p < 1/2 ? xval(d, log(2p)) : xval(d, -log(2(1 - p)))

function testfunc(a)
    (mean(a) - 0.9924999999999999)*2262.443438914029 + 2000
end

function testfunc(a)
    (mean(a) - 0.5569999999999999)*2222 + 1000
end



plot(LinRange(0, 4000.0, 100), iso_splines.splines[4][1].(LinRange(0, 4000.0, 100)))
plot!(LinRange(0, 4000.0, 100), iso_splines.splines[4][2].(LinRange(0, 4000.0, 100)))
plot!(LinRange(0, 4000.0, 100), iso_splines.splines[4][3].(LinRange(0, 4000.0, 100)))
vline!([384.499, 384.499 + 1000.0f0])
plot(LinRange(0, 4000.0, 100), iso_splines[1][1].(LinRange(0, 4000.0, 100)))

test_p = Polynomial(iso_splines[1][1].coeffs[1:4])
plot!(LinRange(0, 1000.0, 100), test_p.(LinRange(0, 1000.0, 100)))

