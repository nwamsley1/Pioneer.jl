x = sort(prosit_detailed[4211121], by = x-> getIntensity(x), rev= true)
cumsum([getIntensity(f) for f in x])./sum([getIntensity(f) for f in x])
plot(cumsum([getIntensity(f) for f in x])./sum([getIntensity(f) for f in x]))

function getTopFragments(precursors::Vector{Vector{LibraryFragment{Float64}}}, id_to_prec_mz::Dict{UInt32, Float64})
    frag_list = Vector{FragmentIon{Float64}}()
    for precursor in precursors
        sort!(precursor, by = x->getIntensity(x), rev = true)
        total = sum([getIntensity(frag) for frag in precursor])
        cummulative_sum = 0
        if !haskey(id_to_prec_mz, getPrecID(precursor[1]))
            continue
        end
        for frag in precursor
            push!(frag_list, FragmentIon(getFragMZ(frag), getPrecID(frag), id_to_prec_mz[getPrecID(frag)], getPrecCharge(frag)))
            cummulative_sum += getIntensity(frag)
            if cummulative_sum > 0.9*total
                break
            end
        end
    end
    return frag_list
end

function getTopFragments(precursors::Vector{Vector{LibraryFragment{Float64}}}, id_to_prec_mz::Dict{UInt32, Float64})
    frag_list = Vector{FragmentIon{Float64}}()
    for precursor in precursors
        #sort!(precursor, by = x->getIntensity(x), rev = true)
        _norm = norm([getIntensity(frag) for frag in precursor])
        if !haskey(id_to_prec_mz, getPrecID(precursor[1]))
            continue
        end
        for frag in precursor
            push!(frag_list, FragmentIon(getFragMZ(frag), getPrecID(frag), id_to_prec_mz[getPrecID(frag)], getIntensity(frag)/_norm, getPrecCharge(frag)))
            #cummulative_sum += getIntensity(frag)
            #if cummulative_sum > 0.9*total
            #    break
            #end
        end
    end
    return sort!(frag_list, by =x->getFragMZ(x))
end



sort!(frag_list_best, by =x->getFragMZ(x))

sort!(prosit_list_simple, by = x->getPrecID(x))

#=id_to_prec_mz = Dict{UInt32, Float64}()
for frag in prosit_list_simple
    id_to_prec_mz[getPrecID(frag)] = getPrecMZ(frag)
end=#

id_to_prec_mz = Dict{UInt32, Float64}()
for prec_bin in prosit_index_all.precursor_bins
    for prec in prec_bin.precs
        id_to_prec_mz[getPrecID(prec)] = getPrecMZ(prec)
    end
end

@time X, H, IDtoROW, time_psms_20 = SearchRAW(MS_TABLE, best_frag_index, prosit_detailed, UInt32(1), min_frag_count = 4, topN = 50, fragment_tolerance = 10.0, scan = 100028)
@time X, H, IDtoROW, time_psms_200 = SearchRAW(MS_TABLE, best_frag_index, prosit_detailed, UInt32(1), min_frag_count = 4, topN = 200, fragment_tolerance = 10.0, scan = 100028)

RTs, iRTs = refinePSMs!(time_psms_20, prosit_precs)
RTs, iRTs = refinePSMs!(time_psms_200, prosit_precs)

X = Matrix(time_psms_20[:,[:hyperscore,:total_ions,:y_ladder,:b_ladder,:intensity_explained,:error,:poisson,:spectral_contrast,:spectrum_peaks,:nmf,:weight,:RT_error]])
X_labels = time_psms_20[:, :decoy]
model = build_forest(X_labels, X, 4, 1000, 0.5, 3)
X = Matrix(time_psms_20[:,[:hyperscore,:total_ions,:y_ladder,:b_ladder,:intensity_explained,:error,:poisson,:spectral_contrast,:spectrum_peaks,:nmf,:weight,:RT_error]])
probs = apply_forest_proba(model, X,[true, false])
time_psms_20[:,:prob] = probs[:,2]

X = Matrix(time_psms_200[:,[:hyperscore,:total_ions,:y_ladder,:b_ladder,:intensity_explained,:error,:poisson,:spectral_contrast,:spectrum_peaks,:nmf,:weight,:RT_error]])
X_labels = time_psms_200[:, :decoy]
model = build_forest(X_labels, X, 4, 1000, 0.5, 3)
X = Matrix(time_psms_200[:,[:hyperscore,:total_ions,:y_ladder,:b_ladder,:intensity_explained,:error,:poisson,:spectral_contrast,:spectrum_peaks,:nmf,:weight,:RT_error]])
probs = apply_forest_proba(model, X,[true, false])
time_psms_200[:,:prob] = probs[:,2]

@time getQvalues!(time_psms_200, time_psms_200[:,:prob], time_psms_200[:,:decoy]);
@time getQvalues!(time_psms_20, time_psms_20[:,:prob], time_psms_20[:,:decoy]);

histogram(time_psms_20[time_psms_20[:,:decoy].==false,:q_values], alpha = 0.5)#, bins = -0.06:0.01:0.0)
histogram!(time_psms_20[time_psms_20[:,:decoy].==true,:q_values], alpha = 0.5)#, b


sum(time_psms_20[time_psms_20[:,:decoy].==true,:q_values].<=0.01)
sum(time_psms_20[time_psms_20[:,:decoy].==false,:q_values].<=0.01)


sum(time_psms_200[time_psms_200[:,:decoy].==true,:q_values].<=0.01)
sum(time_psms_200[time_psms_200[:,:decoy].==false,:q_values].<=0.01)

function testSearch(frag_index, topN::Int, fragment_tolerance::Float64, min_frag_count::Int, lambda::Float64, max_peaks::Int)
    @time PSMs = SearchRAW(MS_TABLE, frag_index, prosit_detailed, UInt32(1), min_frag_count = min_frag_count, topN = topN, fragment_tolerance = fragment_tolerance, lambda = lambda, max_peaks = max_peaks, scan = 100028)
    RTs, iRTs = refinePSMs!(PSMs, prosit_precs)
    rankPSMs!(PSMs, 2);
    @time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
    is_decoy = PSMs[:,:decoy].==true
    is_target = PSMs[:,:decoy].==false
    fdr =  PSMs[:,:q_values].<=0.01
    println("Decoy hits ", length(unique(PSMs[(fdr .& is_decoy),:precursor_idx])))
    println("Target hits ", length(unique(PSMs[(fdr .& is_target),:precursor_idx])))
    return PSMs#unique(PSMs[(fdr .& is_target),:precursor_idx])#PSMs#sum(time_psms_20[time_psms_20[:,:decoy].==true,:q_values].<=0.01), sum(time_psms_20[time_psms_20[:,:decoy].==false,:q_values].<=0.01)
end

@time PSMs_a = testSearch(prosit_index_intensities, 20, 40.0, 5, 1e3, 200); #5529
@time count = testSearch(prosit_index_intensities, 20, 40.0, 5, 1e4, 200); #5529



#dotp*count
function sort!(counter::Counter{I,C,T}, size::Int, topN::Int) where {I,C<:Unsigned,T<:AbstractFloat} 
    return sort!(
                @view(counter.ids[1:size]), 
                by = id -> (getDP(counter, id))*getCount(counter.dotp[id]),
                
                #by = id -> getCount(counter.dotp[id]),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time dotp_count = testSearch(prosit_index_intensities, 20, 40.0, 5, 1e3, 200); #5544
#dotp^2*count
function sort!(counter::Counter{I,C,T}, size::Int, topN::Int) where {I,C<:Unsigned,T<:AbstractFloat} 
    return sort!(
                @view(counter.ids[1:size]), 
                by = id -> (getDP(counter, id)^2)*getCount(counter.dotp[id]),
                
                #by = id -> getCount(counter.dotp[id]),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time dotp2_count = testSearch(prosit_index_intensities, 20, 40.0, 5, 1e3, 200); #5478
#count
function sort!(counter::Counter{I,C,T}, size::Int, topN::Int) where {I,C<:Unsigned,T<:AbstractFloat} 
    return sort!(
                @view(counter.ids[1:size]), 
                by = id -> getCount(counter.dotp[id]),
                
                #by = id -> getCount(counter.dotp[id]),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time count = testSearch(prosit_index_intensities, 20, 40.0, 5, 1e3, 200); #5529
#intensity
function sort!(counter::Counter{I,C,T}, size::Int, topN::Int) where {I,C<:Unsigned,T<:AbstractFloat} 
    return sort!(
                @view(counter.ids[1:size]), 
                by = id -> getObs(counter.dotp[id]),
                
                #by = id -> getCount(counter.dotp[id]),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time intensity = testSearch(prosit_index_intensities, 20, 40.0, 5, 1e3, 200); #4861
#intensity*dotp
function sort!(counter::Counter{I,C,T}, size::Int, topN::Int) where {I,C<:Unsigned,T<:AbstractFloat} 
    return sort!(
                @view(counter.ids[1:size]), 
                by = id -> (getDP(counter, id))*getObs(counter.dotp[id]),
                
                #by = id -> getCount(counter.dotp[id]),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time intensity_dotp = testSearch(prosit_index_intensities, 20, 40.0, 5, 1e3, 200); #5092
#intensity*dotp^2
function sort!(counter::Counter{I,C,T}, size::Int, topN::Int) where {I,C<:Unsigned,T<:AbstractFloat} 
    return sort!(
                @view(counter.ids[1:size]), 
                by = id -> (getDP(counter, id)^2)*getCount(counter.dotp[id]),
                
                #by = id -> getCount(counter.dotp[id]),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time  intensity_dotp2 testSearch(prosit_index_intensities, 20, 40.0, 5, 1e3, 200); # 5523


#sort by Intensity 
@time testSearch(10, 10.0, 4) #(161, 16033)
@time testSearch(20, 10.0, 4) #(158, 15715)
@time testSearch(50, 10.0, 4) #(157, 15551)
@time testSearch(100, 10.0, 4) #(155, 15405)
@time testSearch(200, 10.0, 4) #(151, 15051)

@time testSearch(10, 20.0, 4) #(315, 31228)
@time testSearch(10, 40.0, 4) #(394, 39023)
@time testSearch(10, 60.0, 4) #(352, 34927)

@time testSearch(50, 40.0, 4) #(387, 38491)

#fixed matched peaks

#Intensity
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getCount(counter, id)*getIntensity(counter, id),
                by = id -> getIntensity(counter, id),
                #by = id -> getCount(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time testSearch(best_frag_index, 10, 40.0, 4, 1e0) #(389, 38545)
@time testSearch(best_frag_index, 10, 40.0, 4, 1e1) #(386, 38288)
@time testSearch(best_frag_index, 10, 40.0, 4, 1e2) #(389, 38543)
@time testSearch(best_frag_index, 10, 40.0, 4, 1e3) #(389, 38574)



@time testSearch(best_frag_index, 10, 40.0, 4, 1e4) #(386, 38250)
@time testSearch(prosit_index_all, 10, 40.0, 4, 1e3) #(489, 48481)

@time testSearch(best_frag_index, 50, 40.0, 4, 1e3) #(409, 40573)
@time testSearch(prosit_index_all, 50, 40.0, 4, 1e3) #(555, 54987)

@time testSearch(best_frag_index, 200, 40.0, 4, 1e3) #(393, 38917)
#Count
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getCount(counter, id)*getIntensity(counter, id),
                #by = id -> getIntensity(counter, id),
                by = id -> getCount(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time testSearch(best_frag_index, 10, 40.0, 4,1e3) #(390, 38705)
@time testSearch(best_frag_index, 50, 40.0, 4,1e3) #(408, 40476)
@time testSearch(best_frag_index, 200, 40.0, 4,1e3) #(399, 39526)
#sort by count*intensity
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                by = id -> getCount(counter, id)*getIntensity(counter, id),
                #by = id -> getIntensity(counter, id),
                #by = id -> getCount(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time testSearch(best_frag_index, 10, 40.0, 4,1e3) #(363, 36004)
@time testSearch(best_frag_index, 50, 40.0, 4,1e3) #(384, 38084)
@time testSearch(best_frag_index, 200, 40.0, 4,1e3) #(375, 37182)

@time testSearch(prosit_index_all, 10, 40.0, 4,1e3) #(341, 33761)
@time testSearch(prosit_index_all, 50, 40.0, 4,1e3) #(421, 41799)
@time testSearch(prosit_index_all, 200, 40.0, 4,1e3) #(449, 44517)

#Count
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getCount(counter, id)*getIntensity(counter, id),
                #by = id -> getIntensity(counter, id),
                by = id -> getCount(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time testSearch(prosit_index_all, 10, 40.0, 4,1e3) #(491, 48632)
@time testSearch(prosit_index_all, 50, 40.0, 4,1e3) #(557, 55229)

#@time testSearch(prosit_index_all, 200, 40.0, 4,1e3) #


#Intensity/Count
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                by = id -> getIntensity(counter, id)/getCount(counter, id),
                #by = id -> getIntensity(counter, id),
                #by = id -> getCount(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time testSearch(prosit_index_all, 10, 40.0, 4,1e3) #(31, 3155)
#@time testSearch(prosit_index_all, 200, 40.0, 4,1e3) #(363, 36004)

#Count
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getCount(counter, id)*getIntensity(counter, id),
                #by = id -> getIntensity(counter, id),
                by = id -> getCount(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end

@time testSearch(prosit_index_all, 10, 40.0, 4, 1e3, 200) #(502, 49711)
@time testSearch(prosit_index_all, 10, 40.0, 4, 1e3, 100) #(391, 38790)

@time testSearch(prosit_index_all, 50, 40.0, 4, 1e3, 200) #(545, 53978)
@time testSearch(prosit_index_all, 50, 40.0, 4, 1e3, 100) #(390, 38628)

#Intensity
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getCount(counter, id)*getIntensity(counter, id),
                by = id -> getIntensity(counter, id),
                #by = id -> getCount(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end

@time testSearch(prosit_index_all, 10, 40.0, 4, 1e3, 200) #(500, 49564)
@time testSearch(prosit_index_all, 20, 40.0, 4, 1e3, 200) #(529, 52431)
@time testSearch(prosit_index_all, 50, 40.0, 4, 1e3, 200) #(537, 53234)

function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getCount(counter, id)*getIntensity(counter, id),
                #by = id -> getIntensity(counter, id),
                by = id -> getCount(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end

@time PSMs = testSearch(prosit_index_all, 20, 40.0, 4, 1e3, 200) #(502, 49711)
@time testSearch(prosit_index_all, 10, 40.0, 4, 1e3, 100) #(391, 38790)


#@time testSearch(prosit_index_all, 10, 40.0, 4, 1e3, 100) #(491, 48632)
#end
@time testSearch(prosit_index_intensities, 10, 40.0, 4, 1e3, 200) #(500, 49564)
#Count
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getCount(counter, id)*getIntensity(counter, id),
                #by = id -> getIntensity(counter, id),
                by = id -> getCount(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time testSearch(prosit_index_intensities, 10, 40.0, 4, 1e3, 200) #(500, 49564)

function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getCount(counter, id)*getIntensity(counter, id),
                by = id -> getIntensity(counter, id)/sqrt(getNorm(counter, id)),
                #by = id -> getCount(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end

#dot product 
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getCount(counter, id)*getIntensity(counter, id),
                by = id -> getIntensity(counter, id)/sqrt(getNorm(counter, id)),
                #by = id -> getCount(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time PSMs = testSearch(prosit_index_intensities, 20, 40.0, 4, 1e3, 200) #190, 18844
#count
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getCount(counter, id)*getIntensity(counter, id),
                #by = id -> getIntensity(counter, id)/sqrt(getNorm(counter, id)),
                by = id -> getCount(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time PSMs = testSearch(prosit_index_intensities, 20, 40.0, 4, 1e3, 200) #545, 53957
#count*(dot product)
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getCount(counter, id)*getIntensity(counter, id),
                by = id -> getCount(counter, id)*getIntensity(counter, id)/sqrt(getNorm(counter, id)),
                #by = id -> getCount(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time PSMs = testSearch(prosit_index_intensities, 20, 40.0, 4, 1e3, 200) #505, 50072
#count*(dot^2)
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getIntensity(counter, id)/getCount(counter, id),
                by = id -> getCount(counter, id)*((getIntensity(counter, id)/sqrt(getNorm(counter, id)))^2),
                #by = id -> getCount(counter, id),
                #by = id->getIntensity(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time PSMs = testSearch(prosit_index_intensities, 20, 40.0, 4, 1e3, 200) #310, 30723
#count^2*(dot)
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getIntensity(counter, id)/getCount(counter, id),
                by = id -> (getCount(counter, id)^2)*((getIntensity(counter, id)/sqrt(getNorm(counter, id)))),
                #by = id -> getCount(counter, id),
                #by = id->getIntensity(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time PSMs = testSearch(prosit_index_intensities, 20, 40.0, 4, 1e3, 200) #53619
#intensity
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                #by = id -> getCount(counter, id)*getIntensity(counter, id),
                #by = id -> getCount(counter, id)*getIntensity(counter, id)/sqrt(getNorm(counter, id)),
                #by = id -> getCount(counter, id),
                by = id->getIntensity(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time PSMs = testSearch(prosit_index_intensities, 20, 40.0, 4, 1e3, 200) # 13711
#intensity*count
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                by = id -> getCount(counter, id)*getIntensity(counter, id),
                #by = id -> getCount(counter, id)*getIntensity(counter, id)/sqrt(getNorm(counter, id)),
                #by = id -> getCount(counter, id),
                #by = id->getIntensity(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time PSMs = testSearch(prosit_index_intensities, 20, 40.0, 4, 1e3, 200) #17746
#intensity^2*count
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                by = id -> getCount(counter, id)*(getIntensity(counter, id)^2),
                #by = id -> getCount(counter, id)*getIntensity(counter, id)/sqrt(getNorm(counter, id)),
                #by = id -> getCount(counter, id),
                #by = id->getIntensity(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end
@time PSMs = testSearch(prosit_index_intensities, 20, 40.0, 4, 1e3, 200) #15655
#intensity/count
function sort!(counter::Counter{I,C}, size::Int, topN::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                by = id -> getIntensity(counter, id)/getCount(counter, id),
                #by = id -> getCount(counter, id)*getIntensity(counter, id)/sqrt(getNorm(counter, id)),
                #by = id -> getCount(counter, id),
                #by = id->getIntensity(counter, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )#[1:min(num_precs, end)]
end

######
@time PSMs = testSearch(prosit_index_intensities, 20, 40.0, 4, 1e3, 200) #10875





@time PSMs testSearch(prosit_index_intensities, 10, 40.0, 4, 1e3, 200) #


time_scribe = SearchRAW(MS_TABLE, best_frag_index, prosit_detailed, UInt32(1), min_frag_count = 4, topN = 10, fragment_tolerance = 40.0, scan = 100028)


RTs, iRTs = refinePSMs!(time_psms_20, prosit_precs)

train = randperm(size(time_psms_20)[1])#[1:100000]
X = Matrix(time_psms_20[train,[:hyperscore,:total_ions,:y_ladder,:b_ladder,:intensity_explained,:error,:poisson,:spectral_contrast,:spectrum_peaks,:nmf,:weight,:RT_error,:scribe_score]])
X_labels = time_psms_20[train, :decoy]
model = build_forest(X_labels, X, 4, 1000, 0.05, 3)
X = Matrix(time_psms_20[:,[:hyperscore,:total_ions,:y_ladder,:b_ladder,:intensity_explained,:error,:poisson,:spectral_contrast,:spectrum_peaks,:nmf,:weight,:RT_error,:scribe_score]])
probs = apply_forest_proba(model, X,[true, false])
time_psms_20[:,:prob] = probs[:,2]
@time getQvalues!(time_psms_20, time_psms_20[:,:prob], time_psms_20[:,:decoy]);
time_scribe = time_psms_20

histogram(log.(time_scribe[time_scribe[:,:q_values].<=0.01,:][:,:error]))
histogram!(log.(time_scribe[time_scribe[:,:q_values].>0.01,:][:,:error]))