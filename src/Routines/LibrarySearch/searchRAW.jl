function SearchRAW(
                    spectra::Arrow.Table, 
                    #ptable::PrecursorDatabase,
                    frag_index::FragmentIndex{T},
                    fragment_list::Vector{Vector{LibraryFragment{Float64}}},
                    ms_file_idx::UInt32;
                    isolation_width::Float64 = 4.25,
                    precursor_tolerance::Float64 = 5.0,
                    fragment_tolerance::Float64 = 20.0,
                    topN::Int64 = 20,
                    min_frag_count::Int64 = 4,
                    lambda::Float64 = 1e3,
                    scan_range::Tuple{Int64, Int64} = (0, 0), 
                    max_peaks::Int = 200, 
                    #fragment_match_ppm::U,
                    data_type::Type{T} = Float64
                    ) where {T<:Real}
    
    scored_PSMs = makePSMsDict(XTandem(data_type))
    ms2, MS1, MS1_i = 0, 0, 0
    precs = Counter(UInt32, UInt8, Float32, 9387261) #Prec counter
    times = Dict(:counter => 0.0, :reset => 0.0, :nmf => 0.0, :metrics => 0.0, :match_peaks => 0.0, :build => 0.0, :score => 0.0)
    n = 0
    for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(spectra)))
    #for (i, spectrum) in enumerate(Tables.namedtupleiterator(spectra))

        if spectrum[:msOrder] == 1
            MS1 = spectrum[:masses]
            MS1_i = spectrum[:intensities]
            continue
        else
            ms2 += 1
        end

        if scan_range != (0, 0)
            i < first(scan_range) ? continue : nothing
            i > last(scan_range) ? continue : nothing
        end


        min_intensity = spectrum[:intensities][sortperm(spectrum[:intensities], rev = true)[min(max_peaks, length(spectrum[:intensities]))]]

        #times[:counter] += @elapsed prec_count, match_count = searchScan!(precs,
        prec_count, match_count = searchScan!(precs,
                    frag_index, 
                    min_intensity, spectrum[:masses], spectrum[:intensities], MS1, spectrum[:precursorMZ], 
                    fragment_tolerance, 
                    precursor_tolerance,
                    isolation_width,
                    min_frag_count = min_frag_count, 
                    topN = topN
                    )
        #return precs
        if getSize(precs) <= 1
            #times[:reset] += @elapsed reset!(precs)
            reset!(precs)
            continue
        end
        transitions = selectTransitions(fragment_list, precs, topN)
        reset!(precs)
        if length(transitions) == 0
            continue
        end
        
        #times[:reset] += @elapsed reset!(precs)
       

        #times[:match_peaks] += @elapsed fragmentMatches, fragmentMisses = matchPeaks(transitions, 
        fragmentMatches, fragmentMisses = matchPeaks(transitions, 
                                    spectrum[:masses], 
                                    spectrum[:intensities], 
                                    count_unmatched=true,
                                    δs = zeros(T, (1,)),
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = min_intensity,
                                    ppm = fragment_tolerance
                                    )
        if iszero(length(fragmentMatches))
            continue
        end

        #times[:build] += @elapsed X, H, UNMATCHED, IDtoROW = buildDesignMatrix(fragmentMatches, fragmentMisses, topN)
        X, H, UNMATCHED, IDtoROW = buildDesignMatrix(fragmentMatches, fragmentMisses, topN)
       
        #return X, H, UNMATCHED, IDtoROW, fragmentMatches, fragmentMisses
        #Initialize weights for each precursor template. 
        #Should find a more sophisticated way of doing this. 
        W = reshape([Float32(1000) for x in range(1,H.m)], (1, H.m))
        #W = reshape([Float32(1000) for x in range(1,size(H)[1])], (1, size(H)[1]))
        
        weights = W[1,:]
        #=times[:nmf] += @elapsed weights = (NMF.solve!(NMF.MultUpdate{Float32}(maxiter=50, verbose = false, 
                                                    lambda_w = lambda, 
                                                    tol = 100, #Need a reasonable way to choos lambda?
                                                    update_H = false #Important to keep H constant. 
                                                    ), X, W, H).W[1,:])=#

        #times[:metrics] += @elapsed scribe_score, city_block, matched_ratio, spectral_contrast_matched, spectral_contrast_all = getDistanceMetrics(H, X, UNMATCHED)
        scribe_score, city_block, matched_ratio, spectral_contrast_matched, spectral_contrast_all = getDistanceMetrics(H, X, UNMATCHED)
        
        #For progress and debugging. 

        unscored_PSMs = UnorderedDictionary{UInt32, XTandem{T}}()

        ScoreFragmentMatches!(unscored_PSMs, fragmentMatches)

        #times[:score] += @elapsed Score!(scored_PSMs, unscored_PSMs, 
        Score!(scored_PSMs, unscored_PSMs, 
                length(spectrum[:intensities]), 
                Float64(sum(spectrum[:intensities])), 
                match_count/prec_count, scribe_score, city_block, matched_ratio, spectral_contrast_matched, spectral_contrast_all, weights, IDtoROW,
                scan_idx = Int64(i)
                )
        n += 1
    end

    println("processed $ms2 scans!")
    println("counter: ", times[:counter]/n)
    println("reset: ", times[:reset]/n)
    println("nmf : ", times[:nmf]/n)
    println("metrics : ", times[:metrics]/n)
    println("build : ", times[:build]/n)
    println("score : ", times[:score]/n)
    println("match_peaks : ", times[:match_peaks]/n)
    return DataFrame(scored_PSMs)
end

#=struct runningCount{C<:Unsigned,T<:AbstractFloat}
    count::Base.RefValue{C}
    sum::Base.RefValue{T}
    function runningCount(count::C, sum::T) where {C<:Unsigned, T<:AbstractFloat}
        new{C, T}(Ref(count), Ref(sum))
    end
end=#

#=runningCount(C::DataType, T::DataType) = runningCount(zero(C),zero(T)) 
getCount(rc::runningCount{C,T}) where {C<:Unsigned,T<:AbstractFloat} = rc.count[]
getSum(rc::runningCount{C,T}) where {C<:Unsigned,T<:AbstractFloat} = rc.sum[]=#


mutable struct Counter{I,C<:Unsigned,T<:AbstractFloat}
    ids::Vector{I}
    counts::Vector{C}
    sums::Vector{T}
    size::Int64
    matches::Int64
    function Counter(I::DataType, C::DataType, T::DataType, size::Int) #where {I,C<:Unsigned}
        new{I, C, T}(zeros(I, size), zeros(C,size), zeros(T,size), 1, 0)
    end
end

getCount(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = c.counts[id]
getSums(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = c.sums[id]
function getMatchedRatio(c::Counter{I,C,T}, totals::Vector{Float32}, id::I) where {C,I<:Unsigned,T<:AbstractFloat}
    matched_intensity = getSums(c, id)
    matched_intensity/(totals[id] - matched_intensity)
end

getSize(c::Counter{I,C,T}) where {I,C<:Unsigned,T<:AbstractFloat} = c.size
incSize!(c::Counter{I,C,T}) where {I,C<:Unsigned,T<:AbstractFloat} = c.size += 1
getID(c::Counter{I,C,T}, idx::Int) where {I,C<:Unsigned,T<:AbstractFloat} = c.ids[idx]

function update!(c::Counter{I,C,T}, id::I, pred_intensity::T) where {C,I<:Unsigned,T<:AbstractFloat}
    c.counts[id] += one(C);
    c.sums[id] += pred_intensity;
    return nothing
end

function reset!(c::Counter{I,C,T}, id::I) where {C,I<:Unsigned,T<:AbstractFloat}
    
    c.counts[id] = zero(C);
    c.sums[id] = zero(T);

    return nothing
end
#getObs(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getObs(c.dotp[id])
#getPred(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getPred(c.dotp[id])
#getObsPred(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getObsPred(c.dotp[id])
#getCount(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getCount(c.dotp[id])
#getDP(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getDP(c.[id])

import DataStructures.inc!
function inc!(c::Counter{I,C,T}, id::I, pred_intensity::T) where {I,C<:Unsigned,T<:AbstractFloat} 
    if iszero(c.counts[id])#c.counts[id]<1#iszero(c.counts[id])# == zero(C)
        c.ids[getSize(c)] = id;
        update!(c,id,pred_intensity);
        incSize!(c);
    else
        update!(c,id,pred_intensity);
    end
    return nothing
end

import Base.sort!
function sort!(counter::Counter{I,C,T}, topN::Int) where {I,C<:Unsigned,T<:AbstractFloat}
    sort!(
                @view(counter.ids[1:counter.matches]), 
                by = id -> getMatchedRatio(counter, prosit_totals, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )
    return nothing
end

function reset!(c::Counter{I,C,T}) where {I,C<:Unsigned,T<:AbstractFloat} 
    @turbo for i in 1:(getSize(c) - 1)
        id = c.ids[i]
        c.counts[id] = zero(C)
        c.sums[id] = zero(T)
        #reset!(c, c.ids[i])
    end
    c.size = 1
    c.matches = 0
    return nothing
end

function countFragMatches(c::Counter{I,C,T}, min_count::Int) where {I,C<:Unsigned,T<:AbstractFloat} 
    frag_counts = 0
    for i in 1:(getSize(c) - 1)
        id = c.ids[i]
        frag_count = getCount(c, id)
        frag_counts += frag_count
        if frag_count >= min_count
            if getMatchedRatio(c, prosit_totals, id)>=0.8
                    c.ids[c.matches + 1] = c.ids[i]
                    c.matches += 1
            end
        end
    end
    return frag_counts
end

#for prec in prosit_detailed

#=
big_list = [rand([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]) for x in 1:9000000]
@btime sort(big_list, alg=PartialQuickSort(1:6), rev = true)


@btime sort(big_list[1:10000], alg=PartialQuickSort(1:6), rev = true)


@btime sort(filter(x->x>13,@view(big_list[1:100000])), alg=PartialQuickSort(1:20), rev = true)[1:20]
@btime sort(filter(x->x>0,@view(big_list[1:100000])), alg=PartialQuickSort(1:20), rev = true)[1:20]
@btime sort(filter(x->x>13,big_list[1:100000]), alg=PartialQuickSort(1:20), rev = true)[1:20]
@benchmark @view(sort(filter(x->x>13,@view(big_list[1:100000])), alg=PartialQuickSort(1:20), rev = true)[1:20])#fastest
@benchmark @view(sort(filter(x->x>13,@view(big_list[1:100000])), alg=QuickSort, rev = true)[1:20])

@benchmark sort(filter(x->x>13,@view(big_list[1:100000])), alg=PartialQuickSort(1:20), rev = true)#fastest

@benchmark @view(sort(filter(x->x>13,@view(big_list[1:100000])), rev = true)[1:20])#fastest

sort([5, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 4, 1, 4, 4, 5, 5], alg=PartialQuickSort(1:6), rev = true)

x = rand(100);
k = 50:100;

@time sort!(x[1:10]; alg=PartialQuickSort(1:10), rev = true);
=#