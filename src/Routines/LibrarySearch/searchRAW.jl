function SearchRAW(
                    spectra::Arrow.Table, 
                    #ptable::PrecursorDatabase,
                    frag_index::FragmentIndex{T},
                    fragment_list::Vector{Vector{LibraryFragment{Float64}}},
                    ms_file_idx::UInt32;
                    precursor_tolerance::Float64 = 4.25,
                    fragment_tolerance::Float64 = 20.0,
                    transition_charges::Vector{UInt8} = UInt8[1],
                    transition_isotopes::Vector{UInt8} = UInt8[0],
                    b_start::Int64 = 3,
                    y_start::Int64 = 3,
                    topN::Int64 = 20,
                    min_frag_count::Int64 = 4,
                    #fragment_match_ppm::U,
                    data_type::Type{T} = Float64
                    ) where {T,U<:Real}

    scored_PSMs = makePSMsDict(XTandem(data_type))
    #scored_PSMs = makePSMsDict(FastXTandem(data_type))
    #precursorList needs to be sorted by precursor MZ. 
    #Iterate through rows (spectra) of the .raw data. 
    #i = 0
    ms2 = 0
    min_intensity = Float32(0.0)

    #precs = Dict{UInt32, UInt8}()
    #precs = Accumulator{UInt32, UInt8}()
    precs = Counter(UInt32, UInt8, 9387261)
    match_times = []
    for (i, spectrum) in enumerate(Tables.namedtupleiterator(spectra))
        if spectrum[:msOrder] != 2
            continue
        end
        ms2 += 1
        if ms2 < 100000#50000
            #println("TEST")
            continue
        elseif ms2 > 102000#51000#51000
            continue
        end
        
        fragmentMatches = Vector{FragmentMatch{Float32}}()
        #println(" a scan")
        #precs = Dict{UInt32, UInt8}(zeros(UInt8, pre_aloc_size), zeros(UInt32, pre_aloc_size), zeros(UInt8, pre_aloc_size), 0, 0, 0, pre_aloc_size, 0)
        #precs = Dict{UInt32, UInt8}()
        
        #frag_time = @elapsed 
        frag_time = @elapsed prec_count, match_count = searchScan!(precs,
                    frag_index, 
                    spectrum[:masses], spectrum[:intensities], spectrum[:precursorMZ], 
                    fragment_tolerance, 
                    precursor_tolerance,
                    min_frag_count = min_frag_count, 
                    topN = topN
                    )
        
        #println(length(prec_counts))

        transitions = selectTransitions(fragment_list, pep_id_iterator, topN)
        if precs.size > 1
            frag_time += @elapsed reset!(precs)
            push!(match_times, frag_time)
        else
            push!(match_times, frag_time)
            continue
        end
        #frag_time += @elapsed empty!(precs)
        #push!(match_times, frag_time)
        #frag_time += @elapsed reset!(precs)
        #push!(match_times, frag_time)
        #push!(fragger_times, fragger_time)
        fragmentMatches, fragmentMisses = matchPeaks(transitions, 
                                    spectrum[:masses], 
                                    spectrum[:intensities], 
                                    #δs = params[:δs],
                                    δs = zeros(T, (1,)),#[Float64(0)],
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = min_intensity
                                    )
                                    #push!(fragger_times, fragger_time)
        #println("matches ", length(fragmentMatches))
        #println("misses ", length(fragmentMisses))
        
        X, H, IDtoROW = buildDesignMatrix(fragmentMatches, fragmentMisses, topN)
        #Does this always coincide with there being zero fragmentMatches?
        #Could change to length(fragmentMatches) == 0 ?
        if size(H)[2] == 0
            continue
        end
        #Initialize weights for each precursor template. 
        #Should find a more sophisticated way of doing this. 
        W = reshape([Float32(100) for x in range(1,size(H)[1])], (1, size(H)[1]))

        #Solve NMF. 
        #=nmf_time += @elapsed weights = NMF.solve!(NMF.GreedyCD{Float32}(maxiter=50, verbose = false, 
                                                    lambda_w = 1e3, 
                                                    tol = 1e-6, #Need a reasonable way to choos lambda?
                                                    update_H = false #Important to keep H constant. 
                                                    ), X, W, H).W[1,:]=#

        weights = (NMF.solve!(NMF.ProjectedALS{Float32}(maxiter=50, verbose = false, 
                                                    lambda_w = 1e4, 
                                                    tol = 1e-6, #Need a reasonable way to choos lambda?
                                                    update_H = false #Important to keep H constant. 
                                                    ), X, W, H).W[1,:])
        #weights = W[1,:]                           
        #println(size(H))
        #println(size(X))
        #tH = H'
        
        #nmf_time = @elapsed weights = coef(fit(LassoModel, tH, X[1,:], λ=[1e2], cd_tol=1, criterion=:obj))
        spectral_contrast = getSpectralContrast(H, X)

        #For progress and debugging. 
        if (ms2 % 1000) == 0
            println("ms2: $ms2")
            if ms2 == 8000
                #test_frags = transitions
                #test_matches = fragmentMatches
                #test_misses = fragmentMisses
            end
        end

        unscored_PSMs = UnorderedDictionary{UInt32, XTandem{T}}()

        ScoreFragmentMatches!(unscored_PSMs, fragmentMatches)

        Score!(scored_PSMs, unscored_PSMs, 
                length(spectrum[:intensities]), 
                Float64(sum(spectrum[:intensities])), 
                match_count/prec_count, spectral_contrast, weights, IDtoROW,
                scan_idx = Int64(i)
                )
        #push!(score_times, score)
    
    end

    println("processed $ms2 scans!")
    println("mean matches: ", mean(match_times))
    println("mean matches: ", median(match_times))
    #=println("processed $ms2 scans!")
    println("mean build: ", mean(build_design_times))
    println("mean fragger: ", mean(fragger_times))
    println("mean matches: ", mean(match_times))
    println("mean nmf: ", mean(nmf_times))
    println("median nmf: ", median(nmf_times))
    println("mean s_contrast: ", mean(spectral_contrast_times))
    println("mean score: ", mean(score_times))=#
    return DataFrame(scored_PSMs)#precs#DataFrame(scored_PSMs)# test_frags, test_matches, test_misses#DataFrame(scored_PSMs)
end

#=unction find_nth_largest(array, n)
    # Perform insertion sort in descending order
    for i in 2:length(array)
        j = i
        while j > 1 && array[j] > array[j-1]
            array[j], array[j-1] = array[j-1], array[j]
            j -= 1
        end
    end

    return array[n]
end=#

#prec_ids = UInt32[1, 10, 5, 7, 0, 0, 0, 0, 0, 0]
#prec_counts = UInt8[1, 0, 0, 0, 1, 0, 2, 0, 0, 11]
#Every time we encoutner a precursor, check value to see if it has been assigned. 
#If not, then add one to the indexer and place the precursor key in keys. 
#function sortby(prec_ids::Vector{UInt32}, prec_counts::Vector{UInt8}, max_n::Int)
#    return sort(filter(x->prec_counts[x].>1, @view(prec_ids[1:max_n])), by = x->prec_counts[x], rev = true)
#end

mutable struct Counter{I,C<:Unsigned}
    ids::Vector{I}
    counts::Vector{C}
    size::Int64
    function Counter(I::DataType, C::DataType, size::Int) #where {I,C<:Unsigned}
        new{I, C}(Vector{I}(undef, size), Vector{C}(undef, size), 1)
    end
end


getSize(c::Counter{I,C}) where {I,C<:Unsigned} = c.size
getID(c::Counter{I,C}, idx::Int) where {I,C<:Unsigned} = c.ids[idx]
getCount(c::Counter{I,C}, id::I) where {I,C<:Unsigned} = c.counts[id]
setCount!(c::Counter{I,C}, idx::Int, count::C) where {I,C<:Unsigned} = c.counts[c.ids[idx]] = count


incSize!(c::Counter{I,C}) where {I,C<:Unsigned} = c.size += 1
incCounter!(c::Counter{I,C}, id::I) where {I,C<:Unsigned} = c.counts[id] += one(C)

import DataStructures.inc!
function inc!(c::Counter{I,C}, id::I) where {I,C<:Unsigned} 
    if c.counts[id] == 0
        c.ids[getSize(c)] = id;
        incCounter!(c, id);
        incSize!(c);
    else
        incCounter!(c, id);
    end
end

import Base.sort!
function sort!(counter::Counter{I,C}, size::Int) where {I,C<:Unsigned} 
    return sort!(
                @view(counter.ids[1:size]), 
                by = id -> getCount(counter, id),
                rev = true,
                #alg=PartialQuickSort(1:num_precs)
             )#[1:min(num_precs, end)]
end

function reset!(c::Counter{I,C}) where {I,C<:Unsigned}
    
    @turbo for i in 1:(getSize(c) - 1)
        c.counts[c.ids[i]] = zero(C)
        #setCount(c, idx, zero(C))
    end
    c.size = 1
end

function countFragMatches(c::Counter{I,C}) where {I,C<:Unsigned}
    count = zero(C)
    @turbo for i in 1:(getSize(c) - 1)
        count += c.counts[c.ids[i]]
    end
    return count
end

function countFragMatches(c::Counter{I,C}, min_count::Int) where {I,C<:Unsigned}
    frag_counts = zero(C)
    exceeds_min = 0
    for i in 1:(getSize(c) - 1)
        id = c.ids[i]
        frag_count = c.counts[id]
        frag_counts += frag_count
        if frag_count >= min_count
            c.ids[exceeds_min + 1] = c.ids[i]
            exceeds_min += 1
        end
    end
    return frag_counts, exceeds_min
end

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