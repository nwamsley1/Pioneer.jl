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
                    lambda::Float64 = 1e3,
                    scan::Int = 10000, 
                    max_peaks::Int = 200, 
                    #fragment_match_ppm::U,
                    data_type::Type{T} = Float64
                    ) where {T,U<:Real}
    println("TEST3")
    
    scored_PSMs = makePSMsDict(XTandem(data_type))
    #scored_PSMs = makePSMsDict(FastXTandem(data_type))
    #precursorList needs to be sorted by precursor MZ. 
    #Iterate through rows (spectra) of the .raw data. 
    #i = 0
    ms2 = 0
    min_intensity = Float32(0.0)
    rows = 0
    X, H, IDtoROW = "", "", ""
    precs = Counter(UInt32, UInt8, Float32, 9387261)
    times = Dict(:counter => 0.0, :reset => 0.0, :nmf => 0.0)
    MS1 = 0
    for (i, spectrum) in enumerate(Tables.namedtupleiterator(spectra))
        if spectrum[:msOrder] == 1
            MS1 = spectrum[:masses]
            continue
        end
        ms2 += 1
        #=if i !=  scan
            continue
        end=#
        if ms2 < 73000#50000
            #println("TEST")
            continue
        elseif ms2 > 88000#51000#51000
            continue
        end
        #if (ms2 % 1000) == 0
        if (ms2 % 1000) == 0
            println("$ms2 ", length(scored_PSMs[:entropy]))
            #rows = length(scored_PSMs[:entropy])
        end
        min_intensity = spectrum[:intensities][sortperm(spectrum[:intensities], rev = true)[min(max_peaks, length(spectrum[:intensities]))]]

        #println(min_intensity)

        times[:counter] += @elapsed prec_count, match_count = searchScan!(precs,
                    frag_index, 
                    min_intensity, spectrum[:masses], spectrum[:intensities], MS1, spectrum[:precursorMZ], 
                    fragment_tolerance, 
                    precursor_tolerance,
                    min_frag_count = min_frag_count, 
                    topN = topN
                    )

        if getSize(precs) <= 1
            times[:reset] += @elapsed reset!(precs)
            continue
        end

        transitions = selectTransitions(fragment_list, precs, topN)
        times[:reset] += @elapsed reset!(precs)

        fragmentMatches, fragmentMisses = matchPeaks(transitions, 
                                    spectrum[:masses], 
                                    spectrum[:intensities], 
                                    #δs = params[:δs],
                                    δs = zeros(T, (1,)),#[Float64(0)],
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = min_intensity,
                                    ppm = fragment_tolerance
                                    )
        if iszero(length(fragmentMatches))
            continue
        end

        times[:nmf] += @elapsed X, H, IDtoROW, matched = buildDesignMatrix(fragmentMatches, fragmentMisses, topN)
        #println(norm(X[:,1:matched]))
        #println("matched: $matched")
        ##println(size(X[1,1:matched]))
        #println(size(H[:,1:matched]))
        #Initialize weights for each precursor template. 
        #Should find a more sophisticated way of doing this. 
        W = reshape([Float32(1000) for x in range(1,size(H)[1])], (1, size(H)[1]))
        
        #times[:nmf] += @elapsed weights = (NMF.solve!(NMF.GreedyCD{Float32}(maxiter=50, verbose = false, 
        #                                            lambda_w = lambda, 
        #                                            tol = 1e-4, #Need a reasonable way to choos lambda?
        #                                            update_H = false #Important to keep H constant. 
        #                                            ), X, W, H).W[1,:])
        
        #times[:nmf] += @elapsed weights = (NMF.solve!(NMF.CoordinateDescent{Float32}(maxiter=50, verbose = false, 
        #                                            α = lambda, 
        #                                            regularization=:transformation,
        #                                            l₁ratio=1.0,
        #                                            tol = 1e-4, #Need a reasonable way to choos lambda?
        #                                            update_H = false #Important to keep H constant. 
        #                                            ), X, W, H).W[1,:])
        times[:nmf] += @elapsed weights = (NMF.solve!(NMF.MultUpdate{Float32}(maxiter=50, verbose = false, 
                                                    lambda_w = lambda, 
                                                    tol = 100, #Need a reasonable way to choos lambda?
                                                    update_H = false #Important to keep H constant. 
                                                    ), X[:,1:matched], W, H[:,1:matched]).W[1,:])

        spectral_contrast = getSpectralContrast(H, X)
        scribe_score, city_block, chebyshev, unmatched = getScribeScore(H, X)
        #For progress and debugging. 


        unscored_PSMs = UnorderedDictionary{UInt32, XTandem{T}}()

        ScoreFragmentMatches!(unscored_PSMs, fragmentMatches)

        Score!(scored_PSMs, unscored_PSMs, 
                length(spectrum[:intensities]), 
                Float64(sum(spectrum[:intensities])), 
                match_count/prec_count, spectral_contrast, scribe_score, city_block, chebyshev, unmatched, weights, IDtoROW,
                scan_idx = Int64(i)
                )
        #=scored_PSMs = sort(DataFrame(scored_PSMs),:spectral_contrast)
        transform!( scored_PSMs, AsTable(:) => ByRow(psm -> isDecoy(prosit_precs[psm[:precursor_idx]])) => :decoy)
        transform!( scored_PSMs, AsTable(:) => ByRow(psm -> Float64(getIRT(prosit_precs[psm[:precursor_idx]]))) => :iRT)
        transform!( scored_PSMs, AsTable(:) => ByRow(psm -> Float64(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT)
        scored_PSMs[:,:RT_pred] =  GLM.predict(fit1, ns1( scored_PSMs[:,:iRT]))
        scored_PSMs[:,:RT_error] = abs.( scored_PSMs[:,:RT_pred] .-  scored_PSMs[:,:RT])#abs.(PSMs[:,:RT] .- GLM.predict(fit1, ns1(PSMs[:,:iRT])))
        return fragmentMatches,scored_PSMs
        push!(score_times, score)=#
        #return DataFrame(scored_PSMs), MS1
        #return X, H, W, matched, sort(DataFrame(scored_PSMs),:unmatched)
    end

      println("processed $ms2 scans!")
    println("counter: ", times[:counter])
    println("reset: ", times[:reset])
    println("nmf : ", times[:nmf])
    #=println("mean matches: ", mean(match_times))
    println("mean matches: ", median(match_times))=#
    return DataFrame(scored_PSMs)#X, H, IDtoROW, sort(DataFrame(scored_PSMs),:spectral_contrast)#precs#DataFrame(scored_PSMs)# test_frags, test_matches, test_misses#DataFrame(scored_PSMs)
end

#=for frag in prosit_simple_intensities
    if getPrecID(frag) == 4211123 
        println(getPrecMZ(frag))
    end
end=#
#=
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
    lambda::Float64 = 1e3,
    scan::Int = 10000, 
    max_peaks::Int = 200, 
    #fragment_match_ppm::U,
    data_type::Type{T} = Float64
    ) where {T,U<:Real}
println("TEST3")
scored_PSMs = makePSMsDict(XTandem(data_type))
#scored_PSMs = makePSMsDict(FastXTandem(data_type))
#precursorList needs to be sorted by precursor MZ. 
#Iterate through rows (spectra) of the .raw data. 
#i = 0
ms2 = 0
min_intensity = Float32(0.0)
X, H, IDtoROW = "", "", ""
#precs = Dict{UInt32, UInt8}()
#precs = Accumulator{UInt32, UInt8}()
precs = Counter(UInt32, UInt8, Float32, 9387261)
match_times = []
for (i, spectrum) in enumerate(Tables.namedtupleiterator(spectra))
if spectrum[:msOrder] != 2
continue
end
ms2 += 1
if (ms2 % 1000) == 0
#println("ms2: $ms2")
if ms2 == 8000
#test_frags = transitions
#test_matches = fragmentMatches
#test_misses = fragmentMisses
end
end
if i !=  scan
continue
end
#= if ms2 < 100000#50000
#println("TEST")
continue
elseif ms2 > 110000#51000#51000
continue
end=#

#println("ms2: $ms2")
#if ms2 != scan#100024
#    continue
#end
#println("ms2: $ms2")
fragmentMatches = Vector{FragmentMatch{Float32}}()
#println(" a scan")
#precs = Dict{UInt32, UInt8}(zeros(UInt8, pre_aloc_size), zeros(UInt32, pre_aloc_size), zeros(UInt8, pre_aloc_size), 0, 0, 0, pre_aloc_size, 0)
#precs = Dict{UInt32, UInt8}()

#frag_time = @elapsed 
#Ignore pekas below this intensity
min_intensity = spectrum[:intensities][sortperm(spectrum[:intensities], rev = true)[min(max_peaks, length(spectrum[:intensities]))]]

#println(min_intensity)

frag_time = @elapsed prec_count, match_count = searchScan!(precs,
    frag_index, 
    min_intensity, spectrum[:masses], spectrum[:intensities], spectrum[:precursorMZ], 
    fragment_tolerance, 
    precursor_tolerance,
    min_frag_count = min_frag_count, 
    topN = topN
    )

#println(length(prec_counts))

#println("getSize(precs)", getSize(precs))
if getSize(precs) <= 1
#println("TEST")
continue
end
transitions = selectTransitions(fragment_list, precs, topN)
if precs.size > 1
frag_time += @elapsed reset!(precs)
reset!(precs)
push!(match_times, frag_time)
else
push!(match_times, frag_time)
continue
end
reset!(precs)
#println(length(transitions))
#println(precs.size)
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
                    min_intensity = min_intensity,
                    ppm = fragment_tolerance
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

#=weights = (NMF.solve!(NMF.ProjectedALS{Float32}(maxiter=50, verbose = false, 
                                    lambda_w = lambda, 
                                    tol = 1e-8, #Need a reasonable way to choos lambda?
                                    update_H = false #Important to keep H constant. 
                                    ), X, W, H).W[1,:])=#

weights = (NMF.solve!(NMF.GreedyCD{Float32}(maxiter=50, verbose = false, 
                                    lambda_w = lambda, 
                                    tol = 1e-8, #Need a reasonable way to choos lambda?
                                    update_H = false #Important to keep H constant. 
                                    ), X, W, H).W[1,:])

spectral_contrast = getSpectralContrast(H, X)
scribe_score = getScribeScore(H, X)
#For progress and debugging. 


unscored_PSMs = UnorderedDictionary{UInt32, XTandem{T}}()

ScoreFragmentMatches!(unscored_PSMs, fragmentMatches)

Score!(scored_PSMs, unscored_PSMs, 
length(spectrum[:intensities]), 
Float64(sum(spectrum[:intensities])), 
match_count/prec_count, spectral_contrast, scribe_score, weights, IDtoROW,
scan_idx = Int64(i)
)
scored_PSMs = sort(DataFrame(scored_PSMs),:spectral_contrast)
transform!( scored_PSMs, AsTable(:) => ByRow(psm -> isDecoy(prosit_precs[psm[:precursor_idx]])) => :decoy)
transform!( scored_PSMs, AsTable(:) => ByRow(psm -> Float64(getIRT(prosit_precs[psm[:precursor_idx]]))) => :iRT)
transform!( scored_PSMs, AsTable(:) => ByRow(psm -> Float64(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT)
scored_PSMs[:,:RT_pred] =  GLM.predict(fit1, ns1( scored_PSMs[:,:iRT]))
scored_PSMs[:,:RT_error] = abs.( scored_PSMs[:,:RT_pred] .-  scored_PSMs[:,:RT])#abs.(PSMs[:,:RT] .- GLM.predict(fit1, ns1(PSMs[:,:iRT])))

return fragmentMatches,scored_PSMs
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
return DataFrame(scored_PSMs)#X, H, IDtoROW, sort(DataFrame(scored_PSMs),:spectral_contrast)#precs#DataFrame(scored_PSMs)# test_frags, test_matches, test_misses#DataFrame(scored_PSMs)
end=#
#good_scan 100024
#bad_scan 101357

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
#transform!(test_PSMs_gcd, AsTable(:) => ByRow(psm -> isDecoy(prosit_precs[psm[:precursor_idx]])) => :decoy)
#transform!(test_PSMs_gcd, AsTable(:) => ByRow(psm -> Float64(getIRT(prosit_precs[psm[:precursor_idx]]))) => :iRT)
#transform!(test_PSMs_gcd, AsTable(:) => ByRow(psm -> Float64(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT)
#test_PSMs_gcd[:,:RT_pred] =  GLM.predict(fit1, ns1(test_PSMs_gcd[:,:iRT]))
#test_PSMs_gcd[:,:RT_error] = abs.(test_PSMs_gcd[:,:RT_pred] .- test_PSMs_gcd[:,:RT])#abs.(PSMs[:,:RT] .- GLM.predict(fit1, ns1(PSMs[:,:iRT])))

struct runningDotP{C<:Unsigned,T<:AbstractFloat}
    count::C
    obs::T
    pred::T
    obs_pred::T
end

function update!(dp::runningDotP{C,T}, obs::T, pred::T) where {C<:Unsigned,T<:AbstractFloat}
    return runningDotP{C, T}(
            dp.count + one(C),
            dp.obs + obs*obs,
            dp.pred + pred*pred,
            dp.obs_pred + obs*pred)
end

runningDotP(C::DataType, T::DataType) = runningDotP(zero(C), zero(T), zero(T), zero(T)) 
getCount(dp::runningDotP{C,T}) where {C<:Unsigned,T<:AbstractFloat} = dp.count
getObs(dp::runningDotP{C,T}) where {C<:Unsigned,T<:AbstractFloat} = dp.obs
getPred(dp::runningDotP{C,T}) where {C<:Unsigned,T<:AbstractFloat} = dp.pred
getObsPred(dp::runningDotP{C,T}) where {C<:Unsigned,T<:AbstractFloat} = dp.obs_pred
getDP(dp::runningDotP{C,T}) where {C<:Unsigned,T<:AbstractFloat} = getObsPred(dp)/(sqrt(getPred(dp))*sqrt(getObs(dp)))

mutable struct Counter{I,C<:Unsigned,T<:AbstractFloat}
    ids::Vector{I}
    dotp::Vector{runningDotP{C,T}}
    size::Int64
    matches::Int64
    function Counter(I::DataType, C::DataType, T::DataType, size::Int) #where {I,C<:Unsigned}
        new{I, C, T}(zeros(I, size), Vector{runningDotP{C, T}}([runningDotP(C, T) for x in range(1,size)]), 1, 0)
    end
end

getSize(c::Counter{I,C,T}) where {I,C<:Unsigned,T<:AbstractFloat} = c.size
incSize!(c::Counter{I,C,T}) where {I,C<:Unsigned,T<:AbstractFloat} = c.size += 1
getID(c::Counter{I,C,T}, idx::Int) where {I,C<:Unsigned,T<:AbstractFloat} = c.ids[idx]
getRunningDP(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = c.dotp[id]
#getObs(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getObs(c.dotp[id])
#getPred(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getPred(c.dotp[id])
#getObsPred(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getObsPred(c.dotp[id])
#getCount(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getCount(c.dotp[id])
getDP(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getDP(c.dotp[id])

import DataStructures.inc!
function inc!(c::Counter{I,C,T}, id::I, pred_intensity::T, obs_intensity::T) where {I,C<:Unsigned,T<:AbstractFloat} 
    rdp = c.dotp[id]
    if iszero(getCount(rdp))#c.counts[id]<1#iszero(c.counts[id])# == zero(C)
        c.ids[getSize(c)] = id;
        c.dotp[id] = update!(rdp, obs_intensity, pred_intensity);
        incSize!(c);
    else
        c.dotp[id] = update!(rdp, obs_intensity, pred_intensity);
    end
end

import Base.sort!
function sort!(counter::Counter{I,C,T}, topN::Int) where {I,C<:Unsigned,T<:AbstractFloat} 
    sort!(
                @view(counter.ids[1:counter.matches]), 
                by = id -> getCount(counter.dotp[id]),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )
end

function reset!(c::Counter{I,C,T}) where {I,C<:Unsigned,T<:AbstractFloat} 
    #@turbo  for i in 1:(getSize(c) - 1)
    #println("TESTA")
    for i in 1:(getSize(c) - 1)
        c.dotp[c.ids[i]] = runningDotP{C, T}(
                                            zero(C),
                                            zero(T),
                                            zero(T),
                                            zero(T)
                                            )
    end
    c.size = 1
    c.matches = 0
end

function countFragMatches(c::Counter{I,C,T}, min_count::Int) where {I,C<:Unsigned,T<:AbstractFloat} 
    frag_counts = 0
    for i in 1:(getSize(c) - 1)
        id = c.ids[i]
        frag_count = getCount(getRunningDP(c, id))
        dp = getDP(c, id)
        frag_counts += frag_count
        if frag_count >= min_count
            if dp>=0.65
                c.ids[c.matches + 1] = c.ids[i]
                c.matches += 1
            end
        end
    end
    return frag_counts
end


function countFragMatches(c::Counter{I,C,T}, min_count::Int) where {I,C<:Unsigned,T<:AbstractFloat} 
    frag_counts = 0
    for i in 1:(getSize(c) - 1)
        id = c.ids[i]
        frag_count = getCount(getRunningDP(c, id))
        dp = getDP(c, id)
        frag_counts += frag_count
        if frag_count >= min_count
            if dp>=0.50
                c.ids[c.matches + 1] = c.ids[i]
                c.matches += 1
            end
        end
    end
    return frag_counts
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