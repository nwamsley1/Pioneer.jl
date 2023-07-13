function SearchRAW(
                    spectra::Arrow.Table, 
                    prec_norms::Vector{Float32},
                    #ptable::PrecursorDatabase,
                    frag_index::FragmentIndex{T},
                    fragment_list::Vector{Vector{LibraryFragment{Float64}}},
                    ms_file_idx::UInt32;
                    isolation_width::Float64 = 4.25,
                    precursor_tolerance::Float64 = 5.0,
                    fragment_tolerance::Float64 = 20.0,
                    topN::Int64 = 20,
                    min_frag_count::Int64 = 4,
                    min_matched_ratio::Float32 = Float32(0.8),
                    min_spectral_contrast::Float32 = Float32(0.65),
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
                    prec_norms,
                    frag_index, 
                    min_intensity, spectrum[:masses], spectrum[:intensities], MS1, spectrum[:precursorMZ], 
                    fragment_tolerance, 
                    precursor_tolerance,
                    isolation_width,
                    min_frag_count = min_frag_count, 
                    min_ratio = min_matched_ratio,
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
       
        return X, H, UNMATCHED, IDtoROW, fragmentMatches, fragmentMisses
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
                scan_idx = Int64(i),
                min_spectral_contrast = min_spectral_contrast
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

X, H, UNMATCHED, IDtoROW, fragmentMatches, fragmentMisses = SearchRAW(MS_TABLE, prosit_totals, prosit_index_intensities, prosit_detailed, UInt32(1), 
min_frag_count = 4, 
topN = 200, 
fragment_tolerance = 15.6, 
lambda = 1e5, 
max_peaks = 1000, 
scan_range = (101357, 101357), 
precursor_tolerance = 20.0,
min_spectral_contrast =  Float32(0.65)
)