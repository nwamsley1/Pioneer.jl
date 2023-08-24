function SearchRAW(
                    #Mandatory Args
                    spectra::Arrow.Table, 
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    ion_list::Union{Vector{Vector{LibraryFragment{Float32}}}, Missing},
                    iRT_to_RT_spline::Any,
                    ms_file_idx::UInt32,
                    err_dist::Laplace{Float64},
                    selectIons!::Function,
                    searchScan!::Union{Function, Missing};
                    #keyword args
                    chromatograms::Union{Dict{Symbol, Vector}, Missing} = missing,
                    collect_fmatches = false,
                    expected_matches::Int64 = 100000,
                    frag_ppm_err::Float64 = 0.0,
                    fragment_tolerance::Float64 = 20.0,
                    IonMatchType::DataType = FragmentMatch{Float32},
                    IonTemplateType::DataType = LibraryFragment{Float32},
                    isotope_dict::Union{UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}, Missing} = missing,
                    max_iter::Int = 1000,
                    max_peak_width::Float64 = 2.0,
                    max_peaks::Union{Int64,Bool} = false, 
                    min_frag_count::Int64 = 4,
                    min_matched_ratio::Float32 = Float32(0.8),
                    min_spectral_contrast::Float32 = Float32(0.65),
                    nmf_tol::Float32 = Float32(100.0),
                    precs::Counter{UInt32, UInt8, Float32} = Counter(UInt32, UInt8, Float32, 0),
                    precursor_tolerance::Float64 = 5.0,
                    quadrupole_isolation_width::Float64 = 8.5,
                    regularize::Bool = false,
                    rt_bounds::Tuple{Float64, Float64} = (0.0, 0.0),
                    rt_index::Union{retentionTimeIndex{Float64, Float32}, Vector{Tuple{Float64, UInt32}}, Missing} = missing,
                    rt_tol::Float64 = 30.0,
                    sample_rate::Float64 = 1.0,
                    scan_range::Tuple{Int64, Int64} = (0, 0),
                    scored_PSMs::Union{Dict{Symbol, Vector}, Missing} = missing,
                    spec_order::Set{Int64} = Set(2),
                    topN::Int64 = 20,
                    λ::Float32 = Float32(1e3),
                    γ::Float32 = zero(Float32))

    ##########
    #Initialize 
    weights = Float32[]
    msms_counts = Dict{Int64, Int64}()
    frag_err_idx = 1
    chrom_idx = 1
    prec_idx = 0
    ion_idx = 0
    minimum_rt, maximum_rt = first(rt_bounds), last(rt_bounds)

    ###########
    #Pre-allocate Arrays to save (lots) of time in garbage collection. 
    all_fmatches = Vector{IonMatchType}()
    collect_fmatches ? all_fmatches = [IonMatchType() for x in range(1, expected_matches)] : nothing

    #These are overwritten for every searched spectrum. "expected_matches"
    #is a guess for the largest array size that would be needed for any spectrum. 
    #If the guess is too small, the arrays will simply be increased in size as needed
    #by a pre-determined block-size. 
    ionMatches = [IonMatchType() for _ in range(1, expected_matches)] #IonMatchType is something that inherits from the "Match" class. 
    ionMisses = [IonMatchType() for _ in range(1, expected_matches)]
    ionTemplates = [IonTemplateType() for _ in range(1, expected_matches)] 
    prec_ids = [zero(UInt32) for _ in range(1, expected_matches)]
    H_COLS, H_ROWS, H_VALS = zeros(Int64, expected_matches), zeros(Int64, expected_matches), zeros(Float32, expected_matches)
  
    prec_counts = Int64[]
    println("TESTTEST")
    ##########
    #Iterate through spectra
    for i in range(1, size(spectra[:masses])[1])

        ###########
        #Scan Filtering
        (i%10000) == 0 ? println(i) : nothing
        msn = spectra[:msOrder][i] #An integer 1, 2, 3.. for MS1, MS2, MS3 ...
        msn ∈ spec_order ? nothing : continue #Skip scans outside spec order. (Skips non-MS2 scans is spec_order = Set(2))
        msn ∈ keys(msms_counts) ? msms_counts[msn] += 1 : msms_counts[msn] = 1 #Update counter for each MSN scan type

        (i >= first(scan_range)) & (i <= last(scan_range)) ? nothing : continue #Skip if outside the scan range
        first(rand(1)) <= sample_rate ? nothing : continue #dice-roll. Usefull for random sampling of scans. 

        min_intensity = getMinIntensity(spectra[:intensities][i], max_peaks) #Ignore peaks in the spectrum below this minimum intensity

        #println(typeof(maximum_rt))
        #println(typeof(spectra[:retentionTime][i]))
        iRT_low, iRT_high = getRTWindow(iRT_to_RT_spline(spectra[:retentionTime][i])::Union{Float64,Float32}, maximum_rt, minimum_rt, rt_tol) #Convert RT to expected iRT window

        ##########
        #Ion Template Selection
        #SearchScan! applies a fragment-index based search (MS-Fragger) to quickly identify high-probability candidates to explain the spectrum  
        if !ismissing(searchScan!) | !ismissing(frag_index)
            prec_count, match_count = searchScan!(precs, #counter which keeps track of plausible matches 
                        frag_index, 
                        min_intensity, spectra[:masses][i], spectra[:intensities][i], spectra[:precursorMZ][i], 
                        iRT_low, iRT_high,
                        Float32(fragment_tolerance), 
                        Float32(precursor_tolerance),
                        Float32(quadrupole_isolation_width/2.0),
                        min_frag_count = min_frag_count, 
                        min_ratio = Float32(min_matched_ratio),
                        topN = topN
                        )
            push!(prec_counts, prec_count);
        end
        #selectIons! 
        #Get a sorted list by m/z of ion templates (fills ionTemplates). The spectrum will be searched for matches to these ions only.
        if ismissing(isotope_dict) 
            
            ion_idx, prec_idx = selectIons!(ionTemplates, 
                                               ion_list,
                                               precs,
                                               topN,
                                               prec_ids,
                                               rt_index,
                                               spectra[:retentionTime][i],
                                               Float32(max_peak_width/2.0), #rt_tol
                                               spectra[:precursorMZ][i], #prec_mz
                                               Float32(quadrupole_isolation_width/2.0) #prec_tol
                                               )::Tuple{Int64, Int64}
            
        else
            ion_idx, prec_idx = selectIons!(
                                            ionTemplates,
                                            rt_index,
                                            isotope_dict,
                                            prec_ids,
                                            spectra[:retentionTime][i],
                                            Float32(max_peak_width/2.0)
                                            )
        end
        ion_idx < 2 ? continue : nothing 

        ##########
        #Match sorted list of plausible ions to the observed spectra
        nmatches, nmisses = matchPeaks(ionTemplates, #Search the spectra for these ions 
                                    ion_idx, #search ionTemplates[1:ion_idx]
                                    ionMatches, #Fill with matched ions 
                                    ionMisses, #Fill with unmatched ions 
                                    spectra[:masses][i], 
                                    spectra[:intensities][i], 
                                    count_unmatched=true, #Should we fill "ionMisses"?
                                    δs = [frag_ppm_err], #Mass offsets 
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = min_intensity, #Ignore peaks below this intensity
                                    ppm = fragment_tolerance #Fragment match tolerance in ppm
                                    )

        ##########
        #Spectral Deconvolution and Distance Metrics 
        if nmatches < 2 #Few matches to do not perform de-convolution 
            #reset!(ionMatches, nmatches), reset!(ionMisses, nmisses) #These arrays are pre-allocated so just overwrite to prepare for the next scan 
            IDtoROW = UnorderedDictionary{UInt32, UInt32}()
        else #Spectral deconvolution. Build sparse design/template matrix for nnls regression 
            X, Hs, IDtoROW, last_matched_col = buildDesignMatrix(ionMatches, ionMisses, nmatches, nmisses, H_COLS, H_ROWS, H_VALS)
            #Non-negative least squares coefficients for each precursor template explaining the spectra 
            weights = sparseNMF(Hs, X, λ, γ, regularize, max_iter=max_iter, tol=nmf_tol)[:]
            #Spectral distance metrics between the observed spectrum (X) and the library spectra for each precursor (Hs)
            scores = getDistanceMetrics(X, Hs, last_matched_col)

            ##########
            #Scoring and recording data
            if !ismissing(scored_PSMs)
                unscored_PSMs = UnorderedDictionary{UInt32, XTandem{Float32}}()

                ScoreFragmentMatches!(unscored_PSMs, ionMatches, nmatches, err_dist)
                #println("min_spectral_contrast ", min_spectral_contrast)
                #println("min_frag_count ", min_frag_count)
                #Score unscored_PSMs and write them to scored_PSMs
                Score!(scored_PSMs, 
                        unscored_PSMs, 
                        length(spectra[:intensities][i]), 
                        Float64(sum(spectra[:intensities][i])), 
                        match_count/prec_count, 
                        scores, #Named Tuple of spectrum simmilarity/distance measures 
                        weights, #Coefficients for each precursor in the spectral deconvolution
                        IDtoROW,
                        scan_idx = i,
                        min_spectral_contrast = min_spectral_contrast, #Remove precursors with spectral contrast lower than this ammount
                        min_frag_count = min_frag_count #Remove precursors with fewer fragments 
                        )
            end
        end
        #Add fragment matches to all_fmatches 
        frag_err_idx = collectFragErrs(all_fmatches, ionMatches, nmatches, frag_err_idx, collect_fmatches)
        
        ##########
        #Update Chromatograms 
        if !ismissing(chromatograms)
            frag_counts = counter(UInt32) #Basically a Dictionary that counts the number of matched ions (values) for each precursor (keys)
            for match_idx in range(1, nmatches) #fragmentMatches
                DataStructures.inc!(frag_counts, ionMatches[match_idx].prec_id)
            end
            #Add precursor templates with their weights and retention times to the chromatogram table 
            chrom_idx = fillChroms!(chromatograms, IDtoROW, chrom_idx, prec_ids, prec_idx, frag_counts, weights, spectra[:retentionTime][i])
        end

        ##########
        #Reset pre-allocated arrays 
        reset!(ionTemplates, ion_idx)
        reset!(ionMatches, nmatches), reset!(ionMisses, nmisses)
        fill!(prec_ids, zero(UInt32))

    end



    ############
    #Return Chromatograms and Score/Feature Table
    println("UP HERE")
    if collect_fmatches
        return DataFrame(scored_PSMs), all_fmatches
    else
        if ismissing(chromatograms)
            println("TEST" , length(prec_counts))
            return DataFrame(scored_PSMs), prec_counts
        elseif ismissing(scored_PSMs)
            chromatograms = DataFrame(chromatograms)
            sort!(chromatograms, [:precursor_idx,:rt], alg=QuickSort);
            return groupby(DataFrame(chromatograms), :precursor_idx)
        else
            chromatograms = DataFrame(chromatograms)
            sort!(chromatograms, [:precursor_idx,:rt], alg=QuickSort);
            return DataFrame(scored_PSMs), groupby(DataFrame(chromatograms), :precursor_idx)
        end
    end

end