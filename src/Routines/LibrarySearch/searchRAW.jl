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
    println("max_peak_width $max_peak_width")

    ##########
    #Initialize 
    weights = Float32[]
    msms_counts = Dict{Int64, Int64}()
    frag_err_idx = 1
    chrom_idx = 1
    prec_idx = 0
    ion_idx = 0
    cycle_idx = 0
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
  
    #weights
    precursor_weights = ""
    if ismissing(ion_list)
        precursor_weights = zeros(Float32, maximum(keys(isotope_dict)))
    else
        precursor_weights = zeros(Float32, length(ion_list))
    end

    fragment_intensities = Dictionary{String, Vector{Tuple{Float32, Float32}}}()
    solve_time = 0.0
    index_search_time = 0.0
    prep_time = 0.0
    index_ions_time = 0.0
    ##########
    #Iterate through spectra
    for i in ProgressBar(range(first(scan_range), last(scan_range)))
    #for i in range(1, size(spectra[:masses])[1])

        ###########
        #Scan Filtering
        #(i%10000) == 0 ? println(i) : nothing
        msn = spectra[:msOrder][i] #An integer 1, 2, 3.. for MS1, MS2, MS3 ...
        if msn == 1
            cycle_idx += 1
        end
        msn ∈ spec_order ? nothing : continue #Skip scans outside spec order. (Skips non-MS2 scans is spec_order = Set(2))
        msn ∈ keys(msms_counts) ? msms_counts[msn] += 1 : msms_counts[msn] = 1 #Update counter for each MSN scan type

        (i >= first(scan_range)) & (i <= last(scan_range)) ? nothing : continue #Skip if outside the scan range
        first(rand(1)) <= sample_rate ? nothing : continue #dice-roll. Usefull for random sampling of scans. 

        min_intensity = getMinIntensity(spectra[:intensities][i], max_peaks) #Ignore peaks in the spectrum below this minimum intensity

        iRT_low, iRT_high = getRTWindow(iRT_to_RT_spline(spectra[:retentionTime][i])::Union{Float64,Float32}, maximum_rt, minimum_rt, rt_tol) #Convert RT to expected iRT window

        ##########
        #Ion Template Selection
        #SearchScan! applies a fragment-index based search (MS-Fragger) to quickly identify high-probability candidates to explain the spectrum  
        if !ismissing(searchScan!) | !ismissing(frag_index)
            index_search_time += @elapsed prec_count, match_count = searchScan!(precs, #counter which keeps track of plausible matches 
                        frag_index, 
                        min_intensity, spectra[:masses][i], spectra[:intensities][i], spectra[:precursorMZ][i], 
                        iRT_low, iRT_high,
                        Float32(fragment_tolerance), 
                        Float32(precursor_tolerance),
                        Float32(quadrupole_isolation_width/2.0),
                        min_frag_count = 1,#min_frag_count, 
                        min_ratio = 0.8f0,#Float32(min_matched_ratio),
                        topN = topN
                        )
            
        end
        #selectIons! 
        #Get a sorted list by m/z of ion templates (fills ionTemplates). The spectrum will be searched for matches to these ions only.
        if ismissing(isotope_dict) 
            
            index_ions_time += @elapsed ion_idx, prec_idx = selectIons!(ionTemplates, 
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
        prep_time += @elapsed nmatches, nmisses = matchPeaks(ionTemplates, #Search the spectra for these ions 
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
        #=for m in ionMatches
            if m.prec_id ==9301047
                key = m.ion_type*string(m.frag_index)*"+"*string(m.frag_charge)
                if haskey(fragment_intensities)
                    push!(fragment_intensities[key], (spectra[:retentionTime][i], m.intensity))
                else
                    insert!(fragment_intensities, key, [(spectra[:retentionTime][i], m.intensity)])
                end
            end
        end=#
        ##########
        #Spectral Deconvolution and Distance Metrics 
        if nmatches < 2 #Few matches to do not perform de-convolution 
            IDtoROW = UnorderedDictionary{UInt32, Tuple{UInt32, UInt8}}()
        else #Spectral deconvolution. Build sparse design/template matrix for nnls regression 
            prep_time += @elapsed X, Hs, IDtoROW, last_matched_col = buildDesignMatrix(ionMatches, ionMisses, nmatches, nmisses, H_COLS, H_ROWS, H_VALS)
            #println("Hs.n ", Hs.n)
            #println("Hs.m ", Hs.m)
            #=for match in range(1,nmatches)
                m = ionMatches[match]
                if m.prec_id ==9301047
                    key = m.ion_type*string(m.frag_index)*"+"*string(m.frag_charge)
                    if haskey(fragment_intensities, key)
                        push!(fragment_intensities[key], (spectra[:retentionTime][i], m.intensity))
                    else
                        insert!(fragment_intensities, key, [(spectra[:retentionTime][i], m.intensity)])
                    end
                end
            end=#
            #Initial guess from non-negative least squares. May need to reconsider if this is beneficial
            weights = zeros(eltype(Hs), Hs.n)#sparseNMF(Hs, X, λ, γ, regularize, max_iter=max_iter, tol=nmf_tol)[:]
            for (id, row) in pairs(IDtoROW)
                weights[first(row)] = precursor_weights[id]
            end
            #weights = sparseNMF(Hs, X, λ, γ, regularize, max_iter=max_iter, tol=nmf_tol)[:]
            #weights = sparseNMF(Hs, X, λ, γ, regularize, max_iter=max_iter, tol=Hs.n)[:]
            #println("i $i")
            #println("size(Hs) ", size(Hs))
            #println("size(X) ", size(X))
            solve_time += @elapsed solveHuber!(Hs, Hs*weights .- X, weights, Float32(1000), max_iter_outer = 100, max_iter_inner = 20, tol = Hs.n);

            for (id, row) in pairs(IDtoROW)
                precursor_weights[id] = weights[first(row)]# = precursor_weights[id]
            end
            #weights = sparseNMF(Hs, X, λ, γ, regularize, max_iter=max_iter, tol=nmf_tol)[:]
            #Spectral distance metrics between the observed spectrum (X) and the library spectra for each precursor (Hs)
            scores = getDistanceMetrics(X, Hs, last_matched_col)

            ##########
            #Scoring and recording data
            if !ismissing(scored_PSMs)
                unscored_PSMs = UnorderedDictionary{UInt32, XTandem{Float32}}()

                ScoreFragmentMatches!(unscored_PSMs, ionMatches, nmatches, err_dist)
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
                        min_frag_count = 2#min_frag_count #Remove precursors with fewer fragments 
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
            chrom_idx = fillChroms!(chromatograms, IDtoROW, chrom_idx, i, cycle_idx, prec_ids, prec_idx, frag_counts, weights, spectra[:retentionTime][i])
        end

        ##########
        #Reset pre-allocated arrays 
        reset!(ionTemplates, ion_idx)
        reset!(ionMatches, nmatches), reset!(ionMisses, nmisses)
        fill!(prec_ids, zero(UInt32))

    end

    println("solve_time $solve_time")
    println("prep_time $prep_time")
    println("index_search_time $index_search_time")
    println("index_ions_time $index_ions_time")
    #return fragment_intensities
    ############
    #Return Chromatograms and Score/Feature Table
    #return all_fmatches
    if collect_fmatches
        return DataFrame(scored_PSMs), all_fmatches
    else
        if ismissing(chromatograms)
            #return all_fmatches
            return DataFrame(scored_PSMs)
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
function firstSearch(
    #Mandatory Args
    spectra::Arrow.Table, 
    frag_index::FragmentIndex{Float32},
    ion_list::Vector{Vector{LibraryFragment{Float32}}},
    iRT_to_RT_spline::Any,
    ms_file_idx::UInt32,
    err_dist::Laplace{Float64},
    params::Dict;
    scan_range = (0, 0))

    #=
    main_search_params = (
    expected_matches = 1000000,
    frag_err_dist = frag_err_dist_dict[1],
    frag_tol_quantile = 0.975,
    max_iter = 1000,
    max_peaks = false,
    min_frag_count = 4,
    min_matched_ratio = Float32(0.45),
    min_spectral_contrast = Float32(0.5),
    nmf_tol = Float32(100),
    precursor_tolerance = 5.0,
    quadrupole_isolation_width = 4.25,
    regularize = false,
    rt_tol = 20.0,
    sample_rate = 1.0,
    scan_range = (0, 3000000),
    topN = 100,
    λ = zero(Float32),
    γ = zero(Float32)
    )
    mainLibrarySearch(
           Arrow.Table(MS_TABLE_PATHS[1]),
           prosit_mouse_33NCEcorrected_start1_5ppm_15irt,  
           frags_mouse_detailed_33NCEcorrected_start1, 
           RT_to_iRT_map_dict[1], #RT to iRT map'
           UInt32(1), #MS_FILE_IDX
           frag_err_dist_dict[1],
           main_search_params
       );
    =#
    return SearchRAW(
        spectra, 
        frag_index, ion_list,
        iRT_to_RT_spline,
        ms_file_idx,
        err_dist,
        selectTransitions!,
        searchScan!,
        
        collect_fmatches = true,
        expected_matches = params[:expected_matches],
        frag_ppm_err = params[:frag_ppm_err],
        fragment_tolerance = params[:fragment_tolerance],
        max_iter = params[:max_iter],
        max_peaks = params[:max_peaks],
        min_frag_count = params[:min_frag_count],
        min_matched_ratio = params[:min_matched_ratio],
        min_spectral_contrast = params[:min_spectral_contrast],
        nmf_tol = params[:nmf_tol],
        precs = Counter(UInt32, UInt8, Float32, length(ion_list)),
        precursor_tolerance = params[:precursor_tolerance],
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        regularize = params[:regularize],
        rt_bounds = params[:rt_bounds],
        rt_tol = params[:rt_tol],
        sample_rate = params[:sample_rate],
        scan_range = scan_range,
        scored_PSMs = makePSMsDict(XTandem(Float32)),
        topN = params[:topN],
        λ = params[:λ],
        γ = params[:γ]
    )
end

function mainLibrarySearch(
    #Mandatory Args
    spectra::Arrow.Table, 
    frag_index::FragmentIndex{Float32},
    ion_list::Vector{Vector{LibraryFragment{Float32}}},
    iRT_to_RT_spline::Any,
    ms_file_idx::UInt32,
    err_dist::Laplace{Float64},
    params::Dict;
    scan_range::Tuple{Int64, Int64} = (0, 0))

    frag_ppm_err = err_dist.μ
    fragment_tolerance = quantile(err_dist, params[:frag_tol_quantile])

    #=
    main_search_params = (
    expected_matches = 1000000,
    frag_err_dist = frag_err_dist_dict[1],
    frag_tol_quantile = 0.975,
    max_iter = 1000,
    max_peaks = false,
    min_frag_count = 4,
    min_matched_ratio = Float32(0.45),
    min_spectral_contrast = Float32(0.5),
    nmf_tol = Float32(100),
    precursor_tolerance = 5.0,
    quadrupole_isolation_width = 4.25,
    regularize = false,
    rt_tol = 20.0,
    sample_rate = 1.0,
    scan_range = (0, 3000000),
    topN = 100,
    λ = zero(Float32),
    γ = zero(Float32)
    )
    mainLibrarySearch(
           Arrow.Table(MS_TABLE_PATHS[1]),
           prosit_mouse_33NCEcorrected_start1_5ppm_15irt,  
           frags_mouse_detailed_33NCEcorrected_start1, 
           RT_to_iRT_map_dict[1], #RT to iRT map'
           UInt32(1), #MS_FILE_IDX
           frag_err_dist_dict[1],
           main_search_params
       );
    =#
    return SearchRAW(
        spectra, 
        frag_index, 
        ion_list,
        iRT_to_RT_spline,
        ms_file_idx,
        err_dist,
        selectTransitions!,
        searchScan!,
        
        expected_matches = params[:expected_matches],
        frag_ppm_err = frag_ppm_err,
        fragment_tolerance = fragment_tolerance,
        max_iter = params[:max_iter],
        max_peaks = params[:max_peaks],
        min_frag_count = params[:min_frag_count],
        min_matched_ratio = params[:min_matched_ratio],
        min_spectral_contrast = params[:min_spectral_contrast],
        nmf_tol = params[:nmf_tol],
        precs = Counter(UInt32, UInt8, Float32, length(ion_list)),
        precursor_tolerance = params[:precursor_tolerance],
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        regularize = params[:regularize],
        rt_bounds = params[:rt_bounds],
        rt_tol = params[:rt_tol],
        sample_rate = 1.0,
        scan_range = scan_range,
        scored_PSMs = makePSMsDict(XTandem(Float32)),
        topN = params[:topN],
        λ = params[:λ],
        γ = params[:γ]
    )
end

function integrateMS2(
    #Mandatory Args
    spectra::Arrow.Table,
    ion_list::Vector{Vector{LibraryFragment{Float32}}},
    rt_index::retentionTimeIndex{U, T},
    ms_file_idx::UInt32,
    err_dist::Laplace{Float64},
    params::Dict; 
    N = 600000*10,
    scan_range =  (0, 0)) where {U,T<:AbstractFloat}

    frag_ppm_err = err_dist.μ
    fragment_tolerance = quantile(err_dist, params[:frag_tol_quantile])


    #=integrate_ms2_params = (
    
    expected_matches = 1000000,
    frag_err_dist = frag_err_dist_dict[1],
    frag_tol_quantile = 0.975,
    max_iter = 1000,
    max_peak_width = 2.0,
    max_peaks = false,
    min_frag_count = 4,
    min_matched_ratio = Float32(0.45),
    min_spectral_contrast = Float32(0.5),
    nmf_tol = Float32(100),
    precursor_tolerance = 5.0,
    quadrupole_isolation_width = 4.25,
    regularize = false,
    rt_tol = 20.0,
    sample_rate = 1.0,
    scan_range = (0, 3000000),
    topN = 100,
    λ = zero(Float32),
    γ = zero(Float32)
    )
    
    integrateMS2(
           Arrow.Table(MS_TABLE_PATHS[1]), 
           prosit_mouse_33NCEcorrected_start1_5ppm_15irt,  
           frags_mouse_detailed_33NCEcorrected_start1, 
           rt_index, 
           RT_to_iRT_map_dict[1], #RT to iRT map'
           UInt32(1), #MS_FILE_IDX
           frag_err_dist_dict[1],
           integrate_ms2_params
       )
    =#

    return SearchRAW(
        spectra, 
        missing, 
        ion_list,
        x->x,
        ms_file_idx,
        err_dist,
        selectRTIndexedTransitions!,
        missing,
        
        chromatograms =  Dict(:precursor_idx => zeros(UInt32, N), 
                                :weight => zeros(Float32, N), 
                                :scan_idx => zeros(UInt32, N),
                                :rt => zeros(Float32, N), 
                                :frag_count => zeros(Int64, N),
                                :rank => zeros(UInt8, N),
                                :cycle_idx => zeros(UInt32, N)),
        expected_matches = params[:expected_matches],
        frag_ppm_err = frag_ppm_err,
        fragment_tolerance = fragment_tolerance,
        max_iter = params[:max_iter],
        max_peak_width = params[:max_peak_width],
        max_peaks = params[:max_peaks],
        min_frag_count = params[:min_frag_count],
        min_matched_ratio = params[:min_matched_ratio],
        min_spectral_contrast = params[:min_spectral_contrast],
        nmf_tol = params[:nmf_tol],
        precursor_tolerance = params[:precursor_tolerance],
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        regularize = params[:regularize],
        rt_index = rt_index,
        rt_bounds = params[:rt_bounds],
        rt_tol = params[:rt_tol],
        sample_rate = params[:sample_rate],
        scan_range = scan_range,
        topN = params[:topN],
        λ = params[:λ],
        γ = params[:γ]
    )
end

function integrateMS1(
    #Mandatory Args
    spectra::Arrow.Table, 
    isotope_dict::UnorderedDictionary{UInt32, Vector{Isotope{Float32}}},
    prec_rt_list::Vector{Tuple{T, UInt32}},
    ms_file_idx::UInt32,
    err_dist::Laplace{Float64},
    params::Dict; 
    N = 600000*10,
    scan_range = (0, 0)) where {T<:AbstractFloat}

    frag_ppm_err = err_dist.μ
    fragment_tolerance = quantile(err_dist, params[:frag_tol_quantile])

    #=
    integrate_ms1_params = (
    expected_matches = 1000000,
    frag_err_dist = frag_err_dist_dict[1],
    frag_tol_quantile = 0.975,
    max_iter = 1000,
    max_peak_width = 2.0,
    max_peaks = false,
    min_frag_count = 4,
    min_matched_ratio = Float32(0.45),
    min_spectral_contrast = Float32(0.5),
    nmf_tol = Float32(100),
    precursor_tolerance = 5.0,
    quadrupole_isolation_width = 4.25,
    regularize = false,
    rt_tol = 20.0,
    sample_rate = 1.0,
    scan_range = (0, 3000000),
    topN = 100,
    λ = zero(Float32),
    γ = zero(Float32)

    integrateMS1(
    Arrow.Table(MS_TABLE_PATHS[1]), 
    prosit_mouse_33NCEcorrected_start1_5ppm_15irt,  
    isotopes,
    prec_rt_table, 
    RT_to_iRT_map_dict[1], #RT to iRT map'
    UInt32(1), #MS_FILE_IDX
    frag_err_dist_dict[1],
    integrate_ms1_params)
    )=#
    return SearchRAW(
        spectra, 
        missing, 
        missing, #Not ion list. Instead passing "isotope_dict"
        x->x,
        ms_file_idx,
        err_dist,  #Not really sure how to estimate this yet?
        selectIsotopes!, #Ion Selection Function for MS1 integration 
        missing,
        
        chromatograms =  Dict(:precursor_idx => zeros(UInt32, N), 
        :weight => zeros(Float32, N), 
        :scan_idx => zeros(UInt32, N),
        :rt => zeros(Float32, N), 
        :frag_count => zeros(Int64, N),
        :rank => zeros(UInt8, N),
        :cycle_idx => zeros(UInt32, N)),
                expected_matches = params[:expected_matches],
        frag_ppm_err = frag_ppm_err,
        fragment_tolerance = fragment_tolerance,
        IonTemplateType = Isotope{Float32},
        IonMatchType = PrecursorMatch{Float32},
        isotope_dict = isotope_dict,
        max_iter = params[:max_iter],
        max_peak_width = params[:max_peak_width],
        max_peaks = params[:max_peaks],
        min_frag_count = params[:min_frag_count],
        min_matched_ratio = params[:min_matched_ratio],
        min_spectral_contrast = params[:min_spectral_contrast],
        nmf_tol = params[:nmf_tol],
        precursor_tolerance = params[:precursor_tolerance],
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        regularize = params[:regularize],
        rt_bounds = params[:rt_bounds],
        rt_index = prec_rt_list,
        rt_tol = params[:rt_tol],
        sample_rate = params[:sample_rate],
        scan_range = scan_range,
        spec_order = Set(1),
        topN = params[:topN],
        λ = params[:λ],
        γ = params[:γ]
    )
end


function fillChroms!(chroms::Dict{Symbol, Vector}, id_to_row::UnorderedDictionary{UInt32, Tuple{UInt32, UInt8}}, n::Int64, scan_idx::Int64, cycle_idx::Int64, prec_ids::Vector{UInt32}, prec_idx::Int64, frag_counts::Accumulator{UInt32,Int64}, weights::Vector{T}, retention_time::U; block_size = 100000) where {T,U<:AbstractFloat}
    function inc!(chroms::Dict{Symbol, Vector}, n::Int64, scan_idx::Int64, cycle_idx::Int64, key::UInt32, weight::AbstractFloat, rt::AbstractFloat, frag_count::Int64,rank::UInt8)
        chroms[:precursor_idx][n] = key
        chroms[:weight][n] = weight
        chroms[:rt][n] = rt
        chroms[:scan_idx][n] = scan_idx
        chroms[:frag_count][n] = frag_count
        chroms[:rank][n] = rank
        chroms[:cycle_idx][n] = UInt32(cycle_idx)
    end

    for i in range(1, prec_idx)
        key = prec_ids[i]
        
        if n > length(chroms[:precursor_idx])
            #This block of code is triggered when not enough space has been allocated for the chromatograms
            #So incrase the space allocated by a given "block" size. 
            for key in keys(chroms)
                append!(chroms[key], zeros(eltype(chroms[key]), block_size))
            end
        end
        if haskey(frag_counts, key)
    
            if haskey(id_to_row, key)
                inc!(chroms, n, scan_idx, cycle_idx, key, weights[id_to_row[key][1]], retention_time, frag_counts[key], id_to_row[key][2])
            else
                inc!(chroms, n, scan_idx, cycle_idx, key, Float32(0.0),retention_time, frag_counts[key], zero(UInt8))
            end
        else
            inc!(chroms, n, scan_idx, cycle_idx, key, Float32(0.0), retention_time, 0, zero(UInt8))
        end
        n += 1
    end    
    return n
end

function getRTWindow(irt::U, max_rt::T, min_rt::T, rt_tol::T) where {T,U<:AbstractFloat}
    #(irt < min_rt) ? irt_low = Float32(-Inf) : irt_low = Float32(irt - rt_tol)
    #(irt > max_rt) ? irt_high = Float32(Inf) : irt_high = Float32(irt + rt_tol)
    #return irt_low, irt_high
    return Float32(irt - rt_tol), Float32(irt + rt_tol)
end

function getMinIntensity(intensities::AbstractArray{Union{T, Missing}}, max_peaks::Int) where {T<:AbstractFloat}
    return intensities[sortperm(intensities, rev = true)[min(max_peaks, length(intensities))]]
end

function getMinIntensity(intensities::AbstractArray{Union{T, Missing}}, max_peaks::Bool) where {T<:AbstractFloat}
    return zero(T)
end

function reset!(fms::Vector{M}, last_non_empty::Int64) where {M<:Match}
    for i in range(1, last_non_empty)
        fms[i] = M()
    end
end

#=function reset!(fms::Vector{PrecursorMatch{T}}, last_non_empty::Int64) where {T<:AbstractFloat}
    for i in range(1, last_non_empty)
        fms[i] = PrecursorMatch{T}()
    end
end=#

function reset!(fms::Vector{Isotope{T}}, last_non_empty::Int64) where {T<:AbstractFloat}
    for i in range(1, last_non_empty)
        fms[i] = Isotope{T}()
    end
end

function reset!(lf::Vector{LibraryFragment{T}}, last_non_empty::Int64) where {T<:AbstractFloat}
    for i in range(1, last_non_empty)
        lf[i] = LibraryFragment{T}()
    end
end

function collectFragErrs(all_fmatches::Vector{M}, new_fmatches::Vector{M}, nmatches::Int, n::Int, collect_fmatches::Bool) where {M<:Match}
    if collect_fmatches
        for match in range(1, nmatches)
            if n < length(all_fmatches)
                all_fmatches[n] = new_fmatches[match]
                n += 1
            end
        end
    end
    return n
end

