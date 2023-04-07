function SearchRAW(
                   spectra::Arrow.Table, precursorList::Vector{Precursor}, selectTransitions, 
                   right_precursor_tolerance::Float32,
                   left_precursor_tolerance::Float32,
                   transition_charges::Vector{UInt8},
                   transition_isotopes::Vector{UInt8},
                   b_start::Int64,
                   y_start::Int64,
                   fragment_match_ppm::Float32,
                   ms_file_idx::UInt32,
                   )
    
    scored_PSMs = makePSMsDict(FastXTandem())
    #not goinig to want this for all methods
    allFragmentMatches = Vector{FragmentMatch}()

    #precursorList needs to be sorted by precursor MZ. 
    #Iterate through rows (spectra) of the .raw data. 
    i = 0
    skip = 1
    for spectrum in Tables.namedtupleiterator(spectra)
        i += 1
        if isequal(spectrum[:precursorMZ], missing)
            skip += 1
            continue
        end
        #
        #params = getSpectrumSpecificParams(spectrum, selectParams)

        transitions = selectTransitions(spectrum[:precursorMZ], precursorList, 
                                        right_precursor_tolerance,
                                        left_precursor_tolerance,
                                        transition_charges,
                                        transition_isotopes,
                                        b_start,
                                        y_start,
                                        fragment_match_ppm
                                        )
        #Named tuple for scan 

        fragmentMatches = matchPeaks(transitions, spectrum[:masses], spectrum[:intensities], 
                                    #δs = params[:δs],
                                    δs = [Float64(0)],
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx)
        #sort(getTransitions(selectPrecursors(spectrum)),  
        #                    by = t-> getMZ(t))
        unscored_PSMs = UnorderedDictionary{UInt32, FastXTandem}()

        ScoreFragmentMatches!(unscored_PSMs, fragmentMatches)
        
        #Score!(scored_PSMs, unscored_PSMs, scan_idx = Int64(spectrum[:scanNumber]))
        Score!(scored_PSMs, unscored_PSMs, scan_idx = Int64(i))

        append!(allFragmentMatches, fragmentMatches)
    end
    (scored_PSMs, allFragmentMatches)
end