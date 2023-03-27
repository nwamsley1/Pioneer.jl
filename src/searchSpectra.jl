#=#Need a method that gets the appropriate precursors/transitions for each
#Scan. That might be unique to each method/pipeline. 
using Arrow 
#function CheckRawData(RAW::Arrow.Table)#
#
#
function SearchRAW(spectra::Arrow.Table, precursorList::Vector{Precursor}, selectTransitions, selectParams)
    
    scored_PSMs = makePSMsDict(FastXTandem())

    #not goinig to want this for all methods
    allFragmentMatches = Vector{FragmentMatch}

    #precursorList needs to be sorted by precursor MZ. 
    #Iterate through rows (spectra) of the .raw data. 
    for spectrum in Tables.namedtupleiterator(NRF2_Survey)

        #
        params = getSpectrumSpecificParams(spectrum, selectParams)

        transitions = selectTransitions(spectrum, precursorList, params[:ppm])
        #Named tuple for scan 

        fragmentMatches = matchPeaks(transitions, spectrum[:masses], spectrum[:intensities], δs = params[:δs])
        #sort(getTransitions(selectPrecursors(spectrum)),  
        #                    by = t-> getMZ(t))
        unscored_PSMS = UnorderedDictionary{UInt32, FastXTandem}()

        ScoreFragmentMatches!(unscored_PSMs, fragmentMatches)
        
        Score!(scored_PSMs, unscored_PSMs, scan_idx = spectrum[:scanNumber])

        append!(allFragmentMatches, fragmentMatches)
    end
end

function selectTransitions(spectrum, precursorList, ppm)
lower_bound = 400.0
upper_bound = 401.0
#get all precursors in the tolerance. 
start = searchsortedfirst(test_t, lower_bound,lt=(t,x)->getMZ(t)<x)
stop = searchsortedlast(test_t[start:end], upper_bound,lt=(x,t)->getMZ(t)>x)+start
for precursor in Precursors
#need a getTransitions!() method that appends to a current list. 
#Need arguments for the charge states and isotopes of fragments we want to get. 
#also b_start, y_start etc. 
getTransitions!()
end
=#