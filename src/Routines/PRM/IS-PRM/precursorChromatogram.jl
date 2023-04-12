"""
    PrecursorChromatogram

Type that represents an chromatogram for a `Precursor`

### Fields

- rts::Vector{Float32} -- Retention times for scans targeting this precursor
- last_scan_idx::Vector{Int64} --  Keeps track of the last_scan from which a transition has been added (should change to a Ref in the future)
- transitions::UnorderedDictionary{String, Vector{Float32}} -- Dictionary with keys for each transition, and values for their intensities in each matched scan. Those vectors should all have the same length and correspond to the `rts` field. 
- best_psm::NamedTuple{(:rt, :scan_idx, :name, :mz, :intensity), Tuple{Float32, Int64, Vector{String}, Vector{Float32}, Vector{Float32} -- List of transitions, mz's, and their intensities for the best scoring PSM for the precursor.  

### Examples

- PrecursorChromatogram() = PrecursorChromatogram(Float32(0), 0, Float32[], UnorderedDictionary{String, Vector{Float32}}(), (rt = Float32(0), scan_idx = 0, name = String[], mz = Float32[], intensity = Float32[]))

### GetterMethods

- getRTs(pc::PrecursorChromatogram) = pc.rts
- getLastScanIdx(pc::PrecursorChromatogram) = pc.last_scan_idx
- getTransitions(pc::PrecursorChromatogram) = pc.transitions
- getBestPSM(pc::PrecursorChromatogram) = pc.best_psm

### Methods

- initMatchedPrecursors(bestPSMs::DataFrame) -- Initializes a dictionary container of `PrecursorChromatogram` from a DataFrame of best PSMs. 
"""
struct PairedPrecursorChromatogram
    rts::Dict{String, Vector{Int64}}
    transitions::Dict{String, Dict{String, Vector{Float32}}}
end

PrecursorChromatogram() = PrecursorChromatogram(Float32[], Int64[], UnorderedDictionary{String, Vector{Float32}}(), (rt = Float32(0), scan_idx = 0, name = String[], mz = Float32[], intensity = Float32[]))
getRTs(pc::PrecursorChromatogram) = pc.rts
getLastScanIdx(pc::PrecursorChromatogram) = pc.last_scan_idx
getTransitions(pc::PrecursorChromatogram) = pc.transitions
getBestPSM(pc::PrecursorChromatogram) = pc.best_psm

PairedPrecursorChromatogram

"""
    initPrecursorChromatograms(bestPSMs::DataFrame)

Initializes a dictionary container of `PrecursorChromatogram` from a DataFrame of best PSMs. 

### Input

- `bestPSMs::DataFrame` -- DataFrame of best PSMs. See Examples. Must have fields scan_idx, retentionTime, and precursor_idx


### Output
- Returns a `DataFrame` with one row per unique `precursor_idx`. 

### Notes
    `PSMs` has a key for each field in a type that inherits from PSMs. The values are vectors, each of the same length, for the number of PMSs. 
        For example, 

        Dict{Symbol, Vector} with 6 entries:
        :hyperscore    => [0.000286995, 0.000286995, 0.000286995, ...
        :precursor_idx => UInt32[0x000000f1, 0x00000121, 0x000000a2, ...
        :total_ions    => UInt32[0x00000001, 0x00000001, 0x00000001, ...
        :error         => [0.00274658, 0.0233765, 0.00686646, ...
        :scan_idx      => [1720, 1724, 1731, ...
        :y_ladder      => Int8[1, 1, 1, ...

        This enables easy conversion to a  `DataFrame` of PSMs. 

### Examples 
    julia> bestPSMs = FilterPSMs(scored, testPtable, NRF2_Survey, UInt8(5))
    219x10 DataFrame
    Row │ pep_idx  error       hyperscore  precursor_idx  scan_idx  total_ions  y_ladder  retentionTime  sequence                        proteinNames 
        │ UInt32   Float64     Float64     UInt32         Int64     UInt32      Int8      Float32        String                          String       
    ─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    1 │       1  0.0154419      39.6359            184     15893          10         8        30.3842  YYNAESYEVER[Harg]               ABCB6
    2 │       3  0.0412292      47.5865            760     52286          13         7        59.7062  ALNVLVPIFYR[Harg]               ABCB6
    3 │       4  0.0218811      42.695             436      7588          10         6        23.4598  YVSLPNQNK[Hlys]                 ABHD4
"""
function initPrecursorChromatograms(best_psms::DataFrame, ms_file_idx::UInt32)
    precursor_chromatograms = UnorderedDictionary{UInt32, PrecursorChromatogram}()
    for pep in eachrow(best_psms)
        if pep[:ms_file_idx] != ms_file_idx
            continue
        end
        insert!(precursor_chromatograms, 
                pep[:precursor_idx],
                PrecursorChromatogram(Float32[],
                [0],
                UnorderedDictionary{String, Vector{Float32}}(),
                (rt = pep[:retention_time], scan_idx = pep[:scan_idx], name = Vector{String}(), mz = Float32[], intensity = Float32[])))

    end
    precursor_chromatograms
end

function fillPrecursorChromatograms!(precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram}, fragment_matches::Vector{FragmentMatch}, MS_TABLE::Arrow.Table, rt_tol::Float64, ms_file_idx::UInt32)

    """
        getTransitionName(transition::UnorderedDictionary{String, Vector{Float32}}) 
    """
    function getTransitionName(transition::Transition)
        string(transition.ion_type)*string(transition.ind)*"+"*string(transition.charge)
    end

    function getBestPSMs!(precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram}, fragment_matches::Vector{FragmentMatch}, ms_file_idx::UInt32)

        function addScanToPrecursor!(precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram}, match::FragmentMatch, fragment_name::String, prec_id::UInt32, ms_file_idx::UInt32)
            
            function assignPrecursor!(matched_precursor::PrecursorChromatogram, fragment_name::String)
                if !isassigned(matched_precursor.transitions, fragment_name)
                    insert!(matched_precursor.transitions, fragment_name, Float32[])
                end
            end

            assignPrecursor!(precursor_chromatograms[prec_id], fragment_name);

            if match.scan_idx == precursor_chromatograms[prec_id].best_psm[:scan_idx]
                if match.ms_file_idx == ms_file_idx
                    push!(precursor_chromatograms[prec_id].best_psm[:mz], match.match_mz)
                    push!(precursor_chromatograms[prec_id].best_psm[:intensity], match.intensity)
                    push!(precursor_chromatograms[prec_id].best_psm[:name], fragment_name)
                end
            end

        end

        for match in fragment_matches
            #pep_id = getPepIDFromPrecID(testPtable, match.transition.prec_id)
            prec_id = match.transition.prec_id
            if !isassigned(precursor_chromatograms, prec_id) continue end
            #as in y3+2 or b6+1 (name/length/+charge)
            fragment_name = getTransitionName(match.transition)
            addScanToPrecursor!(precursor_chromatograms, match, fragment_name, prec_id, ms_file_idx)
        end

    end

    function addTransition(precursor_chromatogram::PrecursorChromatogram, match::FragmentMatch, rt::Float32)

        if match.scan_idx>precursor_chromatogram.last_scan_idx[1]
            precursor_chromatogram.last_scan_idx[1]=match.scan_idx
            append!(precursor_chromatogram.rts, rt)
            for key in keys(precursor_chromatogram.transitions)
                push!(precursor_chromatogram.transitions[key],Float32(0))
            end
        end
        precursor_chromatogram.transitions[getTransitionName(match.transition)][end] = match.intensity
    end

    function isInRTWindow(precursor_chromatogram::PrecursorChromatogram, rt, rt_tol::Float64)
        best_rt = precursor_chromatogram.best_psm[:rt]
        return (abs(rt- best_rt)<rt_tol)
    end

    function addTransitionIntensities!(precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram}, fragment_matches::Vector{FragmentMatch}, MS_TABLE::Arrow.Table, rt_tol::Float64)
        for matched_fragment in fragment_matches
            prec_id = matched_fragment.transition.prec_id
            if !isassigned(precursor_chromatograms, prec_id) continue end
            rt = MS_TABLE[:retentionTime][matched_fragment.scan_idx]
            if isInRTWindow(precursor_chromatograms[prec_id], rt, rt_tol) 
                addTransition(precursor_chromatograms[prec_id], matched_fragment, rt) 
                #println("TESTTEST")
            end
        end
    end

    getBestPSMs!(precursor_chromatograms, fragment_matches, ms_file_idx)
    addTransitionIntensities!(precursor_chromatograms, fragment_matches, MS_TABLE, rt_tol)
    return precursor_chromatograms 

end
