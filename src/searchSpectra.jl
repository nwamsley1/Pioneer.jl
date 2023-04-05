#Need a method that gets the appropriate precursors/transitions for each
#Scan. That might be unique to each method/pipeline. 
using Arrow, Table, DataFrames, Tables
using Plots
include("src/precursor.jl")
include("src/matchpeaks.jl")
include("src/getPrecursors.jl")
include("src/PSM_TYPES/PSM.jl")
include("src/PSM_TYPES/FastXTandem.jl")
#function CheckRawData(RAW::Arrow.Table)#
#
# In precursorTable need something that maps prec_id's
# to a more simple precursor representation that is just the string, charge, isotope, and pep_id. 
# then need to get transitions directly from the string w/o going through a `Precursor` intermediate. 
function SearchRAW(spectra::Arrow.Table, precursorList::Vector{Precursor}, selectTransitions, params)
    
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

        transitions = selectTransitions(spectrum[:precursorMZ], precursorList, params)
        #Named tuple for scan 

        fragmentMatches = matchPeaks(transitions, spectrum[:masses], spectrum[:intensities], δs = params[:δs], scan_idx = UInt32(i))
        #sort(getTransitions(selectPrecursors(spectrum)),  
        #                    by = t-> getMZ(t))
        unscored_PSMs = UnorderedDictionary{UInt32, FastXTandem}()

        ScoreFragmentMatches!(unscored_PSMs, fragmentMatches)
        
        #Score!(scored_PSMs, unscored_PSMs, scan_idx = Int64(spectrum[:scanNumber]))
        Score!(scored_PSMs, unscored_PSMs, scan_idx = Int64(i))

        append!(allFragmentMatches, fragmentMatches)
    end
    println("there were this many scans: ", i)
    println("this many scans were skipped: ", skip)
    println("the difference is : ", i-skip)
    (scored_PSMs, allFragmentMatches)
end

"""
    selectTransitionsPRM(window_center::Float32, precursorList::Vector{Precursor}, params)

Given a SORTED vector of `Precursor` in order os ascending MZ, gets the subset of `Precursor` with MZ within the tolerance.
The tolerance is specified based on the `window_center` and `params[:lower_tol]` and `params[:upper_tol]` . Every routine implementing a method `SearchRaw` should
implement a `selectTransitions` method. 

### Input

- `window_center::Float32` -- The isolation window center for an MS2 scan. 
- `precursorList::Vector{Precursor}` -- List of possible `Precursor` 
- `MS_TABLE::Arrow.Table` -- Search parameters. A named tuple that must have fields [:upper_tol] and [:lower_tol]

### Output
- Returns Vector{Precursor} which is a subset of `precursorList` that satisfies the constraints. This output is 
also sorted by MZ just like `precursorLit`

### Notes

### Examples 

"""
function selectTransitionsPRM(window_center::Float32, precursorList::Vector{Precursor}, params)

    function getPrecursors(window_center::Float32, precursorList::Vector{Precursor}, params)
        l_bnd, u_bnd = window_center - params[:lower_tol], window_center + params[:upper_tol]
        start, stop = searchsortedfirst(precursorList, l_bnd,lt=(t,x)->getMZ(t)<x), searchsortedlast(precursorList, u_bnd,lt=(x,t)->getMZ(t)>x)
        return @view(precursorList[start:stop])
    end

    transitions = Vector{Transition}();

    for precursor in getPrecursors(window_center, precursorList, params)
        for charge in params[:transition_charges], isotope in params[:isotopes]
            append!(transitions, getTransitions(precursor, 
                                                charge = charge, 
                                                isotope = isotope, 
                                                b_start = params[:b_start], 
                                                y_start = params[:y_start],
                                                ppm = params[:fragment_ppm]))
        end
    end

    sort!(transitions, by=x->getMZ(x))

    transitions
end

"""
    getPrecursors(fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, n::Int, f_path::String, charges::Vector{UInt8}, isotopes::Vector{UInt8},mass_mods::Dict{String, Float32})

Build a `PrecursorTable` given a file_path to a tab delimited table of protein_name peptide_sequence pairs. Applies fixed and variable modifications. Gets precursors
for each combination of `charges` and `isotopes` supplied. 

### Input

- `fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}` -- Specification of fixed modifications to apply
- `var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}` -- Specification of variable modifications to apply
- `n::Int` -- Will apply all combinations of `n` or fewer variable modifications to each sequence
- `f_path::String` -- File path fo the protein-peptide list. 
- `charges::Vector{UInt8}` -- For each peptide, gets precursors with these charge states 
- `isotopes::Vector{UInt8}` -- For each peptide, gets precursors with these isotopic states
- ` mass_mods::Dict{String, Float32}` -- Specifies mass for each modification in `fixed_mods` and `var_mods` 

### Output
- Returns a `PrecursorTable` struct. See documentation for `PrecursorTable`. Mapps identifiers between precursors,
pepties, peptide groups, and proteins. 

### Notes

### Examples 

"""
function getPrecursors(fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
        var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
        n::Int, 
        f_path::String, charges::Vector{UInt8}, isotopes::Vector{UInt8},
        mass_mods::Dict{String, Float32})
    ptable = PrecursorTable()
    buildPrecursorTable!(ptable, fixed_mods, var_mods, n, f_path)
    addPrecursors!(ptable, charges, isotopes, mass_mods)
    return ptable
end

testPtable = PrecursorTable()
fixed_mods = [(p=r"C", r="C[Carb]"),
(p=r"K$", r="K[Hlys]"),
(p=r"R$", r="R[Harg]")]
#var_mods = [(p=r"(K$)", r="[Hlys]"), (p=r"(R$)", r="[Harg]")]
var_mods = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}()

test_mods::Dict{String, Float32} = 
Dict{String, Float32}(
    "Carb" => Float32(57.021464),
    "Harg" => Float32(10.008269),
    "Hlys" => Float32(8.014199),
)
using Tables
using DataFrames

"""
getBestPSMs(PSMs::Dict{Symbol, Vector}, ptable::PrecursorTable, MS_TABLE::Arrow.Table, min_fragment_count::UInt8)

Converts a dictionary of peptide spectrum matches (PSMs) into a dataframe, adds columns for the peptide id, retention time, sequence, and protein name. Also
retains only the row for each precursor with the highest XTandem hyperscore after filtering for a minimum ion count. 

### Input

- `PSMs::Dict{Symbol, Vector}` -- See `makePSMsDict(::FastXTandem)`. All types that inherit form `PSM` should implement a makePSMsDict. See notes. 
- `ptable::PrecursorTable` -- Mappings between precursor, peptide, peptide group, and protein identifiers. 
- `MS_TABLE::Arrow.Table` -- The RAW mass spec data in memory representation
- `min_fragment_count::UInt8` -- Dict mapping modifications to their masses

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

"""
function getBestPSMs(PSMs::Dict{Symbol, Vector}, ptable::PrecursorTable, MS_TABLE::Arrow.Table, min_fragment_count::UInt8)
    PSMs = DataFrame(PSMs)
    filter!(row -> row.total_ions >= min_fragment_count, PSMs);
    transform!(PSMs, AsTable(:) => ByRow(psm -> getPepIDFromPrecID(ptable, psm[:precursor_idx])) => :pep_idx)
    PSMs = combine(sdf -> sdf[argmax(sdf.hyperscore), :], groupby(PSMs, :pep_idx))
    transform!(PSMs, AsTable(:) => ByRow(psm -> MS_TABLE[:retentionTime][psm[:scan_idx]]) => :retentionTime)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getSeq(ptable.id_to_pep[psm[:pep_idx]])) => :sequence)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getMZ(getSimplePrecursors(ptable)[psm[:precursor_idx]])) => :precursor_mz)
    transform!(PSMs, AsTable(:) => ByRow(psm -> join(sort(collect(getProtNamesFromPepSeq(ptable, psm[:sequence]))), "|")) => :proteinNames)
    sort!(PSMs, [:retentionTime])
    PSMs
end

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
struct PrecursorChromatogram
    rts::Vector{Float32}
    last_scan_idx::Vector{Int64}
    transitions::UnorderedDictionary{String, Vector{Float32}}
    best_psm::NamedTuple{(:rt, :scan_idx, :name, :mz, :intensity), Tuple{Float32, Int64, Vector{String}, Vector{Float32}, Vector{Float32}}}
end

PrecursorChromatogram() = PrecursorChromatogram(Float32[], Int64[], UnorderedDictionary{String, Vector{Float32}}(), (rt = Float32(0), scan_idx = 0, name = String[], mz = Float32[], intensity = Float32[]))
getRTs(pc::PrecursorChromatogram) = pc.rts
getLastScanIdx(pc::PrecursorChromatogram) = pc.last_scan_idx
getTransitions(pc::PrecursorChromatogram) = pc.transitions
getBestPSM(pc::PrecursorChromatogram) = pc.best_psm

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
function initPrecursorChromatograms(bestPSMs::DataFrame)
    precursor_chromatograms = UnorderedDictionary{UInt32, PrecursorChromatogram}()
    for pep in eachrow(bestPSMs)
        insert!(precursor_chromatograms, 
                pep[:precursor_idx],
                PrecursorChromatogram(Float32[],
                [0],
                UnorderedDictionary{String, Vector{Float32}}(),
                (rt = pep[:retentionTime], scan_idx = pep[:scan_idx], name = Vector{String}(), mz = Float32[], intensity = Float32[])))

    end
    precursor_chromatograms
end

NRF2_Survey = Arrow.Table("./data/parquet/Nrf2_SQ_052122_survey.arrow")

testPtable = PrecursorTable()
buildPrecursorTable!(testPtable, fixed_mods, var_mods, 2, "./data/NRF2_SIL.txt")
addPrecursors!(testPtable, UInt8[2, 3, 4], UInt8[0], test_mods)
jldopen("somedata.jld", "w") do file
    addrequire(file, getPrecursors)
    write(file, "x", testPtable)
end
params = (upper_tol = 0.5, lower_tol = 0.5, δs = [0.0], 
          b_start = 3, y_start = 3, 
          prec_charges = UInt8[2,3,4], isotopes = UInt8[0], 
          transition_charges = UInt8[1,2], fragment_ppm = Float32(40),
          min_fragment_count = UInt8(5),
          rtWindow = 0.3)
          
scored, matches =  SearchRAW(NRF2_Survey, testPtable.precursors, selectTransitionsPRM, params)
bestPSMs = getBestPSMs(scored, testPtable, NRF2_Survey, UInt8(5))
precursor_chromatograms = initPrecursorChromatograms(bestPSMs)
fillPrecursorChromatograms!(precursor_chromatograms , matches, NRF2_Survey, 0.15)

testPtable = PrecursorTable()
buildPrecursorTable!(testPtable, fixed_mods, var_mods, 2, "./data/NRF2_SIL.txt")
addPrecursors!(testPtable, UInt8[2, 3, 4], UInt8[0], test_mods)

function ScoreSurveyRun(ptable::PrecursorTable, peptides_path::String, raw_path::String, params)
    MS_TABLE = Arrow.Table(raw_path)
    scored_psms, fragment_matches = SearchRAW(MS_TABLE, getPrecursors(ptable), selectTransitionsPRM, params)
    #Gets the dataframe
    best_psms = FilterPSMs(scored_psms, ptable, MS_TABLE, params[:min_fragment_count])

    ms1_peak_heights = UnorderedDictionary(best_psms[!,:precursor_idx], zeros(Float32, len(best)psms[!,:precursor_idx]))
    
    getMS1PeakHeights!(MS_TABLE[:retentionTime], 
                      MS_TABLE[:masses], 
                      MS_TABLE[:intensities], 
                      MS_TABLE[:msOrder], 
                      MS1_PEAK_HEIGHTS, 
                      best_psms[!,:retentionTime], 
                      best_psms[!,:precursor_idx], 
                      getSimplePrecursors(ptable), 
                      Float32(0.25), 
                      Float32(0.001), 
                      Float32(0.001))
                      
    transform!(bestPSMs, AsTable(:) => ByRow(psm -> ms1_peak_heights[psm[:precursor_idx]]) => :ms1_peak_height)

    precursor_chromatograms = initPrecursorChromatograms(best_psms) |> (best_psms -> fillPrecursorChromatograms(best_psms, fragment_matches, MS_TABLE, params[:rtWindow]))     

    transform!(bestPSMs, AsTable(:) => ByRow(psm -> getBestTransitions(getBestPSM(precursor_chromatograms[psm[:precursor_idx]]))) => :best_transitions)
    transform!(bestPSMs, AsTable(:) => ByRow(psm -> getBestPSM(precursor_chromatograms[psm[:precursor_idx]])[:name][psm[:best_transitions]]) => :transition_names)
    transform!(bestPSMs, AsTable(:) => ByRow(psm -> getBestPSM(precursor_chromatograms[psm[:precursor_idx]])[:mz][psm[:best_transitions]]) => :transition_mzs)
    #select(bestPSMs, [:proteinNames, :sequence,:precursor_mz,:names]) |> CSV.write("./data/method.csv", delim = "\t")
    writeTransitionList(best_psms, joinpath(raw_path, "transition_list.csv"))
    writeIAPIMethod(best_psms, joinpath(raw_path, "iapi_method.csv"))
    #Need to get MS1 heights here. 
    return precursor_chromatograms
end

# open the file for writing
function writeTransitionList(best_psms::DataFrame, f_out::String)
    open(f_out, "w") do io
        # loop over data and write each line
        for row in eachrow(best_psms)

            #data = join(append!([row[:proteinNames]*","*row[:sequence]], row[:names]),",")
            data = append!([row[:proteinNames], 
                            row[:sequence]], #replace(row[:sequence], r"\[(.*?)\]" => "")
                            row[:names])
            write(io, join(data,",")*"\n")
        end
    end
end

function writeIAPIMethod(best_psms::DataFrame, f_out::String)
    open(f_out, "w") do io
        # loop over data and write each line
        write(io, join(["protein_name","sequence","precursor_mz","precursor_intensity","transition_mz"],",")*"\n")
        for row in eachrow(best_psms)

            #data = join(append!([row[:proteinNames]*","*row[:sequence]], row[:names]),",")
            data = append!([row[:proteinNames], row[:sequence], row[:precursor_mz], row[:MS1_PEAK_HEIGHT]], row[:transition_mzs])
            write(io, join(data,",")*"\n")
        end
    end
end

MS1_MAX_HEIGHTS = UnorderedDictionary{UInt32, Float32}()

function getPrecursorsInRTWindow(retentionTimes::Vector{Float32}, l_bnd::Float32, u_bnd::Float32)
    #print("l_bnd ", l_bnd)
    #print("u_bnd ", u_bnd)
    start = searchsortedfirst(retentionTimes, l_bnd ,lt=(t,x)->t<x)
    stop = searchsortedlast(retentionTimes, u_bnd,lt=(x, t)->t>x)
    #println("test inner")
    #println("start1 ", start)
    #println("stop1 ", stop)
    return start, stop
end

searchsortedfirst(bestPSMs[!,:retentionTime], 50 ,lt=(t,x)->t<x)
searchsortedlast(bestPSMs[!,:retentionTime], 50 ,lt=(x, t)->t>x)
MS1_MAX_HEIGHTS = UnorderedDictionary{UInt32, Float32}()
function getMS1Peaks!(MS1::Vector{Union{Missing, Float32}}, INTENSITIES::Vector{Union{Missing, Float32}}, MS1_MAX_HEIGHTS::UnorderedDictionary{UInt32, Float32}, precursor_rts::Vector{Float32}, precursor_idxs::Vector{UInt32}, precursors::UnorderedDictionary{UInt32, SimplePrecursor}, rt::Float32, rt_tol::Float32, left_mz_tol::Float32, right_mz_tol::Float32)
    
    #Get precursors for which the best scan RT is within `rt_tol` of the current scan `rt`
    #precursor_idxs =  @view(bestPSMs[!,:precursor_idx][getPrecursorsInRTWindow(bestPSMs[!,:retentionTime],rt - rt_tol, rt + rt_tol)])
    #println("test")
    start, stop = getPrecursorsInRTWindow(precursor_rts,rt - rt_tol, rt + rt_tol)
    #println("start ", start)
    #println("stop ", stop)
    #stop = 1
    #start = 1
    if (stop-start) >= 0
    #Check to see if the MS1 height for each precursor is greater than the maximum previously observed. 
        for (i, psm_idx) in enumerate(start:stop)

            precursor_idx = precursor_idxs[psm_idx]
            if precursor_idx == UInt32(568)
            end
            if !isassigned(MS1_MAX_HEIGHTS, precursor_idx) #If this precursor has not been encountered before. 
                insert!(MS1_MAX_HEIGHTS, precursor_idx, Float32(0))
            end

            mz = getMZ(precursors[precursor_idx]) #Get the precursor MZ

            idx = binaryGetNearest(MS1, mz, mz-left_mz_tol, mz+right_mz_tol) #Get the peak index of the peak nearest in mz to the precursor. 

            if idx == 0
                continue
            end

            if coalesce(INTENSITIES[idx], Float32(0))>=MS1_MAX_HEIGHTS[precursor_idx] #Replace maximum observed MS1 height for the precursor if appropriate. 
                MS1_MAX_HEIGHTS[precursor_idx] = coalesce(INTENSITIES[idx], Float32(0))
            end
        end
    end
end

function getMS1PeakHeights!(retentionTimes::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}, 
                            masses::Arrow.List{Union{Missing, Vector{Union{Missing, Float32}}}, Int32, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}},
                            intensities::Arrow.List{Union{Missing, Vector{Union{Missing, Float32}}}, Int32, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}},
                            msOrders::Arrow.Primitive{Union{Missing, Int32}, Vector{Int32}},
                            ms1_max_heights::UnorderedDictionary{UInt32, Float32}, 
                            precursor_rts::Vector{Float32}, precursor_idxs::Vector{UInt32}, 
                            precursorsDict::UnorderedDictionary{UInt32, SimplePrecursor}, 
                            rt_tol::Float32, left_mz_tol::Float32, right_mz_tol::Float32)
    #println("tunction")
    #i = 1
    
    for scan_idx in eachindex(retentionTimes)
        if msOrders[scan_idx]!=Int32(1) #Skip non MS1 scans. 
            continue
        end
        #i += 1
        #println("test?")
        @inline getMS1Peaks!(masses[scan_idx], 
                    intensities[scan_idx], 
                    ms1_max_heights, 
                    precursor_rts, precursor_idxs,                          #Precursor retention times and id's (sorted by rt)
                    precursorsDict,                                         #PrecursorTable
                    retentionTimes[scan_idx],                               #RT of current scan
                    rt_tol, left_mz_tol, right_mz_tol #Only consider precursors where the best scan is within the `rt_col` of the current scan
                    );
    end
    #println("i ", i)
end

#need to get RT estimate for each peptide. 
function fillPrecursorChromatograms!(precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram}, fragment_matches::Vector{FragmentMatch}, MS_TABLE::Arrow.Table, rt_tol::Float64)

    """
        getTransitionName(transition::UnorderedDictionary{String, Vector{Float32}}) 
    """
    function getTransitionName(transition::Transition)
        string(transition.ion_type)*string(transition.ind)*"+"*string(transition.charge)
    end

    function getBestPSMs!(precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram}, fragment_matches::Vector{FragmentMatch})

        function addScanToPrecursor!(precursor_chromatograms::UnorderedDictionary{UInt32, PrecursorChromatogram}, match::FragmentMatch, fragment_name::String, prec_id::UInt32)
            
            function assignPrecursor!(matched_precursor::PrecursorChromatogram, fragment_name::String)
                if !isassigned(matched_precursor.transitions, fragment_name)
                    insert!(matched_precursor.transitions, fragment_name, Float32[])
                end
            end

            assignPrecursor!(precursor_chromatograms[prec_id], fragment_name);

            if match.scan_idx == precursor_chromatograms[prec_id].best_psm[:scan_idx]
                push!(precursor_chromatograms[prec_id].best_psm[:mz], match.match_mz)
                push!(precursor_chromatograms[prec_id].best_psm[:intensity], match.intensity)
                push!(precursor_chromatograms[prec_id].best_psm[:name], fragment_name)
            end

        end

        for match in fragment_matches
            #pep_id = getPepIDFromPrecID(testPtable, match.transition.prec_id)
            prec_id = match.transition.prec_id
            if !isassigned(precursor_chromatograms, prec_id) continue end
            #as in y3+2 or b6+1 (name/length/+charge)
            fragment_name = getTransitionName(match.transition)
            addScanToPrecursor!(precursor_chromatograms, match, fragment_name, prec_id)
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

    getBestPSMs!(precursor_chromatograms, fragment_matches)
    addTransitionIntensities!(precursor_chromatograms, fragment_matches, MS_TABLE, rt_tol)
    return precursor_chromatograms 

end

getMatchedPrecursors(matched_precursors, matches, NRF2_Survey, 0.15)



function plotBestSpectra(matched_ions::NamedTuple{(:scan_idx, :name, :mz, :intensity), Tuple{Int64, Vector{String}, Vector{Float32}, Vector{Float32}}}, RAW::Arrow.Table, title::String, out_path::String)

    function plotSpectra!(p::Plots.Plot{Plots.GRBackend}, masses, intensities)
        for (peak_idx, mass) in enumerate(masses)
            plot!(p, [mass, mass], [0, intensities[peak_idx]], legend = false, color = "black")
        end
    end
    #p = plot()
    #plotSpectra(p, NRF2_Survey.masses[2703], NRF2_Survey.intensities[2703])
    #display(p)
    function addFragmentIons!(p::Plots.Plot{Plots.GRBackend}, matched_ions::NamedTuple{(:scan_idx, :name, :mz, :intensity), Tuple{Int64, Vector{String}, Vector{Float32}, Vector{Float32}}})
        i = 1
        for (mz, intensity) in zip(matched_ions[:mz], matched_ions[:intensity])
            plot!(p, [mz, mz], [0, -1*intensity], legend = false, color = "red")
            annotate!(mz + 1, -1*intensity, fontfamily = "helvetica", text(matched_ions[:name][i], :black, :top, :left, 8))
            i += 1
        end
    end

    p = plot(title = title, fontfamily="helvetica")
    plotSpectra!(p, RAW[:masses][matched_ions[:scan_idx]], RAW[:intensities][matched_ions[:scan_idx]])
    addFragmentIons!(p, matched_ions)
    savefig(joinpath(out_path,title*".pdf"))

end

function plotAllBestSpectra(matched_precursors::UnorderedDictionary{UInt32, PrecursorChromatogram}, ptable::PrecursorTable, MS_TABLE::Arrow.Table, out_path::String, fname::String)
    if !isdir(out_path)
        mkpath(out_path)
    end

    @time for key in keys(matched_precursors)
        plotBestSpectra(matched_precursors[key].best_psm, MS_TABLE,
                        ptable.id_to_pep[getPepIDFromPrecID(ptable, key)].sequence,
                        out_path)
    end

    files = filter(x -> isfile(joinpath(out_path, x)) && match(r"\.pdf$", x) != nothing, readdir(out_path))
    merge_pdfs(map(file -> joinpath(out_path,file), files), joinpath(out_path, fname), cleanup=true)

end

plotAllBestSpectra(matched_precursors, testPtable, NRF2_Survey, "./data/figures/spectra/", "MERGED_NRF2_SURVEY_SPECTRA.pdf")

function plotFragmentIonChromatogram(transitions::UnorderedDictionary{String, Vector{Float32}}, rts::Vector{Float32}, title::String, out_path::String)
    p = plot(title = title, fontfamily="helvetica")
    for (color, t) in enumerate(keys(transitions))
        plot!(p, rts, transitions[t], color = color, legend = true, label = t)
        plot!(p, rts, transitions[t], seriestype=:scatter, color = color, label = nothing)
    end
    savefig(joinpath(out_path, title*".pdf"))
end

#testPtable.id_to_pep[1].sequence
function plotAllFragmentIonChromatograms(matched_precursors::UnorderedDictionary{UInt32, PrecursorChromatogram}, ptable::PrecursorTable, out_path::String, fname::String)

    if !isdir(out_path)
        mkpath(out_path)
    end

    @time for key in keys(matched_precursors)
        plotFragmentIonChromatogram(matched_precursors[key].transitions, matched_precursors[key].rts, 
        ptable.id_to_pep[getPepIDFromPrecID(ptable, key)].sequence, out_path)
        # Get all files in the directory that match the pattern
    end
    files = filter(x -> isfile(joinpath(out_path, x)) && match(r"\.pdf$", x) != nothing, readdir(out_path))
    merge_pdfs(map(file -> joinpath(out_path,file), files), joinpath(out_path, fname), cleanup=true)
end
plotAllFragmentIonChromatograms(matched_precursors, testPtable, "./data/figures/", "MERGED_NRF2_SURVEY_CHROMATOGRAMS.pdf")

using Glob
getPepIDFromPrecIDlegend = false
using Glob, Regex
using PDFmerger
# Set the path of the directory and the regular expression pattern
dir_path = "./data/figures/"
pattern = r"\.pdf$"  # Match files that end with .txt

# Get all files in the directory that match the pattern
files = filter(x -> isfile(joinpath(dir_path, x)) && match(pattern, x) != nothing, readdir(dir_path))
@time merge_pdfs(map(file -> joinpath(dir_path,file), files), cleanup=true)
# Print the list of matching files
println(files)