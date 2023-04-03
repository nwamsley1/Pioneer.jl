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

function getPrecursors(fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
        var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
        n::Int, 
        f_path::String, charges::Vector{UInt8}, isotopes::Vector{UInt8},
        mass_mods::Dict{String, Float32})
    ptable = PrecursorTable()
    buildPrecursorTable(ptable, fixed_mods, var_mods, n, f_path)
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

function getPrecursors(fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
                       var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
                       n::Int, 
                       f_path::String, charges::Vector{UInt8}, isotopes::Vector{UInt8},
                       mass_mods::Dict{String, Float32})
    ptable = PrecursorTable()
    buildPrecursorTable(ptable, fixed_mods, var_mods, n, f_path)
    addPrecursors!(ptable, charges, isotopes, mass_mods)
    return ptable
end

testPtable = PrecursorTable()
buildPrecursorTable!(testPtable, fixed_mods, var_mods, 2, "./data/NRF2_SIL.txt")
addPrecursors!(testPtable, UInt8[2, 3, 4], UInt8[0], test_mods)

function ScorePrecursorsPRM(ptable::PrecursorTable, raw_path::String, params)
    MS_TABLE = Arrow.Table(raw_path)
    SearchRaw(MS_TABLE, ptable, selectTransitionsPRM, params)
end

function ProcessPSMs(PSMs::Dict{Symbol, Vector}, ptable::PrecursorTable, MS_TABLE::Arrow.Table, min_ion_count::UInt8)
    PSMs = DataFrame(PSMs)
    filter!(row -> row.total_ions >= min_ion_count, PSMs);
    transform!(PSMs, AsTable(:) => ByRow(psm -> getPepIDFromPrecID(ptable, psm[:precursor_idx])) => :pep_idx)
    PSMs = combine(sdf -> sdf[argmax(sdf.hyperscore), :], groupby(PSMs, :pep_idx))
    transform!(PSMs, AsTable(:) => ByRow(psm -> MS_TABLE[:retentionTime][psm[:scan_idx]]) => :retentionTime)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getSeq(testPtable.id_to_pep[psm[:pep_idx]])) => :sequence)
    transform!(PSMs, AsTable(:) => ByRow(psm -> join(sort(collect(getProtNamesFromPepSeq(ptable, psm[:sequence]))), "|")) => :proteinNames)
    PSMs
end

struct MatchedPrecursor
    rt::Float32
    rts::Vector{Float32}
    last_scan_idx::Vector{Int64}
    transitions::UnorderedDictionary{String, Vector{Float32}}
    best_psm::NamedTuple{(:scan_idx, :name, :mz, :intensity), Tuple{Int64, Vector{String}, Vector{Float32}, Vector{Float32}}}
end

MatchedPrecursor() = MatchedPrecursor(Float32(0), 0, Float32[], UnorderedDictionary{String, Vector{Float32}}(), (scan_idx = 0, name = String[], mz = Float32[], intensity = Float32[]))

function initMatchedPrecursors(bestPSMs::DataFrame)
    matched_precursors = UnorderedDictionary{UInt32, MatchedPrecursor}()
    for pep in eachrow(bestPSMs)
        insert!(matched_precursors, 
                pep[:precursor_idx], 
                MatchedPrecursor(pep[:retentionTime],
                Float32[],
                [0],
                UnorderedDictionary{String, Vector{Float32}}(),
                (scan_idx = pep[:scan_idx], name = Vector{String}(), mz = Float32[], intensity = Float32[])))

    end
    matched_precursors
end

NRF2_Survey = Arrow.Table("./data/parquet/Nrf2_SQ_052122_survey.arrow")
params = (upper_tol = 0.5, lower_tol = 0.5, δs = [0.0], b_start = 3, y_start = 3, prec_charges = UInt8[2,3,4], 
isotopes = UInt8[0], transition_charges = UInt8[1,2], fragment_ppm = Float32(40))
scored, matches =  SearchRAW(NRF2_Survey, testPtable.precursors, selectTransitionsPRM, params)
bestPSMs = ProcessPSMs(scored, testPtable, NRF2_Survey, UInt8(5))
matched_precursors = initMatchedPrecursors(bestPSMs)
getMatchedPrecursors(matched_precursors, matches, NRF2_Survey, 0.15)

#need to get RT estimate for each peptide. 
function getMatchedPrecursors(matched_precursors::UnorderedDictionary{UInt32, MatchedPrecursor}, fragment_matches::Vector{FragmentMatch}, MS_TABLE::Arrow.Table, rt_tol::Float64)

    """
        getTransitionName(transition::UnorderedDictionary{String, Vector{Float32}}) 
    """
    function getTransitionName(transition::Transition)
        string(transition.ion_type)*string(transition.ind)*"+"*string(transition.charge)
    end

    function getBestPSMs!(matched_precursors::UnorderedDictionary{UInt32, MatchedPrecursor}, fragment_matches::Vector{FragmentMatch})

        function addScanToPrecursor!(matched_precursors::UnorderedDictionary{UInt32, MatchedPrecursor}, match::FragmentMatch, fragment_name::String, prec_id::UInt32)
            
            function assignPrecursor!(matched_precursor::MatchedPrecursor, fragment_name::String)
                if !isassigned(matched_precursor.transitions, fragment_name)
                    insert!(matched_precursor.transitions, fragment_name, Float32[])
                end
            end

            assignPrecursor!(matched_precursors[prec_id], fragment_name);

            if match.scan_idx == matched_precursors[prec_id].best_psm[:scan_idx]
                push!(matched_precursors[prec_id].best_psm[:mz], match.match_mz)
                push!(matched_precursors[prec_id].best_psm[:intensity], match.intensity)
                push!(matched_precursors[prec_id].best_psm[:name], fragment_name)
            end

        end

        for match in fragment_matches
            #pep_id = getPepIDFromPrecID(testPtable, match.transition.prec_id)
            prec_id = match.transition.prec_id
            if !isassigned(matched_precursors, prec_id) continue end
            #as in y3+2 or b6+1 (name/length/+charge)
            fragment_name = getTransitionName(match.transition)
            addScanToPrecursor!(matched_precursors, match, fragment_name, prec_id)
        end

    end

    function addTransition(matched_precursor::MatchedPrecursor, match::FragmentMatch, rt::Float32)
        #println("TEST0")
        #println("match ", match)
        #println("matched_precursor", matched_precursor)
        #println("matched_precursor.last_scan_idx[1] ", matched_precursor.last_scan_idx[1])
        if match.scan_idx>matched_precursor.last_scan_idx[1]
            #println("TEST1")
            matched_precursor.last_scan_idx[1]=match.scan_idx
            append!(matched_precursor.rts, rt)
            for key in keys(matched_precursor.transitions)
                #println("TEST");
                push!(matched_precursor.transitions[key],Float32(0))
            end
        end
        matched_precursor.transitions[getTransitionName(match.transition)][end] = match.intensity
    end

    function isInRTWindow(matched_precursor::MatchedPrecursor, rt, rt_tol::Float64)
        best_rt = matched_precursor.rt
        return (abs(rt- best_rt)<rt_tol)
    end

    function addTransitionIntensities!(matched_precursors::UnorderedDictionary{UInt32, MatchedPrecursor}, matches::Vector{FragmentMatch}, MS_TABLE::Arrow.Table, rt_tol::Float64)
        for match in matches
            prec_id = match.transition.prec_id
            if !isassigned(matched_precursors, prec_id) continue end
            rt = MS_TABLE[:retentionTime][match.scan_idx]
            if isInRTWindow(matched_precursors[prec_id], rt, rt_tol) 
                addTransition(matched_precursors[prec_id], match, rt) 
                #println("TESTTEST")
            end
        end
    end

    getBestPSMs!(matched_precursors, fragment_matches)
    addTransitionIntensities!(matched_precursors, fragment_matches, MS_TABLE, rt_tol)
    return matched_precursors 

end

getMatchedPrecursors(matched_precursors, matches, NRF2_Survey, 0.15)

testmatchp(matches, test_dict)


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

function plotAllBestSpectra(matched_precursors::UnorderedDictionary{UInt32, MatchedPrecursor}, ptable::PrecursorTable, MS_TABLE::Arrow.Table, out_path::String, fname::String)
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


#testPtable.id_to_pep[1].sequence
function plotAllFragmentIonChromatograms(matched_precursors::UnorderedDictionary{UInt32, MatchedPrecursor}, ptable::PrecursorTable, out_path::String, fname::String)
    function plotFragmentIonChromatogram(transitions::UnorderedDictionary{String, Vector{Float32}}, rts::Vector{Float32}, title::String, out_path::String)
        p = plot(title = title, fontfamily="helvetica")
        for (color, t) in enumerate(keys(transitions))
            plot!(p, rts, transitions[t], color = color, legend = true, label = t)
            plot!(p, rts, transitions[t], seriestype=:scatter, color = color, label = nothing)
        end
        savefig(joinpath(out_path, title*".pdf"))
    end

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