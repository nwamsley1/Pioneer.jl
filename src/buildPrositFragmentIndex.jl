abstract type FragmentIndexType end

getFragMZ(f::FragmentIndexType) = f.frag_mz
getPrecID(f::FragmentIndexType) = f.prec_id
getPrecCharge(f::FragmentIndexType) = f.prec_charge

import Base.<
import Base.>

<(y::FragmentIndexType, x::T) where {T<:Real} = getFragMZ(y) < x
<(x::T, y::FragmentIndexType) where {T<:Real} = <(y, x)
>(y::FragmentIndexType, x::T) where {T<:Real} = getFragMZ(y) > x
>(x::T, y::FragmentIndexType) where {T<:Real} = >(y, x)

struct FragmentIon{T<:AbstractFloat} <: FragmentIndexType
    frag_mz::T
    prec_id::UInt32
    prec_mz::T
    prec_charge::UInt8
end

getPrecMZ(f::FragmentIon) = f.prec_mz

struct LibraryFragment{T<:AbstractFloat} <: FragmentIndexType
    frag_mz::T
    frag_charge::UInt8
    is_y_ion::Bool
    ion_position::UInt8
    ion_index::UInt8
    intensity::Float32
    prec_charge::UInt8
    prec_id::UInt32
end

getIntensity(f::LibraryFragment) = f.intensity
isyIon(f::LibraryFragment) = f.is_y_ion
getIonIndex(f::LibraryFragment) = f.ion_index
getIonPosition(f::LibraryFragment) = f.ion_position
getFragCharge(f::LibraryFragment) = f.frag_charge

function buildFragmentIndex!(frag_ions::Vector{FragmentIon{T}}, bin_ppm::AbstractFloat; low_frag_mz::AbstractFloat = 150.0, high_frag_mz::AbstractFloat = 1700.0, low_prec_mz::AbstractFloat = 300.0, high_prec_mz::AbstractFloat = 1100.0) where {T<:AbstractFloat}
   
    #The fragment ions are divided into bins of roughtly equal m/z width.
    #That should correspond to roughly half the fragment mass accuracy of the detector?
    frag_index = FragmentIndex(T) 

    function fillPrecursorBin!(frag_index::FragmentIndex{<:AbstractFloat}, frag_ions::Vector{FragmentIon{T}}, frag_bin_idx::UInt32, start::Int, stop::Int, low_prec_mz::AbstractFloat, high_prec_mz::AbstractFloat)
        for ion_index in range(start, stop)
            pep_id = getPrecID(frag_ions[ion_index])
                prec_mz = getPrecMZ(frag_ions[ion_index])#(getPrecMZ(frag_ions[ion_index]) + PROTON*(charge-1))/charge #m/z of the precursor

                if (prec_mz < low_prec_mz) | (prec_mz > high_prec_mz) #Precursor m/z outside the bounds
                    continue
                end
                #Add precursor corresponding to the charge state
                addPrecursorBinItem!(frag_index,
                                     frag_bin_idx,
                                    PrecursorBinItem(pep_id, prec_mz, getPrecCharge(frag_ions[ion_index]))
                                    )
        end
    end

    bin = UInt32(1) #Current fragment bin index
    start = 1 #Fragment index of the first fragment in the current bin

    getPPM(frag_mz::T, ppm::T) = ppm*frag_mz/1e6

    diff = getPPM(getFragMZ(frag_ions[start]), bin_ppm) #ppm tolerance of the current fragment bin

    #Build bins 
    for stop in 2:length(frag_ions)

        #Haven't reached minimum fragment m/z yet
        if getFragMZ(frag_ions[stop]) < low_frag_mz
            start += 1
            continue
        end

        #Doex the difference between the highest and lowest m/z in the bin 
        #enought exceed 10 ppm of the lowest m/z?
        if (getFragMZ(frag_ions[stop]) - getFragMZ(frag_ions[start])) > diff

            #Nedds to be stop - 1 to gaurantee the smallest and largest fragment
            #in the bin differ by less than diff 
            last_frag_in_bin = stop - 1

            #Add a new fragment bin
            addFragmentBin!(frag_index, 
                            FragBin(getFragMZ(frag_ions[start]),
                                    getFragMZ(frag_ions[last_frag_in_bin]),
                                    bin
                                    )
                            )

            addPrecursorBin!(frag_index, 
                                #Preallocate an empty precursor bin of the correct length 
                                PrecursorBin(Vector{PrecursorBinItem{T}}())#undef, (last_frag_in_bin - start + 1)*length(charges)))
                                )
        
            fillPrecursorBin!(frag_index, frag_ions, bin, start, last_frag_in_bin, low_prec_mz, high_prec_mz)

            #Sort the precursor bin by precursor m/z
            sort!(getPrecursors(getPrecursorBin(frag_index, bin)), by = x->getPrecMZ(x));

            #Update counters and ppm tolerance 
            bin += UInt32(1)
            start = stop
            diff = getPPM(getFragMZ(frag_ions[start]), bin_ppm)

            #Maximum fragment m/z reached. Stop adding bins. 
            if getFragMZ(frag_ions[stop]) > high_frag_mz
                break
            end
        end
    end
    return frag_index
end

mutable struct FragmentMatch{T<:AbstractFloat}
    predicted_intensity::T
    intensity::T
    theoretical_mz::T
    match_mz::T
    peak_ind::Int64
    frag_index::UInt8
    frag_charge::UInt8
    frag_isotope::UInt8
    ion_type::Char
    prec_id::UInt32
    count::UInt8
    scan_idx::UInt32
    ms_file_idx::UInt32
end

FragmentMatch() = FragmentMatch(Float64(0), Float64(0), Float64(0), Float64(0), 0, UInt8(0), UInt8(0), UInt8(0),'y', UInt32(0), UInt8(0), UInt32(0), UInt32(0))
getFragMZ(f::FragmentMatch) = f.theoretical_mz
getMatchMZ(f::FragmentMatch) = f.match_mz
getPredictedIntenisty(f::FragmentMatch) = f.predicted_intensity
getIntensity(f::FragmentMatch) = f.intensity

getPeakInd(f::FragmentMatch) = f.peak_ind
getFragInd(f::FragmentMatch) = f.frag_index
getCharge(f::FragmentMatch) = f.frag_charge
getIsotope(f::FragmentMatch) = f.frag_isotope
getIonType(f::FragmentMatch) = f.ion_type

getPrecID(f::FragmentMatch) = f.prec_id
getCount(f::FragmentMatch) = f.count
getScanID(f::FragmentMatch) = f.scan_idx
getMSFileID(f::FragmentMatch) = f.ms_file_idx

function findFirstFragmentBin(frag_index::Vector{FragBin{T}}, frag_min::AbstractFloat, frag_max::AbstractFloat) where {T<:AbstractFloat}
    #Binary Search
    lo, hi = 1, length(frag_index)
    potential_match = nothing
    while lo <= hi

        mid = (lo + hi) ÷ 2

        if (frag_min) <= getHighMZ(frag_index[mid])
            if (frag_max) >= getHighMZ(frag_index[mid]) #Frag tolerance overlaps the upper boundary of the frag bin
                potential_match = mid
            end
            hi = mid - 1
        elseif (frag_max) >= getLowMZ(frag_index[mid]) #Frag tolerance overlaps the lower boundary of the frag bin
            if (frag_min) <= getLowMZ(frag_index[mid])
                potential_match = mid
                #return mid
            end
            lo = mid + 1
        end
    end

    return potential_match#, Int64(getPrecBinID(frag_index[potential_match]))
end

function searchPrecursorBin!(precs::Dictionary{UInt32, UInt8}, precursor_bin::PrecursorBin{T}, window_min::AbstractFloat, window_max::AbstractFloat) where {T<:AbstractFloat}
   
    N = getLength(precursor_bin)
    lo, hi = 1, N

    while lo <= hi
        mid = (lo + hi) ÷ 2
        if getPrecMZ(getPrecursor(precursor_bin, mid)) < window_min
            lo = mid + 1
        else
            hi = mid - 1
        end
    end

    window_start = (lo <= N ? lo : return nothing, nothing)

    if getPrecMZ(getPrecursor(precursor_bin, window_start)) > window_max
        return nothing, nothing
    end

    lo, hi = window_start, N

    while lo <= hi
        mid = (lo + hi) ÷ 2
        if getPrecMZ(getPrecursor(precursor_bin, mid)) > window_max
            hi = mid - 1
        else
            lo = mid + 1
        end
    end

    window_stop = hi

    function addFragmentMatches!(precs::Dictionary{UInt32, UInt8}, precursor_bin::PrecursorBin{T}, window_min::AbstractFloat, window_max::AbstractFloat, start::Int, stop::Int) where {T<:AbstractFloat}
        for precursor_idx in start:stop
            
            precursor = getPrecursor(precursor_bin, precursor_idx)
            prec_mz = getPrecMZ(precursor)
            prec_id = getPrecID(precursor)
            charge = getPrecCharge(precursor)

            prec_mz_min = window_min - 3.0*NEUTRON/charge#upper_tol[charge]
            prec_mz_max = window_max + 1.0*NEUTRON/charge#lower_tol[charge]
            if (prec_mz_min <= prec_mz) & (prec_mz_max >= prec_mz)
                if haskey(precs, prec_id)
                    precs[prec_id] += UInt8(1)
                else
                    insert!(precs, prec_id, UInt8(1))
                end
            end
        end

    end

    addFragmentMatches!(precs, precursor_bin, window_min, window_max, window_start, window_stop)

    return window_start, window_stop

end

#const upper_tol = [(3.0*NEUTRON), (3.0*NEUTRON)/2, (3.0*NEUTRON)/3, (3.0*NEUTRON)/4]
#const lower_tol = [(1*NEUTRON), (1*NEUTRON)/2, (1*NEUTRON)/3, (1*NEUTRON)/4]

function queryFragment!(precs::Dictionary{UInt32, UInt8}, frag_index::FragmentIndex{T}, min_frag_bin::Int64, frag_min::AbstractFloat, frag_max::AbstractFloat, prec_mz::AbstractFloat, prec_tol::AbstractFloat) where {T<:AbstractFloat}
    
    frag_bin = findFirstFragmentBin(getFragBins(frag_index), frag_min, frag_max)
    #No fragment bins contain the fragment m/z
    if (frag_bin === nothing)
        return min_frag_bin
    #This frag bin has already been searched
    elseif frag_bin <= min_frag_bin
        return min_frag_bin
    end

    i = 1
    while (frag_bin < length(getFragBins(frag_index))) #getLowMZ(getFragmentBin(frag_index, frag_bin)) <frag_max
        #Fragment bin matches the fragment ion
        #println(i)
        i += 1
        if (getLowMZ(getFragmentBin(frag_index, frag_bin)) > frag_max)
            return frag_bin
        else
            prec_min = prec_mz - prec_tol - 3.0*NEUTRON
            prec_max = prec_mz + prec_tol + 1.0*NEUTRON
            searchPrecursorBin!(precs, getPrecursorBin(frag_index, UInt32(frag_bin)), prec_min, prec_max)
            frag_bin += 1
        end

    end

    #Only reach this point if frag_bin exceeds length(frag_index)
    return frag_bin - 1
end

function searchScan!(precs::Dictionary{UInt32, UInt8}, f_index::FragmentIndex{T}, massess::Vector{Union{Missing, U}}, intensities::Vector{Union{Missing, U}}, precursor_window::AbstractFloat, ppm::AbstractFloat, width::AbstractFloat; topN::Int = 20, min_frag_count::Int = 3) where {T,U<:AbstractFloat}
    
    getFragTol(mass::U, ppm::AbstractFloat) = mass*(1 - ppm/1e6), mass*(1 + ppm/1e6)

    function filterPrecursorMatches!(precs::Dictionary{UInt32, UInt8}, topN::Int, min_frag_count::Int) where {T<:AbstractFloat}
        #Do not consider peptides wither fewer than 
        match_count = sum(precs)
        prec_count = length(precs)
        filter!(count->(count>=min_frag_count), precs)
        sort!(precs, rev = true)
        #println(precs)
        #Iterator of Peptide ID's for the `topN` scoring peptides
        return Iterators.take(keys(precs), min(topN, length(keys(precs)))), prec_count, match_count
    end

    min_frag_bin = 0

    for (mass, intensity) in zip(massess, intensities)

        mass, intensity = coalesce(mass, 0.0),  coalesce(intensity, 0.0)

        FRAGMIN, FRAGMAX = getFragTol(mass, ppm) 

        min_frag_bin = queryFragment!(precs, f_index, min_frag_bin, FRAGMIN, FRAGMAX, precursor_window, width)
    end 
    #println("PRECS $precs")
    return filterPrecursorMatches!(precs, topN, min_frag_count)
end

function selectTransitions(fragment_list::Vector{Vector{LibraryFragment{T}}}, pep_ids::Base.Iterators.Take{Indices{UInt32}}, ppm::AbstractFloat = 20.0) where {T<:AbstractFloat}
    transitions = Vector{LibraryFragment{T}}()
    for pep_id in pep_ids
        append!(transitions, fragment_list[pep_id])
    end
    return sort!(transitions, by = x->getFragMZ(x))
end

function SearchRAW(
                    spectra::Arrow.Table, 
                    #ptable::PrecursorDatabase,
                    frag_index::FragmentIndex{T},
                    fragment_list::Vector{Vector{LibraryFragment{Float64}}},
                    ms_file_idx::UInt32;
                    precursor_tolerance::Float64 = 0.5,
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
    test_frags = ""
    test_matches = ""
    test_misses = ""
    for (i, spectrum) in enumerate(Tables.namedtupleiterator(spectra))
        if spectrum[:msOrder] != 2
            continue
        end
        ms2 += 1
        fragmentMatches = Vector{FragmentMatch{Float32}}()
        precs = Dictionary{UInt32, UInt8}()
        pep_id_iterator, prec_count, match_count = searchScan!(precs, 
                    frag_index, 
                    spectrum[:masses], spectrum[:intensities], spectrum[:precursorMZ], 
                    fragment_tolerance, 
                    precursor_tolerance,
                    min_frag_count = min_frag_count, 
                    topN = topN
                    )

        transitions = selectTransitions(fragment_list, pep_id_iterator)

        fragmentMatches, fragmentMisses = matchPeaks(transitions, 
                                    spectrum[:masses], 
                                    spectrum[:intensities], 
                                    #δs = params[:δs],
                                    δs = zeros(T, (1,)),#[Float64(0)],
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = min_intensity
                                    )

        X, H, IDtoROW = buildDesignMatrix(fragmentMatches, fragmentMisses, topN)

        #Does this always coincide with there being zero fragmentMatches?
        #Could change to length(fragmentMatches) == 0 ?
        if size(H)[2] == 0
            continue
        end

        #Initialize weights for each precursor template. 
        #Should find a more sophisticated way of doing this. 
        W = reshape([Float32(1000) for x in range(1,size(H)[1])], (1, size(H)[1]))

        #Solve NMF. 
        weights = NMF.solve!(NMF.GreedyCD{Float32}(maxiter=50, verbose = false, 
                                                    lambda_w = 1e4, 
                                                    tol = 1e-6, #Need a reasonable way to choos lambda?
                                                    update_H = false #Important to keep H constant. 
                                                    ), X, W, H).W[1,:]

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
    end
    println("processed $ms2 scans!")
    return DataFrame(scored_PSMs)# test_frags, test_matches, test_misses#DataFrame(scored_PSMs)
end

@time PSMs = SearchRAW(MS_TABLE, prosit_index, prosit_list_detailed, UInt32(1))
@time PSMs = SearchRAW(MS_TABLE, prosit_index, UInt32(1))

@time test_frags, test_matches, test_misses = SearchRAW(MS_TABLE, prosit_index, prosit_list_detailed, UInt32(1))

struct LibraryPrecursor{T<:AbstractFloat}
    iRT::T
    isDecoy::Bool
    charge::UInt8
end

prosit_simple_targets, prosit_detailed_targets, prosit_precs_targets, prec_id = readPrositLib("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_targets.csv")
prosit_simple_decoys, prosit_detailed_decoys, prosit_precs_decoys, prec_id = readPrositLib("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_decoys.csv", isDecoys = true, first_prec_id = prec_id)


@save "/Users/n.t.wamsley/Projects/prosit_simple_targets.jld2"  prosit_simple_targets
@save "/Users/n.t.wamsley/Projects/prosit_detailed_targets.jld2"  prosit_detailed_targets
@save "/Users/n.t.wamsley/Projects/prosit_precs_targets.jld2"  prosit_precs_targets

@save "/Users/n.t.wamsley/Projects/prosit_simple_decoys.jld2"  prosit_simple_decoys
@save "/Users/n.t.wamsley/Projects/prosit_detailed_decoys.jld2"  prosit_detailed_decoys
@save "/Users/n.t.wamsley/Projects/prosit_precs_decoys.jld2"  prosit_precs_decoys


function readPrositLib(prosit_lib_path::String; precision::DataType = Float64, isDecoys::Bool = false, first_prec_id = UInt32(0))
    frag_list = Vector{FragmentIon{precision}}()
    frag_detailed = Vector{Vector{LibraryFragment{precision}}}()
    precursor_list = Vector{LibraryPrecursor}()

    rows = CSV.Rows(prosit_lib_path, reusebuffer=false, select = [:RelativeIntensity, :FragmentMz, :PrecursorMz, :iRT, :Stripped, :ModifiedPeptide,:FragmentNumber,:FragmentCharge,:PrecursorCharge,:FragmentType])
    current_peptide = ""
    current_charge = ""
    prec_id = first_prec_id
    id = UInt32(0)
    ion_position = UInt8(1)
    for (i, row) in enumerate(rows)
        if (row.ModifiedPeptide::PosLenString != current_peptide) | (row.PrecursorCharge::PosLenString != current_charge)
            current_peptide = row.ModifiedPeptide::PosLenString
            current_charge = row.PrecursorCharge::PosLenString
            prec_id += UInt32(1)
            id = UInt32(1)
            ion_position = UInt8(1)
            push!(frag_detailed, Vector{LibraryFragment{precision}}())
            push!(precursor_list, LibraryPrecursor(
                                                    parse(precision, row.iRT::PosLenString),
                                                    isDecoys,
                                                    parse(UInt8, row.PrecursorCharge::PosLenString)
                                                ))
        end

        #Track progress
        if (i % 1_000_000) == 0
            println(i/1_000_000)
        end

        #Exclude y1, y2, b1, and b2 ions. 
        if parse(Int, row.FragmentNumber::PosLenString) < 3
            continue
        end

        push!(frag_list, FragmentIon(parse(precision, row.FragmentMz::PosLenString), 
                                    prec_id, 
                                    parse(precision, row.PrecursorMz::PosLenString), 
                                    parse(UInt8, row.PrecursorCharge::PosLenString)))

        push!(frag_detailed[id], LibraryFragment(parse(precision, row.FragmentMz::PosLenString), 
                                                      parse(UInt8, row.FragmentCharge::PosLenString),
                                                      occursin("y", row.FragmentType::PosLenString),
                                                      parse(UInt8, row.FragmentNumber::PosLenString),
                                                      ion_position,
                                                      parse(Float32, row.RelativeIntensity::PosLenString),
                                                      parse(UInt8, row.PrecursorCharge::PosLenString),
                                                      prec_id,
                                                    )
                                 )
        ion_position += UInt8(1)
    end
    sort!(frag_list, by = x->getFragMZ(x))
    return frag_list, frag_detailed, precursor_list, prec_id
end


using Dictionaries
using CSV, Arrow, Tables, DataFrames, StatsBase
include("src/precursor.jl")
#include("src/buildFragmentIndex.jl")

@time prosit_list_simple, prosit_list_detailed = readPrositLib("/Users/n.t.wamsley/Desktop/myPrositLib.csv")

@save "/Users/n.t.wamsley/Projects/prosit_list_simple.jld2"  prosit_list_simple
@save "/Users/n.t.wamsley/Projects/prosit_list_detailed.jld2"  prosit_list_detailed

@load "/Users/n.t.wamsley/Projects/prosit_list_simple.jld2"  prosit_list_simple
@load "/Users/n.t.wamsley/Projects/prosit_list_detailed.jld2"  prosit_list_detailed

prosit_index = buildFragmentIndex!(prosit_list_simple, 10.0)

X = reshape([Float64(x) for x in range(1, 10)], (1, 10))

W, H = NMF.randinit(X, 2)

W = reshape([1.0 for x in range(1,30)], (1, 30))

H = Matrix([1.0 2 3 4 5 0 0 0 0 0; 0 0 0 0 0 6 7 8 9 10]')

NMF.solve!(NMF.CoordinateDescent{Float64}(maxiter=50, α=0.5, l₁ratio=1.0, update_H = false), X, W, H)


NMF.solve!(NMF.CoordinateDescent{Float64}(maxiter=50, α=0.5, l₁ratio=1.0, update_H = false), test_X, 1000*W, test_H)


NMF.solve!(NMF.CoordinateDescent{Float64}(maxiter=50, α=0.5, l₁ratio=1.0, update_H = false), test_X, 1000*W, test_H)


W = reshape([Float32(1000) for x in range(1,30)], (1, 30))
out = NMF.solve!(NMF.GreedyCD{Float32}(maxiter=50, verbose = false, lambda_w = 1e3, tol = 1e-6, update_H = false), test_X, W, test_H)

function spectralContrast(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
    dot(a, b)/(norm(a)*norm(b))
end

function buildDesignMatrix(matches::Vector{FragmentMatch{Float32}},  misses::Vector{FragmentMatch{Float32}}, seed_size::Int) #where {T<:AbstractFloat}

    #Number of unique matched peaks.
    matched_peaks = length(unique([getPeakInd(x) for x in matches]))
    #Number of rows equals the number of unique matched peaks + the number of expected fragments that 
    #failed to match a peak in the spectrm
    M = (matched_peaks + length(misses))
    #Design matrix. One row for every precursor template. One column for every matched + missed peak. 
    H = zeros(Float32, (seed_size, M))
    #Spectrum/empirical intensities for each peak. Zero by default (for unmatched/missed fragments)
    X = zeros(Float32, (1, M))

    #Maps a precursor id to a row of H. 
    precID_to_row = UnorderedDictionary{UInt32, UInt8}()

    #Current highest row encountered
    prec_row = UInt8(0)
    col = 0
    #Number of unique peaks encountered. 
    last_peak_ind = 0
    for match in matches
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(match))
            prec_row  += UInt8(1)
            insert!(precID_to_row, getPrecID(match), prec_row )
        end

        #If this peak has not been encountered yet, then start filling a new column
        if getPeakInd(match) != last_peak_ind
            col += 1
            last_peak_ind = getPeakInd(match)
        end

        row = precID_to_row[getPrecID(match)]
        H[row, col] = getPredictedIntenisty(match)
        X[1, col] = getIntensity(match)
    end

    for miss in misses
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(miss))
            prec_row  += UInt8(1)
            insert!(precID_to_row, getPrecID(miss), prec_row)
        end
        col = getPeakInd(miss) + matched_peaks
        row = precID_to_row[getPrecID(miss)]
        H[row, col] = getPredictedIntenisty(miss)
    end

    return X, H, precID_to_row
end

function getSpectralContrast(H::Matrix{T}, X::Matrix{T}) where {T<:AbstractFloat}
    function spectralContrast(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
        dot(a, b)/(norm(a)*norm(b))
    end

    N = size(H)[1]
    spectral_contrast = Vector{T}(undef, N)

    for row in range(1, N)
        non_zero = H[N,:].!=0
        spectral_contrast[row] = spectralContrast(H[N, non_zero], X[1, non_zero])
    end

    return spectral_contrast 
end

spectral_contrast = [spectralContrast(H[N,H[N,:].!=0], X[1,H[N,:].!=0]) for N in range(1, topN)]