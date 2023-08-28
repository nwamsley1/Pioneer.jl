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



@time test_table = SearchRAW(MS_TABLE, prosit_index_all, prosit_detailed, UInt32(1))
PSMs = test_table
transform!(PSMs, AsTable(:) => ByRow(psm -> isDecoy(prosit_precs[psm[:precursor_idx]])) => :decoy)
transform!(PSMs, AsTable(:) => ByRow(psm -> getIRT(prosit_precs[psm[:precursor_idx]])) => :iRT)
transform!(PSMs, AsTable(:) => ByRow(psm -> MS_TABLE[:retentionTime][psm[:scan_idx]]) => :RT)
transform!(PSMs, AsTable(:) => ByRow(psm -> psm[:weight] == 0) => :nmf)


best_PSMs = combine(sdf -> sdf[argmax(sdf.hyperscore), :], groupby(PSMs, [:scan_idx])) 
#transform!(best_PSMs, AsTable(:) => ByRow(psm -> MS_TABLE[:retentionTime][psm[:scan_idx]]) => :RT)
plot(best_PSMs[:,:iRT][best_PSMs[:,:decoy].==false], best_PSMs[:,:RT][best_PSMs[:,:decoy].==false], seriestype = :scatter)
plot!(best_PSMs[:,:iRT][best_PSMs[:,:decoy].==true], best_PSMs[:,:RT][best_PSMs[:,:decoy].==true], seriestype = :scatter)

iRTs = best_PSMs[:,:iRT][(best_PSMs[:,:decoy].==false) .& (best_PSMs[:,:spectral_contrast].>=0.7)]
iRTs = reshape(iRTs, (length(iRTs), 1))
RTs = best_PSMs[:,:RT][(best_PSMs[:,:decoy].==false) .& (best_PSMs[:,:spectral_contrast].>=0.7)]
RTs  = [Float64(x) for x in RTs]
iRTs= hcat(iRTs, ones(length(iRTs)))
rlm(iRTs, RTs, TauEstimator{TukeyLoss}(), initial_scale=:mad, maxiter = 200)


plot(iRTs, RTs, seriestype = :scatter)
plot!([0, 150], [18.6381, 150*0.263446 + 18.6381])

transform!(PSMs, AsTable(:) => ByRow(psm -> psm[:iRT]*0.263446 + 18.6381) => :predRT)

transform!(PSMs, AsTable(:) => ByRow(psm -> abs(psm[:RT] - psm[:predRT])) => :RTdiff)
X = Matrix(PSMs[1:1:end,[:hyperscore,:total_ions,:y_ladder,:b_ladder,:intensity_explained,:error,:poisson,:spectral_contrast,:spectrum_peaks,:nmf,:weight,:RTdiff]])'
X_labels = Vector(PSMs[1:1:end, :decoy])
lda = fit(MulticlassLDA, X, X_labels; outdim=1)
Ylda = predict(lda, X)
PSMs[:,:score] = Ylda[1,:]
histogram(Ylda[X_labels.==true], alpha = 0.5, normalize = :pdf)#, bins = -0.06:0.01:0.0)
histogram!(Ylda[X_labels.==false], alpha = 0.5, normalize = :pdf)#, bins = -0.06:0.01:0.0)

plot(PSMs[X_labels.==true,:spectral_contrast], PSMs[X_labels.==true,:RTdiff], alpha = 0.5, seriestype=:scatter)#, bins = -0.06:0.01:0.0)

plot!(PSMs[X_labels.==false,:spectral_contrast], PSMs[X_labels.==false,:RTdiff], alpha = 0.5, seriestype=:scatter)#, bins = -0.06:0.01:0.0)


sum(Ylda[X_labels.==true].<-0.0085)
sum(Ylda[X_labels.==false].<-0.0085)

model = build_forest(X_labels, X', 4, 2000, 0.5, 3)
probs = apply_forest_proba(model, X',[true, false])
PSMs[:,:prob] = probs[:,2]


histogram(PSMs[X_labels.==true,:prob], alpha = 0.5)#, bins = -0.06:0.01:0.0)
histogram!(PSMs[X_labels.==false,:prob], alpha = 0.5)#, bins = -0.06:0.01:0.0)


sum(PSMs[X_labels.==true,:prob].>0.91)
sum(PSMs[X_labels.==false,:prob].>0.91)

unique(PSMs[PSMs[:,:prob].>0.91,:precursor_idx])

struct LibraryPrecursor{T<:AbstractFloat}
    iRT::T
    isDecoy::Bool
    charge::UInt8
end
isDecoy(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.isDecoy
getIRT(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.iRT

prosit_simple_targets, prosit_detailed_targets, prosit_precs_targets, prec_id = readPrositLib("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_targets.csv")
prosit_simple_decoys, prosit_detailed_decoys, prosit_precs_decoys, prec_id = readPrositLib("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_decoys.csv", isDecoys = true, first_prec_id = prec_id)


@save "/Users/n.t.wamsley/Projects/prosit_simple_targets.jld2"  prosit_simple_targets
@save "/Users/n.t.wamsley/Projects/prosit_detailed_targets.jld2"  prosit_detailed_targets
@save "/Users/n.t.wamsley/Projects/prosit_precs_targets.jld2"  prosit_precs_targets

@save "/Users/n.t.wamsley/Projects/prosit_simple_decoys.jld2"  prosit_simple_decoys
@save "/Users/n.t.wamsley/Projects/prosit_detailed_decoys.jld2"  prosit_detailed_decoys
@save "/Users/n.t.wamsley/Projects/prosit_precs_decoys.jld2"  prosit_precs_decoys


@load "/Users/n.t.wamsley/Projects/prosit_simple_decoys.jld2" prosit_simple_decoys
@load "/Users/n.t.wamsley/Projects/prosit_simple_targets.jld2" prosit_simple_targets
max_ = 0
function getMaxID()
max_ = 0
for decoy in prosit_simple_all
    if getPrecID(decoy)>max_
        max_ = getPrecID(decoy)
    end
end

prosit_detailed = append!(prosit_detailed_targets, prosit_detailed_decoys)
@save "/Users/n.t.wamsley/Projects/prosit_detailed.jld2"  prosit_detailed 
prosit_simple_all = append!(prosit_simple_targets, prosit_simple_decoys)
sort!(prosit_simple_all, by = x->getFragMZ(x))
@save "/Users/n.t.wamsley/Projects/prosit_simple_all.jld2"  prosit_simple_all
prosit_precs = append!(prosit_precs_targets, prosit_precs_decoys)
@save "/Users/n.t.wamsley/Projects/prosit_precs.jld2"  prosit_precs

prosit_index_all = buildFragmentIndex!(prosit_simple_all, 10.0)
@save "/Users/n.t.wamsley/Projects/prosit_index_all.jld2"  prosit_index_all

function readPrositLib(prosit_lib_path::String; precision::DataType = Float64, isDecoys::Bool = false, first_prec_id = UInt32(0))
    frag_list = Vector{FragmentIon{precision}}()
    frag_detailed = Vector{Vector{LibraryFragment{precision}}}()
    precursor_list = Vector{LibraryPrecursor}()

    rows = CSV.Rows(prosit_lib_path, reusebuffer=false, select = [:RelativeIntensity, :FragmentMz, :PrecursorMz, :iRT, :Stripped, :ModifiedPeptide,:FragmentNumber,:FragmentCharge,:PrecursorCharge,:FragmentType])
    current_peptide = ""
    current_charge = ""
    prec_id = UInt32(first_prec_id)
    id = UInt32(0)
    ion_position = UInt8(1)
    for (i, row) in enumerate(rows)
        if (row.ModifiedPeptide::PosLenString != current_peptide) | (row.PrecursorCharge::PosLenString != current_charge)
            current_peptide = row.ModifiedPeptide::PosLenString
            current_charge = row.PrecursorCharge::PosLenString
            prec_id += UInt32(1)
            id += UInt32(1)
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
        non_zero = H[row,:].!=0
        spectral_contrast[row] = spectralContrast(H[row, non_zero], X[1, non_zero])
    end

    return spectral_contrast 
end

spectral_contrast = [spectralContrast(H[N,H[N,:].!=0], X[1,H[N,:].!=0]) for N in range(1, topN)]