
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
    prec_intensity::Base.RefValue{T} #Needs to be updated
    prec_rt::T
    prec_charge::UInt8
end

getPrecMZ(f::FragmentIon) = f.prec_mz
getIntensity(f::FragmentIon) = f.prec_intensity[]
getRT(f::FragmentIon) = f.prec_rt

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

struct LibraryPrecursor{T<:AbstractFloat}
    iRT::T
    mz::T
    total_intensity::Base.RefValue{T}
    isDecoy::Bool
    charge::UInt8
end

isDecoy(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.isDecoy
getIRT(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.iRT
getCharge(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.charge
getMz(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.mz
getTotalIntensity(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.total_intensity
addIntensity!(p::LibraryPrecursor{T}, intensity::T) where {T<:AbstractFloat} = p.total_intensity[] += intensity

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
                                                    parse(precision, row.PrecursorMz::PosLenString),
                                                    Ref(zero(precision)),
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

        addIntensity!(precursor_list[id], parse(precision, row.RelativeIntensity))

        push!(frag_list, FragmentIon(parse(precision, row.FragmentMz::PosLenString), 
                                    prec_id, 
                                    parse(precision, row.PrecursorMz::PosLenString), 
                                    Ref(parse(precision, row.RelativeIntensity::PosLenString)),
                                    parse(precision, row.iRT::PosLenString),
                                    parse(UInt8, row.PrecursorCharge::PosLenString)
                                    ))

        push!(frag_detailed[id], LibraryFragment(parse(precision, row.FragmentMz::PosLenString), 
                                                      parse(UInt8, row.FragmentCharge::PosLenString),
                                                      occursin("y", row.FragmentType::PosLenString),
                                                      parse(UInt8, row.FragmentNumber::PosLenString),
                                                      ion_position,
                                                      parse(precision, row.RelativeIntensity::PosLenString),
                                                      parse(UInt8, row.PrecursorCharge::PosLenString),
                                                      prec_id,
                                                    )
                                 )
        ion_position += UInt8(1)
    end
    sort!(frag_list, by = x->getFragMZ(x))
    return frag_list, frag_detailed, precursor_list, prec_id
end

##########
#Example Use
##########
#=
prec_id = 1
target_frag_list, target_frag_detailed, target_precursor_list, prec_id = readPrositLib("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_targets.csv",
                                                                                        precision = Float32,
                                                                                        isDecoys = false)

decoy_frag_list, decoy_frag_detailed, decoy_precursor_list, prec_id = readPrositLib("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_decoys.csv",
                                                                                        precision = Float32,
                                                                                        isDecoys = true,
                                                                                        first_prec_id = prec_id)
#Merge target and decoy lists
frag_list = append!(target_frag_list, decoy_frag_list)
target_frag_list = nothing
decoy_frag_list = nothing
frag_detailed = append!(target_frag_detailed, decoy_frag_detailed)
target_frag_detailed = nothing
decoy_frag_detailed = nothing
precursor_list = append!(target_precursor_list, decoy_precursor_list)
target_precursor_list = nothing
decoy_precursor_list = nothing
#Normalize
function NormalizeIntensities!(frag_list::Vector{FragmentIon{T}}, prec_list::Vector{LibraryPrecursor}) where {T<:AbstractFloat}
    for i in ProgressBar(1:length(frag_list))
        prec_id = frag_list[i].prec_id
        frag_list[i].prec_intensity[] = frag_list[i].prec_intensity[]/prec_list[prec_id].total_intensity[]
    end
end
NormalizeIntensities!(frag_list, precursor_list)
#Save to Disc
@save "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/frag_list.jld2" frag_list
@save "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/frag_detailed.jld2" frag_detailed
@save "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/precursor_list.jld2" precursor_list

#Build Fragment Index 
@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/frag_list.jld2" frag_list
@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/frag_detailed.jld2" frag_detailed
@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/precursor_list.jld2" precursor_list

prosit_index_5ppm_5irt = buildFragmentIndex!(frag_list, Float32(5.0), Float32(5.0))
for i in 1:length(prosit_index_5ppm_5irt.precursor_bins)
    sort!(prosit_index_5ppm_5irt.precursor_bins[i].precs, by = x->getPrecMZ(x))
end
@save "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/prosit_index_5ppm_5irt.jld2" prosit_index_5ppm_5irt
prosit_index_5ppm_5irt = nothing

prosit_index_5ppm_10irt = buildFragmentIndex!(frag_list, Float32(5.0), Float32(10.0))
for i in 1:length(prosit_index_5ppm_10irt.precursor_bins)
    sort!(prosit_index_5ppm_10irt.precursor_bins[i].precs, by = x->getPrecMZ(x))
end
@save "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/prosit_index_5ppm_10irt.jld2" prosit_index_5ppm_10irt
prosit_index_5ppm_10irt = nothing

prosit_index_5ppm_15irt = buildFragmentIndex!(frag_list, Float32(5.0), Float32(15.0))
for i in 1:length(prosit_index_5ppm_15irt.precursor_bins)
    sort!(prosit_index_5ppm_15irt.precursor_bins[i].precs, by = x->getPrecMZ(x))
end
@save "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/prosit_index_5ppm_15irt.jld2" prosit_index_5ppm_15irt
prosit_index_5ppm_15irt = nothing

prosit_index_5ppm_20irt = buildFragmentIndex!(frag_list, Float32(5.0), Float32(20.0))
for i in 1:length(prosit_index_5ppm_20irt.precursor_bins)
    sort!(prosit_index_5ppm_20irt.precursor_bins[i].precs, by = x->getPrecMZ(x))
end
@save "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/prosit_index_5ppm_20irt.jld2" prosit_index_5ppm_20irt
prosit_index_5ppm_20irt = nothing

@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/prosit_index_5ppm_20irt.jld2" prosit_index_5ppm_15irt

prosit_index_5ppm_10irt = buildFragmentIndex!(frag_list, Float32(5.0))
@save "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/prosit_index_5ppm_10irt.jld2" prosit_index_5ppm_10irt
@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/prosit_index_5ppm_10irt.jld2" prosit_index_5ppm_10irt
for i in 1:length(prosit_index_5ppm_5irt.precursor_bins)
    sort!(prosit_index_5ppm_5irt.precursor_bins[i].precs, by = x->getPrecMZ(x))
end
prosit_index_5ppm_10irt = nothing

prosit_index_3ppm_5irt = buildFragmentIndex!(frag_list, Float32(5.0))
@save "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/prosit_index_3ppm_10irt.jld2" prosit_index_3ppm_10irt
prosit_index_3ppm_5irt = nothing

@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/prosit_index_5ppm_15irt.jld2" prosit_index_5ppm_15irt
@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/frag_detailed.jld2" frag_detailed
@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/precursor_list.jld2" precursor_list

for i in 1:length(prosit_index_5ppm_5irt.precursor_bins)
    sort!(prosit_index_5ppm_5irt.precursor_bins[i], by = x->getPrecMZ(x))
end
=#
