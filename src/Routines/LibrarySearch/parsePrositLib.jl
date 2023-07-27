
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
    prec_intensity::T #Needs to be updated
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
    total_intensity::T
    isDecoy::Bool
    charge::UInt8
    pep_id::UInt32
    prot_ids::Vector{UInt32}
    accession_numbers::String
    missed_cleavages::UInt8
    variable_mods::UInt8
    length::UInt8
end

isDecoy(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.isDecoy
getIRT(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.iRT
getCharge(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.charge
getMz(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.mz
getTotalIntensity(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.total_intensity
getPepID(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.pep_id
#addIntensity!(p::LibraryPrecursor{T}, intensity::T) where {T<:AbstractFloat} = p.total_intensity[] += intensity
ArrowTypes.arrowname(::Type{LibraryFragment{Float32}}) = :LibraryFragment
ArrowTypes.JuliaType(::Val{:LibraryFragment}) = LibraryFragment
fixed_mods = [(p=r"C", r="C[Carb]")]
mods_dict = Dict("Carb" => Float64(57.021464),
                 "Ox" => Float64(10.0)
                 )

function parseSequence(seq::String, charge::UInt8, fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, mods_dict::Dict{String, T}) where {T<:AbstractFloat}

    seq = fixedMods(replace(seq, "M(ox)"=>"M[Ox]"), fixed_mods)

    mz = getMZ(Precursor(
                            seq, 
                            charge = charge,
                            mods_dict = mods_dict
                        )
                )

    function getPlainSeq(seq::String)
        pattern = r"\[.*?\]"
        # Remove all occurrences of the pattern from the input string
        return replace(seq, pattern => "")
    end

    function countMissedCleavages(seq::String)
        pattern = r"[KR][^P|$]"
        # Count the number of matches of the pattern in the input string
        return length(collect((eachmatch(pattern, seq))))
    end
    seq = getPlainSeq(seq)

    return Float32(mz), UInt8(length(seq)), UInt8(countMissedCleavages(seq))
end

function parsePrositLib(prosit_csv_path::String, fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, mods_dict::Dict{String, T}) where {T<:AbstractFloat}

    #"/Users/n.t.wamsley/Projects/PROSIT/prosit1/my_prosit/examples/peptidelist.msms"
    prosit_library = CSV.File(prosit_csv_path)
    pattern = r"\((\w+)\+\)"

    ###########
    #Initialize Containers
    #List of lists. The n'th element of the outer list is indexed by the precursor id.
    #Each inner list contains a list of "LibraryFragment"s. These are a detailed specificaion of a framgent ion
    frags_detailed = Vector{Vector{LibraryFragment{Float32}}}(undef, length(prosit_library))

    #Detailed specification for each precursor. 
    precursors = Vector{LibraryPrecursor{Float32}}(undef, length(prosit_library))

    #Used to fuild the framgment index for MSFragger-like searching. 
    frags_simple = [Vector{FragmentIon{Float32}}() for _ in 1:Threads.nthreads()];
    #Loop through rows of prosit library (1-1 correspondence between rows and precursors)
    n = 0
    lk = ReentrantLock()
    Threads.@threads for precursor in ProgressBar(collect(prosit_library))
        lock(lk) do
            n += 1
        end
        #############
        #Parse Column/Precursor
        #############
        matches =  split(precursor[:matched_ions],';')
        masses = parse.(Float32, split(precursor[:masses],';'))
        intensities = parse.(Float32, split(precursor[:intensities],';'))
        intensities = intensities./sum(intensities)
        charge = UInt8(precursor[:Charge])
        pep_id =  UInt32(precursor[:pep_id])
        decoy =  precursor[:decoy]=="TRUE" ? true : false
        iRT = Float32(precursor[:iRT])
        sequence = precursor[:modified_sequence]
        accession_numbers = precursor[:accession_numbers]
        mz, len_AA, missed_cleavages = parseSequence(String(sequence), charge, fixed_mods, mods_dict)

        #pre-allocate library fragments for the n'th precursor
        nth_precursor_frags = Vector{LibraryFragment{Float32}}(undef, length(matches))

        #########
        #Parse each fragment ion 
        for i in eachindex(matches)
            fragment_name = matches[i]
            ion_type = match(pattern, fragment_name)
            frag_charge = ion_type === nothing ? one(UInt8) : parse(UInt8, ion_type.captures[1])

            #Add the n'th fragment 
            nth_precursor_frags[i] = LibraryFragment(
                                                masses[i], #frag_mz
                                                frag_charge, #frag_charge
                                                fragment_name[1] == 'y' ? true : false, #is_y_ion
                                                parse(UInt8, fragment_name[2]), #ion_position
                                                UInt8(i), #ion_index
                                                intensities[i], #intensity
                                                charge, #prec_charge
                                                UInt32(1) #prec_id
            )
            push!(frags_simple[Threads.threadid()],
                    FragmentIon(
                        masses[i], #frag_mz
                        UInt32(n), #prec_id
                        mz, #prec_mz 
                        intensities[i], #prec_intensity
                        iRT, #prec_rt
                        charge, #prec_charge
                    )
            )
        end

        precursors[n] = LibraryPrecursor(
                                            iRT, #iRT
                                            mz, #mz
                                            sum(intensities), #total_intensity
                                            decoy, #isDecoy
                                            charge, #charge
                                            pep_id, #pep_id
                                            UInt32[1], #prot_ids
                                            String(accession_numbers), #accession_numbers
                                            missed_cleavages, #missed_cleavages
                                            UInt8(1), #variable_mods
                                            len_AA, #length
                                        )

        frags_detailed[n] = nth_precursor_frags
    end
    return vcat(frags_simple...), frags_detailed, precursors
end

frags2 = DetailedPrecursor(append!(frags, frags))
table = (col1=frags)
io = IOBuffer()
Arrow.write(io, table)
seekstart(io)
table2 = Arrow.Table(io)
#=
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
=#
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
