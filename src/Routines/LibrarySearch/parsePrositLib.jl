abstract type IonType end


abstract type FragmentIndexType <: IonType end
getMZ(f::FragmentIndexType) = f.frag_mz
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
    sequence::String
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
                 "Ox" => Float64(15.994915)
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

function split_array_into_chunks(arr::Vector, num_chunks::Int)
    len = length(arr)
    chunk_size = div(len, num_chunks)
    
    # Create a list of tuples for each chunk
    chunks = [[start_idx, min(start_idx + chunk_size - 1, len)] for start_idx in 1:chunk_size:len]
    
    # Adjust the last chunk's stopping index to the end of the array
    last_chunk = last(chunks)
    if last_chunk[2] < len
        last_chunk = [last_chunk[1], len]
        chunks[end] = last_chunk
    end
    
    return chunks
end
#CSV.write("outputfile.csv",test_prosit[1000000:1000036,:],delim=',')

function parsePrositLib(prosit_csv_path::String, fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, mods_dict::Dict{String, T}; start_ion::Int = 3) where {T<:AbstractFloat}

    #"/Users/n.t.wamsley/Projects/PROSIT/prosit1/my_prosit/examples/peptidelist.msms"
    prosit_library = CSV.File(prosit_csv_path)
    charge_pattern = r"\((\w+)\+\)"
    index_pattern = r"[yb]([0-9]{1,2})[\(]*"
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
    #lk = ReentrantLock()
    N = split_array_into_chunks(frags_detailed, Threads.nthreads())
    Threads.@threads for precursor in ProgressBar(collect(prosit_library))
        if N[Threads.threadid()][1]>N[Threads.threadid()][2]
            #sleep(rand())
            continue
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
        decoy =  precursor[:decoy]
        iRT = Float32(precursor[:iRT])
        sequence = precursor[:modified_sequence]

        accession_numbers = precursor[:accession_numbers]
        mz, len_AA, missed_cleavages = parseSequence(String(sequence), charge, fixed_mods, mods_dict)

        #pre-allocate library fragments for the n'th precursor
        nth_precursor_frags = Vector{LibraryFragment{Float32}}()#undef, length(matches))
        total_intensity = zero(Float32)
        #########
        #Parse each fragment ion 
        for (i, _match) in enumerate(matches)
            if parse(UInt8, match(index_pattern, _match).captures[1]) < start_ion
                continue
            end
            total_intensity += intensities[i]
        end

        for i in eachindex(matches)
            fragment_name = matches[i]
            #ion_type = match(pattern, fragment_name)
            ion_index = parse(UInt8, match(index_pattern, fragment_name).captures[1]) #for y1 or b12 this is "1" and "12 respectively
            ion_charge = match(charge_pattern, fragment_name)
            frag_charge = ion_charge === nothing ? one(UInt8) : parse(UInt8, ion_charge.captures[1])
            if ion_index < start_ion
                continue
            end
            #Add the n'th fragment 
            push!(nth_precursor_frags,  LibraryFragment(
                                                masses[i], #frag_mz
                                                frag_charge, #frag_charge
                                                fragment_name[1] == 'y' ? true : false, #is_y_ion
                                                ion_index, #ion_position
                                                UInt8(i), #ion_index
                                                intensities[i]/total_intensity, #intensity
                                                charge, #prec_charge
                                                UInt32(N[Threads.threadid()][1]) #prec_id
            ))

            push!(frags_simple[Threads.threadid()],
                    FragmentIon(
                        masses[i], #frag_mz
                        UInt32(N[Threads.threadid()][1]), #prec_id
                        mz, #prec_mz 
                        intensities[i]/total_intensity, #prec_intensity
                        iRT, #prec_rt
                        charge, #prec_charge
                    )
            )

        end
        precursors[N[Threads.threadid()][1]] = LibraryPrecursor(
                                            iRT, #iRT
                                            mz, #mz
                                            total_intensity, #total_intensity
                                            decoy, #isDecoy
                                            charge, #charge
                                            pep_id, #pep_id
                                            UInt32[1], #prot_ids
                                            String(accession_numbers), #accession_numbers
                                            String(sequence),
                                            missed_cleavages, #missed_cleavages
                                            UInt8(1), #variable_mods
                                            len_AA, #length
                                        )

        frags_detailed[N[Threads.threadid()][1]] = nth_precursor_frags

        N[Threads.threadid()][1] += 1

    end
    return vcat(frags_simple...), frags_detailed, precursors
    #return frags_simple, frags_detailed, precursors
end
@time frags_simple, frags_detailed, precursors = parsePrositLib("/Users/n.t.wamsley/Projects/PROSIT/prosit1/my_prosit/prosit_mouse_NCE33_dynamicNCE_073123.csv", fixed_mods, mods_dict);
frags_mouse_simple_33NCEdynamic = frags_simple
frags_mouse_detailed_33NCEdynamic = frags_detailed
precursors_mouse_detailed_33NCEdynamic = precursors
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/frags_mouse_simple_33NCEdynamic.jld2" frags_mouse_simple_33NCEdynamic
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/frags_mouse_detailed_33NCEdynamic.jld2" frags_mouse_detailed_33NCEdynamic
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/precursors_mouse_detailed_33NCEdynamic.jld2" precursors_mouse_detailed_33NCEdynamic
prosit_mouse_33NCEdynamic_5ppm_15irt = buildFragmentIndex!(frags_mouse_simple_33NCEdynamic, Float32(5.0), Float32(15.0))
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/prosit_mouse_33NCEdynamic_5ppm_15irt.jld2" prosit_mouse_33NCEdynamic_5ppm_15irt



@time frags_simple, frags_detailed, precursors = parsePrositLib("/Users/n.t.wamsley/Projects/PROSIT/prosit1/my_prosit/prosit_mouse_NCE33_fixedNCE_073123.csv", fixed_mods, mods_dict);
frags_mouse_simple_33NCEfixed = frags_simple
frags_mouse_detailed_33NCEfixed = frags_detailed
precursors_mouse_detailed_33NCEfixed = precursors
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/frags_mouse_simple_33NCEfixed.jld2" frags_mouse_simple_33NCEfixed
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/frags_mouse_detailed_33NCEfixed.jld2" frags_mouse_detailed_33NCEfixed
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/precursors_mouse_detailed_33NCEfixed.jld2" precursors_mouse_detailed_33NCEfixed
prosit_mouse_33NCEfixed_5ppm_15irt = buildFragmentIndex!(frags_mouse_simple_33NCEfixed, Float32(5.0), Float32(15.0))
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/prosit_mouse_33NCEfixed_5ppm_15irt.jld2" prosit_mouse_33NCEfixed_5ppm_15irt




@time frags_simple, frags_detailed, precursors = parsePrositLib("/Users/n.t.wamsley/Projects/PROSIT/prosit1/my_prosit/prosit_mouse_NCE33_fixedNCE_073123.csv", fixed_mods, mods_dict);
frags_mouse_simple_33NCEfixed = frags_simple
frags_mouse_detailed_33NCEfixed = frags_detailed
precursors_mouse_detailed_33NCEfixed = precursors
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/frags_mouse_simple_33NCEfixed.jld2" frags_mouse_simple_33NCEfixed
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/frags_mouse_detailed_33NCEfixed.jld2" frags_mouse_detailed_33NCEfixed
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/precursors_mouse_detailed_33NCEfixed.jld2" precursors_mouse_detailed_33NCEfixed
prosit_mouse_33NCEfixed_5ppm_15irt = buildFragmentIndex!(frags_mouse_simple_33NCEfixed, Float32(5.0), Float32(15.0))
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_080123/prosit_mouse_33NCEfixed_5ppm_15irt.jld2" prosit_mouse_33NCEfixed_5ppm_15irt



@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/frags_simple_33NCEfixed.jld2" frags_simple_33NCEfixed
@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/frags_detailed_33NCEfixed.jld2" frags_detailed_33NCEfixed
@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/precursors_33NCEfixed.jld2" precursors_33NCEfixed
@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/prosit_33NCEfixed_5ppm_15irt.jld2" prosit_33NCEfixed_5ppm_15irt 

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
    rows = []
    println("TEST")
    for (i, row) in enumerate(rows)
        println(i)
        if (i % 1_000_000) == 0
            println(i/1_000_000)
        end
        #if i < 13000000
        #    continue
        #end
        if current_peptide == "_GGPGSAVSPYPSFNVSSDVAALHK_"
            println("current_peptide $current_peptide")
            println(row)
            rows += [row]
        end
        if (row.ModifiedPeptide::PosLenString != current_peptide) | (row.PrecursorCharge::PosLenString != current_charge)
            current_peptide = row.ModifiedPeptide::PosLenString
            current_charge = row.PrecursorCharge::PosLenString
            prec_id += UInt32(1)
            #if current_peptide == "_GGPGSAVSPYPSFNVSSDVAALHK_"
            #    println("current_peptide $current_peptide")
            #    println(row)
            #end
            id += UInt32(1)
            ion_position = UInt8(1)
            push!(frag_detailed, Vector{LibraryFragment{precision}}())
            #=push!(precursor_list, LibraryPrecursor(
                                                    parse(precision, row.iRT::PosLenString),
                                                    parse(precision, row.PrecursorMz::PosLenString),
                                                    Ref(zero(precision)),
                                                    isDecoys,
                                                    parse(UInt8, row.PrecursorCharge::PosLenString)
                                                ))=#
        end

        #Track progress
        if (i % 1_000_000) == 0
            println(i/1_000_000)
        end

        #Exclude y1, y2, b1, and b2 ions. 
        if parse(Int, row.FragmentNumber::PosLenString) < 3
            continue
        end

        #=addIntensity!(precursor_list[id], parse(precision, row.RelativeIntensity))

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
                                 )=#
        ion_position += UInt8(1)
    end
    sort!(frag_list, by = x->getFragMZ(x))
    return frag_list, frag_detailed, precursor_list, prec_id
end

rows = []
current_peptide = ""
current_charge = ""
for (i, row) in enumerate(CSV.Rows("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_targets.csv"))
    if (i % 1_000_000) == 0
        println(i/1_000_000)
    end
    if i < 13000000
        continue
    end
    current_peptide = row.ModifiedPeptide::PosLenString
    current_charge = row.PrecursorCharge::PosLenString
    if current_peptide == "_GGPGSAVSPYPSFNVSSDVAALHK_"
        if current_charge == "3"
            #println("current_peptide $current_peptide")
            println(row.ModifiedPeptide,",",row.RelativeIntensity,",",row.FragmentMz,",",row.FragmentType,",",",",row.FragmentNumber,",",row.FragmentCharge)
        end
        #rows += [row]
    end
end
#=

@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/frags_simple.jld2" frags_simple
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/frags_detailed.jld2" frags_detailed
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/precursors.jld2" precursors
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/prosit_index_5ppm_15irt.jld2" prosit_index_5ppm_15irt


@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/frags_simple.jld2" frags_simple
@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/frags_detailed.jld2" frags_detailed
@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/precursors.jld2" precursors


prosit_index_5ppm_15irt = buildFragmentIndex!(frags_simple, Float32(5.0), Float32(15.0))
#for i in 1:length(prosit_index_5ppm_10irt.precursor_bins)
#    sort!(prosit_index_5ppm_10irt.precursor_bins[i].precs, by = x->getPrecMZ(x))
#end
@save "/Users/n.t.wamsley/Projects/PROSIT/mouse_072723/prosit_index_5ppm_15irt.jld2" prosit_index_5ppm_15irt

frags2 = DetailedPrecursor(append!(frags, frags))
table = (col1=frags)
io = IOBuffer()
Arrow.write(io, table)
seekstart(io)
table2 = Arrow.Table(io)
=#
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
