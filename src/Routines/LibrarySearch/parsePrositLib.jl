
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

#prosit_index = buildFragmentIndex!(prosit_list_simple, 10.0)
#prosit_simple_decoys, prosit_detailed_decoys, prosit_precs_decoys, prec_id = readPrositLib("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_targets.csv", isDecoys = true, first_prec_id = UInt32(0))
#prosit_simple_decoys, prosit_detailed_decoys, prosit_precs_decoys, prec_id = readPrositLib("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_decoys.csv", isDecoys = true, first_prec_id = prec_id)


#=

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
=#