function getBestPSMs(PSMs::Dict{Symbol, Vector}, ptable::PrecursorDatabase, MS_RT::Dict{UInt32, Vector{Float32}}, min_fragment_count::UInt8)
    PSMs = DataFrame(PSMs)
    filter!(row -> row.total_ions >= min_fragment_count, PSMs);
    transform!(PSMs, AsTable(:) => ByRow(psm -> getPepIDFromPrecID(ptable, psm[:precursor_idx])) => :pep_idx)
    #Charge and isotope state of precursor
    transform!(PSMs, AsTable(:) => ByRow(psm -> getCharge(getPrecursor(ptable,psm[:precursor_idx]))) => :precursor_charge)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getIsotope(getPrecursor(ptable,psm[:precursor_idx]))) => :precursor_isotope)
    PSMs = combine(sdf -> sdf[argmax(sdf.hyperscore), :], groupby(PSMs, [:ms_file_idx, :precursor_idx])) #combine on file and pep_idx, retain only the row with the highest hyperscore in each group
    transform!(PSMs, AsTable(:) => ByRow(psm -> MS_RT[psm[:ms_file_idx]][psm[:scan_idx]]) => :retention_time)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getSeq(ptable.id_to_pep[psm[:pep_idx]])) => :sequence)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getMZ(getPrecursor(ptable,psm[:precursor_idx]))) => :precursor_mz)
    #(collect(getProtIDs(getIDToPepGroup(ptable)[getGroupID(getIDToPep(ptable)[psm[:pep_idx]])])))
    transform!(PSMs, AsTable(:) => ByRow(psm -> join(sort((collect(getProtNamesFromPepSeq(ptable, psm[:sequence])))), "|")) => :protein_names)
    
    sort!(PSMs, [:retention_time])
    PSMs
end