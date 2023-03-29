function buildPrecursorTable!(ptable::PrecursorTable, n::Int, f_path::String)
    open(f_path) do f #"./data/NRF2_SIL.txt"
        pepGroup_id, pep_id, prot_id = UInt32(1), UInt32(1), UInt32(1)
        timetaken = @elapsed for (row, protein_peptide) in enumerate(eachline(f))
            protein, peptide = map(string, split(protein_peptide, "\t")); #Parse input "PROTEIN_NAME\tPEPTIDE_SEQUENCE"
            peptide = fixedMods(peptide, fixed_mods); #Apply fixed modifications
            if !containsProt(ptable, protein) #If the protien hasn't been encountered,
                prot_id += UInt32(1);                  #then add it to the hash table
                addNewProtein!(protein, prot_id, ptable);
            end
            if !containsPepGroup(ptable, peptide) #If the peptide hasn't been encountered,
                pepGroup_id += UInt32(1);                  #then do
                newPeptideGroup!(peptide, pepGroup_id, protein, ptable); #add it to the hash table,
                addPepGroupToProtein!(ptable, protein, peptide); #add a new peptide group to the protein,
                pep_id = applyMods!(var_mods,              #and lastly, apply variable mods and ad them to the peptide hash table
                                    peptide,               #and increase the pep_id for each variable mod applied 
                                    ptable.id_to_pep,
                                    pepGroup_id,
                                    pep_id,
                                    n = 2); #
            else #If this peptide has been encountered before, we don't need to apply the variable modes. Instead,
                addProteinToPepGroup!(ptable, protein, peptide); #Add the current protein to this peptide group
                addPepGroupToProtein!(ptable, protein, peptide); #Add the current peptide group to this protein
            end
        end
        println("Time to build precursor table ", timetaken);
    end
end