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
```
        julia> combined_scored_psms
        Dict{Symbol, Vector} with 7 entries:
          :hyperscore    => [4.78763, 37.7446, 3.17819, 39.4188, 1.79189, 6.57939, 1.79189, 0.000286995, 32.…
          :precursor_idx => UInt32[0x0000013c, 0x000001c6, 0x00000220, 0x000001c6, 0x00000220, 0x000001c6, 0…
          :total_ions    => UInt32[0x00000005, 0x00000009, 0x00000004, 0x0000000a, 0x00000003, 0x00000006, 0…
          :ms_file_idx   => UInt32[0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001, 0…
          :error         => [0.0149231, 0.0289612, 0.0366211, 0.0270996, 0.00613403, 0.0270691, 0.00680542, …
          :scan_idx      => [4, 5, 6, 11, 12, 16, 17, 18, 22, 23  …  52596, 52600, 52608, 52616, 52620, 5262…
          :y_ladder      => Int8[5, 5, 4, 5, 3, 5, 3, 1, 5, 4  …  1, 1, 3, 3, 1, 3, 3, 1, 1, 1]
        
```
        This enables easy conversion to a  `DataFrame` of PSMs. 

### Examples 

"""
function getBestPSMs(PSMs::Dict{Symbol, Vector}, ptable::PrecursorDatabase, MS_RT::Dict{UInt32, Vector{Float32}}, min_fragment_count::UInt8)
    PSMs = DataFrame(PSMs)

    #Remove PSMs where fewer than `min_fragment_count` transitions matched a peak 
    filter!(row -> row.total_ions >= min_fragment_count, PSMs); 

    #Get the peptide_id corresponding to each PSM
    transform!(PSMs, AsTable(:) => ByRow(psm -> getPepIDFromPrecID(ptable, psm[:precursor_idx])) => :pep_idx)

    #Charge and isotope state of precursor
    transform!(PSMs, AsTable(:) => ByRow(psm -> getCharge(getPrecursor(ptable,psm[:precursor_idx]))) => :precursor_charge)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getIsotope(getPrecursor(ptable,psm[:precursor_idx]))) => :precursor_isotope)

    #Combine rows on file and pep_idx, retain only the row with the highest hyperscore in each group. These are the "best scans"
    PSMs = combine(sdf -> sdf[argmax(sdf.hyperscore), :], groupby(PSMs, [:ms_file_idx, :pep_idx])) 

    #Get the retention times, peptide sequence, and precursor mz. 
    transform!(PSMs, AsTable(:) => ByRow(psm -> MS_RT[psm[:ms_file_idx]][psm[:scan_idx]]) => :retention_time)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getSeq(ptable.id_to_pep[psm[:pep_idx]])) => :sequence)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getMZ(getPrecursor(ptable,psm[:precursor_idx]))) => :precursor_mz)

    #Get the protein names for proteins corresponding to the peptide
    transform!(PSMs, AsTable(:) => ByRow(psm -> join(sort((collect(getProtNamesFromPepSeq(ptable, psm[:sequence])))), "|")) => :protein_names)

    #Sort by retention time.
    sort!(PSMs, [:retention_time])
    PSMs
end