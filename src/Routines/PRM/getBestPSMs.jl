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

        Dict{Symbol, Vector} with 6 entries:
        :hyperscore    => [0.000286995, 0.000286995, 0.000286995, ...
        :precursor_idx => UInt32[0x000000f1, 0x00000121, 0x000000a2, ...
        :total_ions    => UInt32[0x00000001, 0x00000001, 0x00000001, ...
        :error         => [0.00274658, 0.0233765, 0.00686646, ...
        :scan_idx      => [1720, 1724, 1731, ...
        :y_ladder      => Int8[1, 1, 1, ...

        This enables easy conversion to a  `DataFrame` of PSMs. 

### Examples 

"""
function getBestPSMs(PSMs::Dict{Symbol, Vector}, ptable::PrecursorTable, MS_TABLE::Dict{UInt32, Arrow.Table}, min_fragment_count::UInt8)
    PSMs = DataFrame(PSMs)
    filter!(row -> row.total_ions >= min_fragment_count, PSMs);
    transform!(PSMs, AsTable(:) => ByRow(psm -> getPepIDFromPrecID(ptable, psm[:precursor_idx])) => :pep_idx)
    #Charge and isotope state of precursor
    transform!(PSMs, AsTable(:) => ByRow(psm -> getCharge(getPrecursor(ptable,psm[:precursor_idx]))) => :precursor_charge)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getIsotope(getPrecursor(ptable,psm[:precursor_idx]))) => :precursor_isotope)
    PSMs = combine(sdf -> sdf[argmax(sdf.hyperscore), :], groupby(PSMs, [:ms_file_idx, :pep_idx])) #combine on file and pep_idx, retain only the row with the highest hyperscore in each group
    transform!(PSMs, AsTable(:) => ByRow(psm -> MS_TABLES[psm[:ms_file_idx]][:retentionTime][psm[:scan_idx]]) => :retention_time)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getSeq(ptable.id_to_pep[psm[:pep_idx]])) => :sequence)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getMZ(getPrecursor(ptable,psm[:precursor_idx]))) => :precursor_mz)
    transform!(PSMs, AsTable(:) => ByRow(psm -> join(sort(collect(getProtNamesFromPepSeq(ptable, psm[:sequence]))), "|")) => :protein_name)
    sort!(PSMs, [:retention_time])
    PSMs
end