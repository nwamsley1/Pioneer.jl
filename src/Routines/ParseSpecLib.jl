# src/ParseSpecLib.jl
"""
    nestedLibrarySort!(spec_lib::BasicEmpiricalLibrary; rt_bin_tol::Float64=0.1)

Sort the library DataFrame in a nested fashion, first by retention time (irt) and then by mass-to-charge ratio (prec_mz) 
within retention time bins of width `rt_bin_tol`.

# Arguments
- `spec_lib::BasicEmpiricalLibrary`: Spectral library containing a DataFrame with `:irt` and `:prec_mz` columns
- `rt_bin_tol::Float64=0.1`: Width of retention time bins within which to sort by prec_mz

# Details
The function performs a two-level sort:
1. First sorts the entire DataFrame by retention time (irt)
2. Then within each retention time bin of width `rt_bin_tol`, sorts by mass-to-charge ratio (prec_mz)

This nested sorting approach maintains retention time as the primary ordering while ensuring 
mass-to-charge ratio ordering within local retention time windows.
"""
function nestedLibrarySort!(spec_lib::BasicEmpiricalLibrary; rt_bin_tol::Float64=0.1)
    df = getDF(spec_lib)  # Get the DataFrame from the library
    
    # Initial sort by retention time
    sort!(df, :irt)
    
    # Initialize indices for RT bin processing
    start_idx = 1
    start_irt = first(df[!, :irt])
    
    # Process each RT bin
    for pid in 1:nrow(df)
        stop_irt = df[pid, :irt]
        if ((stop_irt - start_irt) > rt_bin_tol) && (pid > start_idx)
            # Sort the previous bin by multiple columns
            sort!(@view(df[start_idx:(pid-1), :]), 
                    [:prec_mz, :precursor_idx, :library_intensity], 
                    rev=[false, false, true])

            # Reset start position for next bin
            start_idx = pid
            start_irt = stop_irt
        end
    end
    
    # Sort final bin by prec_mz
    sort!(@view(df[start_idx:end, :]), :prec_mz)
    
    return nothing
end

"""
    ParseSpecLib(params_path::String)

Main function to build a Pioneer formatted spectral library from an empirical librairy.

Parameters:
- params_path: Path to JSON configuration file containing library building parameters

Output:
- Generates a spectral library in the specified output directory
- Creates a detailed log file with timing and performance metrics
- Returns nothing
"""
function ParseSpecLib(path; rt_bin_tol = 1.0)
    #"../Data/SPEC_LIBS/tag6_feb142025/tag6_comet_lib.tsv"
    #csv_file = CSV.File("../Data/SPEC_LIBS/tag6_feb142025/tag6_comet_lib.tsv")
    #seqs = ["YINHSCAPNCVAEVVTFDK","IINEPTAAAIAYGLDR"]
    #return DataFrame([pep_frag for 
    #    pep_frag in csv_file if pep_frag.PeptideSequence ∈ seqs])
    # 
    #println("hello to da world!!")
    #test_lib = BasicEmpiricalLibrary("/Users/nathanwamsley/Desktop/testdf.csv")
    test_lib = BasicEmpiricalLibrary(path)
    nestedLibrarySort!(test_lib, rt_bin_tol = rt_bin_tol)
    return test_lib
end


function parseLib(speclib::BasicEmpericalLibrary, speclib_dir::String)
    speclib_df = getDF(speclib)
    n_precursors = maximum(speclib_df[!,:precursor_idx])
    n_frags = UInt64(size(speclib_df, 1))
    pid_to_fid = zeros(UInt64, n_precursors) #Maps precursor idx to fragment idx range
    prec_idx = one(UInt32) #Precursor idx of current precursor
    unique_base_pep_seqs = Set{String}()
    #Initialize precursor columns
    proteome_identifiers = Vector{String}(undef, n_precursors)
    accession_numbers = Vector{String}(undef, n_precursors)
    sequence = Vector{String}(undef, n_precursors)
    structural_mods = Vector{Union{Missing, String}}(undef, n_precursors)
    isotopic_mods = Vector{Union{Missing, String}}(undef, n_precursors)
    prec_charge = zeros(UInt8, n_precursors)
    collision_energy = zeros(Float32, n_precursors)
    is_decoy = zeros(Bool, n_precursors)
    entrapment_group_id = zeros(UInt8, n_precursors)
    base_pep_id = zeros(UInt32, n_precursors)
    prec_mz = zeros(Float32, n_precursors)
    seq_length = zeros(UInt8, n_precursors)
    missed_cleavages = zeros(UInt8, n_precursors)
    irt = zeros(Float32, n_precursors)
    sulfur_count = zeros(UInt8, n_precursors)
    old_prec_idx = first(speclib_df[!,:precursor_idx])
    pid_to_fid[one(UInt32)] = one(UInt64)
    for frag_idx in ProgressBar(range(one(UInt64), UInt64(n_frags)))
        prec_idx = speclib_df[frag_idx, :precursor_idx]
        #Encountered a new precursor
        if prec_idx .!= old_prec_idx
            old_prec_idx = prec_idx 
            start_idx = frag_idx
            #Index for first fragment corresponding to the precursor
            pid_to_fid[prec_idx] = start_idx
            proteome_identifiers[prec_idx] = getProteomeId(speclib, frag_idx)
            accession_numbers[prec_idx] = getProteinGroupId(speclib, frag_idx)
            sequence[prec_idx] = getSequence(speclib, frag_idx)
            structural_mods[prec_idx] = parseStructuralMods(speclib, frag_idx)
            isotopic_mods[prec_idx] = parseIsotopicMods(speclib, frag_idx)
            prec_charge[prec_idx] = getPrecCharge(speclib, frag_idx)
            collision_energy[prec_idx] = getCollisionEnergy(speclib, frag_idx)
            is_decoy[prec_idx] = getIsDecoy(speclib, frag_idx)
            entrapment_group_id[prec_idx] = getEntrapmentGroupIdx(speclib, frag_idx)
            if getSequence(speclib, frag_idx) ∉ unique_base_pep_seqs
                base_pep_id += one(UInt32)
                push!(unique_base_pep_seqs, getSequence(speclib, frag_idx))
            end
            base_pep_id[prec_idx] = base_pep_id
            prec_mz[prec_idx] = getPrecMz(speclib, frag_idx)
            seq_length[prec_idx] = getLength(speclib, frag_idx)
            missed_cleavages[prec_idx] = getMissedCleavages(speclib, frag_idx) #will need to have specific and en
            irt[prec_idx] = getIrt(speclib, frag_idx)
            sulfur_count[prec_idx] = getPrecursorSulfurCount(speclib, frag_idx)
        end
        #Otherwise, this precursor has already been encountered so ignore
    end
    precursors_df = DataFrame(
        (
            proteome_identifiers = proteome_identifiers,
            accession_numbers = accession_numbers,
            sequence = sequence,
            structural_mods = structural_mods,
            isotopic_mods = isotopic_mods,
            prec_charge = prec_charge,
            collision_energy = collision_energy,
            is_decoy = is_decoy,
            entrapment_group_id = entrapment_group_id,
            base_pep_id = base_pep_id,
            mz = prec_mz,
            length = seq_length,
            missed_cleavages = missed_cleavages,
            irt = irt,
            sulfur_count = sulfur_count,
            isotopic_mods
        )
    )
    Arrow.write(
        joinpath(speclib_dir, "precursors_table.arrow"),
        precursors_df
    )
    Arrow.write(
        joinpath(speclib_dir, "prec_to_frag.arrow"),
        DataFrame(
            (
                start_idx = pid_to_fid,
            )
        )
    )
    #Initialize Fragment columns
    frag_mz = zeros(Float32, n_frags)
    frag_intensity = zeros(Float16, n_frags)
    ion_type = zeros(UInt16, n_frags)
    is_y = zeros(Bool, n_frags)
    is_b = zeros(Bool, n_frags)
    is_p = zeros(Bool, n_frags)
    is_axcz = zeros(Bool, n_frags)
    has_neutral_diff = zeros(Bool, n_frags)
    frag_series_number = zeros(UInt8, n_frags)
    frag_charge = zeros(UInt8, n_frags)
    internal = zeros(Bool, n_frags)
    immonium = zeros(Bool, n_frags)
    internal_ind = Vector{Tuple{UInt8, UInt8}}(undef, n_frags)
    frag_sulfur_count = zeros(UInt8, n_frags)
    for frag_idx in range(one(UInt64), n_frags)
        frag_mz[frag_idx] = getFragMz(speclib, frag_idx)
        frag_intensity[frag_idx] = getFragIntensity(speclib, frag_idx)
        ion_type[frag_idx] = getIonType(speclib, frag_idx)
        is_y[frag_idx] = getIsY(speclib, frag_idx)
        is_b[frag_idx] = getIsB(speclib, frag_idx)
        is_p[frag_idx] = getIsP(speclib, frag_idx)
        is_axcz[frag_idx] = getAXCZ(speclib, frag_idx)
        has_neutral_diff[frag_idx] = getNeutralDiff(speclib, frag_idx)
        frag_series_number[frag_idx] = getFragIndex(speclib, frag_idx)
        frag_charge[frag_idx] = getFragCharge(speclib, frag_idx)
        internal[frag_idx] = isInternal(speclib, frag_idx)
        immonium[frag_idx] = isImmonium(speclib, frag_idx)
        internal_ind[frag_idx] = getInternalInd(speclib, frag_idx)
        frag_sulfur_count[frag_idx] = getFragSulfurCount(speclib, frag_idx)
    end
    Arrow.write(
        joinpath(speclib_dir, "prec_to_frag.arrow"),
        DataFrame(
            (
                start_idx = pid_to_fid,
            )
        )
    )
end

#=
using Random, DataFrames

# Define unique precursors with their properties
precursors = [
    # (ModifiedPeptide, PrecursorCharge, PrecursorMz, RT)
    ("PEPT[+80]IDE", 2, 500.5, 10.5),
    ("SAMPL[+16]E", 3, 600.7, 9.2),
    ("TEST", 2, 800.2, 20.1)
]

# Create multiple fragments for each precursor
fragments = []
for (pep, charge, mz, rt) in precursors
    if pep == "PEPT[+80]IDE"
        push!(fragments, 
            (pep, charge, mz, rt, 1000.0, 200.1, "y", 1, 2),
            (pep, charge, mz, rt, 500.0, 300.2, "b", 1, 3),
            (pep, charge, mz, rt, 100.0, 400.3, "y", 1, 4)
        )
    elseif pep == "SAMPL[+16]E"
        push!(fragments, 
            (pep, charge, mz, rt, 800.0, 250.1, "y", 1, 2),
            (pep, charge, mz, rt, 200.0, 350.2, "b", 1, 3)
        )
    else  # TEST
        push!(fragments,
            (pep, charge, mz, rt, 600.0, 450.3, "y", 1, 4)
        )
    end
end

# Shuffle the fragments
fragments = shuffle(fragments)

# Create DataFrame
test_data = DataFrame(
    PrecursorMz = [f[3] for f in fragments],
    ModifiedPeptide = [f[1] for f in fragments],
    PrecursorCharge = [f[2] for f in fragments],
    Tr_recalibrated = [f[4] for f in fragments],
    LibraryIntensity = [f[5] for f in fragments],
    ProductMz = [f[6] for f in fragments],
    FragmentType = [f[7] for f in fragments],
    FragmentCharge = [f[8] for f in fragments],
    FragmentSeriesNumber = [f[9] for f in fragments]
)
        =#