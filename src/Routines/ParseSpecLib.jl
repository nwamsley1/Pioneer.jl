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
function nestedLibrarySort!(spec_lib::BasicEmpiricalLibrary; rt_bin_tol::AbstractFloat=0.1)
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
    sort!(@view(df[start_idx:end, :]),  
        [:prec_mz, :precursor_idx, :library_intensity], 
        rev=[false, false, true])
    
    return nothing
end
#=
"""
    parse_mods(mod_string::String)::Vector{Tuple{String, Int}}

Parse modification string in the format "[Mod1]Res1[Mod2]Res2" and return a vector of 
(modification, position) tuples.

# Arguments
- `mod_string::String`: String containing modifications in brackets

# Returns
- Vector{Tuple{String, Int}}: Vector of (modification name, position) tuples

# Example
```julia
julia> parse_mods("_[Unimod:35]M[Unimod:4]PEPTIDE_")
2-element Vector{Tuple{String, Int}}:
 ("Unimod:35", 1)
 ("Unimod:4", 2)
"""
function parse_mods(mod_string::Union{String, Missing})::Vector{Tuple{String, Int}} if ismissing(mod_string) || isempty(mod_string) return Tuple{String, Int}[] end
    # Match pattern [ModName]Residue
    mod_pattern = r"\[([^\]]+)\]([A-Z])"

    # Find all matches and their positions
    matches = collect(eachmatch(mod_pattern, mod_string))

    # Convert matches to vector of (mod_name, position) tuples
    mods = Vector{Tuple{String, Int}}(undef, length(matches))
    for (i, m) in enumerate(matches)
        mod_name = m.captures[1]
        position = count(!=('_'), SubString(mod_string, 1, m.offset)) + 1
        mods[i] = (mod_name, position)
    end

    return mods
end
=#

#=
sequence = "PEPMTIDME"

# Create the variable modifications vector
var_mods = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}()
push!(var_mods, (p=r"M", r="Unimod:35"))
=#

# Function to match variable modifications (from your code)
function matchVarMods(sequence::String, var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}})
    var_mod_matches = Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}}()
    for mod in var_mods
        for mod_match in eachmatch(mod[:p], sequence)
            push!(var_mod_matches, (regex_match=mod_match, name=mod[:r]))
        end
    end
    return var_mod_matches
end

"""
    create_precursor_idx!(df::DataFrame)

Creates a UInt32 precursor_idx column that uniquely identifies each unique combination of
sequence, structural_mods, isotopic_mods, and prec_charge.

# Arguments
- `df::DataFrame`: DataFrame containing required columns

# Effects
- Adds or updates 'precursor_idx' column with unique identifiers
"""
function create_precursor_idx!(df::DataFrame)
    # Get unique combinations and assign indices
    unique_combinations = unique(df[!, [:is_decoy, :sequence, :structural_mods, :isotopic_mods, :prec_charge]])
    precursor_dict = Dict(
        tuple(row...) => UInt32(i) 
        for (i, row) in enumerate(eachrow(unique_combinations))
    )
    
    # Create the precursor_idx column
    df[!, :precursor_idx] = UInt32[
        precursor_dict[tuple(row.is_decoy, row.sequence, row.structural_mods, row.isotopic_mods, row.prec_charge)]
        for row in eachrow(df)
    ]
    
    # Reorder columns to put precursor_idx first if it's not already
    if first(names(df)) != "precursor_idx"
        select!(df, :precursor_idx, Not(:precursor_idx))
    end
end

function parseStructralModsFromLib!(speclib::BasicEmpiricalLibrary)
    function _parseStructuralModsFromLib(
        modified_seqs::AbstractVector{<:AbstractString}
    )
        structural_mods = length(modified_seqs) > 0 ? Vector{String}(undef, length(modified_seqs)) : Vector{String}()
        for (i, seq) in enumerate(modified_seqs)
            structural_mods[i] = parseEmpiricalLibraryMods(seq)
        end
        return structural_mods  
    end
    speclib.libdf[!,:structural_mods] = _parseStructuralModsFromLib(speclib.libdf[!,:modified_sequence])
end

"""
    addIsoModsToLib!(speclib::BasicEmpiricalLibrary,
                    mod_key::String,
                    mod_channel::NamedTuple{(:channel, :mass), Tuple{String, Float32}})

Process a spectral library dataframe to add isotopic modifications for a given modification type
and channel. Creates a new channel-specific dataframe for entries containing the target modification.

# Arguments
- `speclib::BasicEmpiricalLibrary`: Spectral library containing precursor information
- `mod_key::String`: Target modification to add isotopic channel for (e.g., "tag6")
- `mod_channel::NamedTuple{(:channel, :mass)}`: Isotopic channel information

# Returns
- DataFrame containing entries that have the target modification, with updated isotopic_mods column

# Example
```julia
# For a library with tag6 modifications:
channel = (channel="d0", mass=0.0f0)
channel_df = addIsoModsToLib!(speclib, "tag6", channel)
```
"""
function addIsoModsToLib(
    speclib::BasicEmpiricalLibrary,
    mod_key::String,
    mod_channel::NamedTuple{(:channel, :mass), Tuple{String, Float32}}
)
    # Create a copy of the library dataframe
    libdf = copy(speclib.libdf)
    
    # Ensure isotopic_mods column exists, initialize if needed
    if !hasproperty(libdf, :isotopic_mods)
        libdf[!, :isotopic_mods] .= ""
    end
    
    # Filter rows that contain the target modification
    has_target_mod = map(row -> contains(row.structural_mods, mod_key), eachrow(libdf))
    channel_df = libdf[has_target_mod, :]
    
    # Process each row to add isotopic modifications
    for row in eachrow(channel_df)
        row.isotopic_mods = addIsoMods(
            row.isotopic_mods,
            row.structural_mods,
            mod_key,
            mod_channel
        )
    end
    
    return channel_df
end

"""
    ParseSpecLib(params_path::String)

Main function to build a Pioneer formatted spectral library from an empirical library.

Parameters:
- params_path: Path to JSON configuration file containing library building parameters

Output:
- Generates a spectral library in the specified output directory
- Creates a detailed log file with timing and performance metrics
- Returns the built spectral library
"""
function ParseSpecLib(params_path::String)

    params = checkParseSpecLibParams(params_path)
    
    # Extract parameters
    lib_path = params["library_params"]["input_lib_path"]
    rt_bin_tol = Float32(params["library_params"]["rt_bin_tol"])
    output_path = params["library_params"]["output_lib_path"]
    generate_decoys = params["library_params"]["generate_decoys"]
    generate_entrapment = params["library_params"]["generate_entrapment"]
    entrapment_groups = params["library_params"]["entrapment_groups"]
    #"../Data/SPEC_LIBS/tag6_feb142025/tag6_comet_lib.tsv"
    #test_lib = ParseSpecLib("/Users/nathanwamsley/Data/SPEC_LIBS/tag6_feb142025/tag6_comet_lib.tsv").libdf;

    # Initialize library
    test_lib = BasicEmpiricalLibrary(lib_path)

    # Process structural modifications
    parseStructralModsFromLib!(test_lib)
    test_lib.libdf[!,:entrapment_group_id] = zeros(UInt8, size(test_lib.libdf, 1))

    # Generate entrapment sequences if requested
    if generate_entrapment
        for i in 1:entrapment_groups
            getShuffledEntrapmentSeqs!(test_lib, UInt8(i))
        end
    end
    
    # Generate reversed decoys if requested
    if generate_decoys
        getRevDecoys!(test_lib)
    end

    # Process modification channels
    mod_channels = Dict{String, Vector{NamedTuple{(:channel, :mass), Tuple{String, Float32}}}}()
    for mod_group in params["isotope_mod_groups"]
        name = mod_group["name"]
        channels = [(channel = ch["channel"], mass = Float32(ch["mass"])) for ch in mod_group["channels"]]
        mod_channels[name] = channels
    end


    # Process modification channels
    mod_channels = Dict{String, Vector{NamedTuple{(:channel, :mass), Tuple{String, Float32}}}}()
    for mod_group in params["isotope_mod_groups"]
        name = mod_group["name"]
        channels = [(channel = ch["channel"], mass = Float32(ch["mass"])) for ch in mod_group["channels"]]
        mod_channels[name] = channels
    end

    # Build isotope modifications dictionary
    iso_mods_dict = Dict{String, Dict{String, Float32}}()
    for mod_group in params["isotope_mod_groups"]
        name = mod_group["name"]
        iso_mods_dict[name] = Dict{String, Float32}()
        for ch in mod_group["channels"]
            iso_mods_dict[name][ch["channel"]] = Float32(ch["mass"])
        end
    end
    
    # Build structural modification masses dictionary
    structural_mod_to_mass = Dict{String, Float32}()
    for (Parwsemass, name) in zip(
        params["fixed_mods"]["mass"],
        params["fixed_mods"]["name"]
    )
        structural_mod_to_mass[name] = Float32(mass)
    end
    println("structural_mod_to_mass $structural_mod_to_mass")
   # Apply isotopic modifications
   channel_dfs = []
   for (mod_key, channels) in pairs(mod_channels)
       for channel in channels
           channel_df = addIsoModsToLib(test_lib, mod_key, channel)
           push!(channel_dfs, channel_df)
       end
   end

    # If channels were processed, update the DataFrame
    if !isempty(channel_dfs)
        channels_df = vcat(channel_dfs...)
        empty!(test_lib.libdf)
        append!(test_lib.libdf, channels_df)
    end

    test_lib.libdf[!,:frag_sulfur_count] = zeros(UInt8, size(test_lib.libdf, 1))
    test_lib.libdf[!,:prec_sulfur_count] = zeros(UInt8, size(test_lib.libdf, 1))
    #Now need to recalculate masses for precursors and fragments with the new modifications 
    # Create precursor indices
    create_precursor_idx!(test_lib.libdf)
    
    # Calculate m/z and sulfur count
    mods_to_sulfur_diff = Dict{String, Int8}()
    for mod_group in get(params, "sulfur_mod_groups", [])
        mods_to_sulfur_diff[mod_group["name"]] = Int8(mod_group["sulfur_count"])
    end

    calculate_mz_and_sulfur_count!(
        test_lib.libdf, 
        structural_mod_to_mass,
        iso_mods_dict,
        mods_to_sulfur_diff
    )

    # Update precursor indices after m/z calculation
    create_precursor_idx!(test_lib.libdf)
    
    # Sort library by retention time and then by precursor_mz within retention time bins
    nestedLibrarySort!(test_lib, rt_bin_tol=rt_bin_tol)

    # Ensure output directory exists
    if !isdir(output_path)
        mkpath(output_path)
    end

    # Parse library to Pioneer format
    parseLib(test_lib, output_path)

    buildPionLib(
        output_path,
        UInt8(params["library_params"]["y_start_index"]),
        UInt8(params["library_params"]["y_start"]),
        UInt8(params["library_params"]["b_start_index"]),
        UInt8(params["library_params"]["b_start"]),
        params["library_params"]["include_p_index"],
        params["library_params"]["include_p"],
        params["library_params"]["include_isotope"],
        params["library_params"]["include_immonium"],
        params["library_params"]["include_internal"],
        params["library_params"]["include_neutral_diff"],
        UInt8(params["library_params"]["max_frag_charge"]),
        UInt8(params["library_params"]["max_frag_rank"]),
        Float32(params["library_params"]["min_frag_intensity"]),
        UInt8.(params["library_params"]["rank_to_score"]),
        FragBoundModel(
            ImmutablePolynomial(zero(Float32)),
            ImmutablePolynomial(Float32(10000.0f0)) 
        ),
        Float32(params["library_params"]["frag_bin_tol_ppm"]),
        Float32(params["library_params"]["rt_bin_tol"]),
        InstrumentSpecificModel(params["library_params"]["instrument_type"])
    )   

    return test_lib
end


function parseLib(speclib::BasicEmpiricalLibrary, speclib_dir::String)
    speclib_df = getDF(speclib)
    n_precursors = maximum(speclib_df[!,:precursor_idx]) #Unique precursors in the data fraqme
    n_frags = UInt64(size(speclib_df, 1)) #Total number of fragments
    pid_to_fid = zeros(UInt64, n_precursors+1) #Maps precursor idx to fragment idx range
    prec_idx = zero(UInt32) #Precursor idx of current precursor
    prev_prec_idx = zero(UInt32) #Previous precursor idx
    #Use to normalize fragment intensities to base peak 
    max_frag_intensities = zeros(Float16, n_precursors)
    new_precursor_idxs = zeros(UInt32, n_frags)

    _base_pep_id_ = one(UInt32) #Base peptide id of current precursor
    base_pep_ids = zeros(UInt32, n_precursors) #Base peptide id of each precursor
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
    prec_mz = zeros(Float32, n_precursors)
    seq_length = zeros(UInt8, n_precursors)
    missed_cleavages = zeros(UInt8, n_precursors)
    irt = zeros(Float32, n_precursors)
    sulfur_count = zeros(UInt8, n_precursors)

    for frag_idx in ProgressBar(range(one(UInt64), UInt64(n_frags)))
        current_prec_idx = speclib_df[frag_idx, :precursor_idx]
        #Encountered a new precursor
        if current_prec_idx != prev_prec_idx 
            prev_prec_idx = current_prec_idx
            prec_idx += one(UInt32)
            #Index for first fragment corresponding to the precursor
            pid_to_fid[prec_idx] = frag_idx
            proteome_identifiers[prec_idx] = getProteomeId(speclib, frag_idx)
            accession_numbers[prec_idx] = getProteinGroupId(speclib, frag_idx)
            sequence[prec_idx] = getSequence(speclib, frag_idx)
            structural_mods[prec_idx] = getStructuralMods(speclib, frag_idx)
            isotopic_mods[prec_idx] = parseIsotopicMods(speclib, frag_idx)
            prec_charge[prec_idx] = getPrecCharge(speclib, frag_idx)
            collision_energy[prec_idx] = getCollisionEnergy(speclib, frag_idx)
            is_decoy[prec_idx] = getIsDecoy(speclib, frag_idx)
            entrapment_group_id[prec_idx] = getEntrapmentGroupIdx(speclib, frag_idx)
            if getSequence(speclib, frag_idx) ∉ unique_base_pep_seqs
                _base_pep_id_ += one(UInt32)
                push!(unique_base_pep_seqs, getSequence(speclib, frag_idx))
            end
            base_pep_ids[prec_idx] = _base_pep_id_
            prec_mz[prec_idx] = getPrecMz(speclib, frag_idx)
            seq_length[prec_idx] = getSeqLength(speclib, frag_idx)
            missed_cleavages[prec_idx] = getMissedCleavages(speclib, frag_idx) #will need to have specific and en
            irt[prec_idx] = getIrt(speclib, frag_idx)
            sulfur_count[prec_idx] = getPrecSulfurCount(speclib, frag_idx)
        end
        new_precursor_idxs[frag_idx] = prec_idx
        if speclib_df[frag_idx, :library_intensity] > max_frag_intensities[prec_idx]
            max_frag_intensities[prec_idx] = speclib_df[frag_idx, :library_intensity]
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
            base_pep_id = base_pep_ids,
            mz = prec_mz,
            length = seq_length,
            missed_cleavages = missed_cleavages,
            irt = irt,
            sulfur_count = sulfur_count
        )
    )
    Arrow.write(
        joinpath(speclib_dir, "precursors_table.arrow"),
        precursors_df
    )
    pid_to_fid[end] = n_frags + 1
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
        max_frag_intensity = max_frag_intensities[new_precursor_idxs[frag_idx]]
        frag_mz[frag_idx] = getFragMz(speclib, frag_idx)
        frag_intensity[frag_idx] = Float16(getFragIntensity(speclib, frag_idx)/max_frag_intensity)
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
        joinpath(speclib_dir, "fragments_table.arrow"),
        DataFrame(
            (
                mz = frag_mz,
                intensity = frag_intensity,
                ion_type = ion_type,
                is_y = is_y,
                is_b = is_b,
                is_p = is_p,
                is_axcz = is_axcz,
                has_neutral_diff = has_neutral_diff,
                fragment_index = frag_series_number,
                charge = frag_charge,
                isotope = zeros(UInt8, n_frags),
                is_internal = internal,
                is_immonium = immonium,
                internal_ind = internal_ind,
                sulfur_count = frag_sulfur_count
            )
        )
    )
end