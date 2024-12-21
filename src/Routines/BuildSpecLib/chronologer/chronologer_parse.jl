# src/chronologer/chronologer_parse.jl

"""
    parse_chronologer_output(
        path_to_precursors::String,
        pion_lib_dir::String,
        mods_to_sulfur_diff::Dict{String, Int8},
        iso_mod_to_mass::Dict{String, Float32},
        isotope_mods_groups::Vector{Any},
        rt_bin_tol::AbstractFloat
    )::String

Parse chronologer output and prepare it for Pioneer library format.

Parameters:
- path_to_precursors: Path to Arrow file containing chronologer results
- pion_lib_dir: Output directory for processed files
- mods_to_sulfur_diff: Maps modification names to sulfur count changes
- iso_mod_to_mass: Maps isotope modification names to mass shifts
- isotope_mods_groups: Configuration for isotope labeling groups
- rt_bin_tol: RT bin tolerance for sorting

Returns:
- String: Path to processed precursors Arrow file
"""
function parse_chronologer_output(
    path_to_precursors::String,
    pion_lib_dir::String,
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32},
    isotope_mods_groups::Vector{Any},
    rt_bin_tol::AbstractFloat
)::String

    function countSulfurs( plain_sequence::AbstractString, 
                            mods::String,
                            mods_to_sulfur_diff::Dict{String, Int8})::Int8
        sulfur_count = zero(UInt8)
        for aa in plain_sequence
            sulfur_count += (aa=='C')|(aa=='M')
        end

        for mod in parseMods(mods)
            if haskey(mods_to_sulfur_diff, getModName(mod.match))
                n_sulfur = mods_to_sulfur_diff[mod_string]
                seq_idx_to_sulfur[getModIndex(mod.match)] += n_sulfur
                sulfur_count += n_sulfur
            end
        end
        return sulfur_count
    end

    function countSulfurs( plain_sequence::AbstractString, 
                            mods::Missing,
                            mods_to_sulfur_diff::Dict{String, Int8})::Int8
        sulfur_count = zero(UInt8)
        for aa in plain_sequence
            sulfur_count += (aa=='C')|(aa=='M')
        end
        return sulfur_count
    end
    # Read chronologer output
    precursors_df = DataFrame(Arrow.Table(path_to_precursors))

    # Rename columns to match Pioneer format
    rename!(precursors_df, Dict(
        :iRT => :irt,
        :upid => :proteome_identifiers
    ))

    # Add sequence length column
    precursors_df[!, :length] = zeros(UInt8, nrow(precursors_df))
    for i in 1:nrow(precursors_df)
        precursors_df[i, :length] = UInt8(length(precursors_df[i, :sequence]))
    end

    # Add sulfur count column
    precursors_df[!, :sulfur_count] = zeros(UInt8, nrow(precursors_df))
    for i in 1:nrow(precursors_df)
        precursors_df[i, :sulfur_count] = countSulfurs(
            precursors_df[i, :sequence],
            precursors_df[i, :mods],
            mods_to_sulfur_diff
        )
    end

    # Initialize isotope modifications column
    precursors_df[!, :isotope_mods] = Vector{Union{Missing, String}}(missing, nrow(precursors_df))

    # Add isotope-modified precursors if specified
    println("Initial precursor count: ", nrow(precursors_df))
    precursors_df = addIsotopeModifiedPrecursors!(
        precursors_df,
        iso_mod_to_mass,
        isotope_mods_groups
    )
    println("Precursor count after isotope mods: ", nrow(precursors_df))

    ########
    #Sort 
    #Need same sort order as fragment index.
    #1) Sort entire data frame in ascending order of irt
    #2) Within irt bins of identical width, sort by mz
    sort!(precursors_df, :irt)
    start_idx, stop_idx = 1, 1
    start_irt, stop_irt = first(precursors_df[!,:irt]), first(precursors_df[!,:irt])

    # Sort by m/z within RT bins
    start_idx = 1
    curr_rt = first(precursors_df[!, :irt])
    
    #Get irt bins and sort by mz within these 
    for pid in range(1, size(precursors_df, 1))
        stop_idx = pid
        stop_irt = precursors_df[stop_idx,:irt]
        if ((stop_irt - start_irt) > rt_bin_tol) & (stop_idx > start_idx)
            stop_idx -= 1
            sort!(@view(precursors_df[start_idx:stop_idx,:]), :mz)
            start_idx = pid
            start_irt = precursors_df[pid,:irt]
        end
    end
    sort!(@view(precursors_df[start_idx:end,:]),:mz)

    # Write processed precursors to Arrow file
    precursors_arrow_path = joinpath(pion_lib_dir, "precursors.arrow")
    Arrow.write(precursors_arrow_path, precursors_df)

    return precursors_arrow_path
end