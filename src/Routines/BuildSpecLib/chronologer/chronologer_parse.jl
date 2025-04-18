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
        :rt => :irt,
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



    @info "Shuffling decoy iRTs..."
    # Set iRT to be equal to the target
    # 1) Extract target rows, rename :irt => :target_irt
    targets_df = filter(row -> row.decoy == false, precursors_df)
    #select!(targets_df, [:irt, :pair_id])

    println("unique pair_ids: ", length(unique(precursors_df.pair_id)) ,"\n")
    duplicates = 0

    ######################################################
    # TODO
    # This checks for duplicate target pairs
    for sub in groupby(targets_df, [:pair_id])
        if nrow(sub) > 1
            #println(sub)
            duplicates += 1
        end
    end
    println("duplicates: ", duplicates)
    ######################################################

    targets_df = combine(groupby(targets_df, [:pair_id])) do sub
        (irt = mean(sub.irt),)  # example aggregator
    end

   
    

    rename!(targets_df, :irt => :target_irt)

    # 2) Left-join in place by :base_pep_id, :mz
    #    This adds a :target_irt column to precursors_df
    leftjoin!(
        precursors_df,
        targets_df,
        on = [:pair_id]
    )

    ######################################################
    # TODO
    # This checks for missing target pairs
    num_missing = 0
    # 3) Update the decoy rows' :irt from :target_irt if available
    for r in eachrow(precursors_df)
        if r.decoy == true && ismissing(r.target_irt)
            #println(r)
            num_missing += 1
        end
    end
    println("num missing: ", num_missing, "\n")
    ######################################################

    #rename!(targets_df, :foobar => :target_irt)

    # 4) Drop the temporary :target_irt column
    #select!(precursors_df, Not(:target_irt))



    precursors_target = subset(precursors_df, [:decoy] => ByRow((d) -> !d); view = true)
    precursors_decoy = subset(precursors_df, [:decoy] => ByRow((d) -> d); view = true)
   
    # 1) Get overall min/max from the target distribution
    target_min_irt = minimum(precursors_target.irt)
    target_max_irt = maximum(precursors_target.irt)

    # 2) Build the target_counts map
    target_counts = build_target_counts(precursors_target; mz_bin_size=2.0, irt_bin_size=Float64(rt_bin_tol), irt_offset=target_min_irt)

    println(target_min_irt, " ", rt_bin_tol)

    # 3) Shift decoys in-place
    shift_decoys!(precursors_decoy, target_counts;
                  mz_bin_size=2.0, irt_bin_size=Float64(rt_bin_tol),
                  irt_bin_window=(-1, 1),
                  target_min_irt=target_min_irt,
                  target_max_irt=target_max_irt)





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












"""
    bin_value(x, binsize)

Returns the integer bin index for `x` using `binsize`.
E.g., with binsize=2.0 and x=4.3, returns floor(Int, 4.3/2.0) = 2.
"""
function bin_value(x, binsize; min_offset=0.0)
    return floor(Int, (x-min_offset) / binsize)
end


"""
    build_target_counts(df; mz_bin_size=2.0, irt_bin_size=1.0)

Bins the `df` DataFrame's columns `:mz` and `:irt`, returning a
Dict{(Int,Int), Int} mapping (mz_bin, irt_bin) => count of how many rows.
"""
function build_target_counts(df; mz_bin_size=2.0, irt_bin_size=1.0, irt_offset=0.0)
    target_counts = Dict{Tuple{Int,Int}, Tuple{Int,Int}}()
    for row in eachrow(df)
        mzb  = bin_value(row.mz,  mz_bin_size)
        irtb = bin_value(row.irt, irt_bin_size; min_offset=irt_offset)
        key = (mzb, irtb)
        new_count = first(get(target_counts, key, (0,0))) + 1
        target_counts[key] = (new_count, new_count) 
    end
    return target_counts
end


"""
    shift_decoys!(
       decoys_df, 
       target_counts; 
       mz_bin_size=2.0, 
       irt_bin_size=1.0,
       irt_bin_window=(-4, -3, 3, 4),
       target_min_irt,
       target_max_irt
    )

For each row in `decoys_df`:
1. Compute its (mz_bin, irt_bin).
2. Gather candidate bins that share the same mz_bin but whose irt_bin is in
   (decoy_irt_bin ± one of irt_bin_window).
3. Summation of all target_counts for these bins => total_count.
   - If total_count > 0, randomly pick one bin weighted by the target_counts.
     * Shift the decoy iRT by (chosen_irt_bin - decoy_irt_bin)*irt_bin_size.
     * Decrement that bin’s count by 1 (but do not let it go below 1).
   - If total_count == 0, fallback:
     * Sample *equally* from the same set of possible bin diffs in `irt_bin_window`.
     * Only consider those diffs that keep final iRT in [target_min_irt, target_max_irt].
     * If none are valid, do nothing or handle as desired.
"""
function shift_decoys!(
    decoys_df::SubDataFrame,
    target_counts::Dict{Tuple{Int,Int},Tuple{Int,Int}}; # value is (original total, remaining total)
    mz_bin_size::Float64 = 2.0,
    irt_bin_size::Float64 = 1.0,
    irt_bin_window::NTuple{2, Int} = (-3, 3),
    target_min_irt::Float32,
    target_max_irt::Float32
)
    

    for i in ProgressBar(shuffle(1:nrow(decoys_df)))
        row = decoys_df[i, :]

        decoy_mz_bin  = bin_value(row.mz,  mz_bin_size)
        decoy_irt_bin = bin_value(row.irt, irt_bin_size; min_offset=target_min_irt)
        #pair_irt_bin = bin_value(row.target_irt, irt_bin_size; min_offset=target_min_irt)

        # Filter to the bins that actually exist in the dictionary
        candidate_keys = filter(k -> (k[1] == decoy_mz_bin) && ((k[2] < decoy_irt_bin) || (k[2] > decoy_irt_bin)) && 
                                        (last(target_counts[k]) > 0),# && (k[2] != pair_irt_bin),
                                keys(target_counts))

        # If we found none in the dict, total_count is effectively 0
        if isempty(candidate_keys)
            # Fallback #1: sampling from original counts if they exist
            candidate_keys = filter(k -> (k[1] == decoy_mz_bin) && ((k[2] < decoy_irt_bin) || (k[2] > decoy_irt_bin)),# && (k[2] != pair_irt_bin),
                                keys(target_counts))

            if isempty(candidate_keys)
                # Fallback #2: sampling from possible_irt_bins
                fallback_shift!(row, irt_bin_size, irt_bin_window, target_min_irt, target_max_irt)
                continue
            else
                # Filter for closest
                #min_dist = minimum([abs(k[2]-decoy_irt_bin) for k in candidate_keys])
                #candidate_keys = filter(k -> abs(k[2]-decoy_irt_bin) == min_dist, candidate_keys)
                # Sum up the original total counts 
                max_distance = 1 + maximum([abs(k[2]-decoy_irt_bin) for k in candidate_keys])
                weights = [first(target_counts[k]) * (max_distance - abs(k[2]-decoy_irt_bin))^2 for k in candidate_keys]
            end
        else
            # Filter for closest
            min_dist = minimum([abs(k[2]-decoy_irt_bin) for k in candidate_keys])
            candidate_keys = filter(k -> abs(k[2]-decoy_irt_bin) == min_dist, candidate_keys)
            # Sum up the counts remaining
            max_distance = 1 + maximum([abs(k[2]-decoy_irt_bin) for k in candidate_keys])
            weights = [last(target_counts[k]) * (max_distance - abs(k[2]-decoy_irt_bin))^2 for k in candidate_keys]
        end

       
        total_count = sum(weights)

        # Weighted random pick
        r = rand() * total_count
        cumulative = 0
        chosen_key = rand(candidate_keys)  # fallback
        for (j, k) in enumerate(candidate_keys)
            cumulative += weights[j]
            if cumulative >= r
                chosen_key = k
                break
            end
        end
        chosen_irt_bin = chosen_key[2]

        # Shift iRT
        bin_diff = chosen_irt_bin - decoy_irt_bin
        row.irt  = row.irt + bin_diff*irt_bin_size

        # Decrement the remaining bin count
        old_total, old_remaining = target_counts[chosen_key]
        target_counts[chosen_key] = (old_total, old_remaining-1)
    end
end

"""
    fallback_shift!(
      row, decoy_irt_bin, irt_bin_size, irt_bin_window,
      target_min_irt, target_max_irt
    )

Helper function that picks an iRT shift *uniformly* from the possible bin offsets
given by `irt_bin_window`, but only among those that keep the final iRT in
[target_min_irt, target_max_irt].
"""
function fallback_shift!(
    row::DataFrameRow,
    irt_bin_size::Float64,
    irt_bin_window::NTuple{2,Int},
    target_min_irt::Float32,
    target_max_irt::Float32
)
    candidate_shifts = Float64[]
    decoy_irt = row.irt

    # For each offset in irt_bin_window:
    for offset in irt_bin_window
        bin_diff = offset
        shifted_irt = decoy_irt + bin_diff * irt_bin_size
        if shifted_irt >= target_min_irt && shifted_irt <= target_max_irt
            push!(candidate_shifts, shifted_irt)
        end
    end

    row.irt = rand(candidate_shifts)
end