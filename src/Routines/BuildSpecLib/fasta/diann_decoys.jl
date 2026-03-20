# DIA-NN-style decoy generation
#
# Instead of shuffling sequences and predicting new spectra/RTs via deep learning,
# DIA-NN generates decoys by:
# 1. Mutating one AA near each terminus (fixed lookup table)
# 2. Shifting all fragment m/z by the resulting mass difference
# 3. Copying fragment intensities and retention times from the target
#
# Reference: DiaNN-1.8/src/diann.cpp, generate_decoy() at line 3590

"""
DIA-NN amino acid mutation table. Each AA maps to a different AA with similar
chemical properties but different mass. Used for decoy fragment m/z generation.
Source: diann.cpp line 1370
"""
const DIANN_MUTATION_TABLE = Dict{Char, Char}(
    'G' => 'L', 'A' => 'L', 'V' => 'L', 'L' => 'V', 'I' => 'V',
    'F' => 'L', 'M' => 'L', 'P' => 'L', 'W' => 'L', 'S' => 'T',
    'C' => 'S', 'T' => 'S', 'Y' => 'S', 'H' => 'S', 'K' => 'L',
    'R' => 'N', 'Q' => 'D', 'E' => 'Q', 'N' => 'E', 'D' => 'Q'
)

"""
    compute_diann_mutation_shifts(sequence, mod_positions) -> (n_shift, c_shift)

Compute the N-terminal and C-terminal mass shifts for DIA-NN-style decoy generation.

`sequence` is the bare amino acid string (no modifications).
`mod_positions` is a Set{Int} of 1-indexed positions that have modifications
(these positions are avoided when selecting the mutation site).

Returns `(n_shift::Float64, c_shift::Float64)` — the mass differences in Da.

Position selection matches DIA-NN exactly (diann.cpp lines 3626-3636, converted to 1-indexed):
- N-terminal: prefers position 2, then min(3, n-1), then 1, fallback 2
- C-terminal: prefers position n-1, then max(1, n-2), then n, fallback n-1
"""
function compute_diann_mutation_shifts(
    sequence::AbstractString,
    mod_positions::Set{Int} = Set{Int}()
)
    n = length(sequence)
    n < 3 && return (0.0, 0.0)  # too short to mutate
    chars = collect(sequence)

    # N-terminal position selection (DIA-NN: 0-indexed pos 1 → min(2,m+1) → 0 → 1)
    # In 1-indexed Julia: pos 2 → min(3, n-1) → 1 → 2
    n_pos = if 2 ∉ mod_positions
        2
    elseif min(3, n - 1) ∉ mod_positions
        min(3, n - 1)
    elseif 1 ∉ mod_positions
        1
    else
        2  # fallback
    end

    # C-terminal position selection (DIA-NN: 0-indexed pos m → max(0,m-1) → m+1 → m)
    # where m = n-2 (0-indexed). In 1-indexed: n-1 → max(1, n-2) → n → n-1
    c_pos = if (n - 1) ∉ mod_positions
        n - 1
    elseif max(1, n - 2) ∉ mod_positions
        max(1, n - 2)
    elseif n ∉ mod_positions
        n
    else
        n - 1  # fallback
    end

    # Compute mass shifts
    n_aa = chars[n_pos]
    c_aa = chars[c_pos]

    n_mutated = get(DIANN_MUTATION_TABLE, n_aa, n_aa)
    c_mutated = get(DIANN_MUTATION_TABLE, c_aa, c_aa)

    n_shift = AA_to_mass[n_mutated] - AA_to_mass[n_aa]
    c_shift = AA_to_mass[c_mutated] - AA_to_mass[c_aa]

    return (n_shift, c_shift)
end

"""
    create_diann_decoy_fragments(target_frags, target_range, n_shift, c_shift, decoy_prec_id)

Create decoy fragments by copying target fragments and shifting m/z values.
Y-ions shift by c_shift/charge, everything else by n_shift/charge.
Intensities are copied unchanged.

Returns a vector of the same fragment type as input.
"""
function create_diann_decoy_fragments(
    target_frags::AbstractVector{F},
    target_range::UnitRange,
    n_shift::Float64,
    c_shift::Float64,
    decoy_prec_id::UInt32
) where {F <: AltimeterFragment{Float32}}
    decoy_frags = Vector{F}(undef, length(target_range))
    for (j, i) in enumerate(target_range)
        frag = target_frags[i]
        charge = max(frag.frag_charge, UInt8(1))
        mz_shift = frag.is_y ? Float32(c_shift / charge) : Float32(n_shift / charge)

        decoy_frags[j] = F(
            decoy_prec_id,
            frag.mz + mz_shift,
            frag.intensity,
            frag.ion_type,
            frag.is_y,
            frag.is_b,
            frag.is_p,
            frag.is_isotope,
            frag.frag_charge,
            frag.ion_position,
            frag.prec_charge,
            frag.rank,
            frag.sulfur_count
        )
    end
    return decoy_frags
end

"""
    apply_diann_decoy_style!(lib_path::String)

Add DIA-NN-style decoys to a target-only Pioneer library.

For each target precursor, creates a decoy by:
1. Mutating terminal AAs (DIA-NN mutation table) to get a new sequence
2. Shifting target fragment m/z by the terminal mass differences
3. Copying target iRT and fragment intensities unchanged
4. Shifting precursor m/z by (n_shift + c_shift) / charge

The resulting library has targets followed by decoys, with proper
pairing (pair_id, partner_precursor_idx) and all required columns.

Modifies `detailed_fragments.jls`, `precursor_to_fragment_indices.jls`,
`precursors_table.arrow`, and rebuilds `partitioned_fragment_index.jls` in place.
"""
function apply_diann_decoy_style!(lib_path::String)
    # Load target-only precursors table
    prec_df = DataFrame(Arrow.Table(joinpath(lib_path, "precursors_table.arrow")))
    n_targets = nrow(prec_df)

    # Load existing fragments and index
    detailed_frags = load_detailed_frags(joinpath(lib_path, "detailed_fragments.jls"))
    pid_to_fid = deserialize_from_jls(joinpath(lib_path, "precursor_to_fragment_indices.jls"))

    @info "DIA-NN decoy generation: $n_targets target precursors, $(length(detailed_frags)) fragments"

    # Build decoy fragments
    F = eltype(detailed_frags)
    all_new_frags = F[]  # will hold target frags + decoy frags
    new_pid_to_fid = UInt64[UInt64(1)]  # start indices, 1-indexed

    # First: copy all target fragments as-is
    for pid in 1:n_targets
        frag_start = Int(pid_to_fid[pid])
        frag_end = Int(pid_to_fid[pid + 1]) - 1
        append!(all_new_frags, @view detailed_frags[frag_start:frag_end])
        push!(new_pid_to_fid, UInt64(length(all_new_frags) + 1))
    end

    # Build decoy precursor rows
    decoy_rows = DataFrame()
    for col in names(prec_df)
        decoy_rows[!, col] = similar(prec_df[!, col], 0)
    end

    n_created = 0
    n_skipped = 0

    for pid in 1:n_targets
        seq = prec_df.sequence[pid]
        if length(seq) < 3
            n_skipped += 1
            continue
        end

        # Compute mutation shifts
        mod_positions = extract_mod_positions(
            hasproperty(prec_df, :structural_mods) ? prec_df.structural_mods[pid] : missing
        )
        n_shift, c_shift = compute_diann_mutation_shifts(seq, mod_positions)

        # Create mutated sequence string
        chars = collect(seq)
        n = length(chars)
        # Apply same position selection as compute_diann_mutation_shifts
        n_pos = if 2 ∉ mod_positions; 2
        elseif min(3, n-1) ∉ mod_positions; min(3, n-1)
        elseif 1 ∉ mod_positions; 1
        else; 2; end

        c_pos = if (n-1) ∉ mod_positions; n-1
        elseif max(1, n-2) ∉ mod_positions; max(1, n-2)
        elseif n ∉ mod_positions; n
        else; n-1; end

        n_aa = chars[n_pos]
        c_aa = chars[c_pos]
        chars[n_pos] = get(DIANN_MUTATION_TABLE, n_aa, n_aa)
        chars[c_pos] = get(DIANN_MUTATION_TABLE, c_aa, c_aa)
        decoy_seq = String(chars)

        # Decoy precursor ID (targets are 1:n_targets, decoys are n_targets+1:2*n_targets)
        decoy_pid = UInt32(n_targets + n_created + 1)

        # Create shifted fragments
        target_frag_start = Int(pid_to_fid[pid])
        target_frag_end = Int(pid_to_fid[pid + 1]) - 1
        target_range = target_frag_start:target_frag_end

        shifted = create_diann_decoy_fragments(
            detailed_frags, target_range,
            n_shift, c_shift, decoy_pid
        )
        append!(all_new_frags, shifted)
        push!(new_pid_to_fid, UInt64(length(all_new_frags) + 1))

        # Build decoy precursor row (copy target row, modify key fields)
        decoy_row = copy(prec_df[pid:pid, :])
        decoy_row.sequence[1] = decoy_seq
        decoy_row.is_decoy[1] = true
        # Shift precursor m/z
        prec_charge = Float32(decoy_row.prec_charge[1])
        decoy_row.mz[1] = prec_df.mz[pid] + Float32((n_shift + c_shift) / prec_charge)
        # Copy target iRT (unchanged)
        decoy_row.irt[1] = prec_df.irt[pid]

        append!(decoy_rows, decoy_row)
        n_created += 1
    end

    @info "DIA-NN decoy generation: created $n_created decoys, skipped $n_skipped (too short)"

    # Combine target + decoy precursor tables
    combined_df = vcat(prec_df, decoy_rows)
    n_total = nrow(combined_df)

    # Set up pairing columns
    combined_df.is_decoy = vcat(falses(n_targets), trues(n_created))

    # pair_id: each target-decoy pair gets the same ID
    pair_ids = Vector{UInt32}(undef, n_total)
    for i in 1:n_targets
        pair_ids[i] = UInt32(i)
    end
    for i in 1:n_created
        pair_ids[n_targets + i] = UInt32(i)  # same pair_id as target
    end
    combined_df.pair_id = pair_ids

    # partner_precursor_idx: target → decoy, decoy → target
    partner_idx = Vector{Union{Missing, UInt32}}(missing, n_total)
    for i in 1:n_created
        target_i = i  # targets that weren't skipped map 1:1
        decoy_i = n_targets + i
        partner_idx[target_i] = UInt32(decoy_i)
        partner_idx[decoy_i] = UInt32(target_i)
    end
    # Handle skipped targets (no decoy partner)
    # They already have missing partner_idx
    combined_df.partner_precursor_idx = partner_idx

    # Save everything
    Arrow.write(joinpath(lib_path, "precursors_table.arrow"), combined_df)
    serialize_to_jls(joinpath(lib_path, "detailed_fragments.jls"), all_new_frags)
    serialize_to_jls(joinpath(lib_path, "precursor_to_fragment_indices.jls"), new_pid_to_fid)

    # Rebuild partitioned fragment indexes
    temp_precursors = SetPrecursors(Arrow.Table(joinpath(lib_path, "precursors_table.arrow")))
    temp_proteins = SetProteins(Arrow.Table(joinpath(lib_path, "proteins_table.arrow")))
    empty_pfi = LocalPartitionedFragmentIndex{Float32}(LocalPartition{Float32}[], Tuple{Float32,Float32}[], 0)
    if F <: SplineDetailedFrag
        spl_knots = if isfile(joinpath(lib_path, "spline_knots.jls"))
            deserialize_from_jls(joinpath(lib_path, "spline_knots.jls"))
        elseif isfile(joinpath(lib_path, "spline_knots.jld2"))
            load(joinpath(lib_path, "spline_knots.jld2"))["spl_knots"]
        else
            error("spline_knots file not found in $lib_path")
        end
        temp_lookup = SplineFragmentLookup(all_new_frags, new_pid_to_fid, Tuple(spl_knots), 3)
        temp_lib = SplineFragmentIndexLibrary(empty_pfi, empty_pfi, temp_precursors, temp_proteins, temp_lookup)
    else
        temp_lookup = StandardFragmentLookup(all_new_frags, new_pid_to_fid)
        temp_lib = FragmentIndexLibrary(empty_pfi, empty_pfi, temp_precursors, temp_proteins, temp_lookup)
    end

    # Read library config for bin tolerances
    config = JSON.parsefile(joinpath(lib_path, "config.json"))
    lib_params = get(config, "library_params", Dict())
    frag_bin_tol_ppm = Float32(get(lib_params, "frag_bin_tol_ppm", 10.0))
    rt_bin_tol = Float32(get(lib_params, "rt_bin_tol", 1.0))

    partitioned_index = build_partitioned_index_from_lib(temp_lib;
        partition_width=5.0f0, frag_bin_tol_ppm=frag_bin_tol_ppm, rt_bin_tol=rt_bin_tol)
    serialize_to_jls(joinpath(lib_path, "partitioned_fragment_index.jls"), partitioned_index)

    presearch_partitioned_index = build_partitioned_index_from_lib(temp_lib;
        partition_width=5.0f0, frag_bin_tol_ppm=frag_bin_tol_ppm, rt_bin_tol=typemax(Float32))
    serialize_to_jls(joinpath(lib_path, "presearch_partitioned_fragment_index.jls"), presearch_partitioned_index)

    @info "DIA-NN decoy generation complete: $n_total total precursors ($n_targets targets + $n_created decoys)"
    return nothing
end

"""
    extract_mod_positions(mods_string) -> Set{Int}

Extract 1-indexed positions of modifications from a modification string.
Returns empty set if no mods or missing.
"""
function extract_mod_positions(mods_string)
    positions = Set{Int}()
    (ismissing(mods_string) || isempty(mods_string)) && return positions

    # Parse mod string format: "pos1:mod1,pos2:mod2,..." or "(pos)mod" etc.
    # Try the common Pioneer format: semicolon or comma separated "position:mod_name"
    for part in split(mods_string, r"[;,]")
        part = strip(part)
        isempty(part) && continue
        m = match(r"^(\d+)", part)
        if m !== nothing
            push!(positions, parse(Int, m.captures[1]))
        end
    end
    return positions
end
