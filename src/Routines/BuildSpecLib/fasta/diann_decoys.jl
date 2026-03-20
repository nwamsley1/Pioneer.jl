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

Returns a Vector{DetailedFrag{Float32}}.
"""
function create_diann_decoy_fragments(
    target_frags::AbstractVector{DetailedFrag{Float32}},
    target_range::UnitRange,
    n_shift::Float64,
    c_shift::Float64,
    decoy_prec_id::UInt32
)
    decoy_frags = Vector{DetailedFrag{Float32}}(undef, length(target_range))
    for (j, i) in enumerate(target_range)
        frag = target_frags[i]
        charge = max(frag.frag_charge, UInt8(1))
        mz_shift = frag.is_y ? Float32(c_shift / charge) : Float32(n_shift / charge)

        decoy_frags[j] = DetailedFrag{Float32}(
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

Post-process a built Pioneer library to convert decoy fragments to DIA-NN style.

For each target-decoy pair:
1. Computes N_shift and C_shift from the DIA-NN mutation of the target sequence
2. Replaces decoy fragments with shifted copies of the target's fragments
3. Copies target iRT to the decoy

Modifies `detailed_fragments.jls`, `precursor_to_fragment_indices.jls`,
`precursors_table.arrow`, and rebuilds `partitioned_fragment_index.jls` in place.
"""
function apply_diann_decoy_style!(lib_path::String)
    # Load precursors table
    prec_df = DataFrame(Arrow.Table(joinpath(lib_path, "precursors_table.arrow")))

    # Load existing fragments and index
    detailed_frags = load_detailed_frags(joinpath(lib_path, "detailed_fragments.jls"))
    pid_to_fid = deserialize_from_jls(joinpath(lib_path, "precursor_to_fragment_indices.jls"))

    n_precs = nrow(prec_df)

    # Build new fragments array and new pid_to_fid
    new_frags = DetailedFrag{Float32}[]
    new_pid_to_fid = Vector{UInt64}(undef, n_precs + 1)
    new_pid_to_fid[1] = UInt64(1)

    n_converted = 0
    n_skipped = 0

    for pid in 1:n_precs
        frag_start = Int(pid_to_fid[pid])
        frag_end = Int(pid_to_fid[pid + 1]) - 1

        if !prec_df.is_decoy[pid]
            # Target: keep existing fragments
            append!(new_frags, @view detailed_frags[frag_start:frag_end])
        else
            # Decoy: find paired target and create shifted fragments
            partner_idx = prec_df.partner_precursor_idx[pid]

            if ismissing(partner_idx) || prec_df.is_decoy[partner_idx]
                # No valid target partner — keep existing fragments as fallback
                append!(new_frags, @view detailed_frags[frag_start:frag_end])
                n_skipped += 1
            else
                # Get target sequence and compute shifts
                target_seq = prec_df.sequence[partner_idx]
                # Get mod positions for the target (positions where mods exist)
                target_mods = prec_df.structural_mods[partner_idx]
                mod_positions = extract_mod_positions(target_mods)

                n_shift, c_shift = compute_diann_mutation_shifts(target_seq, mod_positions)

                # Get target fragment range
                target_frag_start = Int(pid_to_fid[partner_idx])
                target_frag_end = Int(pid_to_fid[partner_idx + 1]) - 1
                target_range = target_frag_start:target_frag_end

                # Create shifted copies of target fragments
                shifted = create_diann_decoy_fragments(
                    detailed_frags, target_range,
                    n_shift, c_shift, UInt32(pid)
                )
                append!(new_frags, shifted)

                # Copy target iRT to decoy
                prec_df.irt[pid] = prec_df.irt[partner_idx]

                n_converted += 1
            end
        end

        new_pid_to_fid[pid + 1] = UInt64(length(new_frags) + 1)
    end

    @info "DIA-NN decoy conversion: $n_converted decoys converted, $n_skipped skipped (no target partner)"

    # Save updated data
    serialize_to_jls(joinpath(lib_path, "detailed_fragments.jls"), new_frags)
    serialize_to_jls(joinpath(lib_path, "precursor_to_fragment_indices.jls"), new_pid_to_fid)
    Arrow.write(joinpath(lib_path, "precursors_table.arrow"), prec_df)

    # Rebuild partitioned fragment indexes
    temp_precursors = SetPrecursors(Arrow.Table(joinpath(lib_path, "precursors_table.arrow")))
    temp_proteins = SetProteins(Arrow.Table(joinpath(lib_path, "proteins_table.arrow")))
    temp_lookup = StandardFragmentLookup(new_frags, new_pid_to_fid)

    empty_pfi = LocalPartitionedFragmentIndex{Float32}(LocalPartition{Float32}[], Tuple{Float32,Float32}[], 0)
    temp_lib = FragmentIndexLibrary(empty_pfi, empty_pfi, temp_precursors, temp_proteins, temp_lookup)

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

    @info "DIA-NN decoy conversion complete. Rebuilt partitioned fragment indexes."
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
