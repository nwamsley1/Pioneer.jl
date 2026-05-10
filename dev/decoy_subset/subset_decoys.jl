#=
subset_decoys.jl

Build a copy of a Pioneer .poin spectral library with a fraction of decoys
removed at random. Targets are kept intact. Pioneer's library_fdr_scale_factor
(n_targets / n_decoys) is computed at search init and applied throughout
get_qvalues!/get_PEP!, so no Pioneer source change is needed.

Usage:
    julia --project=/private/tmp/pioneer-develop-latest \
        dev/decoy_subset/subset_decoys.jl \
        <input_lib.poin> <output_lib.poin> <decoy_fraction> [seed]

Steps:
  1. Filter precursors_table.arrow → all targets + random `fraction` of decoys
  2. Renumber precursor_idx to dense 1..N (precursors_table row index)
  3. Update partner_precursor_idx (set missing for orphaned partners)
  4. Build new detailed_fragments.jls + precursor_to_fragment_indices.jls
     (fragments copied in new precursor order; cumulative offsets)
  5. Build new partitioned_fragment_index.jls + presearch_*.jls via
     build_partitioned_index_from_lib(...)
  6. Copy other files unchanged: spline_knots.jls, frag_name_to_idx.jls,
     ion_annotations.jls, proteins_table.arrow, config.json
=#

using Random, Arrow, DataFrames, Tables, Serialization, Printf

# Pioneer must be loaded so we can access internal types & helpers
using Pioneer
using Pioneer: deserialize_from_jls, serialize_to_jls,
               SetPrecursors, SetProteins,
               SplineFragmentLookup, StandardFragmentLookup,
               SplineFragmentIndexLibrary, FragmentIndexLibrary,
               LocalPartitionedFragmentIndex, LocalPartition,
               build_partitioned_index_from_lib,
               OutputSchemaPolicy,
               SplineCompactFrag, SplineDetailedFrag, CompactFrag, DetailedFrag

function main(args::Vector{String})
    length(args) < 3 && error("usage: subset_decoys.jl <input.poin> <output.poin> <fraction> [seed]")
    input_dir   = args[1]
    output_dir  = args[2]
    fraction    = parse(Float64, args[3])
    seed        = length(args) >= 4 ? parse(Int, args[4]) : 1234
    @assert 0.0 < fraction <= 1.0
    @assert isdir(input_dir) "input library not found: $input_dir"
    isdir(output_dir) && error("output already exists: $output_dir (refusing to overwrite)")

    @info "Subsetting decoys" input_dir output_dir fraction seed
    mkpath(output_dir)

    # ── Step 1: load + filter precursors ────────────────────────────────────
    @info "Loading precursors_table.arrow ..."
    prec_tbl = Arrow.Table(joinpath(input_dir, "precursors_table.arrow"))
    n_old = length(prec_tbl.is_decoy)
    is_decoy = collect(Bool.(prec_tbl.is_decoy))
    n_targets_old = count(!, is_decoy)
    n_decoys_old  = count(is_decoy)
    @info "Original library" n_total=n_old n_targets=n_targets_old n_decoys=n_decoys_old

    # Pick which old indices to keep
    rng = MersenneTwister(seed)
    keep_old_to_new = zeros(UInt32, n_old)
    n_keep = 0
    decoy_keep_threshold = fraction
    for i in 1:n_old
        if !is_decoy[i] || rand(rng) < decoy_keep_threshold
            n_keep += 1
            keep_old_to_new[i] = UInt32(n_keep)
        end
    end
    kept_old_indices = findall(!iszero, keep_old_to_new)  # in old order
    n_decoys_new = count(i -> is_decoy[i], kept_old_indices)
    @info "Filtered" n_keep=length(kept_old_indices) n_targets_kept=(length(kept_old_indices)-n_decoys_new) n_decoys_kept=n_decoys_new
    @info @sprintf("New target:decoy ratio = %.2f (fdr_scale_factor will be %.2f)",
                   (length(kept_old_indices)-n_decoys_new) / max(n_decoys_new, 1),
                   (length(kept_old_indices)-n_decoys_new) / max(n_decoys_new, 1))

    # ── Step 2: build new precursors_table ──────────────────────────────────
    @info "Building new precursors_table.arrow ..."
    prec_df = DataFrame(prec_tbl; copycols=false)
    new_df = prec_df[kept_old_indices, :]

    # Remap partner_precursor_idx
    if hasproperty(new_df, :partner_precursor_idx)
        old_partner = new_df.partner_precursor_idx
        E = nonmissingtype(eltype(old_partner))
        new_partner = Vector{Union{Missing, E}}(missing, length(old_partner))
        for i in 1:length(old_partner)
            p = old_partner[i]
            if !ismissing(p)
                old_p = Int(p)
                if old_p >= 1 && old_p <= n_old
                    new_p = keep_old_to_new[old_p]
                    if new_p != 0
                        new_partner[i] = E(new_p)
                    end
                end
            end
        end
        new_df.partner_precursor_idx = new_partner
    end

    Arrow.write(joinpath(output_dir, "precursors_table.arrow"), new_df)
    new_df = nothing  # release

    # ── Step 3: build new fragments + pid_to_fid ────────────────────────────
    @info "Loading detailed_fragments.jls ..."
    old_frags = deserialize_from_jls(joinpath(input_dir, "detailed_fragments.jls"))
    @info "Loading precursor_to_fragment_indices.jls ..."
    old_pid_to_fid = deserialize_from_jls(joinpath(input_dir, "precursor_to_fragment_indices.jls"))
    @assert length(old_pid_to_fid) == n_old + 1 "old pid_to_fid length=$(length(old_pid_to_fid)) expected $(n_old+1)"

    # Total kept fragments
    total_new = UInt64(0)
    for old_i in kept_old_indices
        total_new += UInt64(old_pid_to_fid[old_i+1] - old_pid_to_fid[old_i])
    end
    @info "Building new fragments vector" n_new_frags=Int(total_new)

    F = eltype(old_frags)
    new_frags = Vector{F}(undef, Int(total_new))
    new_pid_to_fid = Vector{UInt64}(undef, length(kept_old_indices) + 1)
    new_pid_to_fid[1] = UInt64(1)

    cursor = UInt64(1)
    @inbounds for (new_i, old_i) in enumerate(kept_old_indices)
        lo = old_pid_to_fid[old_i]
        hi = old_pid_to_fid[old_i+1] - one(UInt64)
        for fi in lo:hi
            new_frags[Int(cursor)] = old_frags[Int(fi)]
            cursor += one(UInt64)
        end
        new_pid_to_fid[new_i+1] = cursor
    end
    @assert cursor == total_new + 1

    @info "Serializing new fragments + pid_to_fid ..."
    serialize_to_jls(joinpath(output_dir, "detailed_fragments.jls"), new_frags)
    serialize_to_jls(joinpath(output_dir, "precursor_to_fragment_indices.jls"), new_pid_to_fid)
    old_frags = nothing  # release

    # ── Step 4: rebuild fragment indexes ────────────────────────────────────
    @info "Building new SetPrecursors / SetProteins ..."
    temp_precursors = SetPrecursors(Arrow.Table(joinpath(output_dir, "precursors_table.arrow")))
    temp_proteins   = SetProteins(Arrow.Table(joinpath(input_dir, "proteins_table.arrow")))

    is_spline = (F <: SplineCompactFrag) || (F <: SplineDetailedFrag)
    if is_spline
        spl_knots_path = joinpath(input_dir, "spline_knots.jls")
        @assert isfile(spl_knots_path) "spline_knots.jls missing"
        spl_knots = deserialize_from_jls(spl_knots_path)
        temp_lookup = SplineFragmentLookup(new_frags, new_pid_to_fid, Tuple(spl_knots), 3)
        empty_pfi = LocalPartitionedFragmentIndex{Float32}(LocalPartition{Float32}[], Tuple{Float32,Float32}[], 0)
        temp_lib = SplineFragmentIndexLibrary(empty_pfi, empty_pfi, temp_precursors, temp_proteins, temp_lookup, OutputSchemaPolicy())
    else
        temp_lookup = StandardFragmentLookup(new_frags, new_pid_to_fid)
        empty_pfi = LocalPartitionedFragmentIndex{Float32}(LocalPartition{Float32}[], Tuple{Float32,Float32}[], 0)
        temp_lib = FragmentIndexLibrary(empty_pfi, empty_pfi, temp_precursors, temp_proteins, temp_lookup, OutputSchemaPolicy())
    end

    # Match BuildSpecLib.jl defaults: frag_bin_tol_ppm = 0 (mDa mode), rt_bin_tol = 3
    frag_bin_tol_ppm = 0.0f0
    frag_bin_tol_mda = 2.0f0
    rt_bin_tol       = 3.0f0

    @info "Building partitioned_fragment_index.jls (main, rt_bin_tol=$rt_bin_tol) ..."
    partitioned_index = build_partitioned_index_from_lib(temp_lib;
        partition_width=5.0f0,
        frag_bin_tol_ppm=frag_bin_tol_ppm,
        frag_bin_tol_mda=frag_bin_tol_mda,
        rt_bin_tol=rt_bin_tol)
    serialize_to_jls(joinpath(output_dir, "partitioned_fragment_index.jls"), partitioned_index)
    partitioned_index = nothing

    @info "Building presearch_partitioned_fragment_index.jls (rt_bin_tol=∞) ..."
    presearch_index = build_partitioned_index_from_lib(temp_lib;
        partition_width=5.0f0,
        frag_bin_tol_ppm=frag_bin_tol_ppm,
        frag_bin_tol_mda=frag_bin_tol_mda,
        rt_bin_tol=typemax(Float32))
    serialize_to_jls(joinpath(output_dir, "presearch_partitioned_fragment_index.jls"), presearch_index)
    presearch_index = nothing

    # ── Step 5: copy unchanged auxiliary files ──────────────────────────────
    for fname in ["spline_knots.jls", "frag_name_to_idx.jls", "ion_annotations.jls",
                  "proteins_table.arrow", "config.json"]
        src = joinpath(input_dir, fname)
        if isfile(src)
            cp(src, joinpath(output_dir, fname); force=false)
        end
    end

    @info "Done." output_dir
end

main(ARGS)
