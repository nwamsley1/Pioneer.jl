#=
iRT_match_decoys.jl

Build a copy of a Pioneer .poin spectral library where every decoy precursor
has its iRT replaced by its paired target's iRT. All other fields, fragment
intensities, and library tables are unchanged.

Fragment indexes are rebuilt because their RT bins are derived from precursor
iRT (see PartitionedFragmentIndex/build.jl).

Usage:
    julia --project=. dev/decoy_subset/iRT_match_decoys.jl <input.poin> <output.poin>
=#

using Arrow, DataFrames, Tables, Serialization, Printf, Statistics
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
    length(args) < 2 && error("usage: iRT_match_decoys.jl <input.poin> <output.poin>")
    input_dir, output_dir = args[1], args[2]
    @assert isdir(input_dir) "input library not found: $input_dir"
    isdir(output_dir) && error("output already exists: $output_dir (refusing to overwrite)")

    @info "Building iRT-matched-decoy library" input_dir output_dir
    mkpath(output_dir)

    # ── Step 1: load precursors_table, override decoy iRT to target iRT ────
    @info "Loading precursors_table.arrow ..."
    prec_tbl = Arrow.Table(joinpath(input_dir, "precursors_table.arrow"))
    df = DataFrame(prec_tbl; copycols=true)
    n = nrow(df)
    @info "  rows: $n  targets: $(count(!, df.is_decoy))  decoys: $(count(df.is_decoy))"

    # New iRT vector
    new_irt = Float32.(df.irt)
    n_paired = 0
    n_orphan = 0
    @inbounds for i in 1:n
        if df.is_decoy[i] && !ismissing(df.partner_precursor_idx[i])
            target_idx = Int(df.partner_precursor_idx[i])
            if target_idx >= 1 && target_idx <= n
                new_irt[i] = Float32(df.irt[target_idx])
                n_paired += 1
            else
                n_orphan += 1
            end
        end
    end
    @info "  decoys updated to target iRT: $n_paired   orphan decoys: $n_orphan"
    df.irt = new_irt

    # Stats: how much did decoy iRT shift?
    deltas = Float32[]
    @inbounds for i in 1:n
        if df.is_decoy[i] && !ismissing(df.partner_precursor_idx[i])
            target_idx = Int(df.partner_precursor_idx[i])
            if target_idx >= 1 && target_idx <= n
                push!(deltas, new_irt[i] - Float32(prec_tbl.irt[i]))
            end
        end
    end
    if !isempty(deltas)
        abs_d = abs.(deltas)
        @info "  decoy iRT shift  median |Δ|: $(round(quantile(abs_d, 0.5), digits=3))   max |Δ|: $(round(maximum(abs_d), digits=3))"
    end

    Arrow.write(joinpath(output_dir, "precursors_table.arrow"), df)

    # ── Step 2: copy unchanged auxiliary files ────────────────────────────
    @info "Copying detailed_fragments.jls (fragments unchanged) ..."
    cp(joinpath(input_dir, "detailed_fragments.jls"),
       joinpath(output_dir, "detailed_fragments.jls"); force=false)
    @info "Copying precursor_to_fragment_indices.jls ..."
    cp(joinpath(input_dir, "precursor_to_fragment_indices.jls"),
       joinpath(output_dir, "precursor_to_fragment_indices.jls"); force=false)
    for fname in ["spline_knots.jls", "frag_name_to_idx.jls", "ion_annotations.jls",
                   "proteins_table.arrow", "config.json"]
        src = joinpath(input_dir, fname)
        isfile(src) && cp(src, joinpath(output_dir, fname); force=false)
    end

    # ── Step 3: rebuild fragment indexes (RT bins depend on iRT) ──────────
    @info "Rebuilding fragment indexes ..."
    temp_precursors = SetPrecursors(Arrow.Table(joinpath(output_dir, "precursors_table.arrow")))
    temp_proteins   = SetProteins(Arrow.Table(joinpath(output_dir, "proteins_table.arrow")))
    detailed_frags  = deserialize_from_jls(joinpath(output_dir, "detailed_fragments.jls"))
    pid_to_fid      = deserialize_from_jls(joinpath(output_dir, "precursor_to_fragment_indices.jls"))
    F = eltype(detailed_frags)
    is_spline = (F <: SplineCompactFrag) || (F <: SplineDetailedFrag)
    if is_spline
        spl_knots = deserialize_from_jls(joinpath(output_dir, "spline_knots.jls"))
        temp_lookup = SplineFragmentLookup(detailed_frags, pid_to_fid, Tuple(spl_knots), 3)
        empty_pfi = LocalPartitionedFragmentIndex{Float32}(LocalPartition{Float32}[], Tuple{Float32,Float32}[], 0)
        temp_lib = SplineFragmentIndexLibrary(empty_pfi, empty_pfi, temp_precursors, temp_proteins, temp_lookup, OutputSchemaPolicy())
    else
        temp_lookup = StandardFragmentLookup(detailed_frags, pid_to_fid)
        empty_pfi = LocalPartitionedFragmentIndex{Float32}(LocalPartition{Float32}[], Tuple{Float32,Float32}[], 0)
        temp_lib = FragmentIndexLibrary(empty_pfi, empty_pfi, temp_precursors, temp_proteins, temp_lookup, OutputSchemaPolicy())
    end

    # Match BuildSpecLib defaults
    frag_bin_tol_ppm = 0.0f0
    frag_bin_tol_mda = 2.0f0
    rt_bin_tol       = 3.0f0

    @info "  building partitioned_fragment_index.jls (rt_bin_tol=$rt_bin_tol) ..."
    pfi_main = build_partitioned_index_from_lib(temp_lib;
        partition_width=5.0f0, frag_bin_tol_ppm=frag_bin_tol_ppm,
        frag_bin_tol_mda=frag_bin_tol_mda, rt_bin_tol=rt_bin_tol)
    serialize_to_jls(joinpath(output_dir, "partitioned_fragment_index.jls"), pfi_main)
    pfi_main = nothing

    @info "  building presearch_partitioned_fragment_index.jls (rt_bin_tol=∞) ..."
    pfi_pre = build_partitioned_index_from_lib(temp_lib;
        partition_width=5.0f0, frag_bin_tol_ppm=frag_bin_tol_ppm,
        frag_bin_tol_mda=frag_bin_tol_mda, rt_bin_tol=typemax(Float32))
    serialize_to_jls(joinpath(output_dir, "presearch_partitioned_fragment_index.jls"), pfi_pre)
    pfi_pre = nothing

    @info "Done." output_dir
end

main(ARGS)
