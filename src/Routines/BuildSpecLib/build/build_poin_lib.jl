"""
    sort_detailed_fragments_by_mz!(frags, prec_ranges) -> Int

Sort each precursor's fragments in `frags` in place by ascending m/z.
`prec_ranges[k]` is the starting index of precursor `k`; precursor `k`'s
range is `prec_ranges[k] : prec_ranges[k+1] - 1`.

This must run AFTER `build_partitioned_index_from_lib` because the
partitioned index assigns bitmask bits via iteration order (`rank += 1`
in `PartitionedFragmentIndex/build.jl`), which requires the rank-sorted
layout. After the index is built, sorting `frags` by m/z satisfies
`run_fused!`'s pre-condition (`verify_mz_sorted` in fusedScan.jl) without
changing partitioned-index semantics — the index keeps its own array of
`SimpleFrag`s with bitmask bits already encoded.

Returns the number of precursors that needed mutation (the rest were
already m/z-sorted — e.g. single-fragment precursors).
"""
function sort_detailed_fragments_by_mz!(frags::AbstractVector,
                                         prec_ranges::AbstractVector)
    n_prec = length(prec_ranges) - 1
    n_mutated = 0
    @inbounds for pid in 1:n_prec
        lo = Int(prec_ranges[pid])
        hi = Int(prec_ranges[pid + 1]) - 1
        lo > hi && continue
        v = view(frags, lo:hi)
        if !issorted(v, by = f -> f.mz)
            sort!(v, by = f -> f.mz, alg = InsertionSort)
            n_mutated += 1
        end
    end
    return n_mutated
end

function validate_fragment_index_filters(y_start_index::UInt8,
                                         y_start::UInt8,
                                         b_start_index::UInt8,
                                         b_start::UInt8,
                                         include_p_index::Bool,
                                         include_p::Bool)
    y_start_index < y_start && throw(ArgumentError(
        "y_start_index ($y_start_index) cannot be less stringent than y_start ($y_start)"))
    b_start_index < b_start && throw(ArgumentError(
        "b_start_index ($b_start_index) cannot be less stringent than b_start ($b_start)"))
    include_p_index && !include_p && throw(ArgumentError(
        "include_p_index=true requires include_p=true because the index is built from detailed fragments"))
    return nothing
end


"""
    buildPionLib(spec_lib_path::String,
                y_start_index::UInt8,
                y_start::UInt8,
                b_start_index::UInt8,
                b_start::UInt8,
                include_p_index::Bool,
                include_p::Bool,
                include_isotope::Bool,
                include_immonium::Bool,
                include_internal::Bool,
                include_neutral_diff::Bool,
                max_frag_charge::UInt8,
                max_frag_rank::UInt8,
                min_frag_intensity::Float32,
                frag_bounds::FragBoundModel,
                frag_bin_tol_ppm::Float32,
                rt_bin_tol_ppm::Float32,
                model_type::SplineCoefficientModel)

Build a Pioneer spectral library from preprocessed fragment and precursor data.

# Parameters
- `spec_lib_path`: Path to the directory containing input files and where output files will be written
- `y_start_index`: Minimum index for y-ions to include in the fragment index
- `y_start`: Minimum index for y-ions to include in detailed fragments
- `b_start_index`: Minimum index for b-ions to include in the fragment index
- `b_start`: Minimum index for b-ions to include in detailed fragments
- `include_p_index`: Whether to include precursor ions in the fragment index
- `include_p`: Whether to include precursor ions in detailed fragments
- `include_isotope`: Whether to include isotope peaks
- `include_immonium`: Whether to include immonium ions
- `include_internal`: Whether to include internal fragment ions
- `include_neutral_diff`: Whether to include fragments with neutral losses
- `max_frag_charge`: Maximum fragment charge state to include
- `max_frag_rank`: Maximum number of fragments per precursor (ranked by intensity)
- `min_frag_intensity`: Minimum relative intensity threshold for fragments
- `frag_bounds`: Model defining minimum and maximum fragment m/z bounds
- `frag_bin_tol_ppm`: Fragment binning tolerance in parts per million
- `rt_bin_tol_ppm`: Retention time binning tolerance in parts per million
- `model_type`: SplineCoefficientModel describing the Altimeter prediction backend

# Returns
- `nothing`: Results are written to files in `spec_lib_path`

# Input Files
Requires the following files in `spec_lib_path`:
- fragments_table.arrow: Table of fragment data
- prec_to_frag.arrow: Table mapping precursors to fragments
- precursors_table.arrow: Table of precursor data

# Output Files
Creates the following files in `spec_lib_path`:
- f_index_fragments.arrow: Fragment index
- f_index_rt_bins.arrow: Retention time bins
- f_index_fragment_bins.arrow: Fragment m/z bins
- presearch_f_index_fragments.arrow: Presearch fragment index
- presearch_f_index_rt_bins.arrow: Presearch retention time bins
- presearch_f_index_fragment_bins.arrow: Presearch fragment m/z bins
- detailed_fragments.jls: Detailed fragment information
- precursor_to_fragment_indices.jls: Mapping of precursors to fragment indices
"""
function buildPionLib(spec_lib_path::String,
                      y_start_index::UInt8,
                      y_start::UInt8,
                      b_start_index::UInt8,
                      b_start::UInt8,
                      include_p_index::Bool,
                      include_p::Bool,
                      include_isotope::Bool,
                      include_immonium::Bool,
                      include_internal::Bool,
                      include_neutral_diff::Bool,
                      max_frag_charge::UInt8,
                      max_frag_rank::UInt8,
                      length_to_frag_count_multiple::AbstractFloat,
                      min_frag_intensity::Float32,
                      frag_bounds::FragBoundModel,
                      frag_bin_tol_ppm::Float32,
                      rt_bin_tol_ppm::Float32,
                      model_type::SplineCoefficientModel;
                      frag_bin_tol_mda::Float32 = 2.0f0,
                      detailed_frags = nothing,
                      pid_to_fid = nothing,
                      )
    validate_fragment_index_filters(
        y_start_index, y_start,
        b_start_index, b_start,
        include_p_index, include_p,
    )

    # Fused path (production): detailed_frags + pid_to_fid are handed in from
    # build_detailed_frags_from_raw — no fragments_table.arrow round-trip.
    # Disk path (resume with predict_fragments=false, and the unit tests): read
    # fragments_table.arrow + prec_to_frag.arrow and decode via getDetailedFrags.
    if detailed_frags === nothing
        fragments_table, prec_to_frag, precursors_table = nothing, nothing, nothing
        try
            fragments_table = Arrow.Table(joinpath(spec_lib_path,"fragments_table.arrow"));
            prec_to_frag = Arrow.Table(joinpath(spec_lib_path,"prec_to_frag.arrow"));
            precursors_table = Arrow.Table(joinpath(spec_lib_path,"precursors_table.arrow"));
        catch e
            @error "could not find library..."
            return nothing
        end

        #println("Get full fragments list...")
        detailed_frags, pid_to_fid = getDetailedFrags(
        fragments_table[:mz],
        fragments_table[:coefficients],
        fragments_table[:intensity],
        fragments_table[:is_y],
        fragments_table[:is_b],
        fragments_table[:is_p],
        fragments_table[:fragment_index],
        fragments_table[:charge],
        fragments_table[:sulfur_count],
        fragments_table[:ion_type],
        fragments_table[:isotope],
        fragments_table[:is_internal],
        fragments_table[:is_immonium],
        fragments_table[:has_neutral_diff],
        precursors_table[:mz],
        precursors_table[:prec_charge],#:precursor_charge],
        precursors_table[:length],
        prec_to_frag[:start_idx],
        y_start,
        b_start,
        include_p,
        include_isotope,
        include_immonium,
        include_internal,
        include_neutral_diff,
        max_frag_charge,
        frag_bounds,
        max_frag_rank,
        length_to_frag_count_multiple,
        min_frag_intensity,
        model_type
        );

        # fragments_table.arrow + prec_to_frag.arrow are now fully consumed —
        # getDetailedFrags above is their only reader. Release the mmap handles
        # and delete them before serializing the .jls outputs.
        fragments_table = nothing
        prec_to_frag = nothing
        GC.gc()
        safeRm(joinpath(spec_lib_path, "fragments_table.arrow"), nothing; force=true)
        safeRm(joinpath(spec_lib_path, "prec_to_frag.arrow"), nothing; force=true)
    end

    # Build partitioned fragment indexes BEFORE sorting detailed_frags by m/z.
    # See sort_detailed_fragments_by_mz! for why the order matters.
    precursors_arrow = Arrow.Table(joinpath(spec_lib_path, "precursors_table.arrow"))
    temp_precursors = SetPrecursors(precursors_arrow)
    # Load spline knots for SplineFragmentLookup
    spl_knots = if isfile(joinpath(spec_lib_path, "spline_knots.jls"))
        deserialize_from_jls(joinpath(spec_lib_path, "spline_knots.jls"))
    elseif isfile(joinpath(spec_lib_path, "spline_knots.jld2"))
        load(joinpath(spec_lib_path, "spline_knots.jld2"))["spl_knots"]
    else
        error("spline_knots file not found in $spec_lib_path")
    end
    temp_lookup = SplineFragmentLookup(detailed_frags, pid_to_fid, Tuple(spl_knots), 3)
    temp_proteins = SetProteins(Arrow.Table(joinpath(spec_lib_path, "proteins_table.arrow")))
    empty_pfi = LocalPartitionedFragmentIndex{Float32}(LocalPartition{Float32}[], Tuple{Float32,Float32}[], 0)
    temp_lib = SplineFragmentIndexLibrary(empty_pfi, empty_pfi, temp_precursors, temp_proteins, temp_lookup, OutputSchemaPolicy())

    partitioned_index = build_partitioned_index_from_lib(temp_lib;
        partition_width=5.0f0, frag_bin_tol_ppm=frag_bin_tol_ppm, frag_bin_tol_mda=frag_bin_tol_mda,
        rt_bin_tol=rt_bin_tol_ppm,
        y_start_index=y_start_index, b_start_index=b_start_index,
        include_p_index=include_p_index)

    presearch_partitioned_index = build_partitioned_index_from_lib(temp_lib;
        partition_width=5.0f0, frag_bin_tol_ppm=frag_bin_tol_ppm, frag_bin_tol_mda=frag_bin_tol_mda,
        rt_bin_tol=typemax(Float32),
        y_start_index=y_start_index, b_start_index=b_start_index,
        include_p_index=include_p_index)

    # Sort detailed_frags by m/z within each precursor (run_fused! pre-condition).
    sort_detailed_fragments_by_mz!(detailed_frags, pid_to_fid)

    serialize_to_jls(
        joinpath(spec_lib_path, "detailed_fragments.jls"),
        detailed_frags
    )

    serialize_to_jls(
        joinpath(spec_lib_path, "precursor_to_fragment_indices.jls"),
        pid_to_fid
    )

    serialize_to_jls(joinpath(spec_lib_path, "partitioned_fragment_index.jls"), partitioned_index)
    serialize_to_jls(joinpath(spec_lib_path, "presearch_partitioned_fragment_index.jls"), presearch_partitioned_index)

    # Arrow tables can keep Windows file handles alive until GC runs.
    # Release references before BuildSpecLib removes the intermediate Arrow files.
    fragments_table = nothing
    prec_to_frag = nothing
    precursors_table = nothing
    detailed_frags = nothing
    pid_to_fid = nothing
    GC.gc()

    return nothing
end

"""
    buildPionLib(..., model_type::InstrumentAgnosticModel; detailed_frags, pid_to_fid)

Prosit analog of the spline `buildPionLib`. Builds a `StandardFragmentLookup`
(constant scalar intensity, `ConstantType`) instead of a `SplineFragmentLookup`,
so there are no spline knots. Fused path only — `detailed_frags`/`pid_to_fid` are
handed in from `build_detailed_frags_from_raw(::InstrumentAgnosticModel)` (Prosit
resume-from-disk is not supported yet).
"""
function buildPionLib(spec_lib_path::String,
                      y_start_index::UInt8,
                      y_start::UInt8,
                      b_start_index::UInt8,
                      b_start::UInt8,
                      include_p_index::Bool,
                      include_p::Bool,
                      include_isotope::Bool,
                      include_immonium::Bool,
                      include_internal::Bool,
                      include_neutral_diff::Bool,
                      max_frag_charge::UInt8,
                      max_frag_rank::UInt8,
                      length_to_frag_count_multiple::AbstractFloat,
                      min_frag_intensity::Float32,
                      frag_bounds::FragBoundModel,
                      frag_bin_tol_ppm::Float32,
                      rt_bin_tol_ppm::Float32,
                      model_type::InstrumentAgnosticModel;
                      frag_bin_tol_mda::Float32 = 2.0f0,
                      detailed_frags = nothing,
                      pid_to_fid = nothing,
                      )
    validate_fragment_index_filters(
        y_start_index, y_start,
        b_start_index, b_start,
        include_p_index, include_p,
    )

    if detailed_frags === nothing
        error("buildPionLib(::InstrumentAgnosticModel) requires the fused path " *
              "(detailed_frags/pid_to_fid); Prosit resume-from-disk is not supported yet.")
    end

    precursors_arrow = Arrow.Table(joinpath(spec_lib_path, "precursors_table.arrow"))
    temp_precursors = SetPrecursors(precursors_arrow)
    temp_lookup = StandardFragmentLookup(detailed_frags, pid_to_fid)
    temp_proteins = SetProteins(Arrow.Table(joinpath(spec_lib_path, "proteins_table.arrow")))
    empty_pfi = LocalPartitionedFragmentIndex{Float32}(LocalPartition{Float32}[], Tuple{Float32,Float32}[], 0)
    temp_lib = FragmentIndexLibrary(empty_pfi, empty_pfi, temp_precursors, temp_proteins, temp_lookup, OutputSchemaPolicy())

    partitioned_index = build_partitioned_index_from_lib(temp_lib;
        partition_width=5.0f0, frag_bin_tol_ppm=frag_bin_tol_ppm, frag_bin_tol_mda=frag_bin_tol_mda,
        rt_bin_tol=rt_bin_tol_ppm,
        y_start_index=y_start_index, b_start_index=b_start_index,
        include_p_index=include_p_index)

    presearch_partitioned_index = build_partitioned_index_from_lib(temp_lib;
        partition_width=5.0f0, frag_bin_tol_ppm=frag_bin_tol_ppm, frag_bin_tol_mda=frag_bin_tol_mda,
        rt_bin_tol=typemax(Float32),
        y_start_index=y_start_index, b_start_index=b_start_index,
        include_p_index=include_p_index)

    sort_detailed_fragments_by_mz!(detailed_frags, pid_to_fid)

    serialize_to_jls(joinpath(spec_lib_path, "detailed_fragments.jls"), detailed_frags)
    serialize_to_jls(joinpath(spec_lib_path, "precursor_to_fragment_indices.jls"), pid_to_fid)
    serialize_to_jls(joinpath(spec_lib_path, "partitioned_fragment_index.jls"), partitioned_index)
    serialize_to_jls(joinpath(spec_lib_path, "presearch_partitioned_fragment_index.jls"), presearch_partitioned_index)

    detailed_frags = nothing
    pid_to_fid = nothing
    GC.gc()

    return nothing
end

"""
    cleanUpLibrary(spec_lib_path::String)

Remove intermediate files after library building is complete.

# Parameters
- `spec_lib_path::String`: Path to the library directory

# Effects
Removes the following files if they exist:
- raw_fragments.arrow
- fragments_table.arrow
- prec_to_frag.arrow
- precursors.arrow

# Returns
- `nothing`
"""
function cleanUpLibrary(spec_lib_path::String)
    # Debug toggle: keep build intermediates (raw_fragments/fragments_table/
    # prec_to_frag) for disk-usage analysis. Off by default.
    haskey(ENV, "PIONEER_KEEP_BUILD_TEMP") && return nothing
    GC.gc()
    for fname in ["raw_fragments.arrow", "fragments_table.arrow", "prec_to_frag.arrow", "precursors.arrow"]
        fpath = joinpath(spec_lib_path, fname)
        if isfile(fpath)
            try
                safeRm(fpath, nothing; force=true)
            catch e
                @user_warn "Failed to remove temporary file $(fpath): $(sprint(showerror, e))"
            end
        end
    end
    return nothing
end

"""
    fragFilter(
        frag_is_y::Bool,
        frag_is_b::Bool,
        frag_is_p::Bool,
        frag_index::UInt8,
        frag_charge::UInt8,
        frag_isotope::UInt8,
        frag_internal::Bool,
        frag_immonium::Bool,
        frag_neutral_diff::Bool,
        frag_mz::Float32,
        frag_bounds::FragBoundModel,
        prec_mz::Float32,
        y_start::UInt8,
        b_start::UInt8,
        include_p::Bool,
        include_isotope::Bool,
        include_immonium::Bool,
        include_internal::Bool,
        include_neutral_diff::Bool,
        max_frag_charge::UInt8)::Bool

Filter fragments based on type, index, charge, and additional criteria.

# Parameters
- `frag_is_y`: Whether the fragment is a y-ion
- `frag_is_b`: Whether the fragment is a b-ion
- `frag_is_p`: Whether the fragment is a precursor ion
- `frag_index`: Index of the fragment in the peptide sequence
- `frag_charge`: Charge state of the fragment
- `frag_isotope`: Isotope state of the fragment
- `frag_internal`: Whether the fragment is an internal fragment
- `frag_immonium`: Whether the fragment is an immonium ion
- `frag_neutral_diff`: Whether the fragment has neutral losses
- `frag_mz`: m/z value of the fragment
- `frag_bounds`: Model defining valid m/z range based on precursor m/z
- `prec_mz`: m/z value of the precursor
- `y_start`: Minimum index for y-ions to include
- `b_start`: Minimum index for b-ions to include
- `include_p`: Whether to include precursor ions
- `include_isotope`: Whether to include isotope peaks
- `include_immonium`: Whether to include immonium ions
- `include_internal`: Whether to include internal fragment ions
- `include_neutral_diff`: Whether to include fragments with neutral losses
- `max_frag_charge`: Maximum fragment charge state to include

# Returns
- `Bool`: true if fragment passes all filters, false otherwise
"""
function fragFilter(
    frag_is_y::Bool,
    frag_is_b::Bool,
    frag_is_p::Bool,
    frag_index::UInt8,
    frag_charge::UInt8,
    frag_isotope::UInt8,
    frag_internal::Bool,
    frag_immonium::Bool,
    frag_neutral_diff::Bool,
    frag_mz::Float32,
    frag_bounds::FragBoundModel,
    prec_mz::Float32,
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8)
    
    #println("frag_neutral_diff $frag_neutral_diff")
    min_frag_mz, max_frag_mz = frag_bounds(prec_mz)
    if (frag_mz < min_frag_mz) | (frag_mz > max_frag_mz)
        return false
    end
    if frag_is_y
        if frag_index < y_start
            return false
        end
    end
    if frag_is_b
        if frag_index < b_start
            return false
        end
    end
    if frag_is_p
        if !include_p
            return false
        end
    end
    if frag_immonium
        if !include_immonium
            return false
        end
    end
    if frag_internal
        if !include_internal
            return false
        end
    end
    if frag_neutral_diff
        if !include_neutral_diff
            return false
        end
    end
    if !include_isotope
        if !iszero(frag_isotope)
            return false
        end
    end
    if frag_charge > max_frag_charge
        return false
    end
    return true
end

"""
    getSimpleFrags(
        frag_mz::AbstractVector{Float32},
        frag_is_y::AbstractVector{Bool},
        frag_is_b::AbstractVector{Bool},
        frag_is_p::AbstractVector{Bool},
        frag_index::AbstractVector{UInt8},
        frag_charge::AbstractVector{UInt8},
        frag_isotope::AbstractVector{UInt8},
        frag_internal::AbstractVector{Bool},
        frag_immonium::AbstractVector{Bool},
        frag_neutral_diff::AbstractVector{Bool},
        precursor_mz::AbstractVector{Float32},
        precursor_irt::AbstractVector{Float32},
        precursor_charge::AbstractVector{UInt8},
        prec_to_frag_idx::AbstractVector{UInt64},
        y_start::UInt8,
        b_start::UInt8,
        include_p::Bool,
        include_isotope::Bool,
        include_immonium::Bool,
        include_internal::Bool,
        include_neutral_diff::Bool,
        max_frag_charge::UInt8,
        frag_bounds::FragBoundModel,
    )::Vector{SimpleFrag{Float32}}

Extract fragments for the fragment index from raw fragment data.

# Parameters
- `frag_mz`: Fragment m/z values
- `frag_is_y`: Whether each fragment is a y-ion
- `frag_is_b`: Whether each fragment is a b-ion  
- `frag_is_p`: Whether each fragment is a precursor ion
- `frag_index`: Index of each fragment in its peptide sequence
- `frag_charge`: Charge state of each fragment
- `frag_isotope`: Isotope state of each fragment
- `frag_internal`: Whether each fragment is an internal fragment
- `frag_immonium`: Whether each fragment is an immonium ion
- `frag_neutral_diff`: Whether each fragment has neutral losses
- `precursor_mz`: m/z values of precursors
- `precursor_irt`: Retention time values of precursors
- `precursor_charge`: Charge states of precursors
- `prec_to_frag_idx`: Index mapping precursors to their fragments
- `y_start`: Minimum index for y-ions to include
- `b_start`: Minimum index for b-ions to include
- `include_p`: Whether to include precursor ions
- `include_isotope`: Whether to include isotope peaks
- `include_immonium`: Whether to include immonium ions
- `include_internal`: Whether to include internal fragment ions
- `include_neutral_diff`: Whether to include fragments with neutral losses
- `max_frag_charge`: Maximum fragment charge state to include
- `frag_bounds`: Model defining valid m/z range based on precursor m/z

# Returns
- Vector of SimpleFrag objects, containing filtered fragments for the index
"""
function getSimpleFrags(
    frag_mz::AbstractVector{Float32},
    frag_is_y::AbstractVector{Bool},
    frag_is_b::AbstractVector{Bool},
    frag_is_p::AbstractVector{Bool},
    frag_index::AbstractVector{UInt8},
    frag_charge::AbstractVector{UInt8},
    frag_isotope::AbstractVector{UInt8},
    frag_internal::AbstractVector{Bool},
    frag_immonium::AbstractVector{Bool},
    frag_neutral_diff::AbstractVector{Bool},
    precursor_mz::AbstractVector{Float32},
    precursor_irt::AbstractVector{Float32},
    precursor_charge::AbstractVector{UInt8},
    prec_to_frag_idx::AbstractVector{UInt64},
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8,
    frag_bounds::FragBoundModel,
    )
    if (length(prec_to_frag_idx) - 1) != (length(precursor_mz))
        #println("mistake")
    end
    # Score is a per-rank bitmask packed into a UInt8 — 8 ranks max.
    max_rank_index = 8
    #Number of precursors
    n_precursors = UInt32(length(precursor_mz))
    simple_frags = Vector{SimpleFrag{Float32}}(undef, n_precursors*max_rank_index)
    simple_frag_idx = 0
    for pid in range(one(UInt32), n_precursors)
        prec_mz = precursor_mz[pid]
        frag_start_idx, frag_stop_idx = prec_to_frag_idx[pid], prec_to_frag_idx[pid+1] - 1
        rank = 1
        for frag_idx in range(frag_start_idx, frag_stop_idx)
            if fragFilter(
                    frag_is_y[frag_idx],
                    frag_is_b[frag_idx],
                    frag_is_p[frag_idx],
                    frag_index[frag_idx],
                    frag_charge[frag_idx],
                    frag_isotope[frag_idx],
                    frag_internal[frag_idx],
                    frag_immonium[frag_idx],
                    frag_neutral_diff[frag_idx],
                    frag_mz[frag_idx],
                    frag_bounds,
                    prec_mz,
                    y_start,
                    b_start,
                    include_p,
                    include_isotope,
                    include_immonium,
                    include_internal,
                    include_neutral_diff,
                    max_frag_charge)==false
                continue
            end
            simple_frag_idx += 1
            simple_frags[simple_frag_idx] = SimpleFrag(
                frag_mz[frag_idx],
                pid,
                precursor_mz[pid],
                precursor_irt[pid],
                precursor_charge[pid],
                UInt8(1) << UInt8(rank - 1)  # bitmask: rank 1→bit0, rank 2→bit1, ...
            )
            rank += 1
            if rank > max_rank_index  # cap at 8 bits
                break
            end
        end

    end
    return simple_frags[1:simple_frag_idx]
end

"""
    buildFragmentIndex!(
        folder_out::String,
        frag_ions::Vector{SimpleFrag{T}}, 
        frag_bin_tol_ppm::AbstractFloat, 
        rt_bin_tol::AbstractFloat;
        index_name::String = ""
    ) where {T<:AbstractFloat}

Build a hierarchical fragment index from a list of fragments.

# Parameters
- `folder_out`: Path to the directory where index files will be written
- `frag_ions`: Vector of SimpleFrag objects to index
- `frag_bin_tol_ppm`: Fragment binning tolerance in parts per million
- `rt_bin_tol`: Retention time bin width
- `index_name`: Prefix for output file names (default: "")

# Effects
Creates and writes the following files to `folder_out`:
- `index_name`f_index_fragments.arrow: Fragment entries sorted by retention time and m/z
- `index_name`f_index_rt_bins.arrow: Two-level binning structure for retention time
- `index_name`f_index_fragment_bins.arrow: Fragment m/z bins within each RT bin

# Returns
- `nothing`

# Implementation Notes
The indexing is hierarchical:
1. First bins by retention time
2. Within each RT bin, fragments are sorted by m/z
3. Then bins by m/z within each RT bin
This enables efficient fragment lookup during search.
"""
function buildFragmentIndex!(
                            folder_out::String,
                            frag_ions::Vector{SimpleFrag{T}}, 
                            frag_bin_tol_ppm::AbstractFloat, 
                            rt_bin_tol::AbstractFloat;
                            index_name::String = ""
                            ) where {T<:AbstractFloat}

    function buildFragIndex!(
        index_fragments::Vector{IndexFragment{T}}, 
        rt_bins::Vector{FragIndexBin{T}},
        frag_bins::Vector{FragIndexBin{T}},
        frag_ions::Vector{SimpleFrag{T}}, 
        frag_bin_tol_ppm::AbstractFloat,
        rt_bin_tol::AbstractFloat) where {T<:AbstractFloat}

        start_irt, stop_irt = getIRT(first(frag_ions)), getIRT(first(frag_ions))
        start_idx, stop_idx = 1, 1
        rt_bin_idx = 1
        frag_bin_idx = 1
        #Within each iRT bin (defined by `rt_size`) sort the fragments by precursor_mz

        for (i, frag_ion) in enumerate(frag_ions)
            stop_irt = getIRT(frag_ion)
            stop_idx = i
            #diff_mz = stop_fragmz - start_fragmz
            #mean_mz = (stop_fragmz + start_fragmz)/2
            #Is the difference between the first and last fragment greater than the bin_ppm?
            if ((stop_irt - start_irt) > rt_bin_tol) & (stop_idx != start_idx)#(diff_mz/(mean_mz/1e6) > frag_bin_tol_ppm) & (stop_idx != start_idx)
                stop_idx = i - 1 #i - 1 fragment is the last that should be incluced in the bin
                stop_irt = getIRT(frag_ions[stop_idx])
                #Within the fragment bin Sort by IRT
                sort!(@view(frag_ions[start_idx:stop_idx]), by = x->getMZ(x))#Try stable sorting algorithm for now, alg=QuickSort)
                #Add new fragbin
                #Build fragment bins for the retention time bin 
                first_frag_bin_idx = frag_bin_idx
                frag_bin_idx = buildFragBins!(index_fragments,
                            frag_bins,
                            frag_bin_idx,
                            frag_ions,
                            start_idx,stop_idx,
                            frag_bin_tol_ppm)
                rt_bins[rt_bin_idx] = FragIndexBin(start_irt, 
                                                    stop_irt, #important that stop_idx is i - 1 and not i
                                                        UInt32(first_frag_bin_idx),
                                                        UInt32(frag_bin_idx-1)
                                                        ) #-1 is critical
                rt_bin_idx += 1

                start_idx, stop_idx = i, i
                start_irt = getIRT(frag_ion)
            end
        end


        #Last bin is special case 
        if start_idx != length(frag_ions)
            stop_irt =  getIRT(frag_ions[stop_idx])
            sort!(@view(frag_ions[start_idx:stop_idx]), by = x->getMZ(x))# Try stable sorting for now, alg=QuickSort)
            #Add new fragbin
            #Build RT bins for the frag bin
            first_frag_bin_idx = frag_bin_idx
            frag_bin_idx = buildFragBins!(index_fragments,
                        frag_bins,
                        frag_bin_idx,
                        frag_ions,
                        start_idx,stop_idx,frag_bin_tol_ppm)

            rt_bins[rt_bin_idx] = FragIndexBin(start_irt, 
                                                    stop_irt, #important that stop_idx is i - 1 and not i
                                                    UInt32(first_frag_bin_idx),
                                                    UInt32(frag_bin_idx-1)) #-1 is critical
            rt_bin_idx += 1
        else
            first_frag_bin_idx = frag_bin_idx

            frag_bin_idx = buildFragBins!(index_fragments,
                        frag_bins,
                        frag_bin_idx,
                        frag_ions,
                        start_idx,stop_idx,frag_bin_tol_ppm)

            #Add new fragbin
            rt_bins[rt_bin_idx] = FragIndexBin(start_irt, 
                                                    getIRT(frag_ions[stop_idx]), #important that stop_idx is i - 1 and not i
                                                    UInt32(first_frag_bin_idx),
                                                    UInt32(frag_bin_idx-1)) #-1 is critical
            rt_bin_idx += 1
        end

        return frag_bin_idx, rt_bin_idx
    end

    function buildFragBins!(index_fragments::Vector{IndexFragment{T}}, 
                            frag_bins::Vector{FragIndexBin{T}},
                            frag_bin_idx::Int64,
                            frag_ions::Vector{SimpleFrag{T}}, 
                            start::Int64, 
                            stop::Int64, 
                            frag_bin_tol_ppm::AbstractFloat) where {T<:AbstractFloat}
        start_idx, stop_idx = start, start
        start_fragmz, stop_fragmz = getMZ(frag_ions[start]), getMZ(frag_ions[start])
        for i in range(start, stop)
            frag_ion = frag_ions[i]
            stop_fragmz = getMZ(frag_ion)
            stop_idx = i

            diff_mz = stop_fragmz - start_fragmz
            mean_mz = (stop_fragmz + start_fragmz)/2
            if (diff_mz/(mean_mz/1e6) > frag_bin_tol_ppm) & (stop_idx != start_idx)
                stop_idx = i - 1
                stop_fragmz =  getMZ(frag_ions[stop_idx]) #Need to set before sorting 
                sort!(@view(frag_ions[start_idx:stop_idx]), by = x->getPrecMZ(x))# Try stable sorting for now, alg=QuickSort)
                #Add new rt bin
                frag_bins[frag_bin_idx] = FragIndexBin(start_fragmz, 
                                                        stop_fragmz, #important that stop_idx is i - 1 and not i
                                                    UInt32(start_idx),
                                                    UInt32(stop_idx)
                                                )
                frag_bin_idx += 1
                for idx in range(start_idx, stop_idx)
                    index_fragments[idx] = IndexFragment(
                                                    getPrecID(frag_ions[idx]),
                                                    getPrecMZ(frag_ions[idx]),
                                                    getScore(frag_ions[idx]),
                                                    getPrecCharge(frag_ions[idx])
                                                    )
                end
                start_idx, stop_idx = i, i
                start_fragmz = getMZ(frag_ions[stop_idx])
            end
        end

        #Last bin is special case 
        if start_idx != stop
            stop_fragmz = getMZ(frag_ions[stop_idx])
            sort!(@view(frag_ions[start_idx:stop_idx]), by = x->getPrecMZ(x))# Try stable sorting for now, alg=QuickSort), alg=QuickSort)
            #Add new fragbin
            frag_bins[frag_bin_idx] = FragIndexBin(start_fragmz, 
                        stop_fragmz, #important that stop_idx is i - 1 and not i
                        UInt32(start_idx),
                        UInt32(stop_idx)
                    )
            frag_bin_idx += 1
            for idx in range(start_idx, stop_idx)
                index_fragments[idx] = IndexFragment(
                                                getPrecID(frag_ions[idx]),
                                                getPrecMZ(frag_ions[idx]),
                                                getScore(frag_ions[idx]),
                                                getPrecCharge(frag_ions[idx])
                                                )
            end
        else
            frag_bins[frag_bin_idx] = FragIndexBin(start_fragmz, 
                        getMZ(frag_ions[stop]), #important that stop_idx is i - 1 and not i
                        UInt32(start_idx),
                        UInt32(stop_idx)
                    )
            frag_bin_idx += 1
            index_fragments[stop_idx] = IndexFragment(
                getPrecID(frag_ions[stop_idx]),
                getPrecMZ(frag_ions[stop_idx]),
                getScore(frag_ions[stop_idx]),
                getPrecCharge(frag_ions[stop_idx])
                )
        end
        return frag_bin_idx
    end

    sort!(frag_ions, by = x->getIRT(x))# Try stable sorting for now, alg=QuickSort), alg=QuickSort)
    #diff = getPPM(getMZ(frag_ions[start]), bin_ppm) #ppm tolerance of the current fragment bin
    #Get smallest iRT in the library
    index_fragments = Vector{IndexFragment{T}}(undef, length(frag_ions))
    rt_bins = Vector{FragIndexBin{T}}(undef, length(frag_ions))
    frag_bins = Vector{FragIndexBin{T}}(undef, length(frag_ions))
    #println("building fragment index...")
    # Legacy Arrow-based index (not used at search time — LocalPartitionedFragmentIndex
    # from .jls files is used instead). Use fallback ppm if mDa mode is active.
    legacy_ppm = frag_bin_tol_ppm > 0 ? frag_bin_tol_ppm : 2.5f0
    frag_bin_idx, rt_bin_idx = buildFragIndex!(index_fragments,
                    rt_bins,
                    frag_bins,
                    frag_ions,
                    legacy_ppm,
                    rt_bin_tol)
                    
    fragments = (IndexFragment = index_fragments,)
    rt_bins  = (FragIndexBin = rt_bins[1:rt_bin_idx-1],)
    frag_bins = (FragIndexBin = frag_bins[1:frag_bin_idx-1],)
    #println("writing tables...")
    if !isdir(folder_out)
        mkdir(folder_out)
    end
    Arrow.write(joinpath(folder_out, index_name*"f_index_fragments.arrow"), fragments)
    Arrow.write(joinpath(folder_out, index_name*"f_index_rt_bins.arrow"), rt_bins)
    Arrow.write(joinpath(folder_out, index_name*"f_index_fragment_bins.arrow"), frag_bins)

    return 
end

"""
    Helper function to decide if the fragment should be filtered out. Filters based on maximum rank or
        a constant multiple of the peptide length, whichever is least 
"""
function filterFrag(rank::Int64, prec_len::UInt8, max_frag_rank::UInt8, length_to_frag_count_multiple::AbstractFloat)
    return rank > min(max_frag_rank, round((prec_len)*length_to_frag_count_multiple)+1)
end


"""
    getDetailedFrags(
        frag_mz::AbstractVector{Float32},
        frag_coef::AbstractVector{NTuple{N, Float32}},
        frag_intensity::AbstractVector{Float16},
        frag_is_y::AbstractVector{Bool},
        frag_is_b::AbstractVector{Bool},
        frag_is_p::AbstractVector{Bool},
        frag_index::AbstractVector{UInt8},
        frag_charge::AbstractVector{UInt8},
        frag_sulfur_count::AbstractVector{UInt8},
        frag_ion_type::AbstractVector{UInt16},
        frag_isotope::AbstractVector{UInt8},
        frag_internal::AbstractVector{Bool},
        frag_immonium::AbstractVector{Bool},
        frag_neutral_diff::AbstractVector{Bool},
        precursor_mz::AbstractVector{Float32},
        precursor_charge::AbstractVector{UInt8},
        prec_to_frag_idx::AbstractVector{UInt64},
        y_start::UInt8,
        b_start::UInt8,
        include_p::Bool,
        include_isotope::Bool,
        include_immonium::Bool,
        include_internal::Bool,
        include_neutral_diff::Bool,
        max_frag_charge::UInt8,
        frag_bounds::FragBoundModel,
        max_frag_rank::UInt8,
        min_frag_intensity::AbstractFloat,
        koina_model::SplineCoefficientModel
    )::Tuple{Vector{SplineCompactFrag{N, Float32}}, Vector{UInt64}} where {N}

Extract detailed fragment information for peptide scoring, including spline coefficients.

# Parameters
- `frag_mz`: Fragment m/z values
- `frag_coef`: Spline coefficients for each fragment
- `frag_intensity`: Fragment intensities
- `frag_is_y`: Whether each fragment is a y-ion
- `frag_is_b`: Whether each fragment is a b-ion
- `frag_is_p`: Whether each fragment is a precursor ion
- `frag_index`: Index of each fragment in its peptide sequence
- `frag_charge`: Charge state of each fragment
- `frag_sulfur_count`: Sulfur count for each fragment
- `frag_ion_type`: Ion type identifier for each fragment
- `frag_isotope`: Isotope state of each fragment
- `frag_internal`: Whether each fragment is an internal fragment
- `frag_immonium`: Whether each fragment is an immonium ion
- `frag_neutral_diff`: Whether each fragment has neutral losses
- `precursor_mz`: m/z values of precursors
- `precursor_charge`: Charge states of precursors
- `prec_to_frag_idx`: Index mapping precursors to their fragments
- `y_start`: Minimum index for y-ions to include
- `b_start`: Minimum index for b-ions to include
- `include_p`: Whether to include precursor ions
- `include_isotope`: Whether to include isotope peaks
- `include_immonium`: Whether to include immonium ions
- `include_internal`: Whether to include internal fragment ions
- `include_neutral_diff`: Whether to include fragments with neutral losses
- `max_frag_charge`: Maximum fragment charge state to include
- `frag_bounds`: Model defining valid m/z range based on precursor m/z
- `max_frag_rank`: Maximum number of fragments per precursor (ranked by intensity)
- `min_frag_intensity`: Minimum relative intensity threshold for fragments
- `koina_model`: Type of spline coefficient model used

# Returns
- Tuple containing:
  - Vector of SplineDetailedFrag objects with fragment info including spline coefficients
  - Vector mapping precursor indices to fragment indices

# Type Parameters
- `N`: Number of spline coefficients per fragment
"""
function getDetailedFrags(
    frag_mz::AbstractVector{Float32},
    frag_coef::AbstractVector{NTuple{N, Float32}},
    frag_intensity::AbstractVector{Float16},
    frag_is_y::AbstractVector{Bool},
    frag_is_b::AbstractVector{Bool},
    frag_is_p::AbstractVector{Bool},
    frag_index::AbstractVector{UInt8},
    frag_charge::AbstractVector{UInt8},
    frag_sulfur_count::AbstractVector{UInt8},
    frag_ion_type::AbstractVector{UInt16},
    frag_isotope::AbstractVector{UInt8},
    frag_internal::AbstractVector{Bool},
    frag_immonium::AbstractVector{Bool},
    frag_neutral_diff::AbstractVector{Bool},
    precursor_mz::AbstractVector{Float32},
    precursor_charge::AbstractVector{UInt8},
    prec_len::AbstractVector{UInt8},
    prec_to_frag_idx::AbstractVector{UInt64},
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8,
    frag_bounds::FragBoundModel,
    max_frag_rank::UInt8,
    length_to_frag_count_multiple::AbstractFloat,
    min_frag_intensity::AbstractFloat,
    koina_model::SplineCoefficientModel) where {N}

    # As of §9, every fragment present in the input table is already known to
    # have passed all per-fragment metadata filters and the per-precursor
    # rank cap — `filter_fragments!(::SplineCoefficientModel, ::SplineFragFilterCtx)`
    # in fragment_predict.jl drops anything that wouldn't survive. The filter
    # knobs (`y_start`, `b_start`, `include_*`, `max_frag_charge`,
    # `max_frag_rank`, `length_to_frag_count_multiple`, `frag_bounds`,
    # `min_frag_intensity`) and `precursor_mz`/`prec_len` remain in this
    # signature for backward compatibility with existing callers but are no
    # longer consulted here.
    n_precursors = UInt32(length(prec_to_frag_idx) - 1)
    n_frags = UInt64(length(frag_mz))

    n_tuple_p = typeof(first(frag_coef))
    n_tuple_size = length(n_tuple_p.parameters)
    n_tuple_type = eltype(n_tuple_p)
    detailed_frags = Vector{SplineCompactFrag{n_tuple_size, n_tuple_type}}(undef, n_frags)
    prec_to_frag_idx_new = Vector{UInt64}(undef, n_precursors + 1)
    detailed_frag_idx = 1

    for pid in ProgressBar(range(one(UInt32), n_precursors))
        prec_to_frag_idx_new[pid] = UInt64(detailed_frag_idx)
        frag_start_idx, frag_stop_idx = prec_to_frag_idx[pid], prec_to_frag_idx[pid+1] - 1
        rank = 1
        for frag_idx in range(frag_start_idx, frag_stop_idx)
            detailed_frags[detailed_frag_idx] = SplineCompactFrag(
                UInt32(pid),
                frag_mz[frag_idx],
                frag_coef[frag_idx],
                frag_is_y[frag_idx],
                frag_is_b[frag_idx],
                frag_is_p[frag_idx],
                frag_isotope[frag_idx] > 0,
                frag_charge[frag_idx],
                frag_index[frag_idx],
                precursor_charge[pid],
                UInt8(rank),
                frag_sulfur_count[frag_idx],
            )
            detailed_frag_idx += 1
            rank += 1
        end
    end
    prec_to_frag_idx_new[end] = UInt64(detailed_frag_idx)
    return detailed_frags[1:detailed_frag_idx-1], prec_to_frag_idx_new
end
