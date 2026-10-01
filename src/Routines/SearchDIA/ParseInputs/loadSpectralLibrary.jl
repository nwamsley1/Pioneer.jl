# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.


#=
function load_detailed_frags(filename::String)
    jldopen(filename, "r") do file
        data = read(file, "data")
        spline_type_name = eltype(Vector{SplineDetailedFrag{4, Float32}}(undef, 0)).name
        if eltype(data).name != spline_type_name
            return map(x -> DetailedFrag{Float32}(
                x.prec_id,
                x.mz,
                x.intensity,
                x.ion_type,
                x.is_y,
                x.is_b,
                x.is_p,
                x.is_isotope,
                x.frag_charge,
                x.ion_position,
                x.prec_charge,
                x.rank,
                x.sulfur_count
            ), data)
        else
            return map(x -> eltype(data)(
                x.prec_id,
                x.mz,
                x.intensity,
                x.ion_type,
                x.is_y,
                x.is_b,
                x.is_p,
                x.is_isotope,
                x.frag_charge,
                x.ion_position,
                x.prec_charge,
                x.rank,
                x.sulfur_count
            ), data)
        end
    end
end
=#


"""
    median_ms2_isolation_width(ms_path) -> Union{Nothing, Float64}

The median MS2 isolation window width (m/z) of one MS data file; `nothing` without MS2 isolation metadata.
"""
function median_ms2_isolation_width(ms_path::AbstractString)
    spectra = loadMassSpecData(ms_path)
    orders = getMsOrders(spectra); widths = getIsolationWidthMzs(spectra)
    w = Float64[Float64(widths[i]) for i in eachindex(orders)
                if orders[i] == 2 && !ismissing(widths[i]) && isfinite(widths[i]) && widths[i] > 0]
    return isempty(w) ? nothing : median(w)
end

"""
    target_partition_width(window, tims) -> Float64

The precursor partition width (Da) that suits the data: always 10 Da for timsTOF (`tims`) data, whatever its
window width (10 Da was fastest on diaPASEF, and ion mobility spreads each window's precursors further);
otherwise the size class of the MS2 isolation window `window` (m/z): below 3.75 → 2.5 Da, below 7.5 → 5, else 10
(dev_docs/fragment_index/PARTITION_WIDTH_SWEEP.md). `nothing` when the window is unknown and the data not timsTOF.
"""
target_partition_width(window, tims::Bool) =
    tims ? 10.0 : window === nothing ? nothing : window < 3.75 ? 2.5 : window < 7.5 ? 5.0 : 10.0

"""
    choose_fragment_index(lib_dir, ms_paths) -> (main, presearch, width, window, tims)

The fragment index to search with. A library listing several (`fragment_indices.json`) offers one per precursor
partition width; the one nearest (in ratio) to `target_partition_width` of the data is used. The data is the first
MS file: a search's files share an acquisition method. Its window, the median MS2 isolation width, is read only
for non-timsTOF data. A library without the descriptor has its single index under the historical names (`width`
`nothing`); with no target (window unknown) the library's first index is used.
"""
function choose_fragment_index(lib_dir::AbstractString, ms_paths::AbstractVector{<:AbstractString})
    tims = !isempty(ms_paths) && is_tdfs_path(first(ms_paths))
    legacy = (main = "partitioned_fragment_index.jls", presearch = "presearch_partitioned_fragment_index.jls",
              width = nothing, window = nothing, tims = tims)
    desc_path = joinpath(lib_dir, FRAGMENT_INDEX_DESCRIPTOR)
    isfile(desc_path) || return legacy
    entries = JSON.parsefile(desc_path)["indexes"]
    window = isempty(ms_paths) || tims ? nothing : median_ms2_isolation_width(first(ms_paths))
    target = target_partition_width(window, tims)
    e = first(entries)
    if target !== nothing && length(entries) > 1
        e = entries[argmin([abs(log(Float64(x["partition_width_da"]) / target)) for x in entries])]
    end
    return (main = String(e["main"]), presearch = String(e["presearch"]),
            width = Float64(e["partition_width_da"]), window = window, tims = tims)
end

function loadSpectralLibrary(SPEC_LIB_DIR::String,
                             params::PioneerParameters;
                             fragment_index = choose_fragment_index(SPEC_LIB_DIR, String[]))
    # Note: Can't use @user_info here as LoggingSystem is loaded after ParseInputs
    # This message will be captured by the logging system when SearchDIA runs
    spec_lib = Dict{String, Any}()

    # Load detailed fragments and convert to CompactFrag for the search pipeline
    frag_path = joinpath(SPEC_LIB_DIR, "detailed_fragments")
    raw_frags = if isfile(frag_path * ".jls")
        deserialize_from_jls(frag_path * ".jls")
    elseif isfile(frag_path * ".jld2")
        @user_warn "Loading legacy JLD2 format for detailed_fragments. Consider rebuilding library."
        jldopen(frag_path * ".jld2", "r") do file
            read(file, "data")
        end
    else
        error("Fragment file not found: $(frag_path).jls or $(frag_path).jld2")
    end
    # Convert to compact types for the search pipeline.
    detailed_frags = if eltype(raw_frags) <: DetailedFrag
        map(CompactFrag, raw_frags)
    elseif eltype(raw_frags) <: SplineDetailedFrag
        map(SplineCompactFrag, raw_frags)
    elseif eltype(raw_frags) <: CompactFrag || eltype(raw_frags) <: SplineCompactFrag
        raw_frags  # already compact
    else
        raw_frags  # unknown type — pass through
    end

    # Load precursor-to-fragment indices with backwards compatibility
    prec_frag_ranges = if isfile(joinpath(SPEC_LIB_DIR, "precursor_to_fragment_indices.jls"))
        deserialize_from_jls(joinpath(SPEC_LIB_DIR, "precursor_to_fragment_indices.jls"))
    elseif isfile(joinpath(SPEC_LIB_DIR, "precursor_to_fragment_indices.jld2"))
        @user_warn "Loading legacy JLD2 format for precursor_to_fragment_indices. Consider rebuilding library."
        load(joinpath(SPEC_LIB_DIR, "precursor_to_fragment_indices.jld2"))["pid_to_fid"]
    else
        error("precursor_to_fragment_indices file not found in $SPEC_LIB_DIR")
    end

    library_fragment_lookup_table = nothing
    if eltype(detailed_frags) <: SplineCompactFrag || eltype(detailed_frags) <: SplineDetailedFrag
        try
            # Load spline knots with backwards compatibility
            spl_knots = if isfile(joinpath(SPEC_LIB_DIR, "spline_knots.jls"))
                deserialize_from_jls(joinpath(SPEC_LIB_DIR, "spline_knots.jls"))
            elseif isfile(joinpath(SPEC_LIB_DIR, "spline_knots.jld2"))
                @user_warn "Loading legacy JLD2 format for spline_knots. Consider rebuilding library."
                load(joinpath(SPEC_LIB_DIR, "spline_knots.jld2"))["spl_knots"]
            else
                error("spline_knots file not found in $SPEC_LIB_DIR")
            end
            library_fragment_lookup_table = SplineFragmentLookup(
                detailed_frags,
                prec_frag_ranges,
                Tuple(spl_knots)
            )
        catch e
            @user_warn "Could not load spline_knots"
            throw(e)
        end

    else
        library_fragment_lookup_table = StandardFragmentLookup(detailed_frags, prec_frag_ranges)
    end
    #Is this still necessary?
    #last_range = library_fragment_lookup_table.prec_frag_ranges[end] #0x29004baf:(0x29004be8 - 1)
    #last_range = range(first(last_range), last(last_range) - 1)
    #library_fragment_lookup_table.prec_frag_ranges[end] = last_range
    spec_lib["f_det"] = library_fragment_lookup_table

    precursors = Arrow.Table(joinpath(SPEC_LIB_DIR, "precursors_table.arrow"))
    proteins = Arrow.Table(joinpath(SPEC_LIB_DIR, "proteins_table.arrow"))

    # Load the partitioned fragment indexes of the chosen width (choose_fragment_index)
    partitioned_index = deserialize_from_jls(joinpath(SPEC_LIB_DIR, fragment_index.main))
    presearch_partitioned_index = deserialize_from_jls(joinpath(SPEC_LIB_DIR, fragment_index.presearch))

    # Load the BuildSpecLib config (if present) for output policy and for
    # reconstructing variable-modification counts in older precursor tables.
    build_config_path = joinpath(SPEC_LIB_DIR, "config.json")
    build_config = if isfile(build_config_path)
        try
            JSON.parsefile(build_config_path, dicttype=Dict{String,Any})
        catch
            nothing
        end
    else
        nothing
    end
    output_schema_policy = OutputSchemaPolicy(build_config)
    variable_mod_names = _configured_variable_mod_names(build_config)

    if typeof(library_fragment_lookup_table) == Pioneer.StandardFragmentLookup{Float32}
        return FragmentIndexLibrary(
            presearch_partitioned_index,
            partitioned_index,
            SetPrecursors(precursors; variable_mod_names = variable_mod_names),
            SetProteins(proteins),
            spec_lib["f_det"],
            output_schema_policy
        )
    else
        return SplineFragmentIndexLibrary(
            presearch_partitioned_index,
            partitioned_index,
            SetPrecursors(precursors; variable_mod_names = variable_mod_names),
            SetProteins(proteins),
            spec_lib["f_det"],
            output_schema_policy
        )
    end
end
