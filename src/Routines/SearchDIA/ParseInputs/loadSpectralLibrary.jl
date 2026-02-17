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

using JSON

#==========================================================
Batched Library Loading Support
==========================================================#

"""
    load_library_config(config_path::String) -> BatchedLibraryConfig

Load batched library configuration from JSON file.
"""
function load_library_config(config_path::String)::BatchedLibraryConfig
    config_data = JSON.parsefile(config_path)

    batch_ranges = [(Float32(r[1]), Float32(r[2])) for r in config_data["batch_mz_ranges"]]
    batch_min_charges = [Int(c) for c in config_data["batch_min_charges"]]

    return BatchedLibraryConfig(
        Int(config_data["n_batches"]),
        Int(config_data["n_precursors"]),
        Int(config_data["batch_size"]),
        batch_ranges,
        batch_min_charges,
        String(config_data["library_base_path"]),
        Bool(config_data["is_batched"])
    )
end

"""
    save_library_config(config_path::String, config::BatchedLibraryConfig)

Save batched library configuration to JSON file.
"""
function save_library_config(config_path::String, config::BatchedLibraryConfig)
    config_data = Dict(
        "n_batches" => config.n_batches,
        "n_precursors" => config.n_precursors,
        "batch_size" => config.batch_size,
        "batch_mz_ranges" => config.batch_mz_ranges,
        "batch_min_charges" => config.batch_min_charges,
        "library_base_path" => config.library_base_path,
        "is_batched" => config.is_batched
    )
    open(config_path, "w") do f
        JSON.print(f, config_data, 2)
    end
end

"""
    get_batch_path(config::BatchedLibraryConfig, batch_idx::Int) -> String

Get the directory path for a specific batch.
"""
function get_batch_path(config::BatchedLibraryConfig, batch_idx::Int)::String
    return joinpath(config.library_base_path, "batch_$batch_idx")
end

"""
    loadBatchedSpectralLibrary(lib_dir::String, config::BatchedLibraryConfig, params::PioneerParameters)

Load a batched spectral library. Only loads shared data initially; batches are loaded on demand.

No dictionary needed for precursor->batch mapping! The config contains batch_size,
so we can compute batch_idx = (prec_idx - 1) ÷ batch_size + 1 in O(1).
"""
function loadBatchedSpectralLibrary(
    lib_dir::String,
    config::BatchedLibraryConfig,
    params::PioneerParameters
)::BatchedSpectralLibrary

    # Load shared protein data
    proteins = Arrow.Table(joinpath(lib_dir, "proteins_table.arrow"))

    return BatchedSpectralLibrary(
        config,
        0,  # No batch loaded initially
        nothing,
        SetProteins(proteins)
    )
end

"""
    load_batch!(lib::BatchedSpectralLibrary, batch_idx::Int, params::PioneerParameters)

Load a specific batch into memory. Unloads any previously loaded batch first.
"""
function load_batch!(lib::BatchedSpectralLibrary, batch_idx::Int, params::PioneerParameters)
    # Check if already loaded
    if lib.current_batch_idx == batch_idx && lib.current_batch !== nothing
        return nothing
    end

    # Unload current batch
    unload_current_batch!(lib)

    # Load new batch using existing loadSpectralLibrary_internal for the batch directory
    batch_path = get_batch_path(lib.config, batch_idx)
    lib.current_batch = loadSpectralLibrary_internal(batch_path, params)
    lib.current_batch_idx = batch_idx

    @info "Loaded batch $batch_idx of $(lib.config.n_batches)"

    return nothing
end

"""
    unload_current_batch!(lib::BatchedSpectralLibrary)

Unload the current batch to free memory.
"""
function unload_current_batch!(lib::BatchedSpectralLibrary)
    if lib.current_batch !== nothing
        lib.current_batch = nothing
        lib.current_batch_idx = 0
        GC.gc()  # Trigger garbage collection to free memory
    end
    return nothing
end

#==========================================================
Internal Library Loading (Renamed from loadSpectralLibrary)
==========================================================#

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
    loadSpectralLibrary(SPEC_LIB_DIR::String, params::PioneerParameters) -> SpectralLibrary

Load a spectral library from the given directory. Automatically detects whether this is a
batched library (has library_config.json) or a monolithic library.

For batched libraries, returns a BatchedSpectralLibrary that loads batches on demand.
For monolithic libraries, returns a FragmentIndexLibrary or SplineFragmentIndexLibrary.
"""
function loadSpectralLibrary(SPEC_LIB_DIR::String, params::PioneerParameters)
    # Check for batched library
    config_path = joinpath(SPEC_LIB_DIR, "library_config.json")

    if isfile(config_path)
        config = load_library_config(config_path)
        if config.is_batched
            @info "Detected batched spectral library with $(config.n_batches) batches"
            return loadBatchedSpectralLibrary(SPEC_LIB_DIR, config, params)
        end
    end

    # Original monolithic loading
    return loadSpectralLibrary_internal(SPEC_LIB_DIR, params)
end

"""
    loadSpectralLibrary_internal(SPEC_LIB_DIR::String, params::PioneerParameters)

Internal function to load a monolithic (non-batched) spectral library.
This is also used to load individual batches within a batched library.
"""
function loadSpectralLibrary_internal(SPEC_LIB_DIR::String,
                             params::PioneerParameters)
    f_index_fragments = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_fragments.arrow"))
    f_index_rt_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_rt_bins.arrow"))
    f_index_frag_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_fragment_bins.arrow"))


    presearch_f_index_fragments = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_fragments.arrow"))
    presearch_f_index_rt_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_rt_bins.arrow"))
    presearch_f_index_frag_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_fragment_bins.arrow"))


    # Note: Can't use @user_info here as LoggingSystem is loaded after ParseInputs
    # This message will be captured by the logging system when SearchDIA runs
    # println("Loading spectral library from $SPEC_LIB_DIR into main memory...")
    spec_lib = Dict{String, Any}()
    detailed_frags = load_detailed_frags(joinpath(SPEC_LIB_DIR,"detailed_fragments.jld2"))
    prec_frag_ranges = load(joinpath(SPEC_LIB_DIR,"precursor_to_fragment_indices.jld2"))["pid_to_fid"]
    library_fragment_lookup_table = nothing
    if (eltype(detailed_frags).name == (eltype(Vector{SplineDetailedFrag{4, Float32}}(undef, 0)).name))
        try
            spl_knots = load(joinpath(SPEC_LIB_DIR,"spline_knots.jld2"))["spl_knots"]
            library_fragment_lookup_table = SplineFragmentLookup(
                detailed_frags,
                prec_frag_ranges,
                Tuple(spl_knots),
                3
            )
        catch e
            @user_warn "Could not load `spline_knots.jld2`"
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

    f_index = FragmentIndex(
        f_index_frag_bins[:FragIndexBin],
        f_index_rt_bins[:FragIndexBin],
        f_index_fragments[:IndexFragment],
    );
    presearch_f_index = FragmentIndex(
        presearch_f_index_frag_bins[:FragIndexBin],
        presearch_f_index_rt_bins[:FragIndexBin],
        presearch_f_index_fragments[:IndexFragment],
    );
    spec_lib["f_index"] = f_index;
    spec_lib["presearch_f_index"] = presearch_f_index;
    spec_lib["precursors"] = precursors;
    spec_lib["proteins"] = proteins;

    if typeof(library_fragment_lookup_table) == Pioneer.StandardFragmentLookup{Float32}
        return FragmentIndexLibrary(
            spec_lib["presearch_f_index"], 
            spec_lib["f_index"], 
            SetPrecursors(spec_lib["precursors"]), 
            SetProteins(spec_lib["proteins"]),
            spec_lib["f_det"]
        )
    else
        return SplineFragmentIndexLibrary(
            spec_lib["presearch_f_index"], 
            spec_lib["f_index"], 
            SetPrecursors(spec_lib["precursors"]), 
            SetProteins(spec_lib["proteins"]),
            spec_lib["f_det"]
        )
    end
end