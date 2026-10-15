using Pioneer, Arrow, JSON
include("CompactLibraryStorage.jl")
using .CompactLibraryStorage

"""Assert semantic equality, including SIMD padding and the complete rank mask."""
function verify_index(original, restored)
    @assert typeof(original) == typeof(restored)
    @assert original.partition_bounds == restored.partition_bounds
    @assert original.n_partitions == restored.n_partitions
    for (a, b) in zip(original.partitions, restored.partitions)
        @assert a.fragments == b.fragments
        @assert a.rt_bins == b.rt_bins
        @assert a.local_to_global == b.local_to_global
        @assert a.n_local_precs == b.n_local_precs
        @assert a.skip_hints == b.skip_hints
        for field in fieldnames(typeof(a.fragment_bins))
            @assert isequal(getfield(a.fragment_bins, field), getfield(b.fragment_bins, field))
        end
    end
end

"""
    verify_reference(baseline_directory, output_directory)

Compare the reference codecs against a saved develop benchmark library. This
is an offline verification script, not a legacy reader in the new format.
Output must be empty; input files are never modified. No network is used.
"""
function verify_reference(baseline_directory, output_directory)
    mkpath(output_directory)
    isempty(readdir(output_directory)) || error("Output directory must be empty")
    source_table = Arrow.Table(joinpath(baseline_directory, "precursors_table.arrow"))
    precursor_path = joinpath(output_directory, "precursors_table.arrow")
    write_precursors(precursor_path, source_table)
    stored_table = read_precursors(precursor_path)
    for name in filter(!=(:koina_sequence), propertynames(source_table))
        @assert isequal(collect(getproperty(source_table, name)), collect(getproperty(stored_table, name)))
    end
    before, after = Pioneer.SetPrecursors(source_table), Pioneer.SetPrecursors(stored_table)
    @assert before.pid_to_cv_fold == after.pid_to_cv_fold
    a, b = Pioneer.getPrecursorTextIds(before), Pioneer.getPrecursorTextIds(after)
    for field in fieldnames(typeof(a))
        @assert isequal(getfield(a, field), getfield(b, field))
    end

    fragments = Pioneer.deserialize_from_jls(joinpath(baseline_directory, "detailed_fragments.jls"))
    offsets = Pioneer.deserialize_from_jls(joinpath(baseline_directory, "precursor_to_fragment_indices.jls"))
    detail_path = joinpath(output_directory, "detailed_fragments.jls")
    write_details(detail_path, fragments, offsets, stored_table.prec_charge)
    restored, ranges = read_details(detail_path, stored_table.prec_charge)
    @assert isequal(restored, fragments)
    @assert ranges == offsets
    @assert typeof(restored) == typeof(fragments)

    indexes = filter(name -> occursin("fragment_index", name) && endswith(name, ".jls"), readdir(baseline_directory))
    for name in indexes
        original = Pioneer.deserialize_from_jls(joinpath(baseline_directory, name))
        write_index(joinpath(output_directory, name), original)
        verify_index(original, read_index(joinpath(output_directory, name)))
    end
    for name in ("proteins_table.arrow", "config.json", "spline_knots.jls", "fragment_indices.json")
        source = joinpath(baseline_directory, name)
        isfile(source) && cp(source, joinpath(output_directory, name))
    end

    baseline_bytes = sum(filesize(joinpath(baseline_directory, name)) for name in
        vcat(indexes, ["precursors_table.arrow", "detailed_fragments.jls", "precursor_to_fragment_indices.jls"]))
    compact_bytes = sum(filesize(joinpath(output_directory, name)) for name in
        vcat(indexes, ["precursors_table.arrow", "detailed_fragments.jls"]))
    result = Dict(
        "precursors" => length(source_table.sequence), "fragments" => length(fragments),
        "indexes_checked" => length(indexes), "all_values_equal" => true,
        "develop_bytes" => baseline_bytes, "reference_bytes" => compact_bytes,
        "saved_percent" => 100 * (1 - compact_bytes / baseline_bytes)
    )
    write(joinpath(output_directory, "verification.json"), JSON.json(result, 2))
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 || error("Usage: julia --project reference/library_storage/verify_reference.jl BASELINE OUTPUT")
    println(JSON.json(verify_reference(ARGS[1], ARGS[2]), 2))
end
