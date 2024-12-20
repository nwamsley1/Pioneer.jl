
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
                rank_to_score::Vector{UInt8},
                frag_bounds::FragBoundModel,
                frag_bin_tol_ppm::Float32,
                rt_bin_tol_ppm::Float32,
                model_type::KoinaModelType = InstrumentSpecificModel("default"))

Build Pioneer spectral library including fragment indices.

Parameters:
- spec_lib_path: Path to library directory
- y/b_start_index: Minimum index for index-building
- y/b_start: Minimum index for detailed fragments 
- include_*: Flags for fragment type inclusion
- max_frag_charge: Maximum fragment charge to include
- max_frag_rank: Maximum number of fragments per precursor
- min_frag_intensity: Minimum relative intensity threshold
- rank_to_score: Mapping from rank to scoring value
- frag_bounds: Model for fragment m/z bounds
- frag/rt_bin_tol_ppm: Binning tolerances
- model_type: Type of prediction model used
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
                     min_frag_intensity::Float32,
                     rank_to_score::Vector{UInt8},
                     frag_bounds::FragBoundModel,
                     frag_bin_tol_ppm::Float32,
                     rt_bin_tol_ppm::Float32,
                     model_type::KoinaModelType = InstrumentSpecificModel("default"))

    # Try to load required tables
    fragments_table = nothing
    prec_to_frag = nothing
    precursors_table = nothing
    
    try
        fragments_table = Arrow.Table(joinpath(spec_lib_path, "fragments_table.arrow"))
        prec_to_frag = Arrow.Table(joinpath(spec_lib_path, "prec_to_frag.arrow"))
        precursors_table = Arrow.Table(joinpath(spec_lib_path, "precursors_table.arrow"))
    catch e
        @error "Could not find library..." exception=e
        return nothing
    end

    # Get simple fragments for index building
    println("Get index fragments...")
    simple_frags = getSimpleFrags(
        fragments_table[:mz],
        fragments_table[:is_y],
        fragments_table[:is_b],
        fragments_table[:is_p],
        fragments_table[:frag_index],
        fragments_table[:charge],
        fragments_table[:isotope],
        fragments_table[:internal],
        fragments_table[:immonium],
        fragments_table[:has_neutral_diff],
        precursors_table[:mz],
        precursors_table[:irt],
        precursors_table[:prec_charge],
        prec_to_frag[:start_idx],
        y_start_index,
        b_start_index,
        include_p_index,
        include_isotope,
        include_immonium,
        include_internal,
        include_neutral_diff,
        max_frag_charge,
        frag_bounds,
        rank_to_score
    )

    println("Build fragment index...")
    # Main RT-binned index
    sort!(simple_frags, by = x -> x.prec_irt)
    buildFragmentIndex!(
        spec_lib_path,
        simple_frags,
        frag_bin_tol_ppm,
        rt_bin_tol_ppm,
        index_name = ""
    )

    println("Build presearch fragment index...")
    # Pre-search index without RT binning
    sort!(simple_frags, by = x -> x.prec_irt)
    buildFragmentIndex!(
        spec_lib_path,
        simple_frags,
        frag_bin_tol_ppm,
        typemax(Float32),  # No RT binning
        index_name = "presearch_"
    )

    simple_frags = nothing
    GC.gc()

    println("Get full fragments list...")
    if model_type isa SplineCoefficientModel
        detailed_frags, pid_to_fid = getDetailedFrags(
            fragments_table[:mz],
            fragments_table[:coefficients],
            fragments_table[:intensity],
            fragments_table[:is_y],
            fragments_table[:is_b],
            fragments_table[:is_p],
            fragments_table[:frag_index],
            fragments_table[:charge],
            fragments_table[:sulfur_count],
            fragments_table[:ion_type],
            fragments_table[:isotope],
            fragments_table[:internal],
            fragments_table[:immonium],
            fragments_table[:has_neutral_diff],
            precursors_table[:mz],
            precursors_table[:prec_charge],
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
            min_frag_intensity,
            model_type
        )
    else
        detailed_frags, pid_to_fid = getDetailedFrags(
            fragments_table[:mz],
            fragments_table[:intensity],
            fragments_table[:is_y],
            fragments_table[:is_b],
            fragments_table[:is_p],
            fragments_table[:frag_index],
            fragments_table[:charge],
            fragments_table[:sulfur_count],
            fragments_table[:ion_type],
            fragments_table[:isotope],
            fragments_table[:internal],
            fragments_table[:immonium],
            fragments_table[:has_neutral_diff],
            precursors_table[:mz],
            precursors_table[:prec_charge],
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
            min_frag_intensity,
            model_type
        )
    end

    # Save detailed fragments
    save_detailed_frags(
        joinpath(spec_lib_path, "detailed_fragments.jld2"),
        detailed_frags
    )

    # Save fragment indices
    jldsave(
        joinpath(spec_lib_path, "precursor_to_fragment_indices.jld2");
        pid_to_fid
    )

    return nothing
end
