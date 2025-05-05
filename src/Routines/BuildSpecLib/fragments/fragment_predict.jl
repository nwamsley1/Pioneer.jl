# src/fragments/fragment_predict.jl

"""
    predict_fragments(
        peptide_table_path::String,
        frags_out_path::String,
        model_type::KoinaModelType,
        instrument_type::String, 
        max_koina_batches::Int,
        batch_size::Int,
        model_name::String;
        intensity_threshold::Float32 = 0.001f0
    )

General fragment prediction dispatcher that routes to appropriate method based on model type.
"""
function predict_fragments(
    peptide_table_path::String,
    frags_out_path::String,
    model_type::KoinaModelType,
    instrument_type::String,
    max_koina_batches::Int,
    batch_size::Int,
    model_name::String;
    intensity_threshold::Float32 = 0.001f0
)
    # Verify model configuration
    if !haskey(KOINA_URLS, model_name)
        error("Invalid model name: $model_name. Valid options: $(join(keys(KOINA_URLS), ", "))")
    end

    # Load data
    peptides_df = DataFrame(Arrow.Table(peptide_table_path))

    # Process in batches
    koina_pool_size = max_koina_batches * 5
    nprecs = nrow(peptides_df)
    batch_size = min(batch_size, 1000)
    batch_start_idxs = collect(one(UInt32):UInt32(batch_size*koina_pool_size):UInt32(nprecs))

    rm(frags_out_path, force=true)
    
    for start_idx in ProgressBar(batch_start_idxs)
        stop_idx = min(start_idx + batch_size*koina_pool_size - 1, nrow(peptides_df))
        batch_df = peptides_df[start_idx:stop_idx, :]
        
        # Generate predictions for batch
        frags_out = predict_fragments_batch(
            batch_df,
            model_type,
            instrument_type,
            batch_size,
            max_koina_batches,
            start_idx,
            
        )

        # Write or append results
        if start_idx == 1
            # Create file in stream format to support appending
            open(frags_out_path, "w") do io
                Arrow.write(io, frags_out; file=false)  # file=false creates stream format
            end
        else
            Arrow.append(frags_out_path, frags_out)
        end
    end
end

"""
Fragment prediction for instrument-specific models (e.g., UniSpec, AlphaPeptDeep).
"""
function predict_fragments_batch(
    peptides_df::DataFrame,
    model::InstrumentSpecificModel,
    instrument_type::String,
    batch_size::Int,
    concurrent_koina_requests::Int,
    first_prec_idx::UInt32
)::DataFrame
    # Verify instrument compatibility
    if instrument_type ∉ MODEL_CONFIGS[model.name].instruments
        error("Invalid instrument: $instrument_type for model $(model.name). Valid names are ", MODEL_CONFIGS[model.name].instruments)
    end

    # Prepare batches
    json_batches = prepare_koina_batch(
        model,
        peptides_df,
        instrument_type,
        batch_size=batch_size
    )
    # Request predictions
    responses = make_koina_batch_requests(json_batches, KOINA_URLS[model.name]; concurrency=concurrent_koina_requests)
    # Process responses
    batch_dfs = []
    for (i, response) in enumerate(responses)
        batch_result = parse_koina_batch(model, response)
        start_idx = (i-1) * batch_size + 1 + first_prec_idx - 1
        
        # Add precursor indices
        batch_df = batch_result.fragments
        n_precursors_in_batch = UInt32(fld(size( batch_df , 1), batch_result.frags_per_precursor))
        batch_df[!, :precursor_idx] = repeat(start_idx:(start_idx + n_precursors_in_batch - one(UInt32)), 
                                                inner=batch_result.frags_per_precursor)
        # Filter and sort fragments
        filter_fragments!(batch_df, model)
        push!(batch_dfs, batch_df)
    end

    # Combine and filter results
    fragments_df = vcat(batch_dfs...)
    
    sort_fragments!(fragments_df)

    return fragments_df
end

"""
Fragment prediction for instrument-agnostic models (e.g., Prosit).
"""
function predict_fragments_batch(
    peptides_df::DataFrame,
    model::InstrumentAgnosticModel,
    _::String,  # instrument type not used
    batch_size::Int,
    concurrent_koina_requests::Int,
    first_prec_idx::UInt32
)::DataFrame
    # Prepare batches (no instrument type needed)
    json_batches = prepare_koina_batch(
        model,
        peptides_df,
        batch_size=batch_size
    )

    # Request predictions
    responses = make_koina_batch_requests(json_batches, KOINA_URLS[model.name]; concurrency=concurrent_koina_requests)

    # Process responses
    batch_dfs = []
    for (i, response) in enumerate(responses)
        batch_result = parse_koina_batch(model, response)
        start_idx = (i-1) * batch_size + 1 + first_prec_idx - 1
        
        # Add precursor indices
        batch_df = batch_result.fragments
        n_precursors_in_batch = UInt32(fld(size( batch_df , 1), batch_result.frags_per_precursor))
        batch_df[!, :precursor_idx] = repeat(start_idx:(start_idx + n_precursors_in_batch - one(UInt32)), 
                                                inner=batch_result.frags_per_precursor)
        # Filter and sort fragments
        filter_fragments!(batch_df, model)
        push!(batch_dfs, batch_df)
    end

    fragments_df = vcat(batch_dfs...)
    sort_fragments!(fragments_df)

    return fragments_df
end

"""
Fragment prediction for spline coefficient models (e.g., Altimeter).
"""
function predict_fragments_batch(
    peptides_df::DataFrame,
    model::SplineCoefficientModel,
    instrument_type::String,
    batch_size::Int,
    concurrent_koina_requests::Int,
    first_prec_idx::UInt32
)::DataFrame
    # Similar to InstrumentSpecificModel but handles spline coefficients
    json_batches = prepare_koina_batch(
        model,
        peptides_df,
        instrument_type,
        batch_size=batch_size
    )

    responses = make_koina_batch_requests(json_batches, KOINA_URLS[model.name]; concurrency=concurrent_koina_requests)

    batch_dfs = []
    knot_vectors = []
    
    for (i, response) in enumerate(responses)
        batch_result = parse_koina_batch(model, response)
        start_idx = (i-1) * batch_size + 1 + first_prec_idx - 1
        
        batch_df = batch_result.fragments
        n_precursors_in_batch = UInt32(fld(size( batch_df , 1), batch_result.frags_per_precursor))
        batch_df[!, :precursor_idx] = repeat(start_idx:(start_idx + n_precursors_in_batch - one(UInt32)), 
                                                inner=batch_result.frags_per_precursor)
        filter_fragments!(batch_df, model)
        push!(batch_dfs, batch_df)
        push!(knot_vectors, batch_result.extra_data)  # Store knot vectors
    end

    fragments_df = vcat(batch_dfs...)
    
    # Verify knot vectors are consistent
    if !all(k == first(knot_vectors) for k in knot_vectors)
        error("Inconsistent knot vectors across batches")
    end
    
    # Store knot vector with the data
    fragments_df[!, :knot_vector] .= Ref(first(knot_vectors))
    #For altimeter fragments are already sorted 
    #sort_fragments!(f)

    return fragments_df
end

"""
Filter fragments based on intensity and other criteria.
"""
function filter_fragments!(df::DataFrame, model::KoinaModelType)
    # Basic filtering common to all models
    filter!(:intensities => x -> x > 0.001f0, df)  # Remove very low intensity
    filter!(:mz => x -> x > 0, df)  # Remove invalid m/z
    
    # Model-specific filtering
    if model isa InstrumentSpecificModel
        filter!(row -> !occursin('i', row.annotation), df)  # Remove isotope peaks
    end
end

"""
Filter fragments based on intensity and other criteria.
"""
function filter_fragments!(df::DataFrame, model::SplineCoefficientModel)
    # Basic filtering common to all models
    #filter!(:coefficients => x -> x > zero(Float32), df)  # Remove very low intensity
    filter!(:mz => x -> x > 0, df)  # Remove invalid m/z
    
    # Model-specific filtering
    if model isa InstrumentSpecificModel
        filter!(row -> !occursin('i', row.annotation), df)  # Remove isotope peaks
    end
end


"""
Sort fragments by intensity within each precursor group.
"""
function sort_fragments!(df::DataFrame)
    sort!(df, [:precursor_idx, order(:intensities, rev=true)])
end

function sort_fragments!(df::DataFrame)
    sort!(df, [:precursor_idx, order(:intensities, rev=true)])
end
