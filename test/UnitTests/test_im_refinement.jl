using Test, Random, Statistics, DataFrames, Pioneer

@testset "Sequence-held-out ion mobility refinement" begin
    rng = MersenneTwister(42)
    base = [String(rand(rng, collect("ACDEFGHIKLMNPQRSTVWY"), rand(rng, 8:24))) for _ in 1:500]
    sequences = repeat(base; inner = 3)
    charge = repeat([2, 3, 4], 500)
    mods = fill("", length(sequences))
    tokens = Pioneer.im_composition_tokens(sequences, mods)
    pred = 0.7 .+ 0.02length.(sequences) .+ 0.06charge .+ 0.05rand(rng, length(sequences))
    obs = pred .+ 0.02charge .+ 0.025 .* tokens[:, 13] .* (charge .- 2) .+
          0.001randn(rng, length(pred))
    results = Dict(mode => Pioneer.crossfit_im_correction(pred, obs, charge, tokens, sequences;
                   mode = mode, min_charge = 20) for mode in Pioneer.IM_REFINEMENT_MODES)
    err(mode) = median(abs.(results[mode].refined .- obs))
    @test results[:none].refined == pred
    @test err(:charge) < err(:none)
    @test err(:composition) < err(:charge)
    @test err(:charge_composition) < 0.25err(:composition)
    @test err(:auto) < 0.25err(:composition)
    @test all(results[:auto].folds[1:3:end] .== results[:auto].folds[2:3:end])
    @test all(results[:auto].folds[1:3:end] .== results[:auto].folds[3:3:end])
    @test length(unique(results[:auto].folds)) == 5

    # Changing a held-out fold's observations cannot affect its predictions,
    # including automatic model selection on the inner training holdout.
    held = results[:auto].folds .== 1
    changed = copy(obs); changed[held] .+= 10
    again = Pioneer.crossfit_im_correction(pred, changed, charge, tokens, sequences;
                                          mode = :auto, min_charge = 20)
    @test again.refined[held] == results[:auto].refined[held]
    @test again.selected[1] == results[:auto].selected[1]

    identity = Pioneer.crossfit_im_correction(pred, pred, charge, tokens, sequences;
                                              mode = :auto, min_charge = 20)
    @test all(==(:none), values(identity.selected))
    @test identity.refined == pred
    short = Pioneer.crossfit_im_correction(pred, obs, charge, tokens, sequences;
                 mode = :charge_composition, min_charge = 10_000)
    @test short.refined == pred
    no_anchors = Pioneer.crossfit_im_correction(pred, obs, charge, tokens, sequences;
                  mode = :auto, training_mask = falses(length(pred)))
    @test no_anchors.refined == pred
    model = Pioneer.fit_im_correction(pred, obs, charge, tokens; min_charge = 20)
    @test Pioneer.predict_im_correction(model, [1.0], [7], tokens[1:1, :]) == [0.0]
    @test_throws ArgumentError Pioneer.crossfit_im_correction(pred, obs, charge, tokens, sequences; mode = :invalid)
    bad = copy(pred); bad[1] = NaN
    @test_throws ArgumentError Pioneer.fit_im_correction(bad, obs, charge, tokens)

    masked = Pioneer.crossfit_im_correction(pred, obs, charge, tokens, sequences;
                 mode = :composition, training_mask = charge .== 2, min_charge = 20)
    @test masked.refined[charge .!= 2] == pred[charge .!= 2]
end

struct _ImRefinementPrecursors
    sequence::Vector{String}
    structural_mods::Vector{String}
    inv_ion_mobility::Vector{Float32}
end
Pioneer.getSequence(p::_ImRefinementPrecursors) = p.sequence
Pioneer.getStructuralMods(p::_ImRefinementPrecursors) = p.structural_mods
Pioneer.getInvIonMobility(p::_ImRefinementPrecursors) = p.inv_ion_mobility

@testset "IM correction preserves observed mobility and excludes decoys" begin
    rng = MersenneTwister(9)
    sequences = [String(rand(rng, collect("ACDEFGHIKLMNPQRSTVWY"), 12)) for _ in 1:1200]
    pred = Float32.(0.8 .+ 0.4rand(rng, 1200))
    precursors = _ImRefinementPrecursors(sequences, fill("", 1200), pred)
    psms = DataFrame(precursor_idx = UInt32.(1:1200), target = trues(1200),
        charge = fill(UInt8(2), 1200), im_obs = pred .+ 0.04f0,
        im_error = fill(-4f0, 1200), lgbm_prob = fill(0.99f0, 1200))
    psms.target[1001:end] .= false
    psms.lgbm_prob[1001:end] .= 0.1f0
    psms.im_obs[1001:end] .+= 0.2f0
    original_obs = copy(psms.im_obs)
    Pioneer.refine_im_error!(psms, precursors, 0.01; mode = :charge)
    @test psms.im_obs == original_obs
    @test psms.im_error_uncorrected == fill(-4f0, 1200)
    @test maximum(abs, psms.im_error[1:1000]) < 0.001
    @test all(<(-19f0), psms.im_error[1001:end])
    @test psms.im_pred == pred
    @test eltype(psms.im_pred_refined) == Float32
    unchanged = deepcopy(psms)
    Pioneer.refine_im_error!(psms, precursors, 0.01; mode = :none)
    @test psms == unchanged
end

@testset "IM refinement configuration" begin
    mktempdir() do dir
        path = joinpath(dir, "params.json")
        config = Dict("paths" => Dict("ms_data" => dir, "library" => dir, "results" => dir))
        write(path, Pioneer.JSON.json(config))
        params = Pioneer.parse_pioneer_parameters(path)
        @test Pioneer.MainSearchParameters(params).im_refinement === :auto
        for mode in Pioneer.IM_REFINEMENT_MODES
            config["global"] = Dict("im_refinement" => String(mode))
            write(path, Pioneer.JSON.json(config))
            params = Pioneer.parse_pioneer_parameters(path)
            @test Pioneer.MainSearchParameters(params).im_refinement === mode
        end
        config["global"] = Dict("im_refinement" => "invalid")
        write(path, Pioneer.JSON.json(config))
        params = Pioneer.parse_pioneer_parameters(path)
        @test_throws ArgumentError Pioneer.MainSearchParameters(params)
    end
end
