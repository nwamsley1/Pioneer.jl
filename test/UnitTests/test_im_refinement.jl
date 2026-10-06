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
    result = Pioneer.crossfit_im_correction(pred, obs, charge, tokens, sequences;
                                             min_charge = 20)
    @test median(abs.(result.refined .- obs)) < 0.002
    @test median(abs.(result.refined .- obs)) < 0.05median(abs.(pred .- obs))
    @test all(result.folds[1:3:end] .== result.folds[2:3:end])
    @test all(result.folds[1:3:end] .== result.folds[3:3:end])
    @test length(unique(result.folds)) == 5

    # Held-out observations cannot influence that fold's predictions.
    held = result.folds .== 1
    changed = copy(obs); changed[held] .+= 10
    again = Pioneer.crossfit_im_correction(pred, changed, charge, tokens, sequences;
                                          min_charge = 20)
    @test again.refined[held] == result.refined[held]

    short = Pioneer.crossfit_im_correction(pred, obs, charge, tokens, sequences;
                                          min_charge = 10_000)
    @test short.refined == pred
    no_anchors = Pioneer.crossfit_im_correction(pred, obs, charge, tokens, sequences;
                                               training_mask = falses(length(pred)))
    @test no_anchors.refined == pred
    model = Pioneer.fit_im_correction(pred, obs, charge, tokens; min_charge = 20)
    @test Pioneer.predict_im_correction(model, [1.0], [7], tokens[1:1, :]) == [0.0]
    bad = copy(pred); bad[1] = NaN
    @test_throws ArgumentError Pioneer.fit_im_correction(bad, obs, charge, tokens)

    masked = Pioneer.crossfit_im_correction(pred, obs, charge, tokens, sequences;
                 training_mask = charge .== 2, min_charge = 20)
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
    Pioneer.refine_im_error!(psms, precursors, 0.01)
    @test psms.im_obs == original_obs
    @test psms.im_error_uncorrected == fill(-4f0, 1200)
    @test maximum(abs, psms.im_error[1:1000]) < 0.001
    @test all(<(-19f0), psms.im_error[1001:end])
    @test psms.im_pred == pred
    @test eltype(psms.im_pred_refined) == Float32
end

@testset "IM refinement skips empty tables" begin
    @test Pioneer.refine_im_error!(DataFrame(), _ImRefinementPrecursors(String[], String[], Float32[]), 0.01) === nothing
end
