@testset "Model peptide-length limits" begin
    # Altimeter's ceiling is 40, probed against Koina: a 40-mer is accepted and
    # a 41-mer is rejected outright. The catalog builds at exactly 7-40, so this
    # request must pass through untouched.
    @test Pioneer.clamp_digest_length_to_model("altimeter", 7, 40) == (7, 40)

    # ...and one residue past it must be narrowed here rather than surviving to
    # Koina, where a single over-long peptide fails its entire batch and the
    # retry loop then misreports it as a network fault.
    @test Pioneer.clamp_digest_length_to_model("altimeter", 7, 41) == (7, 40)
    @test Pioneer.clamp_digest_length_to_model("altimeter", 7, 50) == (7, 40)

    # No lower bound: a 1-mer is accepted by the model, so a short digest is the
    # user's business.
    @test Pioneer.clamp_digest_length_to_model("altimeter", 5, 30) == (5, 30)

    # An unrecognised name is a pass-through, not an error: prediction_model is
    # validated where it is read, and this function must not become a second
    # place that rejects it.
    @test Pioneer.clamp_digest_length_to_model("not_a_model", 5, 99) == (5, 99)

    # Install a constrained model for the duration of the test. MODEL_CONFIGS is
    # a const binding to a mutable Dict, so this is a temporary insertion.
    key = "test_only_bounded_model"
    Pioneer.MODEL_CONFIGS[key] = (
        annotation_type = Pioneer.UniSpecFragAnnotation("y1^1"),
        model_type = Pioneer.SplineCoefficientModel(key),
        instruments = Set([]),
        fragmentation_type = nothing,
        peptide_length = (min = 7, max = 30),
        supported_mods = Pioneer.BASE_MODS,
        free_cys = false,
    )
    try
        # Fully inside the model's range: untouched.
        @test Pioneer.clamp_digest_length_to_model(key, 8, 25) == (8, 25)

        # Exactly on the boundaries: untouched, and no off-by-one narrowing.
        @test Pioneer.clamp_digest_length_to_model(key, 7, 30) == (7, 30)

        # The case this whole mechanism exists for: a user max above the model's
        # cap. 31-40mers would otherwise reach Koina and be dropped silently.
        @test Pioneer.clamp_digest_length_to_model(key, 7, 40) == (7, 30)

        # Below the model's minimum, raised.
        @test Pioneer.clamp_digest_length_to_model(key, 4, 30) == (7, 30)

        # Both ends out of range.
        @test Pioneer.clamp_digest_length_to_model(key, 4, 40) == (7, 30)

        # Disjoint ranges would yield an empty library, so this must throw rather
        # than return an inverted or empty window.
        @test_throws ErrorException Pioneer.clamp_digest_length_to_model(key, 35, 40)
        @test_throws ErrorException Pioneer.clamp_digest_length_to_model(key, 2, 5)
    finally
        delete!(Pioneer.MODEL_CONFIGS, key)
    end

    @test !haskey(Pioneer.MODEL_CONFIGS, key)
end

@testset "Retention-time model peptide-length limits" begin
    clamp(frag, lo, hi, rt) = Pioneer.clamp_digest_length_to_model(frag, lo, hi; rt_model = rt)

    # Prosit 2024 iRT rejects 31-mers (probed against Koina). Paired with Altimeter,
    # whose own ceiling is 40, the RT model is the binding limit: without it a
    # 7-40 digest reached Koina and failed its batches (issue #267).
    @test clamp("altimeter", 7, 40, "prosit_2024_irt_ptm") == (7, 30)
    @test clamp("altimeter", 7, 30, "prosit_2024_irt_ptm") == (7, 30)

    # Chronologer declares no limit, so Altimeter's 40 still applies.
    @test clamp("altimeter", 7, 40, "chronologer") == (7, 40)
    @test clamp("altimeter", 7, 41, "chronologer") == (7, 40)

    # The iRT model accepts short peptides, so it never raises the minimum.
    @test clamp("altimeter", 1, 20, "prosit_2024_irt_ptm") == (1, 20)

    # Both models cap at 30: same answer as the fragment model alone.
    @test clamp("prosit_2020_hcd", 7, 40, "prosit_2024_irt_ptm") == (7, 30)

    # No RT model given, or an unknown one: the fragment model alone decides.
    @test clamp("altimeter", 7, 40, nothing) == (7, 40)
    @test clamp("altimeter", 7, 40, "not_a_model") == (7, 40)

    # Nothing in the requested window that the RT model accepts.
    @test_throws ErrorException clamp("altimeter", 35, 40, "prosit_2024_irt_ptm")

    # Every RT model declares the field, so `get` never silently skips a typo.
    @test all(haskey(cfg, :peptide_length) for cfg in values(Pioneer.RT_MODEL_CONFIGS))
end
