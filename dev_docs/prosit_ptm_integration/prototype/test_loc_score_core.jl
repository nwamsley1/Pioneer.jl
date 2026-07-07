using Test
include("loc_score_core.jl")

# Observed spectrum from the *true* isomer: every one of its b/y fragments present at intensity 1.
obs_from(true_mods, L; ion_types = (:b, :y)) = Dict{Tuple{Symbol,Int,Int},Float64}(
    (T, i, frag_modcount(true_mods, T, i, L)) => 1.0 for T in ion_types for i in 1:(L - 1))

@testset "frag_modcount" begin
    L = 7
    @test frag_modcount([1], :b, 1, L) == 1      # b1 covers res 1, mod at 1
    @test frag_modcount([1], :b, 6, L) == 1
    @test frag_modcount([3], :b, 2, L) == 0      # b2 covers 1..2, mod at 3 outside
    @test frag_modcount([3], :b, 3, L) == 1
    @test frag_modcount([1], :y, 6, L) == 0      # y6 covers 2..7, mod at 1 outside
    @test frag_modcount([7], :y, 1, L) == 1      # y1 covers 7..7, mod at 7 inside
    @test frag_modcount([1,4], :b, 5, L) == 2    # b5 covers 1..5, mods 1 and 4 inside
end

@testset "PLAN worked example: S A S A S A K (present = A = pS1)" begin
    # acceptors S1,S3,S5 -> targets A=[1], B=[3], C=[5]; shadow decoy A'=[2]
    L = 7
    A, B, C, Ap = [1], [3], [5], [2]
    obs = obs_from(A, L)                          # scan contains A (pS1)

    # target A vs the OTHER targets {B,C} -> distinguished -> 1.0
    @test S_vs_targets_target(A, [A, B, C], L, obs).S == 1.0

    # decoy A' vs ALL targets incl parent {A,B,C} -> killed by parent (b1) -> 0.0
    @test S_vs_targets_decoy(Ap, [A, B, C], L, obs).S == 0.0

    # contrast: dropping the parent (A' vs {B,C}) makes A' ESCAPE to 1.0 (the rejected design)
    @test S_score(Ap, [B, C], L, obs).S == 1.0

    # a mislocalized target lands where the (kept-parent) decoy does: report B when truth is A
    @test S_vs_targets_target(B, [A, B, C], L, obs).S == 0.0
end

@testset "competitors differ at >1 site (general cumulative-content)" begin
    # SAASAASAAR-like, L=10, acceptors 1,4,7. Di-phospho isomers.
    L = 10
    A = [1, 4]; B = [4, 7]            # share pos4, differ at 1 and 7 (symdiff {1,7})
    obs = obs_from(A, L)             # truth = A

    r = s_pair(A, B, L, obs)
    @test r.ndist > 0                # distinguishing ions exist across b (S1 end) and y (S7 end)
    @test r.s == 1.0                 # observed = A -> all distinguishing ions support A

    # B (as the assigned placement) is rejected when truth is A
    @test s_pair(B, A, L, obs).s == 0.0

    # per-site: site 1 (A's differing site) is independently confirmed by b-ions (not joint-only)
    ps = per_site_support(A, B, L, obs)
    @test haskey(ps.site_supA, 1)
    @test ps.site_supA[1] > 0
    @test ps.joint_only[1] == false
end

@testset "per-site joint-only flag" begin
    # A=[1,4] vs B=[7,10] (L=12): symdiff {1,4,7,10}; A's differing sites {1,4}.
    # No single b/y ion isolates 1 from 4 (they move together on the N-term side under +2),
    # so per-site for site 1 or 4 should be flagged joint-only when the only distinguishing
    # ions span both.
    L = 12
    A = [1, 4]; B = [7, 10]
    obs = obs_from(A, L)
    ps = per_site_support(A, B, L, obs)
    # b-ions in [4,7) carry a +2 shift (both 1 and 4 inside) -> span two differing sites -> joint;
    # b-ions in [1,4) isolate site 1 (only 1 inside) -> site 1 independently confirmed.
    @test ps.joint_only[1] == false      # b1..b3 isolate site 1
    # site 4: ions isolating it (containing 4 but not 1) don't exist on the b-side (1<4 always co-included);
    # on the y-side, y-ions containing 4 also contain 7,10 region... it should remain joint-only.
    @test haskey(ps.joint_only, 4)
end

println("loc_score_core tests done.")
