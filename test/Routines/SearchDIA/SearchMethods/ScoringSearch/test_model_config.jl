using Test
using Pioneer

@testset "ScoringSearch model configurations" begin
    configs = Pioneer.create_model_configurations()

    # Expect five models after addition of simplified probit
    @test length(configs) == 5

    names = getfield.(configs, :name)
    @test "ProbitRegression" in names
    @test "ProbitRegressionSimple" in names

    # Intercept should appear only in probit model feature sets
    model_by_name = Dict(name => cfg for (name, cfg) in zip(names, configs))
    @test :intercept in model_by_name["ProbitRegression"].features
    @test :intercept in model_by_name["ProbitRegressionSimple"].features
    @test :intercept ∉ model_by_name["SimpleXGBoost"].features
    @test :intercept ∉ model_by_name["AdvancedXGBoost"].features
    @test :intercept ∉ model_by_name["SuperSimplified"].features
end
