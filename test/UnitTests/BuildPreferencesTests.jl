using UUIDs

@testset "Build preferences" begin
    host_cpu_features_uuid = UUID("3e5b6fbb-0976-4d2c-9146-d79de83f2fb0")
    prefs = Base.get_preferences(host_cpu_features_uuid)

    @test get(prefs, "allow_runtime_invalidation", false) === true
end
