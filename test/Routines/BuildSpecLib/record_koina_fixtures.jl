# Capture live Koina responses for every BuildSpecLib invocation in
# test/Routines/BuildSpecLib/test_build_spec_lib.jl, save them to
# test/Routines/BuildSpecLib/koina_fixtures.json, then commit that
# fixture file to the repo.
#
# After re-running this script, `test_build_spec_lib.jl` will replay
# the captured responses via `Pioneer.FixtureKoinaClient` instead of
# making live HTTP calls — which is what makes the BuildSpecLib leg of
# `runtests.jl` runnable on GitHub Actions / any no-network env.
#
# Re-run requirements:
#   - Network access to https://koina.wilhelmlab.org
#   - Same Pioneer commit you intend to run tests on (to keep the
#     request payloads byte-identical)
#
# Usage:
#   julia --project=. test/Routines/BuildSpecLib/record_koina_fixtures.jl

using Pioneer
using Pioneer: RecordingKoinaClient, with_koina_client, save_koina_fixtures

# test_build_spec_lib.jl assumes a Test-aware namespace with these imports
# (it's normally included from runtests.jl which loads them). Mirror the
# subset it actually touches so this script can run standalone.
using Test
using Arrow, CSV, DataFrames, JSON, HTTP
using DataStructures: OrderedDict

const FIXTURE_PATH = joinpath(@__DIR__, "koina_fixtures.json")

# Run the entire BuildSpecLib testsets under the recorder. Every
# `Pioneer.BuildSpecLib(params_file)` call inside test_build_spec_lib.jl
# flows its Koina requests through the recorder. Test assertions still
# evaluate (ignored if they fail — we just want the captures); the
# point of this script is the capture, not the test outcome.
recorder = RecordingKoinaClient()
try
    with_koina_client(recorder) do
        # Suppress the test-failure exit so a single failing assertion
        # doesn't prevent us from saving the captures we did get.
        try
            include(joinpath(@__DIR__, "test_build_spec_lib.jl"))
        catch e
            @warn "test_build_spec_lib.jl raised during recording (continuing to save captures)" exception=(e, catch_backtrace())
        end
    end
finally
    save_koina_fixtures(recorder, FIXTURE_PATH)
    println("Wrote $(length(recorder.captures)) Koina captures to $FIXTURE_PATH")
end
