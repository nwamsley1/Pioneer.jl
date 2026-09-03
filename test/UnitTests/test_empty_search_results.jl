# Regression tests for searches that find nothing.
#
# A file that matches nothing is a normal outcome — the wrong library, a blank
# run, or a gradient with no usable MS2 scans — and it must not take the rest of
# the search with it. `library_search` used to return a bare `DataFrame()` on
# its two early exits, so such a file produced a frame with NO COLUMNS, and the
# first downstream reader of `:precursor_idx` (`_add_psm_features!`, via
# `prepare_psm_features!`) threw `ArgumentError: column name :precursor_idx not
# found in the data frame`. One dud file aborted a whole multi-file run.
#
# The fix returns a zero-ROW frame carrying the real schema, taken as a
# zero-length view of the same scored-PSM store the non-empty path uses. These
# tests pin the property that makes that work: the schema comes from the
# element type, so an empty store still yields every column.
#
# Run standalone: julia --project=. test/UnitTests/test_empty_search_results.jl
if !@isdefined(Pioneer)
    using Test
    using Pioneer
    using DataFrames
end

@testset "empty search results carry the full schema" begin
    # Both stores `get_scored_psms` can dispatch to.
    for T in (
        Pioneer.MainSearchScoredPSM{Float32, Float16},
        Pioneer.TuningScoredPSM{Float32, Float16},
    )
        store = T[]
        empty_frame = DataFrame(@view(store[1:0]))

        @test nrow(empty_frame) == 0
        # Not a bare DataFrame(): the columns are what the crash was about.
        @test ncol(empty_frame) > 0
        # Every field of the row type, so the schema cannot drift from the
        # non-empty path the way a hand-written column list would.
        @test Set(Symbol.(names(empty_frame))) == Set(fieldnames(T))
    end

    # The four columns `_add_psm_features!` reads off the frame before doing
    # anything else. If any of these ever leaves MainSearchScoredPSM, the
    # feature pass needs revisiting, not just this test.
    main_frame = DataFrame(@view((Pioneer.MainSearchScoredPSM{Float32, Float16}[])[1:0]))
    for col in (:precursor_idx, :scan_idx, :error, :total_ions)
        @test hasproperty(main_frame, col)
    end
end
