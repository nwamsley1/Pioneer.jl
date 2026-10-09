using Test, Arrow, DataFrames

# Three runs, one file each. Run 1 repeats precursor 10 (two passing rows), run 3 has
# no passing rows. Run-level and global q-values are independent columns.
function _mbr_passed_set_fixture(dir)
    rows = [
        # file, pid, qval,   global_qval, weight
        (1, 10, 0.005f0, 0.001f0, 100.0f0),
        (1, 10, 0.001f0, 0.001f0, 120.0f0),
        (1, 11, 0.5f0,   0.002f0, 50.0f0),
        (1, 12, 0.005f0, 0.5f0,   80.0f0),
        (2, 10, 0.5f0,   0.001f0, 90.0f0),
        (2, 11, 0.005f0, 0.002f0, 60.0f0),
        (2, 12, 0.008f0, 0.5f0,   70.0f0),
        (3, 13, 0.5f0,   NaN32,   40.0f0),
        (3, 10, NaN32,   0.001f0, 30.0f0),
    ]
    paths = String[]
    for f in 1:3
        sel = [r for r in rows if r[1] == f]
        path = joinpath(dir, "run$f.arrow")
        Arrow.write(path, DataFrame(
            ms_file_idx = UInt32[r[1] for r in sel],
            precursor_idx = UInt32[r[2] for r in sel],
            qval = Float32[r[3] for r in sel],
            global_qval = Float32[r[4] for r in sel],
            weight = Float32[r[5] for r in sel],
        ))
        push!(paths, path)
    end
    return paths
end

# (receiver, pid) => eligible: passes global q AND does not pass run q in the receiver.
const _MBR_EXPECTED_ELIGIBLE = Dict(
    (1, 10) => false, (1, 11) => true,  (1, 12) => false, (1, 13) => false,
    (2, 10) => true,  (2, 11) => false, (2, 12) => false, (2, 13) => false,
    (3, 10) => true,  (3, 11) => true,  (3, 12) => false, (3, 13) => false,
)
# pid => runs where it passes run q (each counted once)
const _MBR_EXPECTED_PASSING_RUNS = Dict(10 => [1], 11 => [2], 12 => [1, 2])

# The run-level sets are built per receiver file, as compute_postintegration_mbr_features! does.
function _receiver_eligibility(e, path, q)
    tbl = Arrow.Table(path)
    return Pioneer._MBRCounterfactualEligibility(
        e.global_passed,
        Pioneer._mbr_run_passed_by_file(tbl.precursor_idx, tbl.ms_file_idx, tbl.qval, q),
    )
end

@testset "MBR run-level passing sets" begin
    mktempdir() do dir
        paths = _mbr_passed_set_fixture(dir)
        q = 0.01f0
        eligibility = Pioneer.build_mbr_counterfactual_eligibility(paths; q_value_threshold = q)
        for ((r, p), expected) in _MBR_EXPECTED_ELIGIBLE
            e = _receiver_eligibility(eligibility, paths[r], q)
            @test Pioneer._mbr_counterfactual_eligible(e, UInt32(r), UInt32(p)) == expected
        end

        clusters = Pioneer.build_mbr_receiver_run_clusters(paths; q_value_threshold = q)
        # every run with rows is clustered, including run 3 where nothing passes
        @test sort!(collect(keys(clusters.cluster_by_file))) == UInt32[1, 2, 3]
        for (pid, runs) in _MBR_EXPECTED_PASSING_RUNS
            support = get(clusters.support_by_precursor, UInt32(pid), Dict{UInt32, UInt32}())
            # duplicate passing rows in one run count once
            @test sum(values(support); init = UInt32(0)) == length(runs)
        end
    end
end
