using Test
using Pioneer

@testset "scan-local pair competition gates by selected scan" begin
    prec_ids = UInt32[1, 2, 3, 4]
    targets = Bool[true, false, true, false]
    partner_col = Union{Missing, UInt32}[UInt32(2), UInt32(1), UInt32(4), UInt32(3)]
    row_by_pid = Int32[1, 2, 3, 4]
    score_by_pid = Float32[0.9, 0.8, 0.1, 0.95]
    scan_by_pid = UInt32[10, 11, 20, 30]

    keep_margin_1 = fill(true, length(prec_ids))
    Pioneer._mark_paircomp_losers!(
        keep_margin_1,
        prec_ids,
        targets,
        partner_col,
        row_by_pid,
        score_by_pid,
        scan_by_pid,
        1,
    )
    @test keep_margin_1 == Bool[true, false, true, true]

    keep_margin_10 = fill(true, length(prec_ids))
    Pioneer._mark_paircomp_losers!(
        keep_margin_10,
        prec_ids,
        targets,
        partner_col,
        row_by_pid,
        score_by_pid,
        scan_by_pid,
        10,
    )
    @test keep_margin_10 == Bool[true, false, false, true]
end
