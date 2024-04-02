integratePrecursorMS2(
    MS2_CHROMS[50002]
)

test_chrom = MS2_CHROMS[N][!,[:score,:total_ions,:entropy_score,:weight,:RT]]
plot(test_chrom[!,:RT],
     test_chrom[!,:weight],
     seriestype=:scatter)
N += 1
