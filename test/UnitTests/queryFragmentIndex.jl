##########
#findFirstFragmentBin
##########
#Test cases for fragment bin search
@testset "find_first_fragment_bin" begin
    bin_bounds = Float32.(collect(LinRange(0, 10, 11)))
    test_frag_bins = Vector{FragIndexBin}(undef, length(bin_bounds) - 1)
    for i in range(1, length(bin_bounds) - 1)
        test_frag_bins[i] = FragIndexBin(
            bin_bounds[i],
            bin_bounds[i + 1],
            UInt32(i),
            UInt32(i)
        )
    end
    #=
    julia> test_frag_bins
    10-element Vector{FragIndexBin}:
    FragIndexBin{Float32}(0.0f0, 1.0f0, 0x00000001, 0x00000001)
    FragIndexBin{Float32}(1.0f0, 2.0f0, 0x00000002, 0x00000002)
    FragIndexBin{Float32}(2.0f0, 3.0f0, 0x00000003, 0x00000003)
    FragIndexBin{Float32}(3.0f0, 4.0f0, 0x00000004, 0x00000004)
    FragIndexBin{Float32}(4.0f0, 5.0f0, 0x00000005, 0x00000005)
    FragIndexBin{Float32}(5.0f0, 6.0f0, 0x00000006, 0x00000006)
    FragIndexBin{Float32}(6.0f0, 7.0f0, 0x00000007, 0x00000007)
    FragIndexBin{Float32}(7.0f0, 8.0f0, 0x00000008, 0x00000008)
    FragIndexBin{Float32}(8.0f0, 9.0f0, 0x00000009, 0x00000009)
    FragIndexBin{Float32}(9.0f0, 10.0f0, 0x0000000a, 0x0000000a)
    =#

    @test UInt32(2) === findFirstFragmentBin(
        test_frag_bins,
        UInt32(1), #lowest frag bin to search
        UInt32(10), #Highst frag bin to search
        1.1f0, #query
    )

    #Query is greater than any bin, so return max allowed 
    #bin index 
    @test UInt32(10) === findFirstFragmentBin(
        test_frag_bins,
        UInt32(1), #lowest frag bin to search
        UInt32(10), #Highst frag bin to search
        100.1f0, #query
    )
    @test UInt32(9) === findFirstFragmentBin(
        test_frag_bins,
        UInt32(1), #lowest frag bin to search
        UInt32(9), #Highst frag bin to search
        100.1f0, #query
    )
    @test UInt32(5) === findFirstFragmentBin(
        test_frag_bins,
        UInt32(5), #lowest frag bin to search
        UInt32(5), #Highst frag bin to search
        100.1f0, #query
    )

    #Query is less than any bin, so return min allowed 
    #bin index 
    @test UInt32(1) === findFirstFragmentBin(
        test_frag_bins,
        UInt32(1), #lowest frag bin to search
        UInt32(10), #Highst frag bin to search
        0.0f0, #query
    )
    @test UInt32(2) === findFirstFragmentBin(
        test_frag_bins,
        UInt32(2), #lowest frag bin to search
        UInt32(10), #Highst frag bin to search
        0.0f0, #query
    )
    @test UInt32(5) === findFirstFragmentBin(
        test_frag_bins,
        UInt32(5), #lowest frag bin to search
        UInt32(5), #Highst frag bin to search
        0.0f0, #query
    )

    #What happens when the query falls in a gap
    #between fragment bins? Still gives first bin
    #where lower bound exceeds the query. 
    test_frag_bins = Vector{FragIndexBin}([
    FragIndexBin{Float32}(0.0f0, 1.0f0, 0x00000001, 0x00000001)
    FragIndexBin{Float32}(1.0f0, 2.0f0, 0x00000002, 0x00000002)
    #FragIndexBin{Float32}(2.0f0, 3.0f0, 0x00000003, 0x00000003) #missing bin
    FragIndexBin{Float32}(3.0f0, 4.0f0, 0x00000004, 0x00000004)
    FragIndexBin{Float32}(4.0f0, 5.0f0, 0x00000005, 0x00000005)
    FragIndexBin{Float32}(5.0f0, 6.0f0, 0x00000006, 0x00000006)])

    findFirstFragmentBin(
        test_frag_bins,
        UInt32(2), #lowest frag bin to search
        UInt32(5), #Highst frag bin to search
        2.5f0, #query
    )
end
##########
#exponentailFragmentBinSearch
##########
@testset "exponential_fragment_bin_search" begin
    N = 100000
    bin_bounds = Float32.(collect(LinRange(0, N, N + 1)))
    test_frag_bins = Vector{FragIndexBin}(undef, length(bin_bounds) - 1)
    for i in range(1, length(bin_bounds) - 1)
        test_frag_bins[i] = FragIndexBin(
            bin_bounds[i],
            bin_bounds[i + 1],
            UInt32(i),
            UInt32(i)
        )
    end
    #=
    FragIndexBin{Float32}(0.0f0, 1.0f0, 0x00000001, 0x00000001)
    FragIndexBin{Float32}(1.0f0, 2.0f0, 0x00000002, 0x00000002)
    FragIndexBin{Float32}(2.0f0, 3.0f0, 0x00000003, 0x00000003)
    FragIndexBin{Float32}(99997.0f0, 99998.0f0, 0x0001869e, 0x0001869e)
    FragIndexBin{Float32}(99998.0f0, 99999.0f0, 0x0001869f, 0x0001869f)
    FragIndexBin{Float32}(99999.0f0, 100000.0f0, 0x000186a0, 0x000186a0)
    =#

    #Restrictive absolute minimum
    frag_bin_max_idx = UInt32(N÷2)
    frag_mz_max = 30000.0f0
    frag_mz_absolute_min = 10.0f0
    lower_idx, upper_idx = exponentialFragmentBinSearch(
        test_frag_bins, 
        frag_bin_max_idx, #frag_bin_max_idx
        UInt32(1), #lower_bound_guess
        UInt32(1), #upper_bound_guess
        frag_mz_max ,#frag_mz_max
        frag_mz_absolute_min,#frag_mz_abslute_min
        UInt32(2048)#step_si´
    )
    @test lower_idx == UInt32(1)
    @test getHigh(test_frag_bins[lower_idx]) < frag_mz_absolute_min

    frag_bin_max_idx = UInt32(N÷2)
    frag_mz_max = 30000.0f0
    frag_mz_absolute_min = 29000.0f0
    lower_idx, upper_idx = exponentialFragmentBinSearch(
        test_frag_bins, 
        frag_bin_max_idx, #frag_bin_max_idx
        UInt32(1), #lower_bound_guess
        UInt32(1), #upper_bound_guess
        frag_mz_absolute_min,#frag_mz_abslute_min
        frag_mz_max ,#frag_mz_max
        UInt32(2048)#step_si´
    )
    @test lower_idx <= UInt32((upper_idx >> 1))
    @test getHigh(test_frag_bins[lower_idx]) < frag_mz_absolute_min

    frag_bin_max_idx = UInt32(N÷2)
    frag_mz_max = Float32(Inf)
    frag_mz_absolute_min = 29000.0f0
    lower_idx, upper_idx = exponentialFragmentBinSearch(
        test_frag_bins, 
        frag_bin_max_idx, #frag_bin_max_idx
        UInt32(1), #lower_bound_guess
        UInt32(1), #upper_bound_guess
        frag_mz_absolute_min,#frag_mz_abslute_min
        frag_mz_max ,#frag_mz_max
        UInt32(2048)#step_si´
    )
    @test upper_idx == frag_bin_max_idx
    @test lower_idx != upper_idx
    @test getHigh(test_frag_bins[lower_idx]) < frag_mz_absolute_min


    frag_bin_max_idx = UInt32(N÷2)
    frag_mz_max = Float32(Inf)
    frag_mz_absolute_min = Float32(Inf)
    lower_idx, upper_idx = exponentialFragmentBinSearch(
        test_frag_bins, 
        frag_bin_max_idx, #frag_bin_max_idx
        UInt32(1), #lower_bound_guess
        UInt32(1), #upper_bound_guess
        frag_mz_absolute_min,#frag_mz_abslute_min
        frag_mz_max ,#frag_mz_max
        UInt32(2048)#step_si´
    )
    @test upper_idx == frag_bin_max_idx
    #@test lower_idx == upper_idx
    @test getHigh(test_frag_bins[lower_idx]) < frag_mz_absolute_min
end
##########
#searchFragmentBin!
##########
@testset "search_fragment_bin" begin
    n_precursors = 10
    prec_id_to_score = Counter(UInt32, UInt8, n_precursors)
    precursor_mzs = Float32.(collect(LinRange(1, 10, 10)))
    index_fragments = Vector{IndexFragment}([IndexFragment(
        UInt32(i),
        prec_mz,
        one(UInt8),
        one(UInt8)
    ) for (i, prec_mz) in enumerate(precursor_mzs)])

    frag_id_range =  UInt32(1):UInt32(length(index_fragments))
    prec_mz_min, prec_mz_max = 5.9f0, 7.1f0
    searchFragmentBin!(
        prec_id_to_score,
        index_fragments,
        frag_id_range,
        prec_mz_min,
        prec_mz_max
    )
    @test prec_id_to_score.ids[1] == UInt32(6)
    @test prec_id_to_score.ids[2] == UInt32(7)
    @test prec_id_to_score.ids[3] == UInt32(0)
    @test prec_id_to_score.size == 3
    @test sum(prec_id_to_score.counts) == 2
    @test prec_id_to_score.counts[6] == 1
    @test prec_id_to_score.counts[7] == 1

    reset!(prec_id_to_score)
    frag_id_range =  UInt32(8):UInt32(length(index_fragments))
    prec_mz_min, prec_mz_max = 5.9f0, 7.1f0
    searchFragmentBin!(
        prec_id_to_score,
        index_fragments,
        frag_id_range,
        prec_mz_min,
        prec_mz_max
    )

    @test sum(prec_id_to_score.ids)==0
    @test sum(prec_id_to_score.counts)==0

    reset!(prec_id_to_score)
    frag_id_range =  UInt32(1):UInt32(5)
    prec_mz_min, prec_mz_max = 5.9f0, 7.1f0
    searchFragmentBin!(
        prec_id_to_score,
        index_fragments,
        frag_id_range,
        prec_mz_min,
        prec_mz_max
    )

    @test sum(prec_id_to_score.ids)==0
    @test sum(prec_id_to_score.counts)==0

    reset!(prec_id_to_score)
    frag_id_range =  UInt32(6):UInt32(7)
    prec_mz_min, prec_mz_max = 5.9f0, 7.1f0
    searchFragmentBin!(
        prec_id_to_score,
        index_fragments,
        frag_id_range,
        prec_mz_min,
        prec_mz_max
    )

    @test prec_id_to_score.ids[1] == UInt32(6)
    @test prec_id_to_score.ids[2] == UInt32(7)
    @test prec_id_to_score.ids[3] == UInt32(0)
    @test prec_id_to_score.size == 3
    @test sum(prec_id_to_score.counts) == 2
    @test prec_id_to_score.counts[6] == 1
    @test prec_id_to_score.counts[7] == 1

    reset!(prec_id_to_score)
    frag_id_range =  UInt32(6):UInt32(6)
    prec_mz_min, prec_mz_max = 5.9f0, 7.1f0
    searchFragmentBin!(
        prec_id_to_score,
        index_fragments,
        frag_id_range,
        prec_mz_min,
        prec_mz_max
    )

    @test prec_id_to_score.ids[1] == UInt32(6)
    @test prec_id_to_score.ids[2] == UInt32(0)
    @test prec_id_to_score.ids[3] == UInt32(0)
    @test prec_id_to_score.size == 2
    @test sum(prec_id_to_score.counts) == 1
    @test prec_id_to_score.counts[6] == 1
    @test prec_id_to_score.counts[7] == 0
end
##########
#queryFragment!/searchScan!
##########
rt_bins = Vector{FragIndexBin}([
    FragIndexBin(1.0f0, 2.0f0, UInt32(1), UInt32(5)), 
    FragIndexBin(2.0f0, 3.0f0,  UInt32(6), UInt32(10)), 
    ])

fragment_bins = Vector{FragIndexBin}([
#First RT Bin 
FragIndexBin(100.0f0, 100.1f0, UInt32(1), UInt32(3)), 
FragIndexBin(100.1f0, 100.2f0,  UInt32(4), UInt32(4)), 
FragIndexBin(300.1f0, 301.0f0,  UInt32(5), UInt32(5)), 
FragIndexBin(301.0f0, 302.0f0,  UInt32(6), UInt32(6)), 
FragIndexBin(500.0f0, 500.1f0,  UInt32(7), UInt32(8)), 
#Second RT Bin 
FragIndexBin(100.0f0, 100.1f0, UInt32(9), UInt32(10)), 
FragIndexBin(100.1f0, 100.2f0,  UInt32(11), UInt32(11)), 
FragIndexBin(300.1f0, 301.0f0,  UInt32(12), UInt32(13)), 
FragIndexBin(301.0f0, 302.0f0,  UInt32(14), UInt32(14)), 
FragIndexBin(500.0f0, 500.1f0,  UInt32(15), UInt32(20)), 
])

#For simplicity, all fragments are assigned to the same precursor with an identical
#precursor match. So we will verify that the number of fragments matched is correct in each 
#test case/varaint. 
fragments = Vector{IndexFragment}([
#1 Frag bin
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #1
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #2
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #3
#2 Frag bin
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #4
#3 Frag bin
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #5
#4 Frag bin
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #6
#5 Frag bin
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #7
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #8
#Second RT bin
#6 Frag bin
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #9
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #10
#7 Frag bin
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #11
#8 Frag bin
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #12
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #13
#9 Frag bin
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #14
#10 Frag bin
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #15
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #16
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #17
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #18
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #19
IndexFragment(UInt32(1), 300.0f0, one(UInt8), one(UInt8)), #20
])

@testset "query_fragment" begin
    #=
    Need to build a fragment index for toy example.
    Should include 2 RT bins,
    10 fragment bins for each RT bin,
    and 0-10 precursors. 

    Test with a spectrum of 5-10 fragment ions

    Need to test cases that challenge variable, intensity dependent fragment tolerances. 
    =#

    n_precursors = 1
    prec_id_to_score = Counter(UInt32, UInt8, n_precursors)
    #=
    queryFragment!(prec_id_to_score,
                rt_bins[1].last_bin,
                rt_bins[1].first_bin,
                rt_bins[2].last_bin,
                fragment_bins,
                fragments,
                100.05f0,
                100.15f0,
                299.0f0,
                301.0f0)

    @test prec_id_to_score.counts[1] == 4
    =#
    reset!(prec_id_to_score)
    #Query should not match any fragments becuase fragment mass does not match 
    queryFragment!(prec_id_to_score,
                rt_bins[1].last_bin,
                rt_bins[1].first_bin,
                rt_bins[2].last_bin,
                fragment_bins,
                fragments,
                200.0f0,
                200.15f0,
                299.0f0,
                301.0f0)
    @test prec_id_to_score.counts[1] == 0
    reset!(prec_id_to_score)

    #Query should not match any fragments becuase of precursor mass
    queryFragment!(prec_id_to_score,
                rt_bins[1].last_bin,
                rt_bins[1].first_bin,
                rt_bins[2].last_bin,
                fragment_bins,
                fragments,
                100.05f0,
                100.15f0,
                1000.0f0,
                1008.0f0)
    @test prec_id_to_score.counts[1] == 0
    reset!(prec_id_to_score)
end
#########
#Test sscans 
#=
@testset "search_scan" begin
    n_precursors = 1
    prec_id_to_score = Counter(UInt32, UInt8, n_precursors)

    masses = allowmissing(Float32[
        100.1001, 100.11, 301.1, 500.2, 1000.0
    ])
    intensities = allowmissing(Float32[
        100.0f0, 100.0f0, 100.0f0, 100.2f0, 100.0f0
    ])
    mass_err_model = MassErrorModel(
        0.0f0,
        (10.0f0, 10.0f0)
    )

    #Should hit all fragments in the first rt bin. 
    searchScan!(
        prec_id_to_score,
        rt_bins,
        fragment_bins,
        fragments,
        masses, 
        intensities,
        one(Int64),#rt_bin_idx
        1.9f0,#irt_high
        mass_err_model,
        300.0f0, #prec_mz
        4.0f0, #prec_tol
        (1, 0), #isotope_err_bounds
    )
    @test sum(prec_id_to_score.counts)==6
    reset!(prec_id_to_score)

    #Increase tolerance to add three matches 
    mass_err_model = MassErrorModel(
        0.0f0,
        (100.0f0, 100.0f0)
    )

    #Should hit all fragments in the first rt bin. 
    #Should hit all fragments in the first rt bin. 
    searchScan!(
        prec_id_to_score,
        rt_bins,
        fragment_bins,
        fragments,
        masses, 
        intensities,
        one(Int64),#rt_bin_idx
        1.9f0,#irt_high
        mass_err_model,
        300.0f0, #prec_mz
        4.0f0, #prec_tol
        (1, 0), #isotope_err_bounds
    )
    @test sum(prec_id_to_score.counts)==9
    reset!(prec_id_to_score)

    #Increase mass minimum mass tolerance to hig all bins. 
    mass_err_model = MassErrorModel(
        0.0f0,
        (1000.0f0, 1000.0f0)
    )

    searchScan!(
        prec_id_to_score,
        rt_bins,
        fragment_bins,
        fragments,
        masses, 
        intensities,
        one(Int64),#rt_bin_idx
        1.9f0,#irt_high
        mass_err_model,
        300.0f0, #prec_mz
        4.0f0, #prec_tol
        (1, 0), #isotope_err_bounds
    )
    @test sum(prec_id_to_score.counts)==12
    reset!(prec_id_to_score)

    #Increase RT tolerance 
    searchScan!(
        prec_id_to_score,
        rt_bins,
        fragment_bins,
        fragments,
        masses, 
        intensities,
        one(Int64),#rt_bin_idx
        2.1f0,#irt_high
        mass_err_model,
        300.0f0, #prec_mz
        4.0f0, #prec_tol
        (1, 0), #isotope_err_bounds
    )
    @test sum(prec_id_to_score.counts)==27
    reset!(prec_id_to_score)


    #Change precursor mass so there are no hits 
    searchScan!(
        prec_id_to_score,
        rt_bins,
        fragment_bins,
        fragments,
        masses, 
        intensities,
        one(Int64),#rt_bin_idx
        2.1f0,#irt_high
        mass_err_model,
        304.01f0, #prec_mz
        4.0f0, #prec_tol
        (0, 0), #isotope_err_bounds
    )
    @test sum(prec_id_to_score.counts)==0
    reset!(prec_id_to_score)

    #Change precursor mass so there are no hits 
    searchScan!(
        prec_id_to_score,
        rt_bins,
        fragment_bins,
        fragments,
        masses, 
        intensities,
        one(Int64),#rt_bin_idx
        2.1f0,#irt_high
        mass_err_model,
        295.99f0, #prec_mz
        4.0f0, #prec_tol
        (0, 0), #isotope_err_bounds
    )
    @test sum(prec_id_to_score.counts)==0
    reset!(prec_id_to_score)
end
=#


