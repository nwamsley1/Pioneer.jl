@testset "binaryRangeQuery.jl" begin
    using Missings
    arr = allowmissing([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    @test binaryGetNearest(arr, 5, 2, 2) == 5

    arr = allowmissing([10.1, 10.19, 10.2, 10.21, 10.22, 10.23, 10.3])

    #Ties are broken by first item (10.21 beats 10.22)
    @test binaryGetNearest(arr, 10.215, 0.1, 0.1) == 4
    @test binaryGetNearest(arr, 10.215, 0.0001, 0.1) == 5

    #Empty array input
    @test binaryGetNearest(allowmissing(Int[]), 10, 5, 3) == 0

    #query is out of bounds of the array
    @test binaryGetNearest(allowmissing(Int[10]), 5, 5, 3) == 0
    @test binaryGetNearest(allowmissing(Int[10]), 17, 5, 3) == 0
    @test binaryGetNearest(allowmissing(Int[10]), 14, 5, 3) == 1

    #Condition is > or < not >= or <=
    @test binaryGetNearest(allowmissing(Int[10]), 14, 4, 3) == 0
    @test binaryGetNearest(allowmissing(Int[10]), 7, 4, 3) == 0
    @test binaryGetNearest(allowmissing(Int[10]), 7, 4, 4) == 1

end