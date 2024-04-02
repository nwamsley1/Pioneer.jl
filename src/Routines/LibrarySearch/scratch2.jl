test_values = rand(Float32, 100)
sorted_out = zeros(Float32, 100)
sorted_sub_ranges = Vector{UnitRange{UInt32}}([1:10, 11:15, 12:12, 95:100, 85:87])
for sub_range in sorted_sub_ranges
    sort!(@view(test_values[sub_range]))
end
include("src/Routines/LibrarySearch/scratch3.jl")
test_ttree = TournamentTree(
    Vector{TTreeNode{UInt32, Float32}}(undef, 1),
    one(UInt32)
)

buildTTree!(
    test_ttree,
    sorted_sub_ranges,
    UInt32(length( sorted_sub_ranges)),
    test_values
)
n = 1
max_elements = sum([(last(x) - first(x) + 1) for x in sorted_sub_ranges])
while n < max_elements 
removeSmallestElement!(
    test_ttree,
    test_values,
    sorted_out,
    n,
)
n += 1
end