using PioneerAOTCanary

values = sort!(Float32[mod(i * 37, 257) / 256 for i in 1:257])
for target in Float32[-1, 0, 0.25, 0.5, 0.75, 1, 2]
    expected = findfirst(x -> x >= target, values)
    expected_index = isnothing(expected) ? length(values) + 1 : expected
    @assert PioneerAOTCanary.simd_find_first_ge(
        values,
        1,
        length(values),
        target,
    ) == expected_index
end

x = collect(range(-1.0, 1.0; length=257))
y = collect(range(0.5, 1.5; length=257))
@assert PioneerAOTCanary.weighted_sum(x, y) ≈ sum(x .* y)

println("Pioneer AOT kernel workload: PASS")
