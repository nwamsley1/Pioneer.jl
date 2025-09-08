# Test Authoring Guidance

Please avoid fully qualified calls like `Pioneer.function_name(...)` inside tests.

Instead, import the functions you need at the top of each test file and call them unqualified. This keeps tests concise and avoids module coupling.

Example:

```
using Test
using Pioneer: BasicMassSpecData, FilteredMassSpecData,
               getScanHeader, getScanNumber, getRetentionTime,
               getMzArray, getIntensityArray,
               compute_rt_bins, sort_scans_by_peak_density, create_priority_order

@testset "example" begin
    ms = BasicMassSpecData("path.arrow")
    @test getScanNumber(ms, 1) > 0
end
```

If you need many functions, group them and keep the import list explicit. Do not rely on `using Pioneer` alone, and do not prefix calls with `Pioneer.` in the test body.

