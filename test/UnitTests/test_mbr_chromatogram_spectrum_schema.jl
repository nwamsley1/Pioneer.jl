using DataFrames
using Test

@testset "MS2 chromatogram points expose MBR spectrum fields" begin
    chroms = Pioneer.MS2ChromObject[
        Pioneer.MS2ChromObject(
            12.5f0,
            3.0f0,
            UInt32(7),
            UInt32(11),
            100.0f0,
            1.0f0, 2.0f0, 3.0f0, 4.0f0,
            5.0f0, 6.0f0, 7.0f0, 8.0f0,
            8.0f0, 7.0f0, 6.0f0, 5.0f0,
            4.0f0, 3.0f0, 2.0f0, 1.0f0,
        ),
    ]
    df = DataFrame(chroms)

    @test :spectrum_intensity in propertynames(df)
    @test :shadow_frag1_int in propertynames(df)
    @test :shadow_frag8_int in propertynames(df)
    @test :fitted_frag1_int in propertynames(df)
    @test :fitted_frag8_int in propertynames(df)
    @test df.shadow_frag8_int == Float32[8.0]
    @test df.fitted_frag1_int == Float32[8.0]
end
