@testset "chromatogram integration trace ordering" begin
    chromatograms = DataFrame(
        precursor_idx = UInt32[10, 10, 10, 10, 11],
        isotopes_captured = [(Int8(0), Int8(1)), (Int8(0), Int8(1)),
                             (Int8(1), Int8(2)), (Int8(1), Int8(2)),
                             (Int8(0), Int8(1))],
        rt = Float32[2.0, 4.0, 1.0, 3.0, 1.5],
        scan_idx = UInt32[20, 40, 10, 30, 15],
        intensity = Float32[2.0, 4.0, 1.0, 3.0, 5.0],
        precursor_fraction_transmitted = Float32[0.7, 0.7, 0.35, 0.35, 0.8],
    )

    Pioneer.sort_chromatograms_for_integration!(chromatograms)

    target_rows = findall(==(UInt32(10)), chromatograms.precursor_idx)
    @test chromatograms.rt[target_rows] == Float32[1.0, 2.0, 3.0, 4.0]
    @test chromatograms.scan_idx[target_rows] == UInt32[10, 20, 30, 40]
    @test chromatograms.isotopes_captured[target_rows] == [
        (Int8(1), Int8(2)),
        (Int8(0), Int8(1)),
        (Int8(1), Int8(2)),
        (Int8(0), Int8(1)),
    ]
end
