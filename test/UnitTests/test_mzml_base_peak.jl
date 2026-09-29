# convertMzML: each scan's base peak is its most intense stored peak, even when the mzML reports another
# (MS:1000504, MS:1000505): the Huber weight ceiling assumes no intensity in the scan exceeds it.

using Test
using Pioneer: init_spectrum_dict, parseScanDictToScanElement, MzMLCycleIndexTracker, array_base_peak

function _ms1_dict()
    d = init_spectrum_dict()
    d["ms level"] = "1"; d["scan start time"] = "1.5"
    d["scan window lower limit"] = "100"; d["scan window upper limit"] = "1500"
    d
end

@testset "mzML base peak" begin
    mz = Float32[200, 300, 400]; it = Float32[5, 50, 20]
    s = parseScanDictToScanElement(_ms1_dict(), 1, mz, it, true, MzMLCycleIndexTracker())
    @test s.basePeakMz == 300f0 && s.basePeakIntensity == 50f0
    m, i = array_base_peak(Float32[], Float32[])
    @test isnan(m) && i == 0f0
end
