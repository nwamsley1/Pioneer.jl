using Arrow, DataFrames
t = Arrow.Table(expanduser("~/Projects/pioneer-gui-test/rawtest/arrow_out/HELA_MOUSESIL_METHTEST_07312024_06.arrow"))
println(Tables.columnnames(t)); println("scans ", length(t.msOrder), "; MS1 ", count(==(1), t.msOrder), " MS2 ", count(==(2), t.msOrder))
for i in 1:6
    mz = t.mz_array[i]; it = t.intensity_array[i]
    println("scan $i: order=$(t.msOrder[i]) rt=$(t.retentionTime[i]) center=$(t.centerMz[i]) width=$(t.isolationWidthMz[i]) npk=$(length(mz)) eltype=$(eltype(mz)) first=", collect(zip(mz[1:min(3,end)], it[1:min(3,end)])))
end
for c in Tables.columnnames(t)
    c in (:mz_array, :intensity_array) && continue
    println(c, " => ", Tables.getcolumn(t, c)[1:3])
end
println("packetTypes ", unique(t.packetType), " scanHeader sample: ", t.scanHeader[1], " | ", t.scanHeader[2])
