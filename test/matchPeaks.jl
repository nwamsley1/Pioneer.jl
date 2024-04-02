M = 100
test_matches = [FragmentMatch{Float32}() for _ in range(1, M)];
test_misses = [FragmentMatch{Float32}() for _ in range(1, M)];
test_templates = [
    DetailedFrag{Float32}(0x007a1200, 169.09715f0, Float16(0.08875), false, false, 0x01, 0x02, 0x02, 0x0e, 0x00)
    DetailedFrag{Float32}(0x007a1200, 630.8601f0, Float16(0.0887), true, false, 0x02, 0x0c, 0x02, 0x0f, 0x01)
    DetailedFrag{Float32}(0x007a1200, 704.3943f0, Float16(0.05826), true, false, 0x02, 0x0d, 0x02, 0x13, 0x01)
    DetailedFrag{Float32}(0x007a1200, 747.91034f0, Float16(0.06537), true, false, 0x02, 0x0e, 0x02, 0x12, 0x02)
    DetailedFrag{Float32}(0x007a1200, 747.91034f0, Float16(0.06537), true, false, 0x02, 0x0e, 0x02, 0x12, 0x02) #Duplicate
    DetailedFrag{Float32}(0x007a1200, 800.91034f0, Float16(0.06537), true, false, 0x02, 0x0e, 0x02, 0x12, 0x02) #Duplicate
    DetailedFrag{Float32}(0x007a1200, 962.58514f0, Float16(0.0668), true, false, 0x01, 0x09, 0x02, 0x11, 0x00)
    DetailedFrag{Float32}(0x007a1200, 1668.8107f0, Float16(0.0829), false, false, 0x01, 0x10, 0x02, 0x10, 0x03)
];
masses = allowmissing(Float32[
    169.09715, 
    630.8601f0, 
    500.0f0, #nothing should match to this peak. 
    704.3943f0-0.002f0, #Match to closest of peaks with similar mass
    704.3943f0-0.001f0, 
    747.91034f0, #Two fragments should match to this peak 
    962.58514f0, 
    1668.8107f0
]);
intensities = allowmissing(Float32[
    Float32(1e5),
    Float32(1e5),
    Float32(1e5),
    Float32(1e5),
    Float32(1e5),
    Float32(1e5),
    Float32(1e5),
    Float32(1e5)
]);
mass_err_model = MassErrorModel(
    -0.5f0,
    4000.0f0,
    0.0f0
);
matched_idx, unmatched_idx = matchPeaks!(
    test_matches ,
    test_misses,
    test_templates,
    length(test_templates),
    masses,
    intensities,
    0.0f0,
    mass_err_model,
    (5.0f0, 40.0f0),
    one(UInt32),
    one(UInt32)
)
@test matched_idx == 7
@test unmatched_idx == 1
@test abs(test_matches[3].match_mz - masses[4]) > abs(test_matches[3].match_mz - masses[5]) 
reset!(test_matches, matched_idx)
reset!(test_misses, unmatched_idx)
#Tolerance will be very small for this peak so shouldn't match anymore. 
intensities[4]= Float32(1e10)
intensities[5]= Float32(1e10)
matched_idx, unmatched_idx = matchPeaks!(
    test_matches ,
    test_misses,
    test_templates,
    length(test_templates),
    masses,
    intensities,
    0.0f0,
    mass_err_model,
    (0.0f0, 40.0f0), #change minimum
    one(UInt32),
    one(UInt32)
)
@test matched_idx == 6
@test unmatched_idx == 2
reset!(test_matches, matched_idx)
reset!(test_misses, unmatched_idx)
#Tolerance will be very small for this peak so shouldn't match anymore. 
intensities[5]= Float32(1e1)
matched_idx, unmatched_idx = matchPeaks!(
    test_matches ,
    test_misses,
    test_templates,
    length(test_templates),
    masses,
    intensities,
    0.0f0,
    mass_err_model,
    (0.0f0, 40.0f0), #change minimum
    one(UInt32),
    one(UInt32)
)
@test matched_idx == 7
@test unmatched_idx == 1