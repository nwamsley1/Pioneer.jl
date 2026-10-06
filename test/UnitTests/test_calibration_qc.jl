using Test, DataFrames, Random
using Pioneer

struct TestCalibrationContext
    calibration_qc::Pioneer.CalibrationQCState
    output::String
end
Pioneer.getDataOutDir(c::TestCalibrationContext) = c.output
Pioneer.getParsedFileName(::TestCalibrationContext, i::Int) = "run_$i"

@testset "Calibration QC compact records and independent selection" begin
    @test sizeof(Pioneer.CalibrationQCStatus) == 1
    @test isbitstype(Pioneer.CalibrationQCRecord)
    @test sizeof(Pioneer.CalibrationQCRecord) <= 24
    state = Pioneer.CalibrationQCState()
    normal = Pioneer.assess_calibration_qc(:nce,100,(NaN,NaN,0.1,0.1))
    suspicious = Pioneer.assess_calibration_qc(:nce,100,(NaN,NaN,0.1,0.8))
    warning = Pioneer.assess_calibration_qc(:nce,3,(NaN,NaN,0.1,0.8); min_support=50)
    failed = Pioneer.assess_calibration_qc(:nce,3,(NaN,NaN,NaN,NaN); failed=true)
    @test normal.status == Pioneer.QC_NORMAL
    @test suspicious.status == Pioneer.QC_SUSPICIOUS
    @test warning.status == Pioneer.QC_WARNING
    @test failed.status == Pioneer.QC_FAILED
    @test occursin("endpoint_fraction", Pioneer.calibration_qc_reason(:nce, suspicious))
    selected_nce = Int[]; selected_rt = Int[]
    for i in 1:6000
        Pioneer.record_calibration_qc!(state,:nce,i,suspicious)
        Pioneer.record_calibration_qc!(state,:rt,i,normal)
        Pioneer.select_calibration_plot!(state,:nce,i) && push!(selected_nce, i)
        Pioneer.select_calibration_plot!(state,:rt,i) && push!(selected_rt, i)
    end
    @test selected_nce == collect(1:100)
    @test selected_rt == collect(1:50)
    @test length(state.selected[5]) == 100
    @test length(state.selected[1]) == 50
    @test all(isempty, state.pages)
    @test Base.summarysize(state) < 1_000_000
    # A late RT failure gets its own allowance despite the exhausted NCE budget.
    Pioneer.record_calibration_qc!(state,:rt,6001,failed)
    @test Pioneer.select_calibration_plot!(state,:rt,6001)
    @test Pioneer.select_calibration_plot!(state,:rt,6001)
    @test length(state.selected[1]) == 51
    table = DataFrame(file_idx=1:6001)
    Pioneer.add_calibration_qc_columns!(table,state,6001)
    @test table.nce_calibration_qc[6000] == "suspicious"
    @test table.rt_calibration_qc[6001] == "failed"
    @test all(==("not_assessed"),table.ms1_mass_calibration_qc)
    @test !hasproperty(table,:calibration_qc)
    @test !hasproperty(table,:huber_calibration_qc)
    @test Pioneer.assess_calibration_qc(:nce,100,(NaN,NaN,0.5,0.5)).status == Pioneer.QC_NORMAL
    @test Pioneer.assess_calibration_qc(:ms1_mass,100,(NaN,0.1,0.25,0.1)).status == Pioneer.QC_NORMAL
    @test Pioneer.assess_calibration_qc(:quadrupole,100,(0.9,0.1,NaN,1)).status == Pioneer.QC_SUSPICIOUS
end

@testset "Within-file RT and mass error diagnostics" begin
    rts = Float32.(range(0,60,length=1000))
    irts = 2 .* rts
    good = Pioneer.rt_calibration_qc(rts,irts,x->2x,1.5;min_support=300)
    @test good.status == Pioneer.QC_NORMAL
    poor = Pioneer.rt_calibration_qc(rts,irts,x->x,1.5;min_support=300)
    @test poor.status == Pioneer.QC_SUSPICIOUS
    narrow = Pioneer.rt_calibration_qc(rts,irts,x->2x,1.5;min_support=300)
    @test narrow.status == Pioneer.QC_NORMAL
    @test Pioneer.rt_calibration_qc(Float32[],Float32[],Pioneer.IdentityModel(),Inf;min_support=300).status == Pioneer.QC_WARNING
    model = Pioneer.MassErrorModel(0f0,(10f0,10f0))
    samples = [Pioneer.MassErrSample(500f0,500f0,100f0,Float32(i)) for i in 1:100]
    @test Pioneer.ms2_calibration_qc(samples,model).status == Pioneer.QC_NORMAL
    bad = [Pioneer.MassErrSample(500f0,500.1f0,100f0,Float32(i)) for i in 1:100]
    @test Pioneer.ms2_calibration_qc(bad,model).status == Pioneer.QC_SUSPICIOUS
    @test Pioneer.ms2_calibration_qc(samples,model;fallback=true).status == Pioneer.QC_WARNING
    @test Pioneer.ms2_calibration_qc(samples,Pioneer.MassErrorModel(0f0,(0f0,0f0))).status == Pioneer.QC_FAILED
end

@testset "Incremental rendering and streamed PDF" begin
    mktempdir() do dir
        context = TestCalibrationContext(Pioneer.CalibrationQCState(),dir)
        rec = Pioneer.assess_calibration_qc(:nce,100,(NaN,NaN,0.1,0.8))
        Pioneer.record_calibration_qc!(context.calibration_qc,:nce,51,rec)
        mz = Float32.(range(400,900,length=100))
        best_nce = DataFrame(prec_mz=mz,nce=fill(30f0,100),charge=fill(UInt8(2),100))
        model = Pioneer.fit_binned_median_nce(best_nce.prec_mz,best_nce.nce,best_nce.charge,30f0)
        Random.seed!(123)
        expected = rand()
        Random.seed!(123)
        Pioneer.plot_nce_calibration!(context,51,best_nce,Float32.(21:40),model,UInt8[2])
        @test rand() == expected
        pages = context.calibration_qc.pages[5]
        @test length(pages) == 1
        @test isfile(only(pages))
        pdf = joinpath(dir,"qc.pdf")
        Pioneer.PDFGenerator.write_pngs_to_pdf(vcat(pages,pages),pdf)
        bytes = read(pdf)
        @test startswith(String(copy(bytes[1:8])),"%PDF-1.4")
        @test occursin("/Count 2",String(copy(bytes)))
        @test endswith(String(copy(bytes[end-5:end])),"%%EOF\n")
    end
end

@testset "Quadrupole QC uses the fitted observations" begin
    truth = Pioneer.RazoQuadParams(0.8f0,0.9f0,3f0,4f0)
    x0 = Float32.(range(-1.5,1.5,length=200))
    x1 = x0 .+ 0.5f0
    table = DataFrame(x0=x0,x1=x1,yt=[Pioneer.F(truth,a,b) for (a,b) in zip(x0,x1)])
    fitted, initial, metrics = Pioneer.fit_quad_model(table,2.0)
    @test all(isfinite, (fitted.al,fitted.ar,fitted.bl,fitted.br))
    @test 0 <= metrics[1] <= 1
    @test isfinite(metrics[2]) && metrics[2] >= 0
    @test metrics[4] in (0.0,1.0)
end

@testset "Report failure preserves previous output" begin
    mktempdir() do dir
        path = joinpath(dir,"report.pdf")
        write(path,"previous report")
        @test_throws SystemError Pioneer.PDFGenerator.write_pngs_to_pdf([joinpath(dir,"missing.png")],path)
        @test read(path,String) == "previous report"
        @test readdir(dir) == ["report.pdf"]
    end
end

struct CalibrationTestSpectra end
Base.length(::CalibrationTestSpectra) = 2
Pioneer.getMsOrder(::CalibrationTestSpectra, i) = i == 1 ? 1 : 2
Pioneer.getRetentionTime(::CalibrationTestSpectra, i) = Float32(i)
Pioneer.getImScans(::CalibrationTestSpectra) = nothing   # no ion mobility
Pioneer.getMzArray(::CalibrationTestSpectra, i) = Float32[500,500+Pioneer.C13_C12_MASS_DIFF_F32/2,500+Pioneer.C13_C12_MASS_DIFF_F32]
Pioneer.getPeaks!(::Pioneer.PeakDecodeBuffer, s::CalibrationTestSpectra, i::Integer) = (Pioneer.getMzArray(s, i), ones(Float32, 3))
struct CalibrationTestLibrary end
struct CalibrationTestPrecursors end
Pioneer.getSpecLib(::TestCalibrationContext) = CalibrationTestLibrary()
Pioneer.getPrecursors(::CalibrationTestLibrary) = CalibrationTestPrecursors()
Pioneer.getMz(::CalibrationTestPrecursors) = Float32[500]
Pioneer.getCharge(::CalibrationTestPrecursors) = UInt8[2]

@testset "MS1 diagnostic coordinates do not change calibration observations" begin
    context = TestCalibrationContext(Pioneer.CalibrationQCState(),"")
    psms = DataFrame(scan_idx=UInt32[2],precursor_idx=UInt32[1])
    coords = (Float32[],Float32[])
    precursor_ids = UInt32[]
    baseline = Pioneer.collect_ms1_residuals(CalibrationTestSpectra(),psms,context,1)
    observed = Pioneer.collect_ms1_residuals(CalibrationTestSpectra(),psms,context,1;
        qc_coordinates=coords, qc_precursor_ids=precursor_ids)
    @test observed == baseline == zeros(Float32,3)
    @test coords[1] == Pioneer.getMzArray(CalibrationTestSpectra(),1)
    @test coords[2] == ones(Float32,3)
    @test precursor_ids == ones(UInt32,3)
    @test Pioneer.assess_calibration_qc(:nce,100,(0.25,NaN,0.1,0.1)).status == Pioneer.QC_SUSPICIOUS
end

@testset "Robust peptide-supported MS1 calibration QC" begin
    n = 1000
    ids = collect(1:n)
    coords = (Float64.(1:n), Float64.(1:n))
    assess(errors; groups=ids, coordinates=coords, tolerance=10.0) =
        Pioneer.ms1_calibration_diagnostics(errors, coordinates, zeros(length(errors)),
            fill(tolerance,length(errors)); group_ids=groups)
    @test assess(zeros(n)).status == Pioneer.QC_NORMAL

    # A gross outlier, or one biased decile, cannot flag a whole run.
    outlier = zeros(n); outlier[1] = 1000
    @test assess(outlier).status == Pioneer.QC_NORMAL
    isolated = [i <= 100 ? 6.0 : 0.0 for i in 1:n]
    @test assess(isolated).status == Pioneer.QC_NORMAL
    persistent = [i <= 200 ? 6.0 : 0.0 for i in 1:n]
    for errors in (persistent, -persistent)
        record = assess(errors)
        @test record.status == Pioneer.QC_SUSPICIOUS
        @test record.metrics[2] ≈ 0.2f0
        @test occursin("persistent_biased_peptide_fraction",
            Pioneer.calibration_qc_reason(:ms1_mass,record))
    end
    # Either axis can reveal a persistent bias; constant coordinates cannot.
    @test assess(persistent; coordinates=(ones(n),coords[2])).metrics[2] ≈ 0.2f0
    @test assess(persistent; coordinates=(coords[1],ones(n))).metrics[2] ≈ 0.2f0
    constant_coords = (ones(n),ones(n))
    @test assess(persistent; coordinates=constant_coords).status == Pioneer.QC_NORMAL

    # Relative tolerance replaces the absolute 5 ppm trend and 10 ppm MAD limits.
    @test assess(persistent; tolerance=30.0).status == Pioneer.QC_NORMAL
    wide = [isodd(i) ? 8.0 : -8.0 for i in 1:n]
    @test assess(wide; tolerance=30.0).status == Pioneer.QC_NORMAL
    @test assess(wide .* 2).metrics[4] == 1f0
    @test assess(wide .* 2).status == Pioneer.QC_SUSPICIOUS
    @test assess([i <= 100 ? 11.0 : 0.0 for i in 1:n]).status == Pioneer.QC_NORMAL
    @test assess([i <= 101 ? 11.0 : 0.0 for i in 1:n]).status == Pioneer.QC_SUSPICIOUS
    @test assess(fill(2.5,n)).status == Pioneer.QC_NORMAL
    global_bias = assess(fill(3.0,n))
    @test global_bias.metrics[3] ≈ 0.3f0
    @test global_bias.status == Pioneer.QC_SUSPICIOUS
    @test assess(fill(-3.0,n)).status == Pioneer.QC_SUSPICIOUS
    # A median beyond 25% is insufficient if its interval crosses that threshold.
    uncertain = vcat(fill(2.0,490),fill(3.0,510))
    @test assess(uncertain; coordinates=constant_coords).metrics[3] == 0f0

    # Repeated PSMs/isotopes do not manufacture peptide support or coverage weight.
    repeated_coords = (repeat(Float64.(1:10),100),repeat(Float64.(1:10),100))
    repeated = assess(fill(6.0,n); groups=repeat(1:10,100),coordinates=repeated_coords)
    @test repeated.status == Pioneer.QC_WARNING
    @test repeated.n == 10
    mixed_errors = vcat(zeros(n),fill(100.0,1000))
    mixed_ids = vcat(ids,fill(n+1,1000))
    mixed_coords = (Float64.(mixed_ids),Float64.(mixed_ids))
    coverage = assess(mixed_errors; groups=mixed_ids,coordinates=mixed_coords)
    @test coverage.metrics[4] ≈ 1 / (n+1)
    @test coverage.status == Pioneer.QC_NORMAL
    @test assess(zeros(n); tolerance=0).status == Pioneer.QC_FAILED
    @test assess(fill(NaN,n)).status == Pioneer.QC_FAILED
    @test_throws DimensionMismatch assess(zeros(n); groups=1:n-1)

    # QC uses the supplied production model, with no additional fits or random draws.
    errors = repeat(Float64[-2,0,2],n)
    sequence_ids = repeat(ids; inner=3)
    sequence_coords = (Float64.(sequence_ids),Float64.(sequence_ids))
    snapshot = copy(errors)
    fit = Pioneer.fit_ms1_model_from_residuals(errors)
    Random.seed!(123)
    expected = rand()
    Random.seed!(123)
    production = Pioneer.ms1_calibration_qc(errors,sequence_coords,fit[1]; group_ids=sequence_ids)
    @test rand() == expected
    @test production.status == Pioneer.QC_NORMAL
    @test production.n == n
    @test errors == snapshot
    @test Pioneer.fit_ms1_model_from_residuals(errors)[2:3] == fit[2:3]
    biased_errors = errors .+ [i <= 200 ? 6.0 : 0.0 for i in sequence_ids]
    biased = Pioneer.ms1_calibration_qc(biased_errors,sequence_coords,fit[1]; group_ids=sequence_ids)
    @test biased.status == Pioneer.QC_SUSPICIOUS
    @test biased.metrics[2] ≈ 0.2f0
    permutation = collect(length(errors):-1:1)
    reordered = Pioneer.ms1_calibration_qc(biased_errors[permutation],
        (sequence_coords[1][permutation],sequence_coords[2][permutation]),fit[1];
        group_ids=sequence_ids[permutation])
    @test isequal(reordered.metrics,biased.metrics)
    @test Pioneer.ms1_calibration_qc(errors[1:147],
        (sequence_coords[1][1:147],sequence_coords[2][1:147]),fit[1];
        group_ids=sequence_ids[1:147]).status == Pioneer.QC_WARNING
    @test Pioneer.ms1_calibration_qc(errors[1:150],
        (sequence_coords[1][1:150],sequence_coords[2][1:150]),fit[1];
        group_ids=sequence_ids[1:150]).status == Pioneer.QC_NORMAL
    @test Pioneer.ms1_calibration_qc(errors,sequence_coords,nothing;
        group_ids=sequence_ids).status == Pioneer.QC_WARNING
    @test Pioneer.ms1_calibration_qc(zeros(n),coords,Pioneer.MassErrorModel(0f0,(0f0,0f0));
        group_ids=ids).status == Pioneer.QC_FAILED
    @test_throws DimensionMismatch Pioneer.ms1_calibration_qc(errors,sequence_coords,fit[1];
        group_ids=sequence_ids[1:end-1])
end

@testset "Robust peptide-supported MS2 calibration QC" begin
    model = Pioneer.MassErrorModel(0f0, (10f0, 10f0))
    sample(i, error) = Pioneer.MassErrSample(500f0, 500f0 + Float32(error),
        100f0, Float32(i), UInt32(i))
    samples = [sample(i, 0) for i in 1:1000]
    ids = collect(1:1000)
    assess(s; groups=ids) = Pioneer.ms2_calibration_diagnostics(s, fill(model,length(s)); group_ids=groups)
    @test assess(samples).record.status == Pioneer.QC_NORMAL

    # One gross outlier cannot move a peptide median in a well-supported bin.
    outlier = copy(samples)
    outlier[1] = sample(1, 1)
    @test assess(outlier).record.status == Pioneer.QC_NORMAL

    # An isolated biased decile is a local diagnostic, not a whole-run flag.
    isolated = [sample(i, i <= 100 ? 0.0025 : 0) for i in 1:1000]
    @test assess(isolated).record.status == Pioneer.QC_NORMAL
    persistent = [sample(i, i <= 200 ? 0.0025 : 0) for i in 1:1000]
    result = assess(persistent)
    @test result.record.status == Pioneer.QC_SUSPICIOUS
    @test result.record.metrics[2] ≈ 0.2f0
    @test occursin("persistent_biased_peptide_fraction", Pioneer.calibration_qc_reason(:ms2_mass,result.record))
    @test count(b -> b.persistent_bias, result.bins) == 2

    # Centered data can still have unacceptable actual-window coverage.
    broad = [sample(i, isodd(i) ? 0.1 : -0.1) for i in 1:1000]
    @test assess(broad).record.metrics[4] == 1f0
    @test assess(broad).record.status == Pioneer.QC_SUSPICIOUS

    # Fragment replication cannot manufacture independent peptide support.
    repeated = repeat(samples[1:10], 100)
    @test assess(repeated; groups=repeat(1:10,100)).record.status == Pioneer.QC_WARNING
    @test_throws DimensionMismatch assess(samples; groups=1:999)
    sequences = string.(ids)
    @test Pioneer.ms2_calibration_qc(samples,model,sequences).status == Pioneer.QC_NORMAL
    @test Pioneer.ms2_calibration_qc(persistent,model,sequences).status == Pioneer.QC_SUSPICIOUS
    @test Pioneer.ms2_calibration_qc(repeated,model,sequences).status == Pioneer.QC_WARNING
    @test Pioneer.ms2_calibration_qc(samples,model,sequences;fallback=true).status == Pioneer.QC_WARNING
    @test Pioneer.ms2_calibration_qc(samples,model,String[]).status == Pioneer.QC_FAILED
end


@testset "Robust production RT calibration QC" begin
    n = 1000
    rts = collect(range(0,25,length=n))
    ids = collect(1:n)
    assess(errors; tolerance=1.5, groups=ids, rt=rts, model=identity, fallback=false) =
        Pioneer.rt_calibration_qc(rt,rt .+ errors,model,tolerance; group_ids=groups,fallback)
    @test assess(zeros(n)).status == Pioneer.QC_NORMAL
    # Seer-like main population plus sparse large endpoint mismatches: bin means
    # move far past the old 1%-of-span cutoff, but peptide medians remain centered.
    outliers = zeros(n); outliers[1:20] .= -15
    @test assess(outliers).status == Pioneer.QC_NORMAL
    @test assess(outliers).metrics[4] ≈ 0.02f0
    isolated = [i <= 100 ? 0.6 : 0.0 for i in 1:n]
    @test assess(isolated).status == Pioneer.QC_NORMAL
    persistent = [i <= 200 ? 0.6 : 0.0 for i in 1:n]
    for errors in (persistent,-persistent)
        record = assess(errors)
        @test record.status == Pioneer.QC_SUSPICIOUS
        @test record.metrics[2] ≈ 0.2f0
        @test occursin("persistent_biased_peptide_fraction",Pioneer.calibration_qc_reason(:rt,record))
    end
    @test assess(persistent; tolerance=3).status == Pioneer.QC_NORMAL
    @test assess([i <= 100 ? 0.6 : i <= 200 ? -0.6 : 0.0 for i in 1:n]).status == Pioneer.QC_NORMAL
    @test assess(fill(0.375,n)).status == Pioneer.QC_NORMAL
    @test assess(fill(0.4,n)).status == Pioneer.QC_SUSPICIOUS
    @test assess(fill(-0.4,n)).status == Pioneer.QC_SUSPICIOUS
    @test assess(vcat(fill(0.3,490),fill(0.45,510));rt=ones(n)).status == Pioneer.QC_NORMAL
    @test assess([isodd(i) ? 1.0 : -1.0 for i in 1:n]).status == Pioneer.QC_NORMAL
    @test assess([i <= 100 ? 2.0 : 0.0 for i in 1:n]).status == Pioneer.QC_NORMAL
    coverage = assess([i <= 101 ? 2.0 : 0.0 for i in 1:n])
    @test coverage.status == Pioneer.QC_WARNING
    @test occursin("outside_tolerance_fraction",Pioneer.calibration_qc_reason(:rt,coverage))
    @test assess(persistent .+ [i <= 101 ? 2.0 : 0.0 for i in 1:n]).status == Pioneer.QC_SUSPICIOUS
    @test assess(zeros(n);fallback=true).status == Pioneer.QC_WARNING
    @test assess(zeros(n);tolerance=0).status == Pioneer.QC_FAILED
    @test assess(fill(NaN,n)).status == Pioneer.QC_FAILED
    @test assess(zeros(n);model=x->NaN).status == Pioneer.QC_FAILED
    @test assess(zeros(n);tolerance=Inf).status == Pioneer.QC_WARNING
    @test assess(zeros(n);model=Pioneer.IdentityModel()).status == Pioneer.QC_WARNING
    @test assess(zeros(n);groups=repeat(1:10,100)).n == 10
    @test assess(zeros(n);groups=repeat(1:10,100)).status == Pioneer.QC_WARNING
    # Flat model extensions without observations do not imply bad calibration.
    supported_rts = collect(range(3,22,length=n))
    @test assess(zeros(n);rt=supported_rts,model=x->clamp(x,3,22)).status == Pioneer.QC_NORMAL
    @test_throws DimensionMismatch assess(zeros(n);groups=1:n-1)
    Random.seed!(123); expected = rand(); Random.seed!(123)
    assess(persistent)
    @test rand() == expected
    @test !isdefined(Pioneer,:heldout_ms1_calibration_qc)
    @test !isdefined(Pioneer,:heldout_ms2_calibration_qc)
end
