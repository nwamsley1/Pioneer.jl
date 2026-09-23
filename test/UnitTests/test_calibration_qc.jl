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
    @test Pioneer.assess_calibration_qc(:ms1_mass,100,(NaN,10.1,NaN,NaN)).status == Pioneer.QC_SUSPICIOUS
    @test Pioneer.assess_calibration_qc(:quadrupole,100,(0.9,0.1,NaN,1)).status == Pioneer.QC_SUSPICIOUS
end

@testset "Within-file RT and mass error diagnostics" begin
    rts = Float32.(range(0,60,length=1000))
    irts = 2 .* rts
    good = Pioneer.rt_calibration_qc(rts,irts,x->2x,rts;min_support=300)
    @test good.status == Pioneer.QC_NORMAL
    poor = Pioneer.rt_calibration_qc(rts,irts,x->x,rts;min_support=300)
    @test poor.status == Pioneer.QC_SUSPICIOUS
    narrow = Pioneer.rt_calibration_qc(rts,irts,x->2x,Float32[0,180];min_support=300)
    @test narrow.status == Pioneer.QC_SUSPICIOUS
    @test Pioneer.rt_calibration_qc(Float32[],Float32[],Pioneer.IdentityModel(),rts;min_support=300).status == Pioneer.QC_WARNING
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
Pioneer.getMzArray(::CalibrationTestSpectra, i) = Float32[500,500+Pioneer.C13_C12_MASS_DIFF_F32/2,500+Pioneer.C13_C12_MASS_DIFF_F32]
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
    baseline = Pioneer.collect_ms1_residuals(CalibrationTestSpectra(),psms,context,1)
    observed = Pioneer.collect_ms1_residuals(CalibrationTestSpectra(),psms,context,1; qc_coordinates=coords)
    @test observed == baseline == zeros(Float32,3)
    @test coords[1] == Pioneer.getMzArray(CalibrationTestSpectra(),1)
    @test coords[2] == ones(Float32,3)
    @test Pioneer.assess_calibration_qc(:ms1_mass,100,(NaN,1,6,NaN)).status == Pioneer.QC_SUSPICIOUS
    @test Pioneer.assess_calibration_qc(:nce,100,(0.25,NaN,0.1,0.1)).status == Pioneer.QC_SUSPICIOUS
end
