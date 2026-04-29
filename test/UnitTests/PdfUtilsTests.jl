using Test
using Pioneer: PlotSeriesSpec, RtAlignmentPlotSpec, save_multipage_pdf
using Plots

@testset "PDF utils run headless" begin
    output_path = tempname() * ".pdf"
    plot_obj = Plots.plot(1:3, 1:3, title = "headless pdf smoke test")

    save_multipage_pdf([plot_obj], output_path)

    @test isfile(output_path)
    @test filesize(output_path) > 0
end

@testset "PDF utils render plot specs through root Pioneer" begin
    output_path = tempname() * ".pdf"
    plot_spec = RtAlignmentPlotSpec(
        "rt alignment smoke test",
        Float32[1, 2, 3],
        Float32[10, 20, 30],
        Float32[1, 3],
        Float32[10, 30],
        "fit",
        :red,
        :solid,
        2.0f0,
        nothing,
        nothing,
        nothing,
    )

    save_multipage_pdf([plot_spec], output_path)

    @test isfile(output_path)
    @test filesize(output_path) > 0
end
