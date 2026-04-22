using Test
using Pioneer: save_multipage_pdf
using Plots

@testset "PDF utils run headless" begin
    output_path = tempname() * ".pdf"
    plot_obj = Plots.plot(1:3, 1:3, title = "headless pdf smoke test")

    save_multipage_pdf([plot_obj], output_path)

    @test isfile(output_path)
    @test filesize(output_path) > 0
end
