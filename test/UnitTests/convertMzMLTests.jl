using Test
using Arrow
using DataFrames
using EzXML
using Pioneer: build_convert_mzml_output_dir, get_convert_mzml_output_paths,
               init_spectrum_dict, parsePrecursorList!, parseScanDictToScanElement,
               parse_convert_mzml_args, partition_convert_mzml_queue,
               process_mzml_file

function write_test_mzml(path::String; n_spectra::Int, spectrum_body::String = "")
    spectra = join([
        """
        <spectrum id="scan=$(i)">
          $spectrum_body
        </spectrum>
        """
        for i in 1:n_spectra
    ], "\n")
    write(path, """
    <mzML xmlns="http://psi.hupo.org/ms/mzml">
      <run>
        <spectrumList count="$n_spectra">
          $spectra
        </spectrumList>
      </run>
    </mzML>
    """)
    return path
end

function write_test_scan_numbers(path::String, scan_numbers::Vector{Int32})
    Arrow.write(path, DataFrame(scanNumber = scan_numbers))
    return path
end

@testset "convertMzML metadata parsing" begin
    spectrum_dict = Dict{String, String}(
        "ms level" => "2",
        "total ion current" => "",
        "collision energy" => "30.0",
        "spectrum title" => "",
        "scan start time" => "10.5",
        "scan window upper limit" => "900.0",
        "scan window lower limit" => "100.0",
        "isolation window target m/z" => "500.0",
        "isolation window lower offset" => "1.0",
        "isolation window upper offset" => "2.0"
    )
    mz_array = Float32[100.0, 200.0, 300.0]
    intensity_array = Float32[10.0, 25.0, 15.0]

    scan = parseScanDictToScanElement(
        spectrum_dict,
        7,
        mz_array,
        intensity_array,
        true
    )

    @test scan.TIC == Float32(50.0)
    @test scan.centerMz == Float32(500.5)
    @test scan.isolationWidthMz == Float32(3.0)
    @test scan.collisionEnergyField == Float32(30.0)
end

@testset "convertMzML required metadata validation" begin
    spectrum_dict = Dict{String, String}(
        "ms level" => "1",
        "total ion current" => "",
        "collision energy" => "",
        "spectrum title" => "",
        "scan start time" => "",
        "scan window upper limit" => "900.0",
        "scan window lower limit" => "100.0",
        "isolation window target m/z" => "",
        "isolation window lower offset" => "",
        "isolation window upper offset" => ""
    )

    error = try
        parseScanDictToScanElement(
            spectrum_dict,
            3,
            Float32[100.0],
            Float32[1.0],
            true
        )
        nothing
    catch err
        err
    end

    @test error isa ArgumentError
    @test occursin("scan start time", sprint(showerror, error))
    @test occursin("spectrum 3", sprint(showerror, error))
end

@testset "convertMzML activation accession parsing" begin
    precursor_list = EzXML.root(EzXML.parsexml("""
    <precursorList>
      <precursor>
        <activation>
          <cvParam accession="MS:1000133" name="collision-induced dissociation" value=""/>
          <cvParam accession="MS:1000045" name="collision energy" value="27.5"/>
        </activation>
      </precursor>
      <precursor>
        <activation>
          <cvParam accession="MS:1000422" name="beam-type collision-induced dissociation" value=""/>
          <cvParam accession="MS:1000045" name="collision energy" value="31.0"/>
        </activation>
      </precursor>
    </precursorList>
    """))

    spectrum_dict = init_spectrum_dict()
    parsePrecursorList!(spectrum_dict, precursor_list)

    @test spectrum_dict["collision energy"] == "31.0"
end

@testset "convertMzML ignores unsupported activation types" begin
    precursor_list = EzXML.root(EzXML.parsexml("""
    <precursorList>
      <precursor>
        <activation>
          <cvParam accession="MS:9999999" name="unsupported activation" value=""/>
          <cvParam accession="MS:1000045" name="collision energy" value="42.0"/>
        </activation>
      </precursor>
    </precursorList>
    """))

    spectrum_dict = init_spectrum_dict()
    parsePrecursorList!(spectrum_dict, precursor_list)

    @test isempty(spectrum_dict["collision energy"])
end

@testset "convertMzML CLI argument parsing" begin
    options = parse_convert_mzml_args([
        "input_dir",
        "-o", "custom_out",
        "--skip-existing",
        "-n", "4",
        "--include-scan-header",
    ])

    @test options.mzml_path == "input_dir"
    @test options.output_dir == "custom_out"
    @test options.skip_existing
    @test options.concurrent_files == 4
    @test !options.skip_scan_header

    compat_options = parse_convert_mzml_args(["input_dir", "false"])
    @test compat_options.mzml_path == "input_dir"
    @test !compat_options.skip_scan_header
end

@testset "convertMzML output directory and queue planning" begin
    temp_dir = mktempdir()
    mzml_dir = joinpath(temp_dir, "mzml")
    mkpath(mzml_dir)

    mzml_paths = [
        write_test_mzml(joinpath(mzml_dir, "complete.mzML"), n_spectra = 3),
        write_test_mzml(joinpath(mzml_dir, "incomplete.mzML"), n_spectra = 3),
        write_test_mzml(joinpath(mzml_dir, "missing.mzML"), n_spectra = 2),
    ]

    output_dir = build_convert_mzml_output_dir(mzml_dir, "")
    @test output_dir == joinpath(mzml_dir, "arrow_out")

    output_paths = get_convert_mzml_output_paths(output_dir, mzml_paths)
    write_test_scan_numbers(output_paths[1], Int32[1, 2, 3])
    write_test_scan_numbers(output_paths[2], Int32[1, 2])

    queue = partition_convert_mzml_queue(mzml_paths, output_paths; skip_existing = true)

    @test queue.files_to_convert == [2, 3]
    @test queue.skipped_complete_files == 1
    @test queue.reconverted_incomplete_files == 1
    @test queue.missing_output_files == 1
end

@testset "convertMzML skips profile-mode files" begin
    temp_dir = mktempdir()
    mzml_path = write_test_mzml(
        joinpath(temp_dir, "profile.mzML");
        n_spectra = 1,
        spectrum_body = """
        <cvParam accession="MS:1000128" name="profile spectrum" value=""/>
        """
    )
    output_path = joinpath(temp_dir, "profile.arrow")

    converted = process_mzml_file(mzml_path, output_path, true)

    @test !converted
    @test !isfile(output_path)
end
