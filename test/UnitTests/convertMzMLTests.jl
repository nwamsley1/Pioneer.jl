using Test
using EzXML
using Pioneer: init_spectrum_dict, parsePrecursorList!, parseScanDictToScanElement

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
