# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

struct PioneerScanElement
    mz_array::Vector{Union{Missing, Float32}}
    intensity_array::Vector{Union{Missing, Float32}}
    scanHeader::String
    scanNumber::Int32
    basePeakMz::Float32
    basePeakIntensity::Float32
    packetType::Int32
    retentionTime::Float32
    lowMz::Float32
    highMz::Float32
    TIC::Float32
    centerMz::Union{Missing, Float32}
    isolationWidthMz::Union{Missing, Float32}
    collisionEnergyField::Union{Missing, Float32}
    collisionEnergyEvField::Union{Missing, Float32}
    msOrder::UInt8
end
function parseBinaryDataList(binary_data_list::EzXML.Node)
    mz_array, intensity_array = nothing, nothing
    is_mz_array, is_intensity_array = false, false
    #Assumption here is that the intensities and m/z arrays both have the same precision
    precision = zero(Float32)
    for nl in eachelement(binary_data_list)
        is_mz_array, is_intensity_array = false, false
        for sl in eachelement(nl)
            if sl.name =="cvParam"
                if sl["name"] == "m/z array"
                        is_mz_array = true
                elseif sl["name"] == "intensity array"
                        is_intensity_array = true
                end
                if sl["name"] == "64-bit float"
                    precision = zero(Float64)
                end
            elseif sl.name=="binary"
                if is_mz_array == true
                    mz_array = decodeBinaryArray(sl.content, precision)
                elseif is_intensity_array == true
                    intensity_array = decodeBinaryArray(sl.content, precision)
                end
            end
        end
    end
    return mz_array, intensity_array
end

function decodeBinaryArray(encoded_data::String, ::Float32)
    if length(encoded_data) > 0
        # Decode base64
        decoded_data = base64decode(encoded_data)
        # Decompress (since your data appears to be zlib compressed based on the cvParam)
        decompressed_data = transcode(ZlibDecompressor, decoded_data)   
        # Convert to float array
        return collect(reinterpret(Float32, decompressed_data))::Vector{Float32}
    else
        return Vector{Float32}(undef, 0)
    end
end

function decodeBinaryArray(encoded_data::String, ::Float64)
    if length(encoded_data) > 0
        # Decode base64
        decoded_data = base64decode(encoded_data)
        # Decompress (since your data appears to be zlib compressed based on the cvParam)
        decompressed_data = transcode(ZlibDecompressor, decoded_data)   
        # Convert to float array
        return [Float32(x) for x in reinterpret(Float64, decompressed_data)]::Vector{Float32}
    else
        return Vector{Float32}(undef, 0)
    end
end


function parseScanDictToScanElement(
    spectrum_dict::Dict{String, String},
    scan_number::Int64,
    mz_array::Vector{Float32},
    intensity_array::Vector{Float32},
    skip_scan_header::Bool)::PioneerScanElement

    centerMz,isolationWidthMz = missing, missing
    scanHeader = ""
    if !skip_scan_header
        scanHeader = spectrum_dict["spectrum title"]
    end
    scanNumber = Int32(scan_number)
    basePeakMz = parse(Float32, spectrum_dict["base peak m/z"])
    basePeakIntensity = parse(Float32, spectrum_dict["base peak intensity"])
    retentionTime = parse(Float32, spectrum_dict["scan start time"])
    lowMz = parse(Float32, spectrum_dict["scan window lower limit"])
    highMz = parse(Float32, spectrum_dict["scan window upper limit"])
    msOrder = parse(UInt8, spectrum_dict["ms level"])
    if msOrder > 1
        targetMz = parse(Float32, spectrum_dict["isolation window target m/z"])
        lowerOffset = parse(Float32, spectrum_dict["isolation window lower offset"])
        upperOffset = parse(Float32, spectrum_dict["isolation window lower offset"])
        isolationWidthMz = lowerOffset + upperOffset
        centerMz = targetMz + (upperOffset - lowerOffset)/2.0f0
    end
    TIC = parse(Float32, spectrum_dict["total ion current"])
    dissociation = spectrum_dict["beam-type collision-induced dissociation"]
    collisionEnergyField = (dissociation == "" ? missing : parse(Float32, dissociation))
    return PioneerScanElement(
        allowmissing(mz_array), 
        allowmissing(intensity_array),
        scanHeader,
        scanNumber,
        basePeakMz,
        basePeakIntensity,
        zero(Int32),
        retentionTime,
        lowMz,
        highMz,
        TIC,
        centerMz,
        isolationWidthMz,
        collisionEnergyField,
        zero(Float32),
        msOrder
    )
end

function parseScanCvParam!(
    spectrum_dict::Dict{String, String},
    ScanCvParam::EzXML.Node)
    ScanCvParamName = ScanCvParam["name"]::String
    if haskey(spectrum_dict, ScanCvParamName)
        spectrum_dict[ScanCvParamName] =  ScanCvParam["value"]
    end
end        

function parseScanList!(
    spectrum_dict::Dict{String, String},
    scanList::EzXML.Node)

    for scanListElement in eachelement(scanList)
        if scanListElement.name=="scan"
            for scanListSubElement in eachelement(scanListElement)
                if scanListSubElement.name=="cvParam"
                    parseScanCvParam!(spectrum_dict, scanListSubElement)
                elseif scanListSubElement.name == "scanWindowList"
                    for scanWindow in eachelement(scanListSubElement)
                        for scanWindowElement in eachelement(scanWindow)
                            if scanWindowElement.name=="cvParam"
                                parseScanCvParam!(spectrum_dict, scanWindowElement)
                            end
                        end
                    end
                end
            end
        end
    end

    return  nothing
end

function parsePrecursorList!(
    spectrum_dict::Dict{String, String},
    precursorList::EzXML.Node)

    for precursorListElement in eachelement(precursorList)
        if precursorListElement.name =="precursor"
            for precursorElement in eachelement(precursorListElement)
                if precursorElement.name=="isolationWindow"
                    for isolationWindowParam in eachelement(precursorElement)
                        if isolationWindowParam.name=="cvParam"
                            parseScanCvParam!(spectrum_dict, isolationWindowParam)
                        end
                    end
                elseif precursorElement.name == "activation"
                    for activationElement in eachelement(precursorElement)
                        if activationElement.name == "cvParam"
                            parseScanCvParam!(spectrum_dict, activationElement)
                        end
                    end
                end
            end
        end
    end

    return nothing
end

function parseSpectrumElement!(
    spectrum_dict::Dict{String, String},
    spectrumElement::EzXML.Node,
    scanIndex::Int,
    skip_scan_header::Bool)
    mz_array, intensity_array = missing, missing
    for scanElement in eachelement(spectrumElement)
        if scanElement.name=="binaryDataArrayList"
            mz_array, intensity_array = parseBinaryDataList(scanElement)
        elseif scanElement.name=="cvParam"
            parseScanCvParam!(spectrum_dict, scanElement)
        elseif scanElement.name == "scanList"
            parseScanList!(spectrum_dict, scanElement)
        elseif scanElement.name == "precursorList"
            parsePrecursorList!(spectrum_dict, scanElement)
        end
    end

    return parseScanDictToScanElement(spectrum_dict,
                                scanIndex,
                                mz_array,
                                intensity_array,
                                skip_scan_header)

end

function readMzML(
    mzML_path::String,
    skip_scan_header::Bool)
    root_elements = EzXML.root(EzXML.readxml(mzML_path))
    # Create a namespace map
    ns = Dict("ms" => "http://psi.hupo.org/ms/mzml")
    # Use the namespace prefix in the XPath query
    mzML = findfirst("//ms:mzML", root_elements, ns)
    run_element = findfirst("//ms:run", mzML, ns)
    spectrum_list = findfirst("//ms:spectrumList", run_element, ns)

    spectrum_dict = Dict{String, String}(
        "ms level" => "",
        "base peak intensity" => "",
        "base peak m/z" => "",
        "total ion current" => "",
        "spectrum title" => "",
        "scan start time" => "",
        "scan window upper limit" => "",
        "scan window lower limit" =>"",
        "isolation window target m/z" => "",
        "isolation window lower offset" => "",
        "isolation window upper offset" => "",
        "beam-type collision-induced dissociation" => ""
    )

    pairedSpectra = Vector{PioneerScanElement}(undef, length(collect(eachelement(spectrum_list))))
    for (i, spectrum_element) in enumerate(collect(eachelement(spectrum_list)))
        pairedSpectra[i] = parseSpectrumElement!(
                                        spectrum_dict, 
                                        spectrum_element, 
                                        i,
                                        skip_scan_header
        )
    end

    dir, filename = splitdir(mzML_path)
    base_name = splitext(filename)[1]
    arrow_path = joinpath(dir, "$(base_name).arrow")
    Arrow.write(arrow_path, DataFrame(pairedSpectra))
end

"""
    convertMzML(mzml_dir::String; skip_scan_header::Bool=true)

Convert mzML mass spectrometry data files to Arrow IPC format.

Takes either a directory containing mzML files or a path to a single mzML file and converts them to 
Arrow format, preserving scan data including m/z arrays, intensity arrays, and scan metadata.

# Arguments
- `mzml_dir::String`: Path to either a directory containing mzML files or a path to a single mzML file
- `skip_scan_header::Bool=true`: When true, omits scan header information from the output to reduce file size

# Returns
`nothing`

# Output
Creates Arrow (.arrow) files in the same directory as the input mzML files and with the same base filename.

# Examples
```julia
# Convert all mzML files in a directory
convertMzML("path/to/mzml/files")

# Convert a single mzML file
convertMzML("path/to/single/file.mzML")

# Include scan headers in output
convertMzML("path/to/mzml/files", skip_scan_header=false)
```
# Notes

Each mzML file is converted to a corresponding Arrow IPC (.arrow) file in the same directory. This is particularly useful for Sciex data where direct .wiff/.wiff2 conversion is not supported
"""
function convertMzML(
    mzml_dir::String;
    skip_scan_header= true)

    # Clean up any old file handlers in case the program crashed
    GC.gc()
    mzml_dir = expanduser(mzml_dir)

    mzml_paths = missing
    if isdir(mzml_dir)
        mzml_paths = [fpath for fpath in readdir(mzml_dir, join=true) if endswith(fpath, ".mzML")]
    elseif isfile(mzml_dir)
        if endswith(mzml_dir, ".mzML")
            mzml_paths = [mzml_dir]
        end
    end
    if ismissing(mzml_paths)
        throw("Error $mzml_dir is not a directory containing .mzML nor a path to an .mzML file")
    elseif iszero(length(mzml_paths))
        throw("Error $mzml_dir is not a directory containing .mzML nor a path to an .mzML file")
    end

    for mzml_path in ProgressBar(mzml_paths)
        readMzML(mzml_path, skip_scan_header)
    end
    return nothing
end

#=
using Arrow, EzXML, CodecZlib, ProgressBars, Base64, DataFrames
include("src/Routines/mzmlConverter/convertMzMl.jl")
mzml_path = "C:\\Users\\n.t.wamsley\\Desktop\\SCIEX_CONVERT"
convertMzMl(mzml_path; skip_scan_header = true)
=#