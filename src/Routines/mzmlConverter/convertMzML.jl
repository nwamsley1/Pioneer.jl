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

const CONVERT_MZML_APP_NAME = "convertMzML"

Base.@kwdef mutable struct ConvertMzMLOptions
    mzml_path::String = ""
    output_dir::String = ""
    skip_existing::Bool = false
    skip_scan_header::Bool = true
    concurrent_files::Int = 2
    should_show_help::Bool = false
    should_show_version::Bool = false
end

function parse_convert_mzml_args(args::Vector{String})::ConvertMzMLOptions
    options = ConvertMzMLOptions()

    isempty(args) && return ConvertMzMLOptions(; should_show_help = true)

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "-o" || arg == "--output-dir"
            i += 1
            i > length(args) && throw(ArgumentError("Missing value for $arg"))
            options.output_dir = args[i]
        elseif arg == "--skip-existing"
            options.skip_existing = true
        elseif arg == "-n" || arg == "--concurrent-files"
            i += 1
            i > length(args) && throw(ArgumentError("Missing value for $arg"))
            try
                options.concurrent_files = parse(Int, args[i])
            catch
                throw(ArgumentError("Invalid value for $arg"))
            end
        elseif arg == "--skip-header"
            options.skip_scan_header = true
        elseif arg == "--include-scan-header"
            options.skip_scan_header = false
        elseif arg == "--version"
            options.should_show_version = true
        elseif arg == "-h" || arg == "--help"
            options.should_show_help = true
        elseif startswith(arg, "-")
            throw(ArgumentError("Unknown option: $arg"))
        elseif isempty(options.mzml_path)
            options.mzml_path = arg
        elseif lowercase(arg) in ("true", "false")
            # Backward compatibility for the old positional skip_header argument.
            options.skip_scan_header = parse(Bool, lowercase(arg))
        else
            throw(ArgumentError("Unexpected argument: $arg"))
        end

        i += 1
    end

    return options
end

function show_convert_mzml_help(io::IO=stdout)
    println(io, "$(CONVERT_MZML_APP_NAME) $(get_pioneer_version())")
    println(io)
    println(io, "Usage: $(CONVERT_MZML_APP_NAME) MZML_PATH [options]")
    println(io)
    println(io, "Arguments:")
    println(io, "  MZML_PATH                  Path to .mzML file or directory containing .mzML files")
    println(io)
    println(io, "Options:")
    println(io, "  -o, --output-dir <path>   Output directory for .arrow files (default: <input_dir>/arrow_out)")
    println(io, "      --skip-existing       Skip conversion when existing output appears complete")
    println(io, "  -n, --concurrent-files <n>")
    println(io, "                            Number of files to convert at the same time (default: 2)")
    println(io, "      --skip-header         Omit scan header information from output (default)")
    println(io, "      --include-scan-header Include scan header information in output")
    println(io, "      --version             Show version information")
    println(io, "  -h, --help                Show help information")
end

function format_convert_mzml_duration(elapsed_ns::Integer)::String
    elapsed_seconds = elapsed_ns / 1.0e9
    if elapsed_seconds >= 3600
        total_seconds = round(Int, elapsed_seconds)
        hours = total_seconds ÷ 3600
        minutes = (total_seconds % 3600) ÷ 60
        seconds = total_seconds % 60
        return "$(hours)h $(minutes)m $(seconds)s"
    elseif elapsed_seconds >= 60
        total_seconds = round(Int, elapsed_seconds)
        minutes = total_seconds ÷ 60
        seconds = total_seconds % 60
        return "$(minutes)m $(seconds)s"
    elseif elapsed_seconds >= 1
        return string(round(elapsed_seconds; digits = 1), "s")
    else
        return string(round(elapsed_ns / 1.0e6), "ms")
    end
end

# Entry point for PackageCompiler
function main_convertMzML(argv=ARGS)::Cint
    try
        options = parse_convert_mzml_args(String[arg for arg in argv])

        if options.should_show_version
            println("$(CONVERT_MZML_APP_NAME) $(get_pioneer_version())")
            return 0
        end

        if options.should_show_help
            show_convert_mzml_help()
            return 0
        end

        if isempty(options.mzml_path)
            println("Missing required MZML_PATH argument.")
            show_convert_mzml_help()
            return 1
        end

        convertMzML(
            options.mzml_path;
            skip_scan_header = options.skip_scan_header,
            output_dir = options.output_dir,
            skip_existing = options.skip_existing,
            concurrent_files = options.concurrent_files
        )
    catch e
        println(sprint(showerror, e))
        if e isa ArgumentError
            show_convert_mzml_help()
        end
        return 1
    end
    return 0
end

struct PioneerScanElement
    mz_array::Vector{Union{Missing, Float32}}
    intensity_array::Vector{Union{Missing, Float32}}
    scanHeader::String
    scanNumber::Int32
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
    cycle_idx::Int32
end

mutable struct MzMLCycleIndexTracker
    cycle_idx::Int32
    has_first_ms2_window::Bool
    first_ms2_center_mz::Float32
    first_ms2_isolation_width_mz::Float32
end

MzMLCycleIndexTracker() = MzMLCycleIndexTracker(Int32(1), false, 0.0f0, 0.0f0)

function get_cycle_idx!(
    tracker::MzMLCycleIndexTracker,
    ms_order::UInt8,
    center_mz,
    isolation_width_mz,
)::Int32
    ms_order == UInt8(2) || return tracker.cycle_idx

    center = _cycle_idx_value(center_mz)
    width = _cycle_idx_value(isolation_width_mz)
    if !tracker.has_first_ms2_window
        tracker.has_first_ms2_window = true
        tracker.first_ms2_center_mz = center
        tracker.first_ms2_isolation_width_mz = width
        return tracker.cycle_idx
    end

    if abs(center - tracker.first_ms2_center_mz) <= CYCLE_IDX_ISOLATION_WINDOW_TOLERANCE_MZ &&
       abs(width - tracker.first_ms2_isolation_width_mz) <= CYCLE_IDX_ISOLATION_WINDOW_TOLERANCE_MZ
        tracker.cycle_idx += Int32(1)
    end

    return tracker.cycle_idx
end

# PSI-MS activation terms we support when parsing precursor collision energy.
const CID_ACCESSION = "MS:1000133"  # collision-induced dissociation
const BEAM_TYPE_CID_ACCESSION = "MS:1000422"  # beam-type collision-induced dissociation
const COLLISION_ENERGY_ACCESSION = "MS:1000045"  # collision energy
const CENTROID_SPECTRUM_ACCESSION = "MS:1000127"  # centroid spectrum
const PROFILE_SPECTRUM_ACCESSION = "MS:1000128"  # profile spectrum

struct ProfileModeMzMLError <: Exception
    mzml_path::String
    scan_index::Int
end

function Base.showerror(io::IO, err::ProfileModeMzMLError)
    print(
        io,
        "Profile-mode spectrum detected in $(err.mzml_path) at scan $(err.scan_index). " *
        "Pioneer only supports centroided mzML data."
    )
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

function init_spectrum_dict()::Dict{String, String}
    return Dict{String, String}(
        "ms level" => "",
        "total ion current" => "",
        "collision energy" => "",
        "spectrum title" => "",
        "scan start time" => "",
        "scan window upper limit" => "",
        "scan window lower limit" => "",
        "isolation window target m/z" => "",
        "isolation window lower offset" => "",
        "isolation window upper offset" => ""
    )
end

function parse_required_cv_param(::Type{T},
    spectrum_dict::Dict{String, String},
    key::String,
    scan_number::Int64)::T where {T <: Number}

    value = strip(get(spectrum_dict, key, ""))
    isempty(value) && throw(ArgumentError("missing \"$key\" for spectrum $scan_number"))
    return parse(T, value)
end

function parse_optional_cv_param(::Type{T},
    spectrum_dict::Dict{String, String},
    key::String)::Union{Missing, T} where {T <: Number}

    value = strip(get(spectrum_dict, key, ""))
    return isempty(value) ? missing : parse(T, value)
end


function parseScanDictToScanElement(
    spectrum_dict::Dict{String, String},
    scan_number::Int64,
    mz_array::Vector{Float32},
    intensity_array::Vector{Float32},
    skip_scan_header::Bool,
    cycle_tracker::MzMLCycleIndexTracker)::PioneerScanElement

    centerMz,isolationWidthMz = missing, missing
    scanHeader = ""
    if !skip_scan_header
        scanHeader = spectrum_dict["spectrum title"]
    end
    scanNumber = Int32(scan_number)
    retentionTime = parse_required_cv_param(Float32, spectrum_dict, "scan start time", scan_number)
    lowMz = parse_required_cv_param(Float32, spectrum_dict, "scan window lower limit", scan_number)
    highMz = parse_required_cv_param(Float32, spectrum_dict, "scan window upper limit", scan_number)
    msOrder = parse_required_cv_param(UInt8, spectrum_dict, "ms level", scan_number)
    if msOrder > 1
        targetMz = parse_optional_cv_param(Float32, spectrum_dict, "isolation window target m/z")
        lowerOffset = parse_optional_cv_param(Float32, spectrum_dict, "isolation window lower offset")
        upperOffset = parse_optional_cv_param(Float32, spectrum_dict, "isolation window upper offset")
        if !ismissing(targetMz) && !ismissing(lowerOffset) && !ismissing(upperOffset)
            isolationWidthMz = lowerOffset + upperOffset
            centerMz = targetMz + (upperOffset - lowerOffset)/2.0f0
        end
    end
    TIC = coalesce(
        parse_optional_cv_param(Float32, spectrum_dict, "total ion current"),
        sum(intensity_array; init = zero(Float32))
    )
    collisionEnergyField = parse_optional_cv_param(
        Float32,
        spectrum_dict,
        "collision energy"
    )
    cycle_idx = get_cycle_idx!(cycle_tracker, msOrder, centerMz, isolationWidthMz)
    return PioneerScanElement(
        allowmissing(mz_array), 
        allowmissing(intensity_array),
        scanHeader,
        scanNumber,
        zero(Int32),
        retentionTime,
        lowMz,
        highMz,
        TIC,
        centerMz,
        isolationWidthMz,
        collisionEnergyField,
        zero(Float32),
        msOrder,
        cycle_idx
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

function parseActivation!(
    spectrum_dict::Dict{String, String},
    activation::EzXML.Node)

    has_supported_activation = false
    collision_energy = ""

    for activationElement in eachelement(activation)
        if activationElement.name != "cvParam"
            continue
        end

        accession = activationElement["accession"]
        if accession == CID_ACCESSION || accession == BEAM_TYPE_CID_ACCESSION
            has_supported_activation = true
        elseif accession == COLLISION_ENERGY_ACCESSION
            collision_energy = activationElement["value"]
        end
    end

    if has_supported_activation
        spectrum_dict["collision energy"] = collision_energy
    end

    return nothing
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
                    parseActivation!(spectrum_dict, precursorElement)
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
    skip_scan_header::Bool,
    cycle_tracker::MzMLCycleIndexTracker)
    mz_array, intensity_array = missing, missing
    for scanElement in eachelement(spectrumElement)
        if scanElement.name=="binaryDataArrayList"
            mz_array, intensity_array = parseBinaryDataList(scanElement)
        elseif scanElement.name=="cvParam"
            accession = scanElement["accession"]
            if accession == PROFILE_SPECTRUM_ACCESSION
                throw(ProfileModeMzMLError("", scanIndex))
            elseif accession == CENTROID_SPECTRUM_ACCESSION
                continue
            end
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
                                skip_scan_header,
                                cycle_tracker)

end

function get_mzml_input_dir(mzml_path::String)::String
    mzml_path = expanduser(mzml_path)

    if isdir(mzml_path)
        return abspath(mzml_path)
    elseif isfile(mzml_path)
        return dirname(abspath(mzml_path))
    end

    throw(ArgumentError("File or Directory does not exist: $mzml_path"))
end

function get_mzml_paths(mzml_path::String)::Vector{String}
    mzml_path = expanduser(mzml_path)

    if isfile(mzml_path)
        return endswith(lowercase(mzml_path), ".mzml") ? [abspath(mzml_path)] : String[]
    elseif isdir(mzml_path)
        return sort([
            abspath(fpath) for fpath in readdir(mzml_path, join = true)
            if endswith(lowercase(fpath), ".mzml")
        ])
    end

    return String[]
end

function build_convert_mzml_output_dir(input_dir::String, requested_output_dir::String)::String
    output_dir = isempty(strip(requested_output_dir)) ?
                 joinpath(input_dir, "arrow_out") :
                 abspath(expanduser(requested_output_dir))

    if isfile(output_dir)
        throw(ArgumentError("Output path points to an existing file: $output_dir"))
    end

    mkpath(output_dir)
    return output_dir
end

function get_convert_mzml_output_path(output_dir::String, mzml_path::String)::String
    return joinpath(output_dir, "$(splitext(basename(mzml_path))[1]).arrow")
end

function get_convert_mzml_output_paths(output_dir::String, mzml_paths::Vector{String})::Vector{String}
    return [get_convert_mzml_output_path(output_dir, mzml_path) for mzml_path in mzml_paths]
end

function get_mzml_last_scan_number(mzml_path::String)::Int
    root_elements = EzXML.root(EzXML.readxml(mzml_path))
    ns = Dict("ms" => "http://psi.hupo.org/ms/mzml")
    mzML = findfirst("//ms:mzML", root_elements, ns)
    run_element = findfirst("//ms:run", mzML, ns)
    spectrum_list = findfirst("//ms:spectrumList", run_element, ns)

    count_attr = try
        strip(spectrum_list["count"])
    catch
        ""
    end
    if !isempty(count_attr)
        return parse(Int, count_attr)
    end

    return length(collect(eachelement(spectrum_list)))
end

function get_convert_mzml_output_last_scan_number(output_path::String)::Union{Nothing, Int}
    table = Arrow.Table(output_path)
    scan_numbers = table[:scanNumber]
    isempty(scan_numbers) && return nothing
    return Int(scan_numbers[end])
end

function has_complete_existing_output(mzml_path::String, output_path::String)::Bool
    try
        return get_convert_mzml_output_last_scan_number(output_path) == get_mzml_last_scan_number(mzml_path)
    catch
        return false
    end
end

function partition_convert_mzml_queue(
    mzml_paths::Vector{String},
    output_paths::Vector{String};
    skip_existing::Bool)::NamedTuple

    files_to_convert = Int[]
    skipped_complete_files = 0
    reconverted_incomplete_files = 0
    missing_output_files = 0

    for i in eachindex(mzml_paths, output_paths)
        if !skip_existing
            push!(files_to_convert, i)
            continue
        end

        if !isfile(output_paths[i])
            missing_output_files += 1
            push!(files_to_convert, i)
            continue
        end

        if has_complete_existing_output(mzml_paths[i], output_paths[i])
            skipped_complete_files += 1
            continue
        end

        reconverted_incomplete_files += 1
        push!(files_to_convert, i)
    end

    return (
        files_to_convert = files_to_convert,
        skipped_complete_files = skipped_complete_files,
        reconverted_incomplete_files = reconverted_incomplete_files,
        missing_output_files = missing_output_files,
    )
end

function readMzML(
    mzML_path::String,
    output_path::String,
    skip_scan_header::Bool)
    root_elements = EzXML.root(EzXML.readxml(mzML_path))
    # Create a namespace map
    ns = Dict("ms" => "http://psi.hupo.org/ms/mzml")
    # Use the namespace prefix in the XPath query
    mzML = findfirst("//ms:mzML", root_elements, ns)
    run_element = findfirst("//ms:run", mzML, ns)
    spectrum_list = findfirst("//ms:spectrumList", run_element, ns)

    spectrum_elements = collect(eachelement(spectrum_list))
    pairedSpectra = Vector{PioneerScanElement}(undef, length(spectrum_elements))
    cycle_tracker = MzMLCycleIndexTracker()
    for (i, spectrum_element) in enumerate(spectrum_elements)
        try
            pairedSpectra[i] = parseSpectrumElement!(
                                            init_spectrum_dict(),
                                            spectrum_element, 
                                            i,
                                            skip_scan_header,
                                            cycle_tracker
            )
        catch err
            if err isa ProfileModeMzMLError
                throw(ProfileModeMzMLError(mzML_path, err.scan_index))
            end
            rethrow()
        end
    end

    Arrow.write(output_path, DataFrame(pairedSpectra))
end

function process_mzml_file(
    input_path::String,
    output_path::String,
    skip_scan_header::Bool)

    file_label = splitext(basename(input_path))[1]
    println("Starting Conversion For: $file_label")

    start_ns = time_ns()
    try
        readMzML(input_path, output_path, skip_scan_header)
    catch err
        if err isa ProfileModeMzMLError
            println(
                "Warning: skipping $file_label - profile-mode spectra detected. " *
                "Pioneer only handles centroided mzML data."
            )
            return false
        end
        rethrow()
    end
    elapsed_ns = time_ns() - start_ns

    println("Execution Time: $(round(Int, elapsed_ns / 1.0e6)) ms for $file_label")
    return true
end

"""
    convertMzML(mzml_path::String; skip_scan_header::Bool=true, output_dir::String="", skip_existing::Bool=false, concurrent_files::Int=2)

Convert mzML mass spectrometry data files to Arrow IPC format.

Takes either a directory containing mzML files or a path to a single mzML file and converts them to
Arrow format, preserving scan data including m/z arrays, intensity arrays, and scan metadata.

# Arguments
- `mzml_path::String`: Path to either a directory containing mzML files or a path to a single mzML file
- `skip_scan_header::Bool=true`: When true, omits scan header information from the output to reduce file size
- `output_dir::String=""`: Output directory for `.arrow` files. Defaults to `<input_dir>/arrow_out`
- `skip_existing::Bool=false`: When true, skip files whose existing `.arrow` output appears complete
- `concurrent_files::Int=2`: Number of mzML files to convert at the same time

# Returns
`nothing`

# Output
Creates Arrow (.arrow) files in the requested output directory with the same base filename as the input mzML files.

# Requirements
Only centroided mzML is supported. Files containing profile-mode spectra are skipped during conversion.

# Examples
```julia
# Convert all mzML files in a directory
convertMzML("path/to/mzml/files")

# Convert a single mzML file
convertMzML("path/to/single/file.mzML")

# Convert to a custom output directory
convertMzML("path/to/mzml/files", output_dir="path/to/arrow_out")

# Skip files that already have complete Arrow outputs
convertMzML("path/to/mzml/files", skip_existing=true)

# Include scan headers in output
convertMzML("path/to/mzml/files", skip_scan_header=false)
```
# Notes

Each mzML file is converted to a corresponding Arrow IPC (.arrow) file in the output directory. This is particularly useful for Sciex data where direct .wiff/.wiff2 conversion is not supported
"""
function convertMzML(
    mzml_path::String;
    skip_scan_header::Bool = true,
    output_dir::String = "",
    skip_existing::Bool = false,
    concurrent_files::Int = 2)

    # Clean up any old file handlers in case the program crashed
    GC.gc()

    input_dir = get_mzml_input_dir(mzml_path)
    mzml_paths = get_mzml_paths(mzml_path)

    output_dir = build_convert_mzml_output_dir(input_dir, output_dir)
    output_paths = get_convert_mzml_output_paths(output_dir, mzml_paths)
    queue = partition_convert_mzml_queue(mzml_paths, output_paths; skip_existing = skip_existing)

    input_mode = isdir(expanduser(mzml_path)) ? "directory" : "file"
    concurrent_files = max(1, concurrent_files)
    total_start_ns = time_ns()

    println("$(CONVERT_MZML_APP_NAME) $(get_pioneer_version())")
    println("==================================================")
    println("Config: concurrent-files=$concurrent_files  skip-header=$(lowercase(string(skip_scan_header)))")
    println("Config: output=$output_dir")
    println("Config: skip-existing=$(lowercase(string(skip_existing)))")
    println()
    println("Input : $input_mode $(expanduser(mzml_path))")
    println("Queue : discovered=$(length(mzml_paths))  convert=$(length(queue.files_to_convert))")
    if skip_existing
        println(
            "Queue : skipped-complete=$(queue.skipped_complete_files)  " *
            "reconvert-incomplete=$(queue.reconverted_incomplete_files)  " *
            "missing-output=$(queue.missing_output_files)"
        )
    end
    println("==================================================")

    if isempty(mzml_paths)
        println("No .mzML files found to process")
        println("Total conversion time: $(format_convert_mzml_duration(time_ns() - total_start_ns))")
        return nothing
    end

    if isempty(queue.files_to_convert)
        println("No files to convert.")
        println("Total conversion time: $(format_convert_mzml_duration(time_ns() - total_start_ns))")
        return nothing
    end

    if concurrent_files == 1 || length(queue.files_to_convert) == 1
        for file_index in queue.files_to_convert
            process_mzml_file(mzml_paths[file_index], output_paths[file_index], skip_scan_header)
        end
    else
        sem = Base.Semaphore(concurrent_files)
        @sync for file_index in queue.files_to_convert
            Base.acquire(sem)
            Threads.@spawn begin
                try
                    process_mzml_file(
                        mzml_paths[file_index],
                        output_paths[file_index],
                        skip_scan_header
                    )
                finally
                    Base.release(sem)
                end
            end
        end
    end

    println("Total conversion time: $(format_convert_mzml_duration(time_ns() - total_start_ns))")
    return nothing
end

#=
using Arrow, EzXML, CodecZlib, ProgressBars, Base64, DataFrames
include("src/Routines/mzmlConverter/convertMzMl.jl")
mzml_path = "C:\\Users\\n.t.wamsley\\Desktop\\SCIEX_CONVERT"
convertMzMl(mzml_path; skip_scan_header = true)
=#
