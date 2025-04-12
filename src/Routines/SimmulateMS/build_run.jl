"""
    create_cycle_definition(
        cycle_config::Vector;
        ms1_duration_ms::T = T(100.0)
    ) where {T<:AbstractFloat}

Create a cycle definition from a simplified configuration.

Parameters:
- cycle_config: Vector of tuples (ms_order, center_mz, isolation_width, duration_ms)
               For MS1 scans, provide only (MS1(), duration_ms)
               For MS2 scans, provide (MS2(), center_mz, isolation_width, duration_ms)
- ms1_duration_ms: Default duration for MS1 scans (used when duration not specified)

Returns:
- Vector of ScanDefinition objects
"""
function create_cycle_definition(
    cycle_config::Vector;
    ms1_duration_ms::T = Float32(100.0)
) where {T<:AbstractFloat}
    scan_definitions = Vector{ScanDefinition{T}}()
    
    for config in cycle_config
        if first(config) isa MS1
            # MS1 scan
            if length(config) == 1
                # Use default duration
                push!(scan_definitions, ScanDefinition(MS1(), ms1_duration_ms))
            else
                # Use provided duration
                push!(scan_definitions, ScanDefinition(MS1(), T(last(config))))
            end
        elseif first(config) isa MS2
            # MS2 scan
            if length(config) == 4
                # Full definition
                push!(scan_definitions, ScanDefinition(
                    MS2(), T(config[2]), T(config[3]), T(config[4])
                ))
            else
                error("MS2 scan definition requires (MS2(), center_mz, isolation_width, duration_ms)")
            end
        else
            error("Unsupported MS order: $(first(config))")
        end
    end
    
    return scan_definitions
end

"""
    get_ms1_spectra(run::SimulatedMSRun{T}) where {T<:AbstractFloat}

Extract all MS1 spectra from a simulated MS run with their retention times.

Returns:
- Tuple of (retention_times, spectra)
"""
function get_ms1_spectra(run::SimulatedMSRun{T}) where {T<:AbstractFloat}
    ms1_indices = findall(i -> run.spectra[i].ms_order == 0x01, 1:length(run.spectra))
    return (run.rt_points[ms1_indices], run.spectra[ms1_indices])
end

"""
    get_ms2_spectra(run::SimulatedMSRun{T}; prec_mz::Union{T, Nothing}=nothing) where {T<:AbstractFloat}

Extract MS2 spectra from a simulated MS run.

Parameters:
- run: SimulatedMSRun object
- prec_mz: Optional precursor m/z to filter by

Returns:
- Tuple of (retention_times, center_mzs, spectra)
"""
function get_ms2_spectra(run::SimulatedMSRun{T}; prec_mz::Union{T, Nothing}=nothing) where {T<:AbstractFloat}
    # Find all MS2 indices
    ms2_indices = findall(i -> run.spectra[i].ms_order == 0x02, 1:length(run.spectra))
    
    # Filter by precursor m/z if specified
    if !isnothing(prec_mz)
        ms2_indices = filter(i -> abs(run.spectra[i].centerMz - prec_mz) < 0.01, ms2_indices)
    end
    
    # Extract data
    rt_points = run.rt_points[ms2_indices]
    center_mzs = [run.spectra[i].centerMz for i in ms2_indices]
    spectra = run.spectra[ms2_indices]
    
    return (rt_points, center_mzs, spectra)
end

"""
    simulate_ms_run(
        ms_sim::MSSimulation,
        iso_splines::IsotopeSplineModel,
        cycle_definition::Vector{ScanDefinition{T}},
        start_rt::T,
        end_rt::T;
        rt_points::Union{Vector{T}, Nothing} = nothing,
        n_cycles::Union{Int, Nothing} = nothing,
        rt_tolerance::T = T(0.1),
        n_ions_per_precursor::Int = 1000,
        ppm_tolerance::T = T(5.0)
    ) where {T<:AbstractFloat}

Simulate an entire MS run based on a cycle definition.

Parameters:
- ms_sim: MSSimulation object with precursors
- iso_splines: Model for isotope pattern prediction
- cycle_definition: Vector of ScanDefinition objects defining the MS cycle
- start_rt: Starting retention time
- end_rt: Ending retention time
- rt_points: Optional vector of specific RT points to simulate (overrides n_cycles)
- n_cycles: Optional number of cycles to simulate (calculated from durations if not provided)
- rt_tolerance: Retention time window width for each scan
- n_ions_per_precursor: Number of ions to simulate per precursor
- ppm_tolerance: Parts-per-million tolerance for merging peaks

Returns:
- SimulatedMSRun object containing all simulated spectra
"""
function simulate_ms_run(
    ms_sim::MSSimulation,
    iso_splines::IsotopeSplineModel,
    cycle_definition::Vector{ScanDefinition{T}},
    start_rt::T,
    end_rt::T;
    rt_points::Union{Vector{T}, Nothing} = nothing,
    n_cycles::Union{Int, Nothing} = nothing,
    rt_tolerance::T = T(0.1),
    n_ions_per_precursor::Int = 1000,
    ppm_tolerance::T = T(5.0)
) where {T<:AbstractFloat}
    # Calculate total cycle time in milliseconds
    cycle_time_ms = sum(scan.duration_ms for scan in cycle_definition)
    cycle_time_min = cycle_time_ms / 60000.0  # Convert to minutes
    
    # Generate RT points if not provided
    if isnothing(rt_points)
        # Calculate number of cycles if not provided
        if isnothing(n_cycles)
            rt_range = end_rt - start_rt
            n_cycles = ceil(Int, rt_range / cycle_time_min)
        end
        
        # Generate RT points based on scan durations
        rt_points = Vector{T}()
        for cycle in  ProgressBar(collect(range(1, n_cycles)))
            cycle_start_rt = start_rt + (cycle - 1) * cycle_time_min
            
            # Add a point for each scan in the cycle
            cumulative_time = T(0)
            for scan in cycle_definition
                # Position the RT point at the middle of the scan
                scan_rt = cycle_start_rt + (cumulative_time + scan.duration_ms/2) / 60000.0
                push!(rt_points, scan_rt)
                cumulative_time += scan.duration_ms
            end
        end
    end
    
    # Determine which scan definition applies to each RT point
    n_scans_per_cycle = length(cycle_definition)
    scan_indices = [(i-1) % n_scans_per_cycle + 1 for i in 1:length(rt_points)]
    
    # Simulate spectrum at each RT point
    spectra = Vector{SimulatedSpectrum{T}}(undef, length(rt_points))
    
    for i in ProgressBar(collect(1:length(rt_points)))
        rt = rt_points[i]
        scan_def = cycle_definition[scan_indices[i]]
        
        # Get scan parameters
        ms_order = scan_def.ms_order
        
        # Default values for MS1
        center_mz = T(0)
        isolation_width = T(0)
        
        # Use provided values for MS2
        if ms_order isa MS2
            center_mz = scan_def.center_mz
            isolation_width = scan_def.isolation_width
            
            # Skip if missing parameters for MS2
            if ismissing(center_mz) || ismissing(isolation_width)
                @warn "Missing parameters for MS2 scan at RT=$rt"
                continue
            end
        end
        
        # Simulate the spectrum
        spectra[i] = simulate_combined_spectrum(
            ms_sim,
            iso_splines,
            rt,
            rt_tolerance,
            center_mz,
            isolation_width,
            ms_order;
            n_ions_per_precursor = n_ions_per_precursor,
            ppm_tolerance = ppm_tolerance
        )
    end
    
    # Create parameter dictionary
    parameters = Dict{String, Any}(
        "start_rt" => start_rt,
        "end_rt" => end_rt,
        "cycle_time_ms" => cycle_time_ms,
        "n_cycles" => isnothing(n_cycles) ? length(rt_points) ÷ n_scans_per_cycle : n_cycles,
        "rt_tolerance" => rt_tolerance,
        "n_ions_per_precursor" => n_ions_per_precursor,
        "ppm_tolerance" => ppm_tolerance
    )
    
    # Return the simulated MS run
    return SimulatedMSRun{T}(
        rt_points,
        spectra,
        cycle_definition,
        parameters
    )
end
