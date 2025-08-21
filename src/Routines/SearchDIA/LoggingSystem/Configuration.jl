# Configuration.jl
# Configuration management for the logging system

export LoggingConfig, load_config, save_config, merge_config

"""
    LoggingConfig

Configuration for the Pioneer logging system.
"""
Base.@kwdef struct LoggingConfig
    # Console logger settings
    console_level::Symbol = :normal  # :silent, :minimal, :normal, :verbose, :debug
    
    # Simplified logger settings
    simple_level::Symbol = :info     # What to include in simple log
    
    # Full debug logger settings
    full_level::Symbol = :trace      # Capture everything
    
    # Feature flags
    enable_progress::Bool = true     # Progress bar integration
    enable_warnings::Bool = true     # Warning tracking
    
    # File paths (set during initialization)
    simple_log_path::String = ""
    full_log_path::String = ""
    
    # Rotation settings
    rotation_size_mb::Float64 = 100.0  # Rotate when log exceeds this size
    
    # Performance settings
    max_warnings::Int = 10000        # Maximum warnings to track
    buffer_size::Int = 4096          # IO buffer size
    
    # Custom level definitions
    user_info_level::Int = Int(Logging.Info) - 100
    debug_l1_level::Int = Int(Logging.Debug)
    debug_l2_level::Int = Int(Logging.Debug) - 100
    debug_l3_level::Int = Int(Logging.Debug) - 200
    trace_level::Int = Int(Logging.Debug) - 1000
end

"""
    load_config(path::String) -> LoggingConfig

Load logging configuration from a JSON file.
"""
function load_config(path::String)
    if !isfile(path)
        return LoggingConfig()  # Return defaults
    end
    
    config_dict = JSON.parsefile(path)
    
    # Convert to LoggingConfig, handling missing fields
    return LoggingConfig(;
        console_level = Symbol(get(config_dict, :console_level, "normal")),
        simple_level = Symbol(get(config_dict, :simple_level, "info")),
        full_level = Symbol(get(config_dict, :full_level, "trace")),
        enable_progress = get(config_dict, :enable_progress, true),
        enable_warnings = get(config_dict, :enable_warnings, true),
        rotation_size_mb = get(config_dict, :rotation_size_mb, 100.0),
        max_warnings = get(config_dict, :max_warnings, 10000),
        buffer_size = get(config_dict, :buffer_size, 4096)
    )
end

"""
    save_config(config::LoggingConfig, path::String)

Save logging configuration to a JSON file.
"""
function save_config(config::LoggingConfig, path::String)
    config_dict = Dict(
        :console_level => string(config.console_level),
        :simple_level => string(config.simple_level),
        :full_level => string(config.full_level),
        :enable_progress => config.enable_progress,
        :enable_warnings => config.enable_warnings,
        :rotation_size_mb => config.rotation_size_mb,
        :max_warnings => config.max_warnings,
        :buffer_size => config.buffer_size
    )
    
    open(path, "w") do io
        JSON.print(io, config_dict, 2)  # 2 is for indent
    end
end

"""
    merge_config(base::LoggingConfig, overrides::Dict) -> LoggingConfig

Merge configuration with overrides from a dictionary.
"""
function merge_config(base::LoggingConfig, overrides::Dict)
    fields = Dict{Symbol, Any}()
    
    # Get all field values from base
    for field in fieldnames(LoggingConfig)
        fields[field] = getfield(base, field)
    end
    
    # Apply overrides
    for (key, value) in overrides
        if haskey(fields, Symbol(key))
            fields[Symbol(key)] = value
        end
    end
    
    return LoggingConfig(; fields...)
end

"""
    get_config_from_env() -> Dict

Read logging configuration from environment variables.
"""
function get_config_from_env()
    config = Dict{Symbol, Any}()
    
    # Check for environment variables
    if haskey(ENV, "PIONEER_LOG_LEVEL")
        config[:console_level] = Symbol(ENV["PIONEER_LOG_LEVEL"])
    end
    
    if haskey(ENV, "PIONEER_LOG_PROGRESS")
        config[:enable_progress] = parse(Bool, ENV["PIONEER_LOG_PROGRESS"])
    end
    
    if haskey(ENV, "PIONEER_LOG_WARNINGS")
        config[:enable_warnings] = parse(Bool, ENV["PIONEER_LOG_WARNINGS"])
    end
    
    if haskey(ENV, "PIONEER_LOG_ROTATION_MB")
        config[:rotation_size_mb] = parse(Float64, ENV["PIONEER_LOG_ROTATION_MB"])
    end
    
    return config
end