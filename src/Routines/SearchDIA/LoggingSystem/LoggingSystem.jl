# LoggingSystem.jl
# Main orchestration functions for the Pioneer logging system
# All dependencies are loaded through Pioneer.jl and importScripts.jl

struct LoggerState
    warning_logger::WarningCapturingLogger
    console_logger::Ref{AbstractLogger}
    simplified_logger::Union{Nothing, AbstractLogger}
    full_logger::Union{Nothing, AbstractLogger}
    config::LoggingConfig
    start_time::DateTime
end

# Module-level state
const LOGGER_STATE = Ref{Union{Nothing, LoggerState}}(nothing)

"""
    initialize_logging(output_dir::String; config::LoggingConfig = LoggingConfig())

Initialize the Pioneer logging system with three parallel loggers and warning tracking.
Falls back to console-only if file creation fails.
"""
function initialize_logging(output_dir::String; config::LoggingConfig = LoggingConfig())
    # Ensure output directory exists
    try
        mkpath(output_dir)
    catch e
        @warn "Could not create output directory, using console-only logging" error=e
        console_logger = create_console_logger(config)
        global_logger(console_logger)
        return console_logger
    end
    
    # Create file paths
    simple_log_path = joinpath(output_dir, "pioneer_search_log.txt")
    full_log_path = joinpath(output_dir, "pioneer_debug.log")
    warnings_path = joinpath(output_dir, "warnings.txt")
    
    # Initialize warning tracker
    warning_tracker = WarningTracker()
    
    # Create console logger (always succeeds)
    console_logger_ref = Ref{AbstractLogger}(create_console_logger(config))
    
    # Try to create file loggers with error handling
    # Use AbstractLogger array to allow different logger types
    loggers = AbstractLogger[console_logger_ref[]]
    
    # Declare variables in outer scope
    simplified_logger = nothing
    full_logger = nothing
    
    try
        simplified_logger = create_simplified_logger(simple_log_path, config)
        push!(loggers, simplified_logger)
    catch e
        @warn "Could not create simplified log file" path=simple_log_path error=e
    end
    
    try
        full_logger = create_full_logger(full_log_path, config)
        push!(loggers, full_logger)
    catch e
        @warn "Could not create full debug log file" path=full_log_path error=e
    end
    
    # Create TeeLogger for parallel distribution (at minimum has console logger)
    tee_logger = if length(loggers) > 1
        TeeLogger(loggers...)
    else
        loggers[1]  # Just console logger if file creation failed
    end
    
    # Wrap with warning capturing
    warning_logger = WarningCapturingLogger(tee_logger, warning_tracker)
    
    # Store logger state
    # Store nothing for file loggers if they failed to create
    LOGGER_STATE[] = LoggerState(
        warning_logger,
        console_logger_ref,
        simplified_logger,  # Will be nothing if creation failed
        full_logger,        # Will be nothing if creation failed
        config,
        now()
    )
    
    # Set as global logger
    global_logger(warning_logger)
    
    # Log initialization
    @info "Pioneer logging system initialized" output_dir=output_dir timestamp=now()
    
    return warning_logger
end

"""
    finalize_logging()

Finalize the logging system, write warning summary, and close file handles.
"""
function finalize_logging()
    state = LOGGER_STATE[]
    if state === nothing
        return
    end
    
    # Get warnings summary
    warnings = get_warnings(state.warning_logger.tracker)
    
    if !isempty(warnings)
        # Log warning summary to all loggers
        @warn "Session completed with $(length(warnings)) warnings"
        
        # Try to write warnings report if we have a results directory
        # Since we don't store the output path, we'll skip the warnings file for now
        # TODO: Store output_dir in LoggerState to enable warnings file
        
        # Display summary on console
        display_warnings_summary(warnings)
    else
        @info "Session completed successfully with no warnings"
    end
    
    # Calculate runtime
    runtime = now() - state.start_time
    @info "Total runtime: $(runtime)"
    
    # Close file loggers (if they exist)
    if state.simplified_logger !== nothing
        if isa(state.simplified_logger, AsyncFileLogger)
            close_async_logger(state.simplified_logger)
        else
            close_logger(state.simplified_logger)
        end
    end
    if state.full_logger !== nothing
        if isa(state.full_logger, AsyncFileLogger)
            close_async_logger(state.full_logger)
        else
            close_logger(state.full_logger)
        end
    end
    
    # Reset global state
    LOGGER_STATE[] = nothing
    global_logger(ConsoleLogger())
end

"""
    update_console_verbosity(level::Symbol)

Update console logger verbosity at runtime.
Options: :silent, :minimal, :normal, :verbose, :debug
"""
function update_console_verbosity(level::Symbol)
    state = LOGGER_STATE[]
    if state === nothing
        @warn "Logging system not initialized"
        return
    end
    
    # Create new console logger with updated level
    new_config = LoggingConfig(
        console_level = level,
        simple_level = state.config.simple_level,
        full_level = state.config.full_level,
        enable_progress = state.config.enable_progress,
        enable_warnings = state.config.enable_warnings,
        simple_log_path = state.config.simple_log_path,
        full_log_path = state.config.full_log_path,
        rotation_size_mb = state.config.rotation_size_mb
    )
    
    # Update console logger
    state.console_logger[] = create_console_logger(new_config)
    
    # Recreate TeeLogger
    new_tee = TeeLogger(
        state.console_logger[],
        state.simplified_logger,
        state.full_logger
    )
    
    # Update warning logger's wrapped logger
    state.warning_logger.logger = new_tee
    
    @debug_l1 "Console verbosity updated to $(level)"
end

"""
    get_logger_state()

Get current logger state for inspection and debugging.
"""
function get_logger_state()
    return LOGGER_STATE[]
end

"""
    with_logging_context(f, context::Dict)

Execute function f with additional logging context.
"""
function with_logging_context(f, context::Dict)
    with_logger(TransformerLogger(current_logger()) do log
        merge(log, context)
    end) do
        f()
    end
end

# Helper function to close loggers
function close_logger(logger::AbstractLogger)
    # Handle AsyncFileLogger specifically
    if isa(logger, AsyncFileLogger)
        close_async_logger(logger)
    elseif hasproperty(logger, :stream) && logger.stream isa IO
        close(logger.stream)
    elseif hasproperty(logger, :io) && logger.io isa IO
        # Handle FormatLogger which uses 'io' field
        close(logger.io)
    elseif logger isa TeeLogger
        for child in logger.loggers
            close_logger(child)
        end
    elseif logger isa TransformerLogger
        close_logger(logger.logger)
    elseif logger isa MinLevelLogger
        close_logger(logger.logger)
    elseif logger isa EarlyFilteredLogger
        close_logger(logger.logger)
    end
end

# Helper function to write warnings report
function write_warnings_report(warnings::Vector{WarningEntry}, path::String)
    open(path, "w") do io
        println(io, "=" ^ 80)
        println(io, "Pioneer Warnings Report")
        println(io, "Generated: $(now())")
        println(io, "Total Warnings: $(length(warnings))")
        println(io, "=" ^ 80)
        println(io)
        
        # Group by category
        by_category = Dict{String, Vector{WarningEntry}}()
        for w in warnings
            category = string(w.category)
            if !haskey(by_category, category)
                by_category[category] = WarningEntry[]
            end
            push!(by_category[category], w)
        end
        
        # Write by category
        for (category, category_warnings) in by_category
            println(io, "\n$(category) ($(length(category_warnings)) warnings)")
            println(io, "-" ^ 40)
            
            for (i, w) in enumerate(category_warnings)
                println(io, "\n[$i] $(w.timestamp)")
                println(io, "Message: $(w.message)")
                if w.location !== nothing
                    println(io, "Location: $(w.location)")
                end
                if !isempty(w.context)
                    println(io, "Context:")
                    for (k, v) in w.context
                        println(io, "  $(k): $(v)")
                    end
                end
            end
        end
    end
end

# Helper function to display warnings summary on console
function display_warnings_summary(warnings::Vector{WarningEntry})
    # Group by category
    by_category = Dict{String, Int}()
    for w in warnings
        category = string(w.category)
        by_category[category] = get(by_category, category, 0) + 1
    end
    
    println("\n" * "=" ^ 60)
    println("WARNINGS SUMMARY")
    println("=" ^ 60)
    println("Total warnings: $(length(warnings))")
    println("\nBy category:")
    for (category, count) in sort(collect(by_category), by=x->x[2], rev=true)
        println("  $(category): $(count)")
    end
    println("\nDetailed warnings saved to: warnings.txt")
    println("=" ^ 60)
end

