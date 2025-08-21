# WarningCapturingLogger.jl
# Custom logger that captures warnings while forwarding all messages

export WarningCapturingLogger

"""
    WarningCapturingLogger <: AbstractLogger

A logger that captures warning messages while forwarding all messages to wrapped logger(s).
"""
mutable struct WarningCapturingLogger <: AbstractLogger
    logger::AbstractLogger
    tracker::WarningTracker
end

# Implement AbstractLogger interface
Logging.min_enabled_level(logger::WarningCapturingLogger) = Logging.min_enabled_level(logger.logger)
Logging.shouldlog(logger::WarningCapturingLogger, args...) = Logging.shouldlog(logger.logger, args...)
Logging.catch_exceptions(logger::WarningCapturingLogger) = Logging.catch_exceptions(logger.logger)

function Logging.handle_message(logger::WarningCapturingLogger, level, message, _module, group, id,
                                filepath, line; kwargs...)
    # Capture warnings
    if level == Logging.Warn
        location = filepath === nothing ? nothing : "$(filepath):$(line)"
        context = Dict{Symbol, Any}(kwargs...)
        
        # Extract category if provided in kwargs
        category = get(context, :category, nothing)
        delete!(context, :category)  # Remove from context to avoid duplication
        
        add_warning!(logger.tracker, string(message); 
                    category=category, 
                    location=location, 
                    context=context)
    end
    
    # Forward to wrapped logger
    Logging.handle_message(logger.logger, level, message, _module, group, id, 
                          filepath, line; kwargs...)
end