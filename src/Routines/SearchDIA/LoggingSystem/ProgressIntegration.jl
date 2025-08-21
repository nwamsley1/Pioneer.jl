# ProgressIntegration.jl
# Integration with ProgressBars.jl for progress bars

using ProgressBars
using Logging

export ProgressLogger, update_progress!, finish_progress!
export with_progress_bar

"""
    ProgressLogger

Wrapper for managing progress bars with logging integration.
"""
mutable struct ProgressLogger
    progress::Union{ProgressBar, Nothing}
    name::String
    active::Bool
    suppress_logs::Bool
end

# Global progress state
const ACTIVE_PROGRESS = Ref{Union{ProgressLogger, Nothing}}(nothing)

"""
    ProgressLogger(total::Int, name::String; suppress_logs::Bool = true)

Create a new progress logger.
"""
function ProgressLogger(total::Int, name::String; suppress_logs::Bool = true)
    # Log start of operation
    @info "Starting: $(name)" progress_step=true total=total
    
    # Create progress bar using ProgressBars syntax
    progress = ProgressBar(total=total)
    
    logger = ProgressLogger(progress, name, true, suppress_logs)
    
    # Set as active if no other progress is running
    if ACTIVE_PROGRESS[] === nothing
        ACTIVE_PROGRESS[] = logger
    end
    
    return logger
end

"""
    update_progress!(logger::ProgressLogger, increment::Int = 1)

Update the progress bar.
"""
function update_progress!(logger::ProgressLogger, increment::Int = 1)
    if logger.progress !== nothing && logger.active
        for _ in 1:increment
            update(logger.progress)
        end
    end
end

"""
    finish_progress!(logger::ProgressLogger)

Complete and close the progress bar.
"""
function finish_progress!(logger::ProgressLogger)
    if logger.progress !== nothing && logger.active
        # ProgressBars doesn't have a finish method, just mark as inactive
        logger.active = false
        
        # Clear active progress if this was it
        if ACTIVE_PROGRESS[] === logger
            ACTIVE_PROGRESS[] = nothing
        end
        
        # Log completion
        @info "Completed: $(logger.name)" progress_step=true
    end
end

"""
    with_progress_bar(f, total::Int, name::String; kwargs...)

Execute function f with a progress bar.
"""
function with_progress_bar(f, total::Int, name::String; suppress_logs::Bool = true)
    logger = ProgressLogger(total, name; suppress_logs=suppress_logs)
    
    # Set context for log suppression
    if suppress_logs
        with_logger(TransformerLogger(current_logger()) do log
            if ACTIVE_PROGRESS[] !== nothing
                # Add flag to suppress logs during progress
                return merge(log, (kwargs = merge(log.kwargs, (progress_active = true,)),))
            end
            return log
        end) do
            try
                result = f(logger)
                return result
            finally
                finish_progress!(logger)
            end
        end
    else
        try
            result = f(logger)
            return result
        finally
            finish_progress!(logger)
        end
    end
end

"""
    is_progress_active() -> Bool

Check if a progress bar is currently active.
"""
function is_progress_active()
    return ACTIVE_PROGRESS[] !== nothing && ACTIVE_PROGRESS[].active
end

"""
    temporarily_hide_progress(f)

Temporarily hide the progress bar while executing function f.
"""
function temporarily_hide_progress(f)
    active = ACTIVE_PROGRESS[]
    if active !== nothing && active.active
        # Clear the line
        print("\r" * " "^80 * "\r")
        
        result = f()
        
        # With ProgressBars, we can't easily redraw, so just continue
        
        return result
    else
        return f()
    end
end

# Integration with existing ProgressBars usage

"""
    create_progress_bar(total::Int, desc::String = "Processing") -> ProgressBar

Create a progress bar that integrates with the logging system.
"""
function create_progress_bar(total::Int, desc::String = "Processing")
    @info "Starting: $(desc)" progress_step=true total=total
    
    progress = ProgressBar(total=total)
    
    # Store reference for log suppression
    if ACTIVE_PROGRESS[] === nothing
        ACTIVE_PROGRESS[] = ProgressLogger(progress, desc, true, true)
    end
    
    return progress
end

"""
    update_progress_bar!(p::ProgressBar, increment::Int = 1)

Update a progress bar created with create_progress_bar.
"""
function update_progress_bar!(p::ProgressBar, increment::Int = 1)
    for _ in 1:increment
        update(p)
    end
end

"""
    finish_progress_bar!(p::ProgressBar)

Finish a progress bar created with create_progress_bar.
"""
function finish_progress_bar!(p::ProgressBar)
    # ProgressBars doesn't have a finish method
    
    # Clear active progress if this was it
    if ACTIVE_PROGRESS[] !== nothing && ACTIVE_PROGRESS[].progress === p
        @info "Completed: $(ACTIVE_PROGRESS[].name)" progress_step=true
        ACTIVE_PROGRESS[] = nothing
    end
end