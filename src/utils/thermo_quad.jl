struct TQuad{T<:AbstractFloat}
    h::T
    w::T
    l::T
    r::T
    x_l::T
    x_r::T
    σ_l::T
    σ_r::T
    function TQuad(h::T, w::T, l::T, r::T, c::T) where{T<:AbstractFloat}
        x_l, x_r = c - w/2, c + w/2
        σ_l, σ_r = l/sqrt(2*log(2)), r/sqrt(2*log(2))
        new{T}(
            h, w, l, r,
            x_l, x_r, 
            1/(2*σ_l^2), 1/(2*σ_r^2)
        )
    end
end

function (qtf::TQuad)(x::T) where {T<:AbstractFloat}
    if x < qtf.x_l
        xdiff = -one(T)*(x - qtf.x_l)^2
        return qtf.h*exp(xdiff/qtf.σ_l)
    elseif x > qtf.x_r
        xdiff = -one(T)*(x - qtf.x_r)^2
        return qtf.h*exp(xdiff/qtf.σ_r)       
    else
        return qtf.h
    end
end

test_tq = TQuad(1.0, 2.0, 0.8, 0.8, 0.0)

plot_bins = LinRange(-3, 3, 100)
plot!(
    plot_bins,
    test_tq.(plot_bins)
)

function solve_tquad_log(yt::Vector{Float32}, x0::Vector{Float32}, x1::Vector{Float32}, 
                                    h::Float32, l::Float32, r::Float32, c::Float32, w::Float32)
    m = length(yt)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, b)
    @variable(model, residuals[1:m])
    @constraint(model, residuals == log2.(TQuad(h, w, l, r, c).(x1)) .- log2.(TQuad(h, w, l, r, c).(x0)) .- yt) 
    @constraint(model, h >= 0)
    @constraint(model, l >= 0)
    @constraint(model, c >= 0)
    @constraint(model, w >= 0)
    @constraint(model, r >= 0)
    @objective(model, Min, sum(residuals.^2))
    optimize!(model)
    return TQuad(JuMP.value.(h), JuMP.value.(w), JuMP.value.(l), JuMP.value.(r), JuMP.value.(c))
end


function solve_tquad_log(yt::Vector{Float32}, x0::Vector{Float32}, x1::Vector{Float32}, 
                        h::Float32, l::Float32, r::Float32, c::Float32, w::Float32)
    m = length(yt)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, 0<=h<=1)
    @variable(model, l>=0)
    @variable(model, r>=0)
    @variable(model, c>=0)
    @variable(model, w>=0)
    @variable(model, residuals[1:m])
    @variable(model, tquad0[1:m])
    @variable(model, tquad1[1:m])
    @constraint(model, residuals = (log2.(tquad0) - log2.(tquad1)) - yt)
    #@constraint(model, residuals == (TQuad(JuMP.value.(h), JuMP.value.(w), JuMP.value.(l), JuMP.value.(r), JuMP.value.(c)).(x1)
    #)./(TQuad(JuMP.value.(h), JuMP.value.(w), JuMP.value.(l), JuMP.value.(r), JuMP.value.(c)).(x0)) .- exp2(yt))

    @objective(model, Min, sum(residuals.^2))
    optimize!(model)
    return TQuad(JuMP.value.(h), JuMP.value.(w), JuMP.value.(l), JuMP.value.(r), JuMP.value.(c))
end


test_tquad_solved = solve_tquad_log(yt, x0, x1, 1.0f0, 2.0f0, 0.8f0, 0.8f0, 0.0f0)

qtf = QuadTransmission(0.0f0, 5.0f0)
plot(
    plot_bins,
    qtf.(0.0f0, 1.0f0, Float32.(plot_bins)))
plot!(
    plot_bins,
    test_tquad_solved.(plot_bins)
)

function fit_tquad(x0::Vector{Float32}, x1::Vector{Float32}, y::Vector{Float32};
                  initial_h::Float32=1.0f0, initial_w::Float32=1.0f0,
                  initial_l::Float32=1.0f0, initial_r::Float32=1.0f0,
                  initial_c::Float32=0.0f0)
    
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    
    # Variables to optimize
    @variable(model, 1 >= h >= 0)    # height must be positive
    @variable(model, w >= 0)    # width must be positive
    @variable(model, l >= 0)    # left spread must be positive
    @variable(model, r >= 0)    # right spread must be positive
    @variable(model, c)         # center can be anywhere
    
    # Set initial values
    set_start_value(h, initial_h)
    set_start_value(w, initial_w)
    set_start_value(l, initial_l)
    set_start_value(r, initial_r)
    set_start_value(c, initial_c)
    
    # Residuals for each data point
    @variable(model, residuals[1:length(x)])
    
    # Define the TQuad function within the constraints
    #Can always assume x1[i] > x0[i]
    expr = for i in 1:length(x) 
            op_ifelse(
            op_and(op_strictly_less_than(x0[i], c - w/2), op_strictly_less_than(x1[i], c - w/2)),
                log2(h * exp(-(x0[i] - (c - w/2))^2 / (2 * (l/sqrt(2*log(2)))^2))
                ) - log2(h * exp(-(x1[i] - (c - w/2))^2 / (2 * (l/sqrt(2*log(2)))^2))) - y[i],
            op_ifelse(
                op_and(op_strictly_less_than(x0[i], c - w/2), op_strictly_greater_than(x1[i], c + w/2)),
                    log2(h * exp(-(x0[i] - (c - w/2))^2 / (2 * (l/sqrt(2*log(2)))^2))
                    ) - log2(h * exp(-(x1[i] - (c + w/2))^2 / (2 * (r/sqrt(2*log(2)))^2))) - y[i],
            op_ifelse(
                op_strictly_less_than(x0[i], c - w/2), #Implies x1[i] > c - w/2 and x1[i] < c + w/2
                    log2(h * exp(-(x0[i] - (c - w/2))^2 / (2 * (l/sqrt(2*log(2)))^2))
                    ) - log2(h) - y[i],         
            op_ifelse(
                op_and(op_strictly_less_than(x0[i], c + w/2), op_strictly_less_than(x1[i], c + w/2)),
                    log2(h) - log2(h),      
            op_ifelse(
                op_strictly_less_than(x0[i], c + w/2), #Implies x1[i] > c + w/2
                    log2(h
                    ) - log2(h * exp(-(x1[i] - (c + w/2))^2 / (2 * (r/sqrt(2*log(2)))^2))),
                log2(h * exp(-(x0[i] - (c + w/2))^2 / (2 * (r/sqrt(2*log(2)))^2))
                ) - log2(h * exp(-(x1[i] - (c + w/2))^2 / (2 * (r/sqrt(2*log(2)))^2)))
            )))),
            #op_ifelse(
            #    op_strictly_greater_than(x0[i], c + w/2), #Implies x1[i] > c + w/2
            log2(h * exp(-(x0[i] - (c + w/2))^2 / (2 * (r/sqrt(2*log(2)))^2))
            ) - log2(h * exp(-(x1[i] - (c + w/2))^2 / (2 * (r/sqrt(2*log(2)))^2)))
            )
    end

    @constraint(model, 
    residuals == expr .- yt)

    #=
    for i in 1:length(x)
        @NLconstraint(model, residuals[i] == 
            op_ifelse(
                op_and(op_strictly_less_than(x0[i], c - w/2), op_strictly_less_than(x1[i], c - w/2)),
                    log2(h * exp(-(x0[i] - (c - w/2))^2 / (2 * (l/sqrt(2*log(2)))^2))
                    ) - log2(h * exp(-(x1[i] - (c - w/2))^2 / (2 * (l/sqrt(2*log(2)))^2))) - y[i],
            op_ifelse(
                op_and(op_strictly_less_than(x0[i], c - w/2), op_strictly_greater_than(x1[i], c + w/2)),
                    log2(h * exp(-(x0[i] - (c - w/2))^2 / (2 * (l/sqrt(2*log(2)))^2))
                    ) - log2(h * exp(-(x1[i] - (c + w/2))^2 / (2 * (r/sqrt(2*log(2)))^2))) - y[i],
            op_ifelse(
                op_strictly_less_than(x0[i], c - w/2), #Implies x1[i] > c - w/2 and x1[i] < c + w/2
                    log2(h * exp(-(x0[i] - (c - w/2))^2 / (2 * (l/sqrt(2*log(2)))^2))
                    ) - log2(h) - y[i],         
            op_ifelse(
                op_and(op_strictly_less_than(x0[i], c + w/2), op_strictly_less_than(x1[i], c + w/2)),
                    log2(h) - log2(h) - y[i],      
            op_ifelse(
                op_strictly_less_than(x0[i], c + w/2), #Implies x1[i] > c + w/2
                    log2(h
                    ) - log2(h * exp(-(x1[i] - (c + w/2))^2 / (2 * (r/sqrt(2*log(2)))^2))) - y[i],
                log2(h * exp(-(x0[i] - (c + w/2))^2 / (2 * (r/sqrt(2*log(2)))^2))
                ) - log2(h * exp(-(x1[i] - (c + w/2))^2 / (2 * (r/sqrt(2*log(2)))^2))) - y[i]
            ))))),
            #op_ifelse(
            #    op_strictly_greater_than(x0[i], c + w/2), #Implies x1[i] > c + w/2
            log2(h * exp(-(x0[i] - (c + w/2))^2 / (2 * (r/sqrt(2*log(2)))^2))
            ) - log2(h * exp(-(x1[i] - (c + w/2))^2 / (2 * (r/sqrt(2*log(2)))^2))) - y[i]
            )
    end
    =#
    # Minimize sum of squared residuals
    @objective(model, Min, sum(residuals.^2))
    # Optimize
    optimize!(model)
    
    # Return optimized parameters
    return (
        h = value(h),
        w = value(w),
        l = value(l),
        r = value(r),
        c = value(c)
    )
end

out_tquad_params = fit_tquad(x0, x1, yt)

f(x) = x^2
∇f(x) = 2x
∇²f(x) = 2
model = Model();
@operator(model, op_f, 1, f, ∇f, ∇²f)
@variable(model, x)
@objective(model, Min, op_f(x))