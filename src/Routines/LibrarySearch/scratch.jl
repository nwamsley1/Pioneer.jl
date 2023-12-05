abstract type Params end

struct HuberParams{T<:AbstractFloat} <: Params
    σ::T
    tᵣ::T
    τ::T
    H::T
end

import Base.-
function -(a::HuberParams{T}, b::HuberParams{T}) where {T<:AbstractFloat}
    return HuberParams(
                        a.σ - b.σ, 
                        a.tᵣ - b.tᵣ, 
                        a.τ - b.τ,
                        a.H - b.H)
end

import Base.+
function +(a::HuberParams{T}, b::HuberParams{T}) where {T<:AbstractFloat}
    return HuberParams(
                        a.σ + b.σ, 
                        a.tᵣ + b.tᵣ, 
                        a.τ + b.τ,
                        a.H + b.H)
end

function +(a::HuberParams{T}, d::R) where {T<:AbstractFloat,R<:Real}
    return HuberParams(
                        a.σ+d, 
                        a.tᵣ+d, 
                        a.τ+d,
                        a.H+d)
end

import Base.*
function *(a::HuberParams{T}, b::HuberParams{T}) where {T<:AbstractFloat}
    return HuberParams(
                        a.σ*b.σ, 
                        a.tᵣ*b.tᵣ, 
                        a.τ*b.τ,
                        a.H*b.H)
end

function *(a::HuberParams{T}, d::R) where {T<:AbstractFloat,R<:Real}
    return HuberParams(
                        a.σ*d, 
                        a.tᵣ*d, 
                        a.τ*d,
                        a.H*d)
end

*(d::R, a::HuberParams{T}) where {T<:AbstractFloat,R<:Real} = *(a, d)

import Base./
function /(a::HuberParams{T}, b::HuberParams{T}) where {T<:AbstractFloat}
    return HuberParams(
                        a.σ/b.σ, 
                        a.tᵣ/b.tᵣ, 
                        a.τ/b.τ,
                        a.H/b.H)
end

import Base./
function /(a::HuberParams{T}, d::R) where {T<:AbstractFloat,R<:Real}
    return HuberParams(
                        a.σ/d, 
                        a.tᵣ/d, 
                        a.τ/d,
                        a.H/d)
end

function norm(a::HuberParams{T}) where {T<:AbstractFloat}
    return sqrt(a.σ^2 + a.tᵣ^2 + a.τ^2 + a.H^2)
end

import Base.sqrt
function sqrt(a::HuberParams{T}) where {T<:AbstractFloat}
    return HuberParams(
        sqrt(a.σ),sqrt(a.tᵣ),sqrt(a.τ),sqrt(a.H))#sqrt(a.σ^2 + a.tᵣ^2 + a.τ^2 + a.H^2)
end

mutable struct GD_state{P<:Params, T<:Real, I,J<:Integer}
    params::P
    t::Vector{T}
    y::Vector{T}
    mask::BitVector
    n::I
    max_index::J
end
#σ=p[1], tᵣ=p[2], τ=p[3], H = p[4]
function F!(state::GD_state{HuberParams{T}, U, I, J}) where {T,U<:AbstractFloat, I,J<:Integer}
    #σ=p[1], tᵣ=p[2], τ=p[3], H = p[4]
    #Given parameters in 'p'
    #Evaluate EGH function at eath time point tᵢ and store them in pre-allocated array 'f'. 
    #for (i, tᵢ) in enumerate(x)
    for i in range(1, state.max_index)
        #state.mask[i] ? continue : nothing
        tᵢ = state.t[i]
        d = 2*state.params.σ + state.params.τ*(tᵢ - state.params.tᵣ)
        if real(d) > 0
            state.y[i] = state.params.H*exp((-(tᵢ - state.params.tᵣ)^2)/d)
        else
            state.y[i] = zero(T)
        end
    end
end

function F!(params::HuberParams{T}, state::GD_state{HuberParams{T}, U, I, J}) where {T,U<:AbstractFloat, I,J<:Integer}

    #Evaluate EGH function at eath time point tᵢ and store them in `state.y` 
    for i in range(1, state.max_index)
        state.mask[i] ? continue : nothing #If data point is masked, skip
        tᵢ = state.t[i]
        d = 2*params.σ + params.τ*(tᵢ - params.tᵣ)
        if real(d) > 0
            state.y[i] = params.H*exp((-(tᵢ - params.tᵣ)^2)/d)
        else
            state.y[i] = zero(T)
        end
    end
end

function Jacobian(state::GD_state{HuberParams{T}, U, I, J}, y::Vector{T}, δ::T) where {T,U<:AbstractFloat, I,J<:Integer}

    #Initialize parameters
    H, τ, σ, tᵣ = state.params.H, state.params.τ, state.params.σ, state.params.tᵣ
    J_σ, J_tᵣ, J_τ, J_H = 0.0, 0.0, 0.0, 0.0

    for i in range(1, state.max_index)
        state.mask[i] ? continue : nothing #If data point is masked, skip

        #Compute repeating terms in the partial derivatives 
        tᵢ = state.t[i]
        DT = tᵢ - tᵣ
        T2 = DT^2
        D = 2*σ + DT*τ
        T2D2 = (T2/D^2)
        EXP = exp(-1.0*T2/D)
        N = H*EXP - y[i]
        Denom = sqrt(1 + (N/δ)^2)
        Common = N*EXP/Denom
        if 2*σ + τ*(DT) <= 0.0
            continue
        end

        J_σ += -2.0*H*Common*(-1.0*T2D2)
        J_tᵣ += H*(2*DT/D + τ*T2D2)*Common
        J_τ += -1.0*H*Common*DT*(-1.0*T2D2)
        J_H += Common

    end

    return HuberParams(J_σ, J_tᵣ, J_τ, J_H)
end

function reset!(state::GD_state{P,T,I,J}) where {P<:Params, T<:Real, I,J<:Integer} 
    for i in range(1, state.max_index)
        state.t[i], state.y[i] = zero(T), zero(T)
    end
    state.max_index = 0
    state.n = 0
    return 
end

function updateParams(params::HuberParams{T}, lower::HuberParams{T}, upper::HuberParams{T}) where {T<:AbstractFloat}
     return HuberParams(
         max(min(params.σ, upper.σ), lower.σ),
         max(min(params.tᵣ, upper.tᵣ), lower.tᵣ),
         max(min(params.τ, upper.τ), lower.τ),
         max(min(params.H, upper.H), lower.H)
     )
end

function GD(state::GD_state{P,T,I,J}, lower::P, upper::P; 
            tol::Float64 = 1e-3, 
            max_iter::Int64 = 1000, 
            δ::Float64 = 1000.0,
            α::Float64 = 0.1,
            β1 = 0.9,
            β2 = 0.999,
            ϵ = 1e-8) where {P<:Params, T<:Real, I,J<:Integer} 
    
    #Iteration counter
    state.n = 1

    #Initialize moment estimates 
    m1 = HuberParams(0.0, 0.0, 0.0, 0.0)
    m2 = HuberParams(0.0, 0.0, 0.0, 0.0)
    while state.n <= max_iter
        t = state.n #Time step
        F!(state) #Evaluate function with current parameters
        grad = Jacobian(state, δ) #Compute the gradient 
        m1 = β1*m1 + (1 - β1)*grad #First moment estimate
        m2 = β2*m2 + (1 - β2)*grad*grad #second moment estimate
        αₜ = α*sqrt(1 - β2^t)/(1 - β1^t) #Bias correction 
        diff = αₜ*m1/(sqrt(m2) + ϵ)
        #Get new parameters
        state.params = updateParams(state.params - diff, lower, upper)
        #If the euclidian norm of the update has fallen below a throeshold
        if norm(diff) <= tol
            return 
        end

        state.n += 1
    end    
end
######
lower = HuberParams(0.001, -Inf64, -1.0, 0.05);
upper = HuberParams(1.0, Inf64, 1.0, 2.0);


N = 8
t = collect(LinRange(-2.0, 5.0, N))
y = zeros(Float64, length(t))
TRUE_STATE = GD_state(
    HuberParams(1.0, 0.0, 0.5, 1.0),
    t,
    y,
    falses(length(t)),
    0, 
    length(t)
)
F!(TRUE_STATE)
data = max.(0.0, TRUE_STATE.y .+ 0.0*randn(length(TRUE_STATE.y)))
data[5] += 1.0
plot(TRUE_STATE.t, data, seriestype=:scatter)
######
t = collect(LinRange(-2.0, 5.0, 500))
y = zeros(Float64, length(t))
TRUE_STATE = GD_state(
    HuberParams(1.0, 0.0, 0.5, 1.0),
    t,
    y,
    falses(length(t)),
    0, 
    length(t)
)
F!(TRUE_STATE)
plot!(TRUE_STATE.t, TRUE_STATE.y)
##########
mask = falses(N)
STATE = GD_state(#HuberParams(1.0, 0.0, 1.0, 10000.0),
                HuberParams(0.1, 0.5, 0.001, 0.6),
                collect(LinRange(-2.0, 5.0, N)),
                zeros(Float64, N),
                mask,
                0,
                N)
GD(STATE, data, 
        lower, 
        upper, 
        tol = 5e-3, 
        max_iter = 100, 
        δ = 1e-2,
        α=0.1,
        β1 = 0.9,
        β2 = 0.999,
        ϵ = 1e-8)
##########
STATE_ =GD_state(STATE.params, 
            collect(LinRange(-2.0, 5.0, 500)),
            zeros(Float64, 500),
            falses(length(t)),
            0,
            500)
F!(STATE_)
plot!(STATE_.t, STATE_.y)
STATE.n
STATE_ =GD_state( HuberParams(0.1, 0.5, 0.001, 0.6), 
            collect(LinRange(-2.0, 5.0, 500)),
            zeros(Float64, 500),
            falses(length(t)),
            0,
            500)
F!(STATE_
    )
plot!(STATE_.t, STATE_.y)




H_ = []
loss_square = []
loss_huber = []
δ = 100.0
for i in LinRange(-10*1000, 20*1000, 1000)
    STATE.params = HuberParams(0.1, 0.5, 0.001, 6000.0 + i)
    F!(STATE)
    push!(loss_square, sum((data .- STATE.y).^2))
    push!(loss_huber, sum((δ^2).*(sqrt.((1 .+ ((data .- STATE.y)/δ).^2)) .- 1)))
    push!(H_, 6000.0 + i)
end
plot(H_, loss_square)
plot(H_, loss_huber)


H_ = []
loss_square = []
loss_huber = []
δ = 100.0
for i in LinRange(-1, 2, 1000)
    STATE.params = HuberParams(0.1, 0.5 + i, 0.001, 10000.0)
    F!(STATE)
    push!(loss_square, sum((data .- STATE.y).^2))
    push!(loss_huber, sum((δ^2).*(sqrt.((1 .+ ((data .- STATE.y)/δ).^2)) .- 1)))
    push!(H_, 6000.0 + i)
end
plot(H_, loss_square)
plot(H_, loss_huber)


H_ = []
x = []
y = []
loss_square = []
loss_huber = []
δ = 100.0
for i in LinRange(-1, 3, 100)
    for j in LinRange(-10*1000, 20*1000, 100)
    STATE.params = HuberParams(0.1+ i, 0.5, 0.001, 6000.0 + j)
    F!(STATE)
    push!(loss_square, sum((data .- STATE.y).^2))
    push!(loss_huber, sum((δ^2).*(sqrt.((1 .+ ((data .- STATE.y)/δ).^2)) .- 1)))
    push!(x,  0.1 + i)
    push!(y, 6000.0 + j)
    end
end

surface(x, y, loss_huber)
layout = Layout(

    title="Mt Bruno Elevation", autosize=false,

    width=500, height=500,

    margin=attr(l=65, r=50, b=65, t=90)

)
trace1 =  PlotlyJS.surface(z = loss_huber, x = x, y = y)
PlotlyJS.plot(trace1, layout)


plot(H_, loss_square)
plot(H_, loss_huber)

reset!(STATE)

@time begin 
    for i in range(1, 1000)
        STATE.params = HuberParams(0.01, 0.5, 0.001, 10000.0)
        STATE.n = 0
        GD(STATE, data, tol = 1e-4, max_iter = 10000)
    end
end


t = collect(LinRange(-2.0, 5.0, 5))
y = zeros(Float64, length(t))
EGH_inplace(y, t, [1, 0.0, 1, 10000])
plot(t, y)
data = copy(y .+ 100*randn(length(y)))
plot!(t, max.(0.0, data), seriestype=:scatter)
#σ=p[1], tᵣ=p[2], τ=p[3], H = p[4]
test = HuberParams(1.0, 0.0, 1.0, 10000.0)

H = HuberParams(0.001, 0.5, 0.01, maximum(data))

t2 = collect(LinRange(-2.0, 5.0, 500))
y = zeros(Float64, length(t2))
EGH_inplace(y, t, [1, 0.0, 1, 10000])

for i in range(1, 30)
    EGH_inplace(y, t2, H)
    #println("H $H")
    plot!(t2, y, show = true)
    H_new = Jacobian(H, t, data, 100.0)
    #println("H_new ", H_new/(1000000*i))
    #H = H - H_new/(100000000*i)
    H = H - H_new/(norm(H_new)*i)
    sleep(0.2)
end



function updateGrad(p::QuadParams{T}, x::Vector{T}, r::Vector{T}) where {T<:Real}
    return  QuadParams(
                        sum((x.^3).*sign.(r)),
                        sum(((-1.0).*x).*sign.(r)),
                        sum((-1.0).*sign.(r))
                        )
end

mutable struct SGD_state{T <: Params, I <: Integer, U,V <: Real} # contains information regarding one iterattion sequence
    
    θ::T # iterate x_n
    θ_best::T #best parameters so far
    ∇θ::T # one gradient ∇f(x_n)
    f_best::V 
    γ::U # stepsize
    n::I # iteration counter
    t::Vector{V} #Independent variable
    y::Vector{V} #Dependent variable
    r::Vector{V} #Residuals
 
end

function gradientStep(state::SGD_state{QuadParams{T}}) where {T<:Real}
    θ = (state.θ.θ1 - state.γ*state.∇θ.θ1,
        state.θ.θ2 - state.γ*state.∇θ.θ2,
        state.θ.θ3 - state.γ*state.∇θ.θ3)
    return QuadParams(θ[1], θ[2], θ[3]), sqrt(sum(θ.^2))
end

function SGD(state::SGD_state, data::T, max_iter::Int64) where {T <: AbstractVecOrMat{<: Real}}
    θ_best = state.θ
    evalInplace(state.y, state.t, state.θ)
    state.r = state.y .- data
    state.∇θ = updateGrad(state.θ, state.t, state.r)
    state.f_best = sum(abs.(state.r))
    γ = state.γ
    i = 1
    last_L = Inf64
    while i < max_iter
        new_state, norm = gradientStep(state)#state.θ - state.γ*state.∇θ
        evalInplace(state.y, state.t, new_state)
        state.r = state.y .- data
        state.∇θ = updateGrad(state.θ, state.t, state.r)
        #if (maxNorm(state.θ, new_state) < 0.001) & (i > 100)
        #    break
        #end
        sum_of_abs_residuals = sum(abs.(state.r))
        if sum_of_abs_residuals < state.f_best
            state.f_best = sum_of_abs_residuals
            #println("last_L $last_L")
            #println("sum_of_abs_residuals $sum_of_abs_residuals")
            #if abs(abs((last_L - sum_of_abs_residuals))/sum_of_abs_residuals) < 0.001
            #    break
            #end
            state.θ = new_state
            θ_best = state.θ
        end
  
        last_L = sum_of_abs_residuals
        #println(γ/norm)
        #println("norm $norm")
        #println("∇θ ", state.∇θ)
        state.γ = (γ/norm)/(i + 1)
        state.n = i
        i += 1
    end
    state.θ = θ_best
    println("number of iterations to convergence ", state.n)
    return state#, L, states
end

x = collect(LinRange(-10, 10, 150))
ŷ = zeros(Float64, length(x))
evalInplace(ŷ, x, QuadParams(5.0, 100.0, 0.0))
y = copy(ŷ) .+ 100*randn(length(ŷ))
plot(x, y, seriestype=:scatter)
plot!(x, ŷ)
state = SGD(
    SGD_state(
    QuadParams(5.0, 1.0, 1.0),
    QuadParams(5.0, 1.0, 1.0),
    QuadParams(5.0, 1.0, 1.0),
    Float64(Inf),
    0.01/length(x),
    1,
    x, 
    zeros(Float64, length(x)),
    zeros(Float64, length(x))
    ),
    y,
    200000
);
plot!(x, state.y)




plot(x, y, seriestype=:scatter)
plot!(x, ŷ)
println(sum(abs.(y .- ŷ)))
θ = (1.0, 1.0, 0.0)
quadfunc!(x, ŷ, θ)
∇θ = zeros(Float64, 3)
for i in range(1, 100)
    quadfunc!(x, ŷ, θ)
    residuals = ŷ .- y
    signs = sign.(residuals)
    ∇θ = ∇quadfunc(x, θ, signs)
    println(sum(abs.(residuals)))
    new_theta = θ .- (0.001).*∇θ
    println(sum(abs.(new_theta .- θ)))
    θ = new_theta
end
println(θ)
plot!(x, ŷ)

quadfunc(x::T)




t = collect(LinRange(-2.0, 5.0, 100))
y = zeros(Float64, length(t))
EGH_inplace(y, t, [1, 0.0, 1, 10000])
plot(t, y)
plot!(t, max.(0.0, y .+ 1000*randn(length(y))), seriestype=:scatter)


@variables tᵣ, σ, τ, H, t, δ, y

EGH = H*exp( (-(t - tᵣ)^2)/(2*σ + τ*(t - tᵣ)))
D_tᵣ = Differential(tᵣ)
D_σ = Differential(σ)
D_H = Differential(H)
D_τ = Differential(τ)

expand_derivatives(D_tᵣ(EGH))
latexify(expand_derivatives(D_τ(EGH)))
latexify(expand_derivatives(D_H(EGH)))
latexify(expand_derivatives(D_σ(EGH)))

HUBER_EGH = (δ^2)*(sqrt(1 + ((EGH - y)/δ)^2) - 1)

latexify(HUBER_EGH)
latexify(expand_derivatives(D_tᵣ(HUBER_EGH)))
latexify(expand_derivatives(D_τ(HUBER_EGH)))
latexify(expand_derivatives(D_H(HUBER_EGH)))
latexify(expand_derivatives(D_σ(HUBER_EGH)))

latexify(Symbolics.jacobian([HUBER_EGH], [tᵣ, τ, σ, H]))

function GD(problem::GD_problem, x::T, y::T) where {T <: AbstractVecOrMat{<: Real}}
    state = GD_state(
        problem.x0,
        problem.∇f(problem.x0),
        problem.γ,
        1)

    problem.∇f.(data) problem.f.(state.x0, data) .- data
    for i in range(1, 10)
        new_state = state.x - state.γ*state.∇f_x
        println(new_state .- state.x)
        state.x = new_state
        state.γ = 1/(i + 1)
        state.∇f_x = problem.∇f(state.x, data)
    end
    return state
end

GD(GD_problem(quadfunc, ∇quadfunc, Float64[1, 1], 0.1))

function gradientDesc(x, y, learn_rate, conv_threshold, n, max_iter)
    β = rand(Float64, 1)[1]
    α = rand(Float64, 1)[1]
    ŷ = α .+ β .* x
    MSE = sum((y .- ŷ).^2)/n
    converged = false
    iterations = 0

    while converged == false
        # Implement the gradient descent algorithm
        β_new = β - learn_rate*((1/n)*(sum((ŷ .- y) .* x)))
        α_new = α - learn_rate*((1/n)*(sum(ŷ .- y)))
        α = α_new
        β = β_new
        ŷ = β.*x .+ α
        MSE_new = sum((y.-ŷ).^2)/n
        # decide on whether it is converged or not
        if (MSE - MSE_new) <= conv_threshold
            converged = true
            println("Optimal intercept: $α; Optimal slope: $β")
        end
        iterations += 1
        if iterations > max_iter
            converged = true
            println("Optimal intercept: $α; Optimal slope: $β")
        end
    end
end