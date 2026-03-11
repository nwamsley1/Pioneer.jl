# GLM.jl comparison test — run separately:  julia test_glm_comparison.jl
# Requires: GLM, DataFrames, Distributions

using Random, Printf, GLM, DataFrames, Distributions

include("SparseArray.jl")
include("spectralPoissonRegression.jl")

const NR_ITER   = 25
const BS_ITER   = 100
const OUT_ITER  = 2000
const NR_ACC    = 1f-6
const BS_ACC    = 1f-6
const REL_CONV  = 1f-4

function run_poisson(A::Matrix{Float32}, y::Vector{Float32}; x0::Union{Nothing,Vector{Float32}}=nothing)
    sa = build_SparseArray(A, y)
    m, n = size(A)
    μ  = zeros(Float32, m)
    yy = zeros(Float32, m)
    w  = x0 === nothing ? ones(Float32, n) : copy(x0)
    initObserved!(yy, sa)
    initMu!(μ, sa, w)
    solvePoisson!(sa, μ, yy, w, NR_ITER, BS_ITER, OUT_ITER, NR_ACC, BS_ACC, REL_CONV)
    return w, μ, yy, sa
end

# Simple Poisson sampler
function rand_poisson_manual(λ::Float64)
    λ <= 0 && return 0
    L = exp(-λ)
    k = 0
    p = 1.0
    while true
        k += 1
        p *= rand()
        p < L && return k - 1
    end
end

println("── GLM.jl comparison (identity link, no intercept) ──")
Random.seed!(123)
m, n = 30, 3
A = Float32.(abs.(randn(m, n)) .+ 0.5)
x_true = Float32[4, 8, 2]
mu_true = A * x_true
y = Float32.([rand_poisson_manual(Float64(mu_true[i])) for i in 1:m])

# Our solver
w_ours, _, _, _ = run_poisson(A, y; x0=ones(Float32, n))

# GLM.jl with identity link, no intercept
df = DataFrame(A, :auto)
df.y = Float64.(y)
col_syms = Symbol.(names(df)[1:n])
f = term(:y) ~ sum(term.(col_syms)) + ConstantTerm(0)
glm_model = glm(f, df, Poisson(), IdentityLink())
w_glm = Float32.(coef(glm_model))

@printf("  true = %s\n", string(round.(x_true; digits=3)))
@printf("  ours = %s\n", string(round.(w_ours[1:n]; digits=3)))
@printf("  GLM  = %s\n", string(round.(w_glm; digits=3)))
max_diff = maximum(abs.(w_ours[1:n] .- w_glm))
@printf("  max |diff| = %.4f\n", max_diff)
passed = max_diff < 0.5
println("  ", passed ? "PASS" : "FAIL")
