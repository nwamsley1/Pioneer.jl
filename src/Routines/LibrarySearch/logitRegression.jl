

function model_predict(psms::DataFrame, model_fit::Any, column_names::Vector{Symbol})
    coefs = coef(model_fit)
    psms[!,:prob] = zeros(Float16, size(psms, 1))
    function addColumn(probs::Vector{Float16}, vals::Vector{R}, x::Float64) where {R<:Real}
        Threads.@threads for i in range(1, size(psms, 1))
            probs[i] += vals[i]*x
        end
    end
    for (i, name) in enumerate(column_names)
        x = coefs[i+1]
        addColumn(psms[!,:prob], psms[!,name], x)
    end
    intercept = coefs[1]::Float64
    Threads.@threads for i in range(1, size(psms, 1))
        psms[i,:prob] += Float16(intercept)
    end

    Threads.@threads for i in range(1, size(psms, 1))
        psms[i,:prob] = exp( psms[i,:prob] )/(1 + exp( psms[i,:prob] ))
    end

end

#=

X = Matrix(PSMs[!,:])
logistic = LogisticRegression(0.0, fit_intercept = false)

solver = MLJLinearModels.NewtonCG(
    optim_options = Optim.Options(time_limit = 20),
    newtoncg_options = (eta = 0.2,)
)
@time theta = MLJLinearModels.fit(logistic, X, y, solver = solver)

solver = MLJLinearModels.LBFGS(
    optim_options = Optim.Options(time_limit = 20),
    lbfgs_options = (linesearch = Optim.LineSearches.HagerZhang(),))


@time theta = MLJLinearModels.fit(logistic, X, y, solver = solver)

PSMs[!,:intercept] = ones(Float16, size(PSMs, 1))
X = Matrix(PSMs[!,[:intercept,:spectral_contrast,:scribe,:city_block,:entropy_score,:iRT_error,:missed_cleavage,:Mox,:charge,:TIC,:total_ions,:err_norm,:spectrum_peak_count]])
y = PSMs[!,:target]

logistic = LogisticRegression(fit_intercept = false)
solver = MLJLinearModels.NewtonCG(
    optim_options = Optim.Options(time_limit = 20),
    newtoncg_options = (eta =0.0,)
)
@time theta = MLJLinearModels.fit(logistic, X, y, solver = solver)
PSMs[!,:prob] = X*theta
getQvalues!(PSMs, allowmissing(PSMs[!,:prob]),  allowmissing(PSMs[:,:decoy]));
println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)))
println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))
  
@time theta = MLJLinearModels.fit(logistic, X, y, solver = solver)



PSMs[!,:prob] = X*theta[2:end]



function generate_dataset(n_samples = 100, n_features = 10; shift = 0.0)
    X = randn(n_samples, n_features)
    w = randn(n_features)
    y = sign.(X * w)
    X .+= 0.8 * randn(n_samples, n_features) # add noise
    X .+= shift # shift the points in the feature space
    X = hcat(X, ones(n_samples, 1))
    return X, y
end

function softplus(model, t, u)
    z = @variable(model, [1:2], lower_bound = 0.0)
    @constraint(model, sum(z) <= 1.0)
    @constraint(model, [u - t, 1, z[1]] in MOI.ExponentialCone())
    @constraint(model, [-t, 1, z[2]] in MOI.ExponentialCone())
end

function build_logit_model(X, y, λ)
    n, p = size(X)
    model = Model()
    @variable(model, θ[1:p])
    @variable(model, t[1:n])
    for i in 1:n
        u = -(X[i, :]' * θ) * y[i]
        softplus(model, t[i], u)
    end
    # Add ℓ2 regularization
    @variable(model, 0.0 <= reg)
    @constraint(model, [reg; θ] in SecondOrderCone())
    # Define objective
    @objective(model, Min, sum(t) + λ * reg)
    return model
end

n, p = 200000, 10
X, y = generate_dataset(n, p; shift = 10.0);

λ = 10.0
model = build_logit_model(X, y, λ)
set_optimizer(model, SCS.Optimizer)
set_silent(model)
@time JuMP.optimize!(model)

θ♯ = JuMP.value.(model[:θ])


#=
y = Float16.(PSMs[!,:target])
PSMs[!,[:spectral_contrast,:scribe,:entropy_score,:iRT_error,:missed_cleavage,:Mox,:charge,:TIC,:total_ions,:err_norm,:spectrum_peak_count]]

@time X = Matrix(PSMs[!,[:spectral_contrast,:scribe,:entropy_score,:iRT_error,:missed_cleavage,:Mox,:charge,:TIC,:total_ions,:err_norm,:spectrum_peak_count]])

theta_guess = ones(Float64, size(X, 2))

@time theta, flag, history = gradient_descent_probit(y, X, theta_guess, max_iter=1000, tol=0.00001);

hessian_log_likelihood(y, X, theta_guess)
=#
logistic = LogisticRegression(0.5)
theta = fit(logistic, X, y)



PSMs[!,:intercept] = ones(Float16, size(PSMs, 1))
X = Matrix(PSMs[!,[:intercept,:spectral_contrast,:scribe,:city_block,:entropy_score,:iRT_error,:missed_cleavage,:Mox,:charge,:TIC,:total_ions,:err_norm,:spectrum_peak_count]])
Y = PSMs[!,:target]
n, p = size(X)
beta = Variable(p)
problem = minimize(logisticloss(-Y .* (X * beta)))
solve!(problem, SCS.Optimizer; silent_solver = true)

=#