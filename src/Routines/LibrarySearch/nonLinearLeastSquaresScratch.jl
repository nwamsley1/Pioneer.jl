function rosenbrock_f!(out, x)
    out[1] = 1 - x[1]
    out[2] = 100 * (x[2]-x[1]^2)
end
optimize!(LeastSquaresProblem(x = zeros(2), f! = rosenbrock_f!, output_length = 2, autodiff = :central), Dogleg())
   

function EGH_inplace(F::Vector{T}, x::Vector{T}, p::Vector{T}) where {T<:AbstractFloat}
    #σ=p[1], tᵣ=p[2], τ=p[3], H = p[4]
    #Given parameters in 'p'
    #Evaluate EGH function at eath time point tᵢ and store them in pre-allocated array 'f'. 
     for (i, tᵢ) in enumerate(x)
        d = 2*p[1] + p[3]*(tᵢ - p[2])
        if real(d) > 0
            F[i] = p[4]*exp((-(tᵢ - p[2])^2)/d)
        else
            F[i] = zero(T)
        end
    end
end

function JEGH_inplace(J::Matrix{T}, x::Vector{T}, p::Vector{T}) where {T<:AbstractFloat}
    for (i, tᵢ) in enumerate(x)
        δt = tᵢ - p[2]
        d = 2*p[1] + p[3]*δt
        q = δt/d
        f = exp((-δt^2)/(d))
        #f = exp((-δt^2)/(d))
        if d > 0
            J[i,1] = (2*(q)^2)*p[4]*f
            J[i,2] = p[4]*(2*q - (p[3]*(q^2)))*f
            J[i,3] = ((δt)*(q^2))*p[4]*f
            J[i,4] = f
        end
   end
end
p0 = Float32[0.004885812921411583, 38.0, 0.0, 620000.0]
FIT = LsqFit.curve_fit(EGH_inplace, JEGH_inplace, huber_loss[:,:rt], huber_loss[:,:weight], huber_loss[:,:rank]./4, p0; inplace = true, show_trace = true)