# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
    RazoQuadParams

Piecewise model for quadrupole transmission based on general gaussian

### Fields

- al::T -- Distance from center to half-max of LHS of function
- al::T -- Distance from center to half-max of RHS of function
- bl::T -- Controls slope at half-max of LHS of function 
- br::T -- Controls slope at half-max of RHS of function 

### Examples

rqm = RazoQuadParams(1.0, 2.0, 5.0, 20.0)
plot(LinRange(-3, 3, 100), rqm.(LinRange(-3, 3, 100)))
"""
struct RazoQuadParams{T <: AbstractFloat}
    al::T
    ar::T
    bl::T
    br::T
end

"""
    RazoQuadParams

Evaluate a RazoQuadParams model at x

### Examples

rqm = RazoQuadParams(1.0, 2.0, 5.0, 20.0)
plot(LinRange(-3, 3, 100), rqm.(LinRange(-3, 3, 100)))
"""
function (rqm::RazoQuadParams)(x::T) where {T<:AbstractFloat}
    if x < 0
        return 1/(1 + abs(x/rqm.al)^(2*rqm.bl))
    else
        return 1/(1 + abs(x/rqm.ar)^(2*rqm.br))    
    end
end

import Base: +  # Import + to extend it
function +(a::RazoQuadParams, b::RazoQuadParams)
    return RazoQuadParams(
        a.al + b.al,
        a.ar + b.ar,
        a.bl + b.bl,
        a.br + b.br
    )
end

import Base: -  # Import + to extend it
function -(a::RazoQuadParams, b::RazoQuadParams)
    return RazoQuadParams(
        a.al - b.al,
        a.ar - b.ar,
        a.bl - b.bl,
        a.br - b.br
    )
end



import Base: /  # Import + to extend it
function /(a::RazoQuadParams, b::RazoQuadParams)
    return RazoQuadParams(
        a.al/b.al,
        a.ar/b.ar,
        a.bl/b.bl,
        a.br/b.br
    )
end

import Base: abs  # Import + to extend it
function abs(a::RazoQuadParams)
    return RazoQuadParams(
        abs(a.al),
        abs(a.ar),
        abs(a.bl),
        abs(a.br),
    )
end

import Base: maximum  # Import + to extend it
function maximum(a::RazoQuadParams)
    _max_ = a.al
    if a.ar > _max_
        _max_ = a.ar
    end
    if a.bl > _max_
        _max_ = a.bl
    end
    if a.br > _max_
        _max_ = a.br
    end
    return _max_
end

function RazoQuadParams(params::AbstractVector{T}) where {T<:AbstractFloat}
    return RazoQuadParams(
                        params[1],
                        params[2],
                        params[3],
                        params[4]
        )
end

function inBounds(rqm::RazoQuadParams{T},
                      u_bounds::RazoQuadParams{T},
                      l_bounds::RazoQuadParams{T}) where {T<:AbstractFloat}
    return ((rqm.al > u_bounds.al
    ) |     (rqm.ar > u_bounds.ar
    ) |     (rqm.bl > u_bounds.bl
    ) |     (rqm.br > u_bounds.br
    ) |     (rqm.al < l_bounds.al
    ) |     (rqm.ar < l_bounds.ar
    ) |     (rqm.bl < l_bounds.bl
    ) |     (rqm.br < l_bounds.br)) == false
end

struct RazoQuadModel{T<:AbstractFloat} <: QuadTransmissionModel
    params::RazoQuadParams{T}
end

struct RazoQuadFunction{T<:AbstractFloat} <: QuadTransmissionFunction 
    min_mz::T
    max_mz::T
    center_mz::T
    params::RazoQuadParams{T}
end

function getQuadTransmissionFunction(rqm::RazoQuadModel{T}, centerMz::T, isolationWidthMz::T) where {T<:AbstractFloat}
    al, ar = rqm.params.al, rqm.params.ar
    RazoQuadFunction(
        centerMz - al,
        centerMz + ar,
        centerMz,
        rqm.params
    )
end

function getQuadTransmissionBounds(rqm::RazoQuadModel{T}, centerMz::T, isolationWidthMz::T) where {T<:AbstractFloat}
        al, ar = rqm.params.al, rqm.params.ar
        return T(centerMz - al), T(centerMz + ar)
        #return T(centerMz - isolationWidthMz/2), T(centerMz + isolationWidthMz/2)
end

function getPrecMinBound(rqm::RazoQuadFunction{T}) where {T<:AbstractFloat}
    return rqm.min_mz
end

function getPrecMaxBound(rqm::RazoQuadFunction{T}) where {T<:AbstractFloat}
    return rqm.max_mz
end

function (rqf::RazoQuadFunction{T})(ionMz::U) where {T,U<:AbstractFloat}
    rqf.params(ionMz-rqf.center_mz)
end

##############
#Model Fitting


#=
"""
    
    function getRazoQuadJacobian() 

Returns a piecewise function to evaluate the jacobian for the function 

F = log(f(x0)/f(x1)) with respect to the for parameters al, ar, bl, br (see `RazoQuadParams`)


### Input

### Output

- Dict{Int64, Any}
Output is a dictionary with keys -1, 0, 1 for different parts of the domain. 

* -1 Corresponds to x0 < 0 and x1 < 0
* 0 Corresponds to x0 < 0 and x1 > 0
* 1 Corresponds to x0 > 0 and x1 > 0

### Notes

- The derivative is very complicated. This function uses the "Symbolics.jl" package to get functions for the derivatives. 
- The derivatives will be NaN if either x0 or x1 is exactly equal to zero. 

"""
function getRazoQuadJacobian()
    #From Symbolics.jl
    @variables al, ar, bl, br, x0, x1
    jacobian_dict = Dict{Int64, Any}()
    #x0 < 0, x1 < 0
    z = log(1 + abs((x1)/al)^(2*bl)) - log(1 + abs(x0/al)^(2*bl))
    jacobian_dict[-1] = eval(build_function([simplify(x) for x in Symbolics.jacobian([z], [al, ar, bl, br])],
        [x0, x1, al, ar, bl, br])[1])
    #x0 < 0, x1 > 0
    z = log(1 + abs((x1)/ar)^(2*br)) - log(1 + abs(x0/al)^(2*bl))
    jacobian_dict[0] = eval(build_function([simplify(x) for x in Symbolics.jacobian([z], [al, ar, bl, br])],
        [x0, x1, al, ar, bl, br])[1])
    #x0 > 0, x1 > 0
    z = log(1 + abs((x1)/ar)^(2*br)) - log(1 + abs(x0/ar)^(2*br))
    jacobian_dict[1] = eval(build_function([simplify(x) for x in Symbolics.jacobian([z], [al, ar, bl, br])],
        [x0, x1, al, ar, bl, br])[1])
    return jacobian_dict
end

"""
    function evalJ(Jsymb::Dict{Int64, Any}, params::RazoQuadParams, x0::T, x1::T) where {T<:AbstractFloat} 

Given a piecewise function to evaluate the jacobian of F, evaluates the derivatives given `params`, `x0`, and `x1`
- f -- Razo Quad Transmission function 
- F = log(f(x0)/f(x1)) -- with respect to the for parameters al, ar, bl, br (see `RazoQuadParams`)

### Input

- `Jsymb::Dict{Int64, Any}}`: -- See `getRazoQuadJacobian`. Pieciwise function to evaluate the jacobian given parameters and inputs
- `params::RazoQuadParams` -- parameters for a Razo Quad Transmission function. See `RazoQuadParams`
- `x0::T` -- point to evaluate F(x0, x1) at. 
- `x1::T` -- point to evaluate F(x0, x1) at. 

### Output

- SMatrix{1, 4, Float64, 4} -- The parameters of a `RazoQuadParams` model in this order: al, ar, bl, br

"""
function evalJ(Jsymb::Dict{Int64, Any}, params::RazoQuadParams, x0::T, x1::T)::SMatrix{1, 4, Float64, 4} where {T<:AbstractFloat}
    params = SA[x0, x1, params.al, params.ar, params.bl, params.br]
    if (x0 < 0) & (x1 < 0)
        return Jsymb[-1](params)
    elseif (x0 < 0) & (x1 > 0)
        return Jsymb[0](params)
    else
        return Jsymb[1](params)
    end
end
=#

"""
    
    function getRazoQuadJacobian() 

Returns a piecewise function to evaluate the jacobian for the function 

F = log(f(x0)/f(x1)) with respect to the for parameters al, ar, bl, br (see `RazoQuadParams`)


### Input

### Output

- Dict{Int64, Any}
Output is a dictionary with keys -1, 0, 1 for different parts of the domain. 

* -1 Corresponds to x0 < 0 and x1 < 0
* 0 Corresponds to x0 < 0 and x1 > 0
* 1 Corresponds to x0 > 0 and x1 > 0

### Notes

- The derivative is very complicated. This function uses the "Symbolics.jl" package to get functions for the derivatives. See the commented out function above 

"""
function razoJ(x0::U, x1::U, rqm::RazoQuadParams{T}) where {T,U<:AbstractFloat}
    al, ar, bl, br = rqm.al, rqm.ar, rqm.bl, rqm.br
    #The (abs(x0 / al)^(2bl)) has been pulled out of the log function in these expressions, so can evaluate to NaN
    #Will still get a correct answer if we replace zero values with something small 
    if iszero(x0) 
        x0 = U(1e-5)
    end
    if iszero(x1)
        x1 = U(1e-5)
    end

    if (x0 < 0) & (x1 < 0)
        return SMatrix{1, 4, T, 4}([
            ((2//1)*bl*x0*ifelse(signbit(x0 / al), -1, 1)*(abs(x0 / al)^(-1 + 2bl)) - (2//1)*bl*x1*(abs(x1 / al)^(-1 + 2bl))*ifelse(signbit(x1 / al), -1, 1) + (2//1)*bl*x0*ifelse(signbit(x0 / al), -1, 1)*(abs(x1 / al)^(2bl))*(abs(x0 / al)^(-1 + 2bl)) - (2//1)*bl*x1*(abs(x1 / al)^(-1 + 2bl))*ifelse(signbit(x1 / al), -1, 1)*(abs(x0 / al)^(2bl))) / ((al^2)*(1 + abs(x1 / al)^(2bl))*(1 + abs(x0 / al)^(2bl))),
            0,
            (2(abs(x1 / al)^(2bl))*log(abs(x1 / al)) - 2log(abs(x0 / al))*(abs(x0 / al)^(2bl)) - 2(abs(x1 / al)^(2bl))*log(abs(x0 / al))*(abs(x0 / al)^(2bl)) + 2(abs(x1 / al)^(2bl))*log(abs(x1 / al))*(abs(x0 / al)^(2bl))) / ((1 + abs(x1 / al)^(2bl))*(1 + abs(x0 / al)^(2bl))),
            0,
        ])
    elseif (x0 < 0) & (x1 > 0)
        return SMatrix{1, 4, T, 4}([
            (2bl*x0*ifelse(signbit(x0 / al), -1, 1)*(abs(x0 / al)^(-1 + 2bl))) / ((al^2)*(1 + abs(x0 / al)^(2bl))),
            (-2br*x1*(abs(x1 / ar)^(-1 + 2br))*ifelse(signbit(x1 / ar), -1, 1)) / ((ar^2)*(1 + abs(x1 / ar)^(2br))),
            (-2log(abs(x0 / al))*(abs(x0 / al)^(2bl))) / (1 + abs(x0 / al)^(2bl)),
            (2log(abs(x1 / ar))*(abs(x1 / ar)^(2br))) / (1 + abs(x1 / ar)^(2br))
        ])
    else
        return SMatrix{1, 4, T, 4}([
            0,
            ((2//1)*br*x0*ifelse(signbit(x0 / ar), -1, 1)*(abs(x0 / ar)^(-1 + 2br)) - (2//1)*br*x1*(abs(x1 / ar)^(-1 + 2br))*ifelse(signbit(x1 / ar), -1, 1) + (2//1)*br*x0*ifelse(signbit(x0 / ar), -1, 1)*(abs(x0 / ar)^(-1 + 2br))*(abs(x1 / ar)^(2br)) - (2//1)*br*x1*(abs(x0 / ar)^(2br))*(abs(x1 / ar)^(-1 + 2br))*ifelse(signbit(x1 / ar), -1, 1)) / ((ar^2)*(1 + abs(x0 / ar)^(2br))*(1 + abs(x1 / ar)^(2br))),
            0,
            (-2(abs(x0 / ar)^(2br))*log(abs(x0 / ar)) + 2log(abs(x1 / ar))*(abs(x1 / ar)^(2br)) + 2(abs(x0 / ar)^(2br))*log(abs(x1 / ar))*(abs(x1 / ar)^(2br)) - 2(abs(x0 / ar)^(2br))*log(abs(x0 / ar))*(abs(x1 / ar)^(2br))) / ((1 + abs(x0 / ar)^(2br))*(1 + abs(x1 / ar)^(2br)))
        ])
    end
end


"""
    function  F(rqm::RazoQuadParams, x0::T, x1::T) where {T<:AbstractFloat}

Computes log ratio of a Razo Quad Transmission model (see `RazoQuadParams`) at points x0 and x1
- f -- Razo Quad Transmission function 
- F = log(f(x0)/f(x1)) -- with respect to the for parameters al, ar, bl, br (see `RazoQuadParams`)

### Input
- `rqm::RazoQuadParams` -- parameters for a Razo Quad Transmission function. See `RazoQuadParams`
- `x0::T` -- point to evaluate F(x0, x1) at. 
- `x1::T` -- point to evaluate F(x0, x1) at. 

### Output

- AbstractFloat 

"""
function F(rqm::RazoQuadParams, x0::T, x1::T) where {T<:AbstractFloat}
    al, ar, bl, br = rqm.al, rqm.ar, rqm.bl, rqm.br
    δ = x1 - x0
    if x1 < 0
        return Float32(log(1 + abs(x1/al)^(2*bl)) - log(1 + abs(x0/al)^(2*bl)))
    elseif (x0 < 0) & (x0 + δ > 0)
         return Float32(log(1 + abs(x1/ar)^(2*br)) - log(1 + abs(x0/al)^(2*bl)))
    else
         return Float32(log(1 + abs(x1/ar)^(2*br)) - log(1 + abs(x0/ar)^(2*br)))
    end 
end


"""
    fitRazoQuadModel(
        al0::T, ar0::T, bl0::T, br0::T,
        al_bounds::Tuple{T, T},
        ar_bounds::Tuple{T, T},
        bl_bounds::Tuple{T, T},
        br_bounds::Tuple{T, T},
        x0_dat::Vector{T}, x1_dat::Vector{T}, yt_dat::Vector{T}
        ;
            λ0 = 1e-2, #Initial L-M damping coeficient 
            ϵ1 = 1e-4, #Convergence in the gradient
            ϵ2 = 1e-4, #Convergence in the coeficients
            ϵ3 = 1e-4, #Conergence in squared error
            lup = 9,
            ldown = 11,
            max_iter = 100000, #Maximum iterations
        )

Fits a RazoQuadModel from data. The data are the measurement of ratios 

- See `Pioneer.jl/data/references/levenberg_marquardt.pdf` for a guid to understand the algorithm. The Wikipedia on Levenberg-Marqaurdt is also helpful. 

### Input

- `al0, ar0, etc.`: -- Initial guess parameters for a Razo Quadrupole Transmission Model (see `RazoQuadParams`)
- `al_bounds::Tuple{T, T}, etc.` -- These tuples are lower and upper bounds respectively for each parameters. Will not estimate parameter values outside of these boundary constraints
- `Jsymb::Dict{Int64, Any}}`: -- See `getRazoQuadJacobian`. Pieciwise function to evaluate the jacobian given parameters and inputs
- `x0_dat::Vector{T}`, `x1_dat::Vector{T}` -- Equal length vectors that give the x0 and x1 for each data/measurment. These are m/z offset values relative to the quadrupole window center
- `yt_dat::Vector{T}` -- Equal in length to x0. This is the log of the observed ratio between two precursor isotopes with m/z x0 and x1 respectively. 
- `λ0 = 1e-2` -- Initial L-M damping coeficient 
- `ϵ1 = 1e-4` -- Convergence in the gradient
- `ϵ2 = 1e-4` -- Convergence in the coeficients
- `ϵ3 = 1e-4` -- Convergence in squared error
- `lup = 9` -- Factor to divide λ by if step improved the objective
- `ldown = 11` --Factor to multiply by λ by if step worsened the objective
- `max_iter = 100000` -- Maximum number of updates before terminating optimization 
### Output

- `RazoQuadParams` -- parameters of the fitted model

### Notes

### Algorithm 

Uses the Levenbert-Marquardt algorithm 

### Examples 

"""
function fitRazoQuadModel(
    al0::T, ar0::T, bl0::T, br0::T,
    al_bounds::Tuple{T, T},
    ar_bounds::Tuple{T, T},
    bl_bounds::Tuple{T, T},
    br_bounds::Tuple{T, T},
    x0_dat::Vector{T}, x1_dat::Vector{T}, yt_dat::Vector{T}
    ;
        λ0 = 1e-2, #Initial L-M damping coeficient 
        ϵ1 = 1e-4, #Convergence in the gradient
        ϵ2 = 1e-4, #Convergence in the coeficients
        ϵ3 = 1e-4, #Conergence in squared error
        lup = 9,
        ldown = 11,
        max_iter = 100000, #Maximum iterations
    ) where {T<:AbstractFloat}
    #Initial values and parameter bounds 
    rqm = RazoQuadParams(al0, ar0, bl0, br0)
    l_bounds = RazoQuadParams(first(al_bounds), first(ar_bounds), first(bl_bounds), first(br_bounds))
    u_bounds = RazoQuadParams(last(al_bounds), last(ar_bounds), last(bl_bounds), last(br_bounds))
    #Initialize JacobianMat
    M, N = length(x0_dat), 4
    J = zeros(T, (M, N));
    Jy = zeros(T, (N, 1))
    Jt = J'
    JtJ = zeros(T, (N, N))
    X = zeros(T, (N, N))
    y_diff = zeros(T, (M, 1))
    F_vec = zeros(T, M)
    Δ = zeros(T, N)
    λ = λ0
    for n in range(1, max_iter)
        #Update F and Jacobian 
        for i in range(1, M)
            F_vec[i] = F(rqm, x0_dat[i], x1_dat[i])
            J[i,:] = razoJ(x0_dat[i], x1_dat[i], rqm)
        end
        #Get diff 
        for i in range(1, M)
            y_diff[i] = yt_dat[i] - F_vec[i]
        end
        SE_old = zero(T)
        for i in range(1, M)
            SE_old += y_diff[i]^2
        end
        mul!(JtJ, Jt, J)
        #Jy = Jt*y_diff
        mul!(Jy, Jt, y_diff)
        X .= JtJ
        𝐼 = (λ)*Diagonal(JtJ)
        X .+= 𝐼 
        Δ .= X\Jy
        rqm_Δ = RazoQuadParams(Δ)
        rqm_new = rqm + rqm_Δ
        #If updated guess was out of bounds, increase damping coeficient 
        if !inBounds(rqm_new, u_bounds, l_bounds)
            λ *= 10
            continue
        end
        
        SE_new = zero(T)
        for i in range(1, M) #Compute F_vec with new coefficients 
            F_vec[i] = F(rqm_new, x0_dat[i], x1_dat[i])
            y_diff[i] = yt_dat[i] - F_vec[i]
            SE_new += y_diff[i]^2
        end
        
        #Whether to raise or lower the damping coefficient 
        ph = (SE_old - SE_new)/abs((Δ'*(𝐼*Δ + Jy))[1])
        #println("ph $ph")
        if ph > 0.1
            rqm = rqm_new
            λ = max(λ/ldown, 1)
        else
            λ = min(λ*lup, 1e4)
            if λ>1e3
                rqm = rqm_new
            end
            continue
        end
        
        if n > 10
            if maximum(abs.(Jy)) < ϵ1
                return rqm
            elseif maximum(abs(rqm_Δ/rqm)) < ϵ2
                return rqm        
            elseif SE_new/(M - N) < ϵ3
                return rqm
            end
        end
    end
    @warn "exceeded iters"
    return rqm
end

"""
    fitRazoQuadModel(
    quad_window_width::T,
    x0_dat::Vector{T}, 
    x1_dat::Vector{T}, 
    yt_dat::Vector{T};
        λ0 = 1e-2, #Initial L-M damping coeficient 
        ϵ1 = 1e-5, #Convergence in the gradient
        ϵ2 = 1e-4, #Convergence in the coeficients
        ϵ3 = 1e-5, #Conergence in squared error
        lup = 9,
        ldown = 11,
        max_iter = 100000, #Maximum iterations
    )

Helper function for `fitRazoQuadModel`. Automatically determins upper and lower bounds and initial guess given the `quad_window_width`. 
"""
function fitRazoQuadModel(
    quad_window_width::R,
    x0_dat::Vector{T}, 
    x1_dat::Vector{T}, 
    yt_dat::Vector{T};
        λ0 = 1e1, #Initial L-M damping coeficient 
        ϵ1 = 1e-6, #Convergence in the gradient
        ϵ2 = 1e-5, #Convergence in the coeficients
        ϵ3 = 1e-6, #Conergence in squared error
        lup = 9,
        ldown = 11,
        max_iter = 100000, #Maximum iterations
    ) where {T<:AbstractFloat, R<:Real}
    al0 = T(quad_window_width)/2
    ar0 = T(quad_window_width)/2
    bl0, br0 = T(5.0), T(5.0)
    al_bounds = (T(0.2), T(quad_window_width))
    ar_bounds = (T(0.2), T(quad_window_width))
    bl_bounds = (T(1e-3), T(100.0))
    br_bounds = (T(1e-3), T(100.0))
    return fitRazoQuadModel(
                al0, ar0, bl0, br0,
                al_bounds, ar_bounds,
                bl_bounds, br_bounds,
                x0_dat, x1_dat, yt_dat;
                λ0 =  λ0, #Initial L-M damping coeficient 
                ϵ1 = ϵ1, #Convergence in the gradient
                ϵ2 = ϵ2, #Convergence in the coeficients
                ϵ3 = ϵ3, #Conergence in squared error
                lup = lup,
                ldown = ldown,
                max_iter = max_iter #Maximum iterations
                )
end

function simmulateQuad(
    rqm::RazoQuadParams{T},
    mz_offset_range::Tuple{U, U},
    center_mz_range::Tuple{U,U},
    n::Int64) where {T,U<:AbstractFloat}

    x0 = Float32[]
    x1 = Float32[]
    yt = Float32[]
    for i in range(1, n)
        #Random window center mz 
        center_mz = Float32(rand(
            Uniform(first(center_mz_range), last(center_mz_range))
        ))
        #Random offset of mz from the center 
        offset = Float32(rand(
            Uniform(first(mz_offset_range), last(mz_offset_range))
            ))

        #Precursor charge state, m0 and m0 m/z, and m0 mass 
        prec_charge_state = rand([2])
        m0_mz = center_mz + offset
        m1_mz =  m0_mz + NEUTRON/prec_charge_state
        mono_mass = Float32(m0_mz*prec_charge_state) #Approximate not accounting for proton 


        s_count = 0 #Assume zero sulfurs.
        #Simmulate true abundances of mono and m+1 precursor ions 
        ζ0 = 100000
        δ = iso_splines(s_count, 0, mono_mass)/iso_splines(s_count, 1, mono_mass)
        ζ1 =  ζ0/δ
        #Observed m0 and m+1 abundances are modified by the quad transmission funciton 
        y0_obs = ζ0*rqm(Float32(m0_mz - center_mz))
        y1_obs = ζ1*rqm(Float32(m1_mz - center_mz))
        #Get offsets 
        push!(x1,m1_mz - center_mz) 
        push!(x0, m0_mz - center_mz)
        #Get observed data
        push!(yt, log.(y0_obs/(y1_obs*δ)))
    end
    return DataFrame((x0 = Float64.(x0), x1 = Float64.(x1), yt = Float64.(yt)))
end