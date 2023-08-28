first_search_params = Dict(
    :collect_frag_errs => true,
    :expected_matches => 1000000,
    :frag_ppm_err => 0.0,
    :fragment_tolerance => 30.0,
    :max_iter => 1000,
    :max_peaks => false,
    :min_frag_count => 7,
    :min_matched_ratio => Float32(0.6),
    :min_spectral_contrast => Float32(0.95),
    :nmf_tol => Float32(100),
    :precursor_tolerance => 5.0,
    :quadrupole_isolation_width => 8.5,
    :regularize => false,
    :rt_bounds => (0.0, 200.0),
    :rt_tol => 200.0,
    :sample_rate => 0.01,
    :topN => 5,
    :λ => zero(Float32),
    :γ => zero(Float32)
)

main_search_params = Dict(
    :expected_matches => 1000000,
    :frag_tol_quantile => 0.975,
    :max_iter => 1000,
    :max_peaks => false,
    :min_frag_count => 4,
    :min_matched_ratio => Float32(0.45),
    :min_spectral_contrast => Float32(0.5),
    :nmf_tol => Float32(100),
    :precursor_tolerance => 5.0,
    :quadrupole_isolation_width => 8.5,
    :regularize => false,
    :rt_bounds =>(-20.0, 200.0),
    :rt_tol => 20.0,
    :topN => 1000,
    :λ => zero(Float32),
    :γ => zero(Float32)
)

integrate_ms2_params = Dict(
    :expected_matches => 1000000,
    :frag_tol_quantile => 0.975,
    :max_iter => 1000,
    :max_peak_width => 2.0,
    :max_peaks => false,
    :min_frag_count => 4,
    :min_matched_ratio => Float32(0.45),
    :min_spectral_contrast => Float32(0.5),
    :nmf_tol => Float32(100),
    :precursor_tolerance => 5.0,
    :quadrupole_isolation_width => 8.5,
    :regularize => false,
    :rt_bounds => (0.0, 200.0),
    :rt_tol => 20.0,
    :sample_rate => 1.0,
    :topN => 1000,
    :λ => zero(Float32),
    :γ => zero(Float32)
)
integrate_ms1_params = Dict(
        :expected_matches => 1000000,
        :frag_tol_quantile => 0.975,
        :max_iter => 1000,
        :max_peak_width => 2.0,
        :max_peaks => false,
        :min_frag_count => 4,
        :min_matched_ratio => Float32(0.45),
        :min_spectral_contrast => Float32(0.5),
        :nmf_tol => Float32(100),
        :precursor_tolerance => 5.0,
        :quadrupole_isolation_width => 8.5,
        :regularize => false,
        :rt_tol => 20.0,
        :rt_bounds => (0.0, 200.0),
        :sample_rate => 1.0,
        :topN => 100,
        :λ => zero(Float32),
        :γ => zero(Float32)
)

N = 10000
z_data = zeros(Float64, (N, N))
x = collect(range(0, stop=2, length=N))
y = collect(range(0, stop=2, length=N))
A = [
    1 0 0 0 0 ;
    1 1 1 0 0 ;
    0 0 1 1 1;
]

function f(X::Tuple{T, T, T}, A::Matrix{Int64}; α::T = 0.9, β::T = 0.01) where {T<:AbstractFloat}
    loss = 0.0
    @turbo for col in 1:size(A)[2]
        ∑X = 0
        for row in 1:size(A)[1]
            #println(row)
            ∑X += X[row]*A[row, col]
        end
        loss += (1- ∑X)^2
    end
    return loss
end

N = 100 - 1
x = collect(range(0, stop=2, length=N))
A = [
    1 0 0 0 0 ;
    1 1 1 0 0 ;
    0 0 1 1 1;
]
function getMin(x::Vector{T}, A::Matrix{Int64}, f::Function) where {T<:AbstractFloat}
    min_idx = (0.0, 0.0, 0.0)
    min = Inf
    for i in x
        for j in x
            for k in x
                z = f((i, j, k), A)
                if z < min
                    min_idx = (i, j, k)
                    min = z
                end
            end
        end
    end
    return min_idx, min
end

A = [
    0 0 0 0 0 0;
    1 1 1 1 1 0;
    0 0 1 1 1 1;
    1 1 0 0 0 0
]
A'\[1, 1, 1, 1, 1, 1]

A = [
    1 0 0 0 0 0 1;
    1 1 1 1 1 0 0;
    0 0 1 1 1 1 0;
    1 1 0 0 0 0 0;
]

A = Float64[
    1 0 0 0 0 0 1;
    1 1 1 1 1 0 0;
    0 0 1 1 1 1 0;
    1 1 0 0 0 0 0;
]
A'\[1, 1, 1, 1, 1, 1, 1]
sparseNMF(sparse(A'), ones(Float64, 7), 0.0, 0.0, false, tol = 1e-6)
A'.==1.0

A = Float64[
    1 0 0 0 0 0 1;
    1 1 1 1 1 0 0;
    0 0 1 1 0 0 0;
    1 1 0 0 0 0 0;
]
A'\[1, 1, 1, 1, 1, 1, 1]
sparseNMF(sparse(A'), ones(Float64, 7), 0.0, 0.0, false, tol = 1e-6)
A'.==1.0

A = Float64[
    1 0 1 0 0 0 1;
    1 1 1 1 1 0 0;
    0 0 1 1 0 0 0;
]
A'\[1, 1, 1, 1, 1, 1, 1]
sparseNMF(sparse(A'), ones(Float64, 7), 0.0, 0.0, false, tol = 1e-6)
A'.==1.0




A = Float64[
    1 0 0 0 0 0 1;
    1 0 1 1 1 0 0;
    0 0 1 1 0 1 0;
    1 1 0 0 0 0 0;
]
A'\[1, 1, 1, 1, 1, 1, 1]
sparseNMF(sparse(A'), ones(Float64, 7), 0.0, 0.0, false, tol = 1e-6)
A'.==1.0

A = Float64[
    1 0 0 0 0 0 1;
    1 0 1 1 1 0 0;
    0 0 1 1 1 1 0;
    1 1 0 0 0 0 0;
]
A'\[1, 1, 1, 1, 1, 1, 1]
sparseNMF(sparse(A'), ones(Float64, 7), 0.0, 0.0, false, tol = 1e-6)
A'.==1.0


A = Float64[
    1 1 0;
    1 0 1;
    0 1 1;
]
A'\[1, 1, 1]
sparseNMF(sparse(A'), ones(Float64, size(A)[2]), 0.0, 0.0, false, tol = 1e-6)
A'.==1.0







A'\[-1, -1, -1, -1, -1, -1]

plot(surface(z=z_data, x=x, y=y))
xy_min = argmin(z_data)
println(x[xy_min[1]])
println(y[xy_min[2]])