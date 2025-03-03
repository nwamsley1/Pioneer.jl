Hs = load("/Users/nathanwamsley/Desktop/test_hs.jld2")["Hs"]
N = Hs.n_vals
H = Matrix(sparse(Hs.rowval[1:N], Hs.colval[1:N], Hs.nzval[1:N]))


using LinearAlgebra
using Statistics

# 1. Calculate effective rank using Singular Value Decomposition (SVD)
svd_decomp = svd(H)
singular_values = svd_decomp.S

# Visualize the decay of singular values
using Plots
plot(singular_values, yscale=:log10, 
     title="Singular Values (log scale)", 
     xlabel="Index", ylabel="Value",
     marker=:circle)

# Calculate effective rank using a threshold
# Common approach: count values above some fraction of maximum
threshold = 1e-10 * maximum(singular_values)
effective_rank = count(s -> s > threshold, singular_values)
println("Effective rank with threshold: ", effective_rank)


# 2. Check for similar or identical columns
n_cols = size(H, 2)
similar_columns = []

# Compute column similarity using correlation
cor_matrix = cor(H)

# Find highly correlated columns (close to 1 or -1)
similarity_threshold = 0.99
for i in 1:n_cols
    for j in (i+1):n_cols
        if abs(cor_matrix[i, j]) > similarity_threshold
            push!(similar_columns, (i, j, cor_matrix[i, j]))
        end
    end
end

# Check for nearly identical columns using L2 norm
norm_threshold = 1e-10
for i in 1:n_cols
    for j in (i+1):n_cols
        if norm(H[:, i] - H[:, j]) < norm_threshold
            println("Columns $i and $j appear to be identical")
        elseif norm(H[:, i] + H[:, j]) < norm_threshold
            println("Columns $i and $j appear to be negatives of each other")
        end
    end
end

# Report results
if isempty(similar_columns)
    println("No highly similar columns found.")
else
    println("Similar column pairs (i, j, correlation):")
    for (i, j, corr) in similar_columns
        println("  Columns $i and $j: correlation = $corr")
    end
end



#What is happening with the similair columns. Most likely it is different plexes
#of the same precursor 
sub_H = H[:,[103, 104, 105, 106, 112]]
# Method 3: Using findall for more explicit indexing
nonzero_indices = findall(i -> !all(iszero, sub_H[i,:]), 1:size(sub_H,1))
sub_H_nonzero = sub_H[nonzero_indices, :]

sub_H = H[:,[569, 570, 571]]
# Method 3: Using findall for more explicit indexing
nonzero_indices = findall(i -> !all(iszero, sub_H[i,:]), 1:size(sub_H,1))
sub_H_nonzero = sub_H[nonzero_indices, :]
