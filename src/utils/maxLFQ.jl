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
    getS(peptides::AbstractVector{String}, peptides_dict::Dict{String, Int64}, experiments::AbstractVector{UInt32}, experiments_dict::Dict{UInt32, Int64}, abundance::AbstractVector{Union{T, Missing}}, M::Int, N::Int) where {T<:Real}

Get MxN matrix of intensities of each peptide in each experiment. Rows stand for experiments and columns for peptides. If the j'th peptide 
was not seen in the i'th experiment, then S[i, j] == missing. 

### Input 
- `peptides::AbstractVector{String}` -- MxN List of peptides for the protein. 
- `peptides_dict::Dict{String, Int64}` -- Maps N peptides by name to column numbers in S 
- `experiments::AbstractVector{UInt32}` -- MxN List of experiments 
- `experiments_dict::Dict{UInt32, Int64` -- Maps M experiment IDs by name to rows numbers in S 
- `abundance::AbstractVector{Union{T, Missing}}` -- MxN list of peptide abundances
- `M::Int` -- Number of experiments/rows in S
- `N::Int` -- Number of peptides/columns in S


### Examples 

"""
function getS(
    peptides::AbstractVector{UInt32},
    peptides_dict::Dict{UInt32, Int64},
    experiments::AbstractVector{<:Integer},
    experiments_dict::Dict{E, Int64},
    abundance::AbstractVector{<:Union{Missing, Real}},
    M::Int,
    N::Int
) where {E<:Integer}
    T = Base.nonmissingtype(eltype(abundance))
    S = Matrix{Union{Missing, T}}(missing, M, N)
    for i in eachindex(peptides)
        if !ismissing(abundance[i])
            if abundance[i] > 0.0
                S[peptides_dict[peptides[i]], experiments_dict[experiments[i]]] = abundance[i]
            else
                S[peptides_dict[peptides[i]], experiments_dict[experiments[i]]] = missing
            end
        end
    end
    return S
end

"""
    get_log2_intensity_matrix(S::AbstractMatrix{Union{Missing, T}}) where {T<:Real}

Convert the linear-intensity peptide-by-run matrix `S` into log2 space while
preserving `missing` entries.
"""
function get_log2_intensity_matrix(S::AbstractMatrix{Union{Missing, T}}) where {T<:Real}
    X = Matrix{Union{Missing, Float64}}(missing, size(S)...)
    for j in axes(S, 2)
        for i in axes(S, 1)
            if !ismissing(S[i, j])
                X[i, j] = log2(Float64(S[i, j]))
            end
        end
    end
    return X
end

"""
    get_connected_components(X::AbstractMatrix{Union{Missing, T}}) where {T<:Real}

Return a connected-component label for each run column in `X`. Two runs are
connected when they share at least one quantified peptide.
"""
function get_connected_components(X::AbstractMatrix{Union{Missing, T}}) where {T<:Real}
    n_runs = size(X, 2)
    labels = zeros(Int, n_runs)
    component_id = 0

    for seed_idx in 1:n_runs
        if labels[seed_idx] != 0
            continue
        end

        component_id += 1
        stack = Int[seed_idx]
        labels[seed_idx] = component_id

        while !isempty(stack)
            run_idx = pop!(stack)
            for peptide_idx in axes(X, 1)
                if ismissing(X[peptide_idx, run_idx])
                    continue
                end
                for neighbor_idx in 1:n_runs
                    if !ismissing(X[peptide_idx, neighbor_idx]) && labels[neighbor_idx] == 0
                        labels[neighbor_idx] = component_id
                        push!(stack, neighbor_idx)
                    end
                end
            end
        end
    end

    return labels
end

"""
    get_component_sizes(labels::AbstractVector{Int})

Return the size of the connected component for each element in `labels`.
"""
function get_component_sizes(labels::AbstractVector{Int})
    counts = Dict{Int, Int}()
    for label in labels
        counts[label] = get(counts, label, 0) + 1
    end
    return [counts[label] for label in labels]
end

"""
    get_component_priority(component_indices, run_priorities)

Return the summed per-run priority for the provided component indices, ignoring
`missing` priorities.
"""
function get_component_priority(
    component_indices::AbstractVector{Int},
    run_priorities::AbstractVector{Union{Missing, T}}
) where {T<:Real}
    total = 0.0
    for run_idx in component_indices
        priority = run_priorities[run_idx]
        if !ismissing(priority)
            total += Float64(priority)
        end
    end
    return total
end

"""
    get_best_component(component_labels, run_priorities)

Select the component with the most runs. Break ties by the highest summed
per-run priority.
"""
function get_best_component(
    component_labels::AbstractVector{Int},
    run_priorities::AbstractVector{Union{Missing, T}}
) where {T<:Real}
    best_component_id = 0
    best_size = -1
    best_priority = -Inf

    for component_id in 1:maximum(component_labels)
        component_indices = findall(==(component_id), component_labels)
        component_size = length(component_indices)
        component_priority = get_component_priority(component_indices, run_priorities)

        if component_size > best_size ||
           (component_size == best_size && component_priority > best_priority)
            best_component_id = component_id
            best_size = component_size
            best_priority = component_priority
        end
    end

    return best_component_id
end

"""
    getB(X::AbstractMatrix{Union{Missing, T}}) where {T<:Real}

Build the MaxLFQ normal-equation right-hand side using only valid run pairs that
share at least one quantified peptide. `X` must already be in log2 space.

* Yu F, Haynes SE, Nesvizhskii AI. IonQuant Enables Accurate and Sensitive Label-Free Quantification With FDR-Controlled Match-Between-Runs. Mol Cell Proteomics. 2021;20:100077. doi: 10.1016/j.mcpro.2021.100077. Epub 2021 Apr 2. PMID: 33813065; PMCID: PMC8131922.
* https://rdrr.io/cran/iq/man/maxLFQ.html
* Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep;13(9):2513-26. doi: 10.1074/mcp.M113.031591. Epub 2014 Jun 17. PMID: 24942700; PMCID: PMC4159666.

"""
function getB(X::AbstractMatrix{Union{Missing, T}}) where {T<:Real}
    n_runs = size(X, 2)
    Atb = zeros(Float64, n_runs)
    valid_pairs = 0

    for left_idx in 1:(n_runs - 1)
        for right_idx in (left_idx + 1):n_runs
            ratios = Float64[]
            for peptide_idx in axes(X, 1)
                left_val = X[peptide_idx, left_idx]
                right_val = X[peptide_idx, right_idx]
                if !ismissing(left_val) && !ismissing(right_val)
                    push!(ratios, -Float64(left_val) + Float64(right_val))
                end
            end

            if isempty(ratios)
                continue
            end

            ratio_median = median(ratios)
            Atb[left_idx] -= ratio_median
            Atb[right_idx] += ratio_median
            valid_pairs += 1
        end
    end

    return Atb, valid_pairs
end

"""
    getA(X::AbstractMatrix{Union{Missing, T}}) where {T<:Real}

Build the MaxLFQ normal-equation left-hand side using only valid run pairs that
share at least one quantified peptide. `X` must already be in log2 space.

    * Yu F, Haynes SE, Nesvizhskii AI. IonQuant Enables Accurate and Sensitive Label-Free Quantification With FDR-Controlled Match-Between-Runs. Mol Cell Proteomics. 2021;20:100077. doi: 10.1016/j.mcpro.2021.100077. Epub 2021 Apr 2. PMID: 33813065; PMCID: PMC8131922.
    * https://rdrr.io/cran/iq/man/maxLFQ.html
    * Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep;13(9):2513-26. doi: 10.1074/mcp.M113.031591. Epub 2014 Jun 17. PMID: 24942700; PMCID: PMC4159666.
    
"""
function getA(X::AbstractMatrix{Union{Missing, T}}) where {T<:Real}
    n_runs = size(X, 2)
    AtA = zeros(Float64, n_runs, n_runs)

    for left_idx in 1:(n_runs - 1)
        for right_idx in (left_idx + 1):n_runs
            has_overlap = false
            for peptide_idx in axes(X, 1)
                if !ismissing(X[peptide_idx, left_idx]) && !ismissing(X[peptide_idx, right_idx])
                    has_overlap = true
                    break
                end
            end

            if has_overlap
                AtA[left_idx, right_idx] = -1.0
                AtA[right_idx, left_idx] = -1.0
                AtA[left_idx, left_idx] += 1.0
                AtA[right_idx, right_idx] += 1.0
            end
        end
    end

    return AtA
end

"""
    solve_maxlfq_component(X::AbstractMatrix{Union{Missing, T}}) where {T<:Real}

Solve the MaxLFQ system for a single connected run component in log2 space.
"""
function solve_maxlfq_component(X::AbstractMatrix{Union{Missing, T}}) where {T<:Real}
    n_runs = size(X, 2)
    n_runs == 0 && return Union{Missing, Float32}[]

    observed_values = skipmissing(vec(X))
    if isempty(observed_values)
        return Vector{Union{Missing, Float32}}(missing, n_runs)
    end

    AtA = getA(X)
    Atb, valid_pairs = getB(X)
    if valid_pairs == 0
        return Vector{Union{Missing, Float32}}(missing, n_runs)
    end

    system_matrix = Matrix{Float64}(undef, n_runs + 1, n_runs + 1)
    system_matrix[1:n_runs, 1:n_runs] = 2.0 .* AtA
    system_matrix[1:n_runs, end] .= 1.0
    system_matrix[end, 1:n_runs] .= 1.0
    system_matrix[end, end] = 0.0

    rhs = Vector{Float64}(undef, n_runs + 1)
    rhs[1:n_runs] = 2.0 .* Atb
    rhs[end] = mean(observed_values) * n_runs

    solution = system_matrix \ rhs
    estimates = Vector{Union{Missing, Float32}}(missing, n_runs)
    for run_idx in 1:n_runs
        if isfinite(solution[run_idx])
            estimates[run_idx] = solution[run_idx]
        end
    end

    return estimates
end

"""
    solve_maxlfq(
        X::AbstractMatrix{Union{Missing, T}},
        run_priorities::AbstractVector{Union{Missing, F}}
    ) where {T<:Real, F<:Real}

Solve MaxLFQ in log2 space using only the single best connected component,
leaving all other runs as `missing`. The best component is the one with the most
runs, with ties broken by the highest summed per-run priority.
"""
function solve_maxlfq(
    X::AbstractMatrix{Union{Missing, T}},
    run_priorities::AbstractVector{Union{Missing, F}}
) where {T<:Real, F<:Real}
    n_runs = size(X, 2)
    estimates = Vector{Union{Missing, Float32}}(missing, n_runs)
    n_runs == 0 && return estimates, Int[]
    @assert length(run_priorities) == n_runs "run priorities must align with MaxLFQ run columns"

    component_labels = get_connected_components(X)
    component_sizes = get_component_sizes(component_labels)
    best_component_id = get_best_component(component_labels, run_priorities)

    component_indices = findall(==(best_component_id), component_labels)
    if length(component_indices) > 1
        component_estimates = solve_maxlfq_component(@view X[:, component_indices])
        for (local_idx, global_idx) in enumerate(component_indices)
            estimates[global_idx] = component_estimates[local_idx]
        end
    end

    if maximum(component_sizes; init = 0) < n_runs
        @debug "MaxLFQ protein graph disconnected" component_sizes=component_sizes best_component_id=best_component_id
    end

    return estimates, component_labels
end

"""
    get_linear_abundance(log2_abundance::AbstractVector{Union{Missing, T}}) where {T<:Real}

Convert MaxLFQ log2 abundances back to linear space while preserving a nullable
column type for downstream Arrow writes.
"""
function get_linear_abundance(log2_abundance::AbstractVector{Union{Missing, T}}) where {T<:Real}
    abundance = Vector{Union{Missing, Float32}}(missing, length(log2_abundance))
    for idx in eachindex(log2_abundance)
        if !ismissing(log2_abundance[idx])
            abundance[idx] = exp2(log2_abundance[idx])
        end
    end
    return abundance
end

"""
    get_total_peak_area(S::AbstractMatrix{Union{Missing, T}}) where {T<:Real}

Sum the quantified precursor peak areas for each run column in `S`, preserving
`missing` when a run has no quantified precursors.
"""
function get_total_peak_area(S::AbstractMatrix{Union{Missing, T}}) where {T<:Real}
    total_peak_area = Vector{Union{Missing, Float32}}(missing, size(S, 2))
    for run_idx in axes(S, 2)
        total = 0.0
        has_signal = false
        for precursor_idx in axes(S, 1)
            value = S[precursor_idx, run_idx]
            if !ismissing(value)
                total += Float64(value)
                has_signal = true
            end
        end
        if has_signal
            total_peak_area[run_idx] = total
        end
    end
    return total_peak_area
end


"""
    getProtAbundance(protein::String, peptides::AbstractVector{String}, experiments::AbstractVector{UInt32}, abundance::AbstractVector{Union{T, Missing}},
                            protein_out::Vector{String}, peptides_out::Vector{String}, experiments_out::Vector{UInt32}, log2_abundance_out::Vector{Float64}) where {T <: Real}

Estimates protein-level abundances from peptide level abundances using the MaxLFQ algorithm 
* Yu F, Haynes SE, Nesvizhskii AI. IonQuant Enables Accurate and Sensitive Label-Free Quantification With FDR-Controlled Match-Between-Runs. Mol Cell Proteomics. 2021;20:100077. doi: 10.1016/j.mcpro.2021.100077. Epub 2021 Apr 2. PMID: 33813065; PMCID: PMC8131922.
* https://rdrr.io/cran/iq/man/maxLFQ.html
* Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep;13(9):2513-26. doi: 10.1074/mcp.M113.031591. Epub 2014 Jun 17. PMID: 24942700; PMCID: PMC4159666.

### Input
`peptides`, `experiments`, and `abundance` each have length `N`
- `protein::String` -- Protien Name
- `peptides::AbstractVector{String}` -- Peptide names. Same length as `experiments` and `abundance`
- `experiments::AbstractVector{UInt32}` -- MS file id's. Same length as `peptides` and `abundance`
- `abundance::AbstractVector{Union{T, Missing}}` -- Abundances of the peptides. Same lenght as `peptides` and `experiments`
- `protein_out::Vector{String}` -- The protein name. Repeated `N` times for each experiment
- `peptides_out::Vector{String}` -- List of lists of detected peptides for each of the `N` experiments. [PEPA;PEPB PEPA;PEPC ...]
- `experiments_out::Vector{UInt32}` -- List of expeiments 
- `log2_abundance_out::Vector{Float64}` -- List of log2 estimated protein abundances for each experiment. 

### Examples. 
julia> using DataFrames

julia> prot = DataFrame(Dict(
           :peptides => ["A","A","A","B","B","B","C","C","C","D","D","D"],
           :protein => append!(split(repeat("A",9), ""), ["B","B","B"]),
           :file_idx => UInt32[1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3],
           :abundance => [10, 20, 40, 1, 2, 4, 100, 200, missing, 1000, 2000, 3000],
       ))

12x4 DataFrame
 Row │ abundance  file_idx  peptides  protein   
     │ Int64?     UInt32    String    SubStrin… 
─────┼──────────────────────────────────────────
   1 │        10         1  A         A
   2 │        20         2  A         A
   3 │        40         3  A         A
   4 │         1         1  B         A
   5 │         2         2  B         A
   6 │         4         3  B         A
   7 │       100         1  C         A
   8 │       200         2  C         A
   9 │   missing         3  C         A
  10 │      1000         1  D         B
  11 │      2000         2  D         B
  12 │      3000         3  D         B

julia> function testLFQ(prot)
           out = Dict(
               :protein => String[],
               :peptides => String[],
               :log2_abundance => Float64[],
               :experiments => UInt32[],
           )

           for (protein, data) in pairs(groupby(prot, :protein))
               getProtAbundance(string(protein[:protein]), 
                                   collect(data[!,:peptides]), 
                                   collect(data[!,:file_idx]), 
                                   collect(data[!,:abundance]),
                                   out[:protein],
                                   out[:peptides],
                                   out[:experiments],
                                   out[:log2_abundance]
                               )
           end
           out
       end
testLFQ (generic function with 1 method)

julia> out = testLFQ(prot)
Dict{Symbol, Vector} with 4 entries:
  :protein        => ["A", "A", "A", "B", "B", "B"]
  :peptides       => ["A;B;C", "A;B;C", "A;B", "D", "D", "D"]
  :log2_abundance => [3.15526, 4.15526, 5.15526, 9.96578, 10.9658, 11.5507]
  :experiments    => UInt32[0x00000001, 0x00000002, 0x00000003, 0x00000001, 0x00000002, 0x00000003]

"""
function getProtAbundance(protein::String, 
                            row_idx::Int64,
                            target::Bool,
                            entrap_id::UInt8,
                            species::String,
                            peptides::AbstractVector{UInt32}, 
                            experiments::AbstractVector{<:Integer},
                            use_for_quant::AbstractVector{Bool},
                            abundance::AbstractVector{<:Union{Missing, Real}},
                            global_qvals::AbstractVector{<:Union{Missing, Real}},
                            qvals::AbstractVector{<:Union{Missing, Real}},
                            peps::AbstractVector{<:Union{Missing, Real}},
                            pg_scores::AbstractVector{<:Union{Missing, Real}},
                            global_pg_scores::AbstractVector{<:Union{Missing, Real}},
                            target_out::Vector{Union{Missing, Bool}},
                            entrap_id_out::Vector{Union{Missing, UInt8}},
                            species_out::Vector{Union{Missing, String}},
                            protein_out::Vector{Union{Missing, String}}, 
                            peptides_out::Vector{Union{Missing, Vector{Union{Missing, UInt32}}}},
                            experiments_out::Vector{Union{Missing, UInt32}}, 
                            log2_abundance_out::Vector{Union{Missing, Float32}},
                            global_qval_out::Vector{Union{Missing, Float32}},
                            qval_out::Vector{Union{Missing, Float32}},
                            pep_out::Vector{Union{Missing, Float32}},
                            pg_score_out::Vector{Union{Missing, Float32}},
                            global_pg_score_out::Vector{Union{Missing, Float32}},
                            total_peak_area_out::Vector{Union{Missing, Float32}})

    unique_experiments = unique(experiments)
    unique_peptides = unique(peptides)

    N = length(unique_experiments)
    M = length(unique_peptides)

    peptides_dict = Dict(zip(unique_peptides, 1:M))
    experiments_dict = Dict(zip(unique_experiments, 1:N))

    pep_map = Dict{UInt32, Union{Missing, Float32}}()
    score_map = Dict{UInt32, Union{Missing, Float32}}()
    global_score_map = Dict{UInt32, Union{Missing, Float32}}()
    qval_map = Dict{UInt32, Union{Missing, Float32}}()
    global_qval_map = Dict{UInt32, Union{Missing, Float32}}()
    for j in eachindex(experiments)
        exp = UInt32(experiments[j])
        if !haskey(pep_map, exp)
            pep_map[exp] = peps[j]
            score_map[exp] = pg_scores[j]
            global_score_map[exp] = global_pg_scores[j]
            qval_map[exp] = qvals[j]
            global_qval_map[exp] = global_qvals[j]
        end
    end

    run_pg_scores = Union{Missing, Float32}[get(score_map, UInt32(exp), missing) for exp in unique_experiments]

    #Appends the results to the inputs `protein_out`, `peptides_out`, `experiments_out` and `log2_abundance_out`
    function appendResults!(row_idx::Int64,
                            N::Int64,
                            target::Bool,
                            entrap_id::UInt8,
                            species::String,
                            protein::String,
                            unique_peptides::Vector{UInt32}, 
                            log2_abundances::AbstractVector{Union{Missing, A}},
                            global_qvals::AbstractVector{<:Union{Missing, Real}},
                            qvals::AbstractVector{<:Union{Missing, Real}},
                            peps::AbstractVector{<:Union{Missing, Real}},
                            pg_scores::AbstractVector{<:Union{Missing, Real}},
                            global_pg_scores::AbstractVector{<:Union{Missing, Real}},
                            target_out::Vector{Union{Missing, Bool}},
                            entrap_id_out::Vector{Union{Missing, UInt8}},
                            species_out::Vector{Union{Missing, String}}, 
                            protein_out::Vector{Union{Missing, String}}, 
                            peptides_out::Vector{Union{Missing, Vector{Union{Missing, UInt32}}}}, 
                            log2_abundance_out::Vector{Union{Missing, Float32}}, 
                            global_qval_out::Vector{Union{Missing, Float32}}, 
                            qval_out::Vector{Union{Missing, Float32}}, 
                            pep_out::Vector{Union{Missing, Float32}},
                            pg_score_out::Vector{Union{Missing, Float32}},
                            global_pg_score_out::Vector{Union{Missing, Float32}},
                            total_peak_area_out::Vector{Union{Missing, Float32}},
                            experiments_out::Vector{Union{Missing, I}}, 
                            S::Matrix{Union{Missing, S_T}},
                            total_peak_area::AbstractVector{Union{Missing, Float32}}) where {A<:Real,S_T<:Real,I<:Integer}
        
        function appendPeptides!(peptides_out::Vector{Union{Missing, Vector{Union{Missing, UInt32}}}}, 
                                row_idx::Int64,
                                unique_peptides::Vector{UInt32}, 
                                S::Matrix{Union{Missing,S_T}})
            #Each column of S corresponds to and experiment and each row corresponds to a peptide
            #Need to get each non-missing peptide for each experiment. Concatenate the non-missing peptides
            #For and experiment with a semi-colon. See example in the getProtAbundance docstring 
            for j in eachindex(eachcol(S))
                #Each row in S is for a peptide 
                sample_peptides = Vector{Union{Missing, UInt32}}(undef, size(S, 1))
                for i in eachindex(@view(S[:,j]))
                    if !ismissing(S[i,j])
                        sample_peptides[i] =  unique_peptides[i]
                    else
                        sample_peptides[i] = missing
                    end
                end
                peptides_out[row_idx + j - 1] = sample_peptides
            end
        end


        for i in range(0, N-1)
            exp = UInt32(unique_experiments[i + 1])
            log2_abundance_out[row_idx + i] = log2_abundances[i + 1]
            global_qval_out[row_idx + i] = global_qval_map[exp]
            qval_out[row_idx + i] = qval_map[exp]
            pg_score_out[row_idx + i] = score_map[exp]
            global_pg_score_out[row_idx + i] = global_score_map[exp]
            experiments_out[row_idx + i] = exp
            protein_out[row_idx + i] = protein
            pep_out[row_idx + i] = pep_map[exp]
            species_out[row_idx + i] = species
            target_out[row_idx + i] = target
            entrap_id_out[row_idx + i] = entrap_id
            total_peak_area_out[row_idx + i] = total_peak_area[i + 1]
        end
        appendPeptides!(peptides_out, 
                        row_idx,
                        unique_peptides, 
                        S)
    end

    #ixj matrix where rows are for experiments and columns are for peptides. Each entry is the abundance of the peptide
    #in the given experiment, or missing if peptide j was not seen in experiment i. 
    S = getS(peptides, peptides_dict, experiments, experiments_dict, abundance, M, N)
    X = get_log2_intensity_matrix(S)
    total_peak_area = get_total_peak_area(S)
    log2_abundances, component_labels = solve_maxlfq(X, run_pg_scores)
    quantified_runs = findall(x -> !ismissing(x), total_peak_area)
    if length(quantified_runs) == 1
        run_idx = only(quantified_runs)
        log2_abundances[run_idx] = log2(total_peak_area[run_idx])
    end
    
    # Debug: Check if all abundances are missing/NaN/Inf
    n_valid_abundances = sum(!ismissing(x) && isfinite(x) for x in log2_abundances)
    if n_valid_abundances == 0
        @debug "MaxLFQ produced all invalid abundances for protein" protein=protein n_peptides=M n_experiments=N
    elseif maximum(component_labels; init = 0) > 1
        @debug "MaxLFQ protein solved with disconnected components" protein=protein component_labels=component_labels
    end
    
    appendResults!(
                   row_idx,
                   N,
                   target,
                   entrap_id,
                   species, 
                   protein, 
                   unique_peptides, 
                   log2_abundances,
                   global_qvals,
                   qvals,
                   peps,
                   pg_scores,
                   global_pg_scores,
                   target_out,
                   entrap_id_out,
                   species_out,
                   protein_out, 
                   peptides_out, 
                   log2_abundance_out,
                   global_qval_out,
                   qval_out,
                   pep_out,
                   pg_score_out,
                   global_pg_score_out,
                   total_peak_area_out,
                   experiments_out,
                   S,
                   total_peak_area)

end

"""
    build_accession_to_species(precursors)

Build a `Dict{String, String}` mapping each individual UniProt accession (single
token, never `;`-joined) to the species/proteome label of the FASTA it was read
from. Both `accession_numbers` and `proteome_identifiers` on a precursor are
`;`-joined parallel arrays — same order, same length — so splitting on `;` and
zipping recovers per-accession species.
"""
function build_accession_to_species(precursors)
    accs = getAccessionNumbers(precursors)
    proteomes = getProteomeIdentifiers(precursors)
    accession_to_species = Dict{String, String}()
    # 6.8M precursors map to only ~tens of thousands of distinct accession strings,
    # so memoize processed strings and skip duplicates wholesale — avoids re-splitting
    # and re-interning for the ~99.5% of rows that repeat. Bit-identical: a repeated
    # accession string yields the same tokens, and get!/first-wins means re-processing
    # cannot change the dict. `haskey` on a SubString is allocation-free, so tokens are
    # interned only on genuine insert. `seen` is updated only after the length guard
    # passes, so a row first seen with a mismatched proteome string can still be
    # processed later with a matching one.
    seen = Set{String}()
    @inbounds for i in eachindex(accs, proteomes)
        a_str = accs[i]
        (a_str in seen) && continue
        p_str = proteomes[i]
        a_toks = split(a_str, ';')
        p_toks = split(p_str, ';')
        length(a_toks) == length(p_toks) || continue
        push!(seen, String(a_str))
        for k in eachindex(a_toks)
            ak = a_toks[k]
            isempty(ak) && continue
            haskey(accession_to_species, ak) && continue
            accession_to_species[String(ak)] = String(p_toks[k])
        end
    end
    return accession_to_species
end

# FileReference-based implementation with TransformPipeline preprocessing
function LFQ(prot_ref,  # PSMFileReference - using Any to avoid dependency issues
             protein_quant_path::String,
            quant_col::Symbol,
            file_id_to_parsed_name::Vector{String},
            precursor_sequences::AbstractVector{<:AbstractString},
            precursor_structural_mods::AbstractVector,
            precursor_isotopic_mods::AbstractVector,
            q_value_threshold::Float32,
            accession_to_species::Dict{String, String};
            output_schema_policy::OutputSchemaPolicy = OutputSchemaPolicy(),
            batch_size = 100000,
            writer_ref::Base.RefValue{Union{Nothing, Arrow.Writer}} = Ref{Union{Nothing, Arrow.Writer}}(nothing))
    
    # Use eager DataFrame loading (allows editing for filtering)
    prot = DataFrame(Tables.columntable(Arrow.Table(file_path(prot_ref))))

    # Filter out rows with missing inferred_protein_group values
    filter!(:inferred_protein_group => x -> !ismissing(x), prot)

    # Create pipeline operations for batch-wise application
    preprocessing_pipeline = TransformPipeline() |>
        filter_by_multiple_thresholds([
            (:pg_qval, q_value_threshold),
            (:global_pg_qval, q_value_threshold)
        ]) |>
        filter_rows(row -> row.use_for_protein_quant; desc="filter_for_protein_quant")
    
    # Main processing logic (inlined from original LFQ function)
    nfiles = length(unique(prot[!, :ms_file_idx]))
    batch_start_idx, batch_end_idx = 1, min(batch_size, size(prot, 1))
    @assert issorted(prot[!, :inferred_protein_group], rev=true) "LFQ input must be sorted by :inferred_protein_group (descending)"

    while batch_start_idx <= size(prot, 1)
        last_prot_idx = prot[batch_end_idx, :inferred_protein_group]
        while batch_end_idx < size(prot, 1)
            if prot[batch_end_idx+1, :inferred_protein_group] != last_prot_idx
                break
            end
            batch_end_idx += 1
        end
        subdf = prot[range(batch_start_idx, batch_end_idx), :]
        batch_start_idx = batch_end_idx + 1
        batch_end_idx = min(batch_start_idx + batch_size, size(prot, 1))
        
        # Apply pipeline operations to this batch
        initial_rows = nrow(subdf)
        for (desc, op) in preprocessing_pipeline.operations
            subdf = op(subdf)
            @debug "Pipeline operation" operation=desc rows_before=initial_rows rows_after=nrow(subdf)
            initial_rows = nrow(subdf)
        end
        
        if nrow(subdf) == 0
            continue  # Skip empty batches
        end
        
        # Continue with original LFQ logic
        #Exclude precursors with mods that impact quantitation
        #filter!(x->!occursin("M,Unimod:35", coalesce(x.structural_mods, "")), subdf)
        # Group by protein-group identity. :species is intentionally NOT a key
        # — peptides shared across FASTAs (e.g. a peptide matching both the
        # contaminant and human copies of an accession) carry a multi-species
        # string, while sibling peptides matching only one source carry the
        # single-species string. Keying on :species would split a single
        # inferred protein group into per-species rows. Aggregate the species
        # union below instead.
        gpsms = groupby(
            subdf,
            [:target, :entrapment_group_id, :inferred_protein_group]
        )
        ngroups = length(gpsms)
        nrows = nfiles*ngroups
        
        # Pre-allocate the batch with missing values
        out = Dict(
            :target => Vector{Union{Missing, Bool}}(missing, nrows),
            :entrap_id => Vector{Union{Missing, UInt8}}(missing, nrows),
            :species => Vector{Union{Missing, String}}(missing, nrows),
            :protein => Vector{Union{Missing, String}}(missing, nrows),
            :peptides => Vector{Union{Missing, Vector{Union{Missing, UInt32}}}}(missing, nrows),
            :log2_abundance => Vector{Union{Missing, Float32}}(missing, nrows),
            :experiments => Vector{Union{Missing, UInt32}}(missing, nrows),
            :global_qval => Vector{Union{Missing, Float32}}(missing, nrows),
            :qval => Vector{Union{Missing, Float32}}(missing, nrows),
            :pg_pep => Vector{Union{Missing, Float32}}(missing, nrows),
            :pg_score => Vector{Union{Missing, Float32}}(missing, nrows),
            :global_pg_score => Vector{Union{Missing, Float32}}(missing, nrows),
            :total_peak_area => Vector{Union{Missing, Float32}}(missing, nrows),
        )

        for (group_idx, (protein, data)) in enumerate(pairs(gpsms))
            # Species is derived from the accessions inside the inferred protein
            # group, not from the per-PSM :species column. Per-PSM species can
            # include FASTA sources of *other* proteins that share peptides with
            # this group (e.g. a yeast protein whose peptide is also found in a
            # human homolog), which would mislabel the group's organism. The
            # inferred_protein_group is the parsimonious accession set; its
            # species union is the species union of those accessions only.
            species_set = Set{String}()
            for acc in split(String(protein[:inferred_protein_group]), ';')
                isempty(acc) && continue
                sp = get(accession_to_species, String(acc), "")
                isempty(sp) && continue
                push!(species_set, sp)
            end
            species_agg = join(sort!(collect(species_set)), ';')

            getProtAbundance(protein[:inferred_protein_group],
                                (group_idx*nfiles) - nfiles + 1,
                                protein[:target],
                                protein[:entrapment_group_id],
                                species_agg,
                                data[!,:precursor_idx], 
                                data[!,:ms_file_idx], 
                                data[!,:use_for_protein_quant],
                                data[!,quant_col],
                                data[!,:global_pg_qval],
                                data[!,:pg_qval],
                                data[!,:pg_pep],
                                data[!,:pg_score],
                                data[!,:global_pg_score],
                                out[:target],
                                out[:entrap_id],
                                out[:species],
                                out[:protein],
                                out[:peptides],
                                out[:experiments],
                                out[:log2_abundance],
                                out[:global_qval],
                                out[:qval],
                                out[:pg_pep],
                                out[:pg_score],
                                out[:global_pg_score],
                                out[:total_peak_area]
                            )
        end
        out = DataFrame(out)
        n_precursors, n_modified_peptides, n_peptides = countProteinSupport(
            out[!,:peptides],
            precursor_sequences,
            precursor_structural_mods,
            precursor_isotopic_mods
        )
        out[!,:n_precursors] = n_precursors
        out[!,:n_modified_peptides] = n_modified_peptides
        out[!,:n_peptides] = n_peptides
        out[!,:file_name] = Vector{Union{Missing, String}}(undef, size(out, 1))
        for i in range(1, size(out, 1))
            if ismissing(out[i,:experiments])
                out[i,:file_name] = missing
            else
                out[i,:file_name] = file_id_to_parsed_name[out[i,:experiments]]
            end
        end
        filter!(x->(!ismissing(x.n_peptides)), out);#&(x.n_peptides>=min_peptides), out);
        out[!,:abundance] = get_linear_abundance(out[!,:log2_abundance])
        output_cols = enabled_output_columns(output_schema_policy, :protein_groups, Symbol[
            :file_name,
            :target,
            :entrap_id,
            :species,
            :protein,
            :peptides,
            :n_precursors,
            :n_modified_peptides,
            :n_peptides,
            :global_qval,
            :qval,
            :pg_pep,
            :pg_score,
            :global_pg_score,
            :total_peak_area,
            :abundance,
        ])

        if isnothing(writer_ref[])
            isfile(protein_quant_path) && rm(protein_quant_path)
            writer_ref[] = open(Arrow.Writer, protein_quant_path; file=true)
        end
        Arrow.write(
            writer_ref[],
            select!(
                out,
                output_cols)
        )
    end
    
    # Check if we processed all rows
    if batch_start_idx <= size(prot, 1)
        @user_warn "Not all rows processed!" last_processed_idx=batch_end_idx total_rows=size(prot, 1) unprocessed_rows=size(prot, 1)-batch_end_idx
        throw(ErrorException("Not all rows processed! last_processed_idx=$batch_end_idx, total_rows=$(size(prot, 1)), unprocessed_rows=$(size(prot, 1)-batch_end_idx)"))
    end
    
    return nothing
end


"""
    LFQ_chunked(chunk_refs, protein_quant_path, quant_col, file_id_to_parsed_name, q_value_threshold; batch_size)

Process a vector of chunk FileReferences independently for MaxLFQ protein quantification.
Each chunk is loaded one at a time (~1 GB), processed, and results are appended to a single
output Arrow file.  Memory is bounded by the chunk size rather than total dataset size.
"""
function LFQ_chunked(
    chunk_refs::Vector{<:Any},  # Vector of FileReferences
    protein_quant_path::String,
    quant_col::Symbol,
    file_id_to_parsed_name::Vector{String},
    precursor_sequences::AbstractVector{<:AbstractString},
    precursor_structural_mods::AbstractVector,
    precursor_isotopic_mods::AbstractVector,
    q_value_threshold::Float32,
    accession_to_species::Dict{String, String};
    output_schema_policy::OutputSchemaPolicy = OutputSchemaPolicy(),
    batch_size::Int = 100000
)
    n_chunks = length(chunk_refs)
    # Skip progress bar when there's only one chunk — ProgressBars shows
    # "Inf:Inf, InfGs/it" on n=1 (rate divide-by-zero).
    pbar = n_chunks > 1 ? ProgressBar(total=n_chunks) : nothing
    pbar !== nothing && set_description(pbar, "MaxLFQ chunks:")
    writer_ref = Ref{Union{Nothing, Arrow.Writer}}(nothing)
    try
        for chunk_ref in chunk_refs
            LFQ(chunk_ref, protein_quant_path, quant_col,
                file_id_to_parsed_name,
                precursor_sequences,
                precursor_structural_mods,
                precursor_isotopic_mods,
                q_value_threshold,
                accession_to_species;
                output_schema_policy=output_schema_policy,
                batch_size=batch_size, writer_ref=writer_ref)
            pbar !== nothing && update(pbar)
        end
    finally
        if !isnothing(writer_ref[])
            close(writer_ref[])
        end
    end
    return nothing
end

function countProteinSupport(
    peptides::Vector{Union{Missing, Vector{Union{Missing, UInt32}}}},
    precursor_sequences::AbstractVector{<:AbstractString},
    precursor_structural_mods::AbstractVector,
    precursor_isotopic_mods::AbstractVector
)
    n_rows = length(peptides)
    n_precursors = zeros(Union{Missing, UInt32}, n_rows)
    n_modified_peptides = zeros(Union{Missing, UInt32}, n_rows)
    n_peptides = zeros(Union{Missing, UInt32}, n_rows)

    for (i, pep) in enumerate(peptides)
        if ismissing(pep)
            n_precursors[i] = missing
            n_modified_peptides[i] = missing
            n_peptides[i] = missing
            continue
        end

        n = 0
        observed_sequences = Set{String}()
        observed_modified = Set{Tuple{String, String, String}}()
        for p in pep
            if !ismissing(p)
                pid = Int(p)
                n += 1
                sequence = String(precursor_sequences[pid])
                structural_mods = ismissing(precursor_structural_mods[pid]) ? "" : String(precursor_structural_mods[pid])
                isotopic_mods = ismissing(precursor_isotopic_mods[pid]) ? "" : String(precursor_isotopic_mods[pid])
                push!(observed_sequences, sequence)
                push!(observed_modified, (sequence, structural_mods, isotopic_mods))
            end
        end
        n_precursors[i] = UInt32(n)
        n_modified_peptides[i] = UInt32(length(observed_modified))
        n_peptides[i] = UInt32(length(observed_sequences))
    end
    return n_precursors, n_modified_peptides, n_peptides
end
