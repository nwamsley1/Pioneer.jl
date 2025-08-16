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
function getS(peptides::AbstractVector{UInt32}, peptides_dict::Dict{UInt32, Int64}, experiments::AbstractVector{UInt16}, experiments_dict::Dict{UInt16, Int64}, abundance::AbstractVector{Union{T, Missing}}, M::Int, N::Int) where {T<:Real}
    #Initialize
    S = Array{Union{Missing,T}}(undef, (M, N))
    for i in eachindex(peptides)
            if !ismissing(abundance[i]) #Abundance of the the peptide 
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
    getB(S::Matrix{Union{Missing, T}}, N::Int, M::Int) where {T<:Real}

Column vector of sum of median peptide ratios. See the following references. 

* Yu F, Haynes SE, Nesvizhskii AI. IonQuant Enables Accurate and Sensitive Label-Free Quantification With FDR-Controlled Match-Between-Runs. Mol Cell Proteomics. 2021;20:100077. doi: 10.1016/j.mcpro.2021.100077. Epub 2021 Apr 2. PMID: 33813065; PMCID: PMC8131922.
* https://rdrr.io/cran/iq/man/maxLFQ.html
* Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep;13(9):2513-26. doi: 10.1074/mcp.M113.031591. Epub 2014 Jun 17. PMID: 24942700; PMCID: PMC4159666.

"""
function getB(S::Matrix{Union{Missing, T}}, N::Int, M::Int) where {T<:Real}
    B = zeros(Float32, N + 1)
    for i in 1:N
        for j in 1:N
            if j != i
                #Ratios of peptides commonly identified bewteen experiment i and j. 
                r_i_j = skipmissing(-log2.( @view(S[:,j]) ) .+ log2.(@view(S[:,i])))
                #r_i_j could be emptry if no peptides were commonly identified between experiment i and j. 
                if !isempty(r_i_j)
                    #median(r_i_j) is the median of all ratios between peptides commonly quantified between experiment
                    #i and experiment j. 
                    B[i] += median(r_i_j)
                end
            end
        end
    end
    B.*=2
    norm = 0.0

    #Relative abundances are correct without this step, but the absolute values should make sense. 
    for row in eachrow(S)
        norm += sum(log2.(skipmissing(row)))*(length(row)/(length(row) - sum(ismissing.(row))))
    end
    B[end] = norm/M #For scaling
    B
end

"""
    getA(N::Int)

Design matrix for maxLFQ. See the following

    * Yu F, Haynes SE, Nesvizhskii AI. IonQuant Enables Accurate and Sensitive Label-Free Quantification With FDR-Controlled Match-Between-Runs. Mol Cell Proteomics. 2021;20:100077. doi: 10.1016/j.mcpro.2021.100077. Epub 2021 Apr 2. PMID: 33813065; PMCID: PMC8131922.
    * https://rdrr.io/cran/iq/man/maxLFQ.html
    * Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep;13(9):2513-26. doi: 10.1074/mcp.M113.031591. Epub 2014 Jun 17. PMID: 24942700; PMCID: PMC4159666.
    
"""
function getA(N::Int)
    A = ones(Float32, N+1, N+1)
    for i in 1:(N)
        for j in 1:N
            if i == j
                A[i, j] = 2*(N - 1)
            else
                A[i, j] = -2
            end
        end
    end
    A[end, end] = 0

    A
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
                            experiments::AbstractVector{UInt16}, 
                            use_for_quant::AbstractVector{Bool},
                            abundance::AbstractVector{Union{T, Missing}},
                            global_qvals::AbstractVector{Union{F,Missing}},
                            qvals::AbstractVector{Union{F,Missing}},
                            peps::AbstractVector{Union{F,Missing}},
                            pg_scores::AbstractVector{Union{F,Missing}},
                            global_pg_scores::AbstractVector{Union{F,Missing}},
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
                            global_pg_score_out::Vector{Union{Missing, Float32}}) where {T <: Real, F <: Real}

    unique_experiments = unique(experiments)
    unique_peptides = unique(peptides)

    N = length(unique_experiments)
    M = length(unique_peptides)

    peptides_dict = Dict(zip(unique_peptides, 1:M))
    experiments_dict = Dict(zip(unique_experiments, 1:N))

    #Appends the results to the inputs `protein_out`, `peptides_out`, `experiments_out` and `log2_abundance_out`
    function appendResults!(row_idx::Int64,
                            N::Int64,
                            target::Bool,
                            entrap_id::UInt8,
                            species::String,
                            protein::String,
                            unique_peptides::Vector{UInt32}, 
                            log2_abundances::Vector{T},
                            global_qvals::AbstractVector{Union{Float32,Missing}},
                            qvals::AbstractVector{Union{Float32,Missing}},
                            peps::AbstractVector{Union{Float32,Missing}},
                            pg_scores::AbstractVector{Union{Float32,Missing}},
                            global_pg_scores::AbstractVector{Union{Float32,Missing}},
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
                            experiments_out::Vector{Union{Missing, I}}, 
                            S::Matrix{Union{Missing, T}}) where {T<:Real,I<:Integer}
        
        function appendPeptides!(peptides_out::Vector{Union{Missing, Vector{Union{Missing, UInt32}}}}, 
                                row_idx::Int64,
                                unique_peptides::Vector{UInt32}, 
                                S::Matrix{Union{Missing,T}})
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


        # Map experiment → metrics using the first occurrence for each experiment
        pep_map = Dict{UInt32, Float32}()
        score_map = Dict{UInt32, Float32}()
        global_score_map = Dict{UInt32, Float32}()
        qval_map = Dict{UInt32, Float32}()
        global_qval_map = Dict{UInt32, Float32}()
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
        end
        appendPeptides!(peptides_out, 
                        row_idx,
                        unique_peptides, 
                        S)
    end

    #ixj matrix where rows are for experiments and columns are for peptides. Each entry is the abundance of the peptide
    #in the given experiment, or missing if peptide j was not seen in experiment i. 
    S = allowmissing(getS(peptides, peptides_dict, experiments, experiments_dict, abundance, M, N))

    #Column vector. The response matrix in Ax=B
    B = getB(S, N, M)

    #Design matrix. See references in docstring. 
    A = getA(N)

    #Solve linear system to get log-2 abundances 
    log2_abundances = (A\B)[1:(end - 1)]
    
    # Debug: Check if all abundances are missing/NaN/Inf
    n_valid_abundances = sum(!ismissing(x) && isfinite(x) for x in log2_abundances)
    if n_valid_abundances == 0
        @debug "MaxLFQ produced all invalid abundances for protein" protein=protein n_peptides=M n_experiments=N
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
                   experiments_out,
                   S)

end

# FileReference-based implementation with TransformPipeline preprocessing
function LFQ(prot_ref,  # PSMFileReference - using Any to avoid dependency issues
             protein_quant_path::String,
            quant_col::Symbol,
            file_id_to_parsed_name::Vector{String},
            q_value_threshold::Float32;
            batch_size = 100000)
    
    # Always use lazy DataFrame loading (memory efficient)
    prot = DataFrame(Arrow.Table(file_path(prot_ref)))
    
    # Log initial data state
    @info "LFQ input data" total_rows=nrow(prot) unique_proteins=length(unique(prot.inferred_protein_group))
    
    # Check for missing values in key columns
    n_missing_pg_qval = sum(ismissing.(prot.pg_qval))
    n_missing_global_qval = sum(ismissing.(prot.qlobal_pg_qval))
    n_not_for_quant = sum(.!prot.use_for_protein_quant)
    @info "Missing values check" missing_pg_qval=n_missing_pg_qval missing_global_qval=n_missing_global_qval not_for_quant=n_not_for_quant
    
    # Create pipeline operations for batch-wise application
    preprocessing_pipeline = TransformPipeline() |>
        filter_by_multiple_thresholds([
            (:pg_qval, q_value_threshold),
            (:qlobal_pg_qval, q_value_threshold)
        ]) |>
        filter_rows(row -> row.use_for_protein_quant; desc="filter_for_protein_quant")
    
    # Main processing logic (inlined from original LFQ function)
    batch_start_idx, batch_end_idx = 1, min(batch_size, size(prot, 1))
    n_writes = 0
    is_prot_sorted = issorted(prot, :inferred_protein_group, rev = true)
    @info "Is prot sorted? $is_prot_sorted"

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
        
        # Log batch status after filtering
        if nrow(subdf) == 0
            @warn "Batch filtered to 0 rows" batch_start=batch_start_idx batch_end=batch_end_idx
            continue  # Skip empty batches
        end
        
        # Continue with original LFQ logic
        #Exclude precursors with mods that impact quantitation
        #filter!(x->!occursin("M,Unimod:35", coalesce(x.structural_mods, "")), subdf)
        gpsms = groupby(
            subdf,
            [:target, :entrapment_group_id, :species, :inferred_protein_group]
        )
        ngroups = length(gpsms)
        nfiles = length(unique(prot[!,:ms_file_idx]))
        nrows = nfiles*ngroups
        
        # Pre-allocate the batch
        out = Dict(
            :target => Vector{Union{Missing, Bool}}(undef, nrows),
            :entrap_id => Vector{Union{Missing, UInt8}}(undef, nrows),
            :species => Vector{Union{Missing, String}}(undef, nrows),
            :protein => Vector{Union{Missing, String}}(undef, nrows),
            :peptides => Vector{Union{Missing, Vector{Union{Missing, UInt32}}}}(undef, nrows),
            :log2_abundance => zeros(Union{Missing, Float32}, nrows),
            :experiments => zeros(Union{Missing, UInt32}, nrows),
            :global_qval => zeros(Union{Missing, Float32}, nrows),
            :qval => zeros(Union{Missing, Float32}, nrows),
            :pg_pep => zeros(Union{Missing, Float32}, nrows),
            :pg_score => zeros(Union{Missing, Float32}, nrows),
            :global_pg_score => zeros(Union{Missing, Float32}, nrows),
        )
        for i in range(1, nrows)
            out[:target][i] = missing
            out[:entrap_id][i] = missing
            out[:species][i] = missing
            out[:protein][i] = missing
            out[:peptides][i] = missing
            out[:experiments][i] = missing
            out[:log2_abundance][i] = missing
            out[:global_qval][i] = missing
            out[:qval][i] = missing
            out[:pg_pep][i] = missing
            out[:pg_score][i] = missing
            out[:global_pg_score][i] = missing
        end

        for (group_idx, (protein, data)) in enumerate(pairs(gpsms))
            getProtAbundance(protein[:inferred_protein_group], 
                                (group_idx*nfiles) - nfiles + 1,
                                protein[:target],
                                protein[:entrapment_group_id],
                                protein[:species],
                                data[!,:precursor_idx], 
                                data[!,:ms_file_idx], 
                                data[!,:use_for_protein_quant],
                                data[!,quant_col],
                                data[!,:qlobal_pg_qval],
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
                                out[:global_pg_score]
                            )
        end
        out = DataFrame(out)
        out[!,:n_peptides] = countPeptides(out[!,:peptides]);
        out[!,:file_name] = Vector{Union{Missing, String}}(undef, size(out, 1))
        for i in range(1, size(out, 1))
            if ismissing(out[i,:experiments])
                out[i,:file_name] = missing
            else
                out[i,:file_name] = file_id_to_parsed_name[out[i,:experiments]]
            end
        end
        filter!(x->(!ismissing(x.n_peptides)), out);#&(x.n_peptides>=min_peptides), out);
        out[!,:abundance] = exp2.(out[!,:log2_abundance])
        
        # Write results
        if iszero(n_writes)
            if isfile(protein_quant_path)
                rm(protein_quant_path)
            end
            open(protein_quant_path, "w") do io
                Arrow.write(io, 
                select!(
                    out,
                    [:file_name,
                    :target,
                    :entrap_id,
                    :species,:protein,:peptides,:n_peptides,:global_qval,:qval,:pg_pep,:pg_score,:global_pg_score,:abundance]);
                file=false)  # file=false creates stream format
            end
        else
            Arrow.append(
                protein_quant_path,
                select!(
                    out,
                    [:file_name,
                    :target,
                    :entrap_id,
                    :species,:protein,:peptides,:n_peptides,:global_qval,:qval,:pg_pep,:pg_score,:global_pg_score,:abundance])
            )
        end
        n_writes += 1
    end
    
    # Check if we processed all rows
    if batch_start_idx <= size(prot, 1)
        @warn "Not all rows processed!" last_processed_idx=batch_end_idx total_rows=size(prot, 1) unprocessed_rows=size(prot, 1)-batch_end_idx
        throw(ErrorException("Not all rows processed! last_processed_idx=$batch_end_idx, total_rows=$(size(prot, 1)), unprocessed_rows=$(size(prot, 1)-batch_end_idx)"))
    end
    
    return nothing
end


function countPeptides(peptides::Vector{Union{Missing, Vector{Union{Missing, UInt32}}}})
    
    #Number of identified peptides for the protein group in a sample
    n_peptides = zeros(Union{Missing, UInt32}, length(peptides))

    for (i, pep) in enumerate(peptides)
        if ismissing(pep)
            n_peptides[i] = missing
            continue
        end

        n = 0
        for p in pep
            if !ismissing(p)
                n += 1
            end
        end
        n_peptides[i] = n
    end
    return n_peptides
end