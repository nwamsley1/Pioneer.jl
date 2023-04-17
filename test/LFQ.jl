prot = DataFrame(Dict(
    :peptide => ["A","A","A","B","B","B","C","C","C","D","D","D"],
    :protein => append!(split(repeat("A",9), ""), ["B","B","B"]),
    :file_idx => UInt32[1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3],
    :abundance => [10, 20, 40, 1, 2, 4, 100, 200, missing, 1000, 2000, 3000],
))
r=sample(1:size(prot,1), size(prot,1),  replace=false)
@time prot = prot[r,:]
using StatsBase
out = Dict(
    :protein => String[],
    :peptides => String[],
    :log2_abundance => Float64[],
    :experiments => UInt32[],
)

for (protein, data) in pairs(groupby(prot, :protein))
    println(typeof(protein[:protein]))
    getProtAbundance(string(protein[:protein]), 
                        collect(data[!,:peptide]), 
                        collect(data[!,:file_idx]), 
                        collect(data[!,:abundance]),
                        out[:protein],
                        out[:peptides],
                        out[:experiments],
                        out[:log2_abundance]
                    )
    #println(protein[:parent])
end
quant = select(best_psms[(best_psms.isotope .== "light"),:], [:ms_file_idx, :sequence, :protein_names, :par, :isotope, :dev_ratio])

for (protein, data) in pairs(groupby(quant, :protein_names))
    getProtAbundance(string(protein[:protein_names]), 
                        collect(data[!,:sequence]), 
                        collect(data[!,:ms_file_idx]), 
                        (collect(data[!,:par])),
                        out[:protein],
                        out[:peptides],
                        out[:experiments],
                        out[:log2_abundance]
                    )
    #println(protein[:parent])
end
