

#=
file_to_condition = ["A","A","A","B","B","B"]
tout[!,:condition] = [file_to_condition[x] for x in tout[!,:experiments]]

log2fc = zeros(Union{Missing, Float32}, size(gtout, 1))
cv_a = zeros(Union{Missing, Float32}, size(gtout, 1))
cv_b = zeros(Union{Missing, Float32}, size(gtout, 1))
species = Vector{Union{Missing, String}}(undef, size(gtout, 1))
accession_numbers = Vector{Union{Missing, String}}(undef, size(gtout, 1))
for (idx, (key, value)) in ProgressBar(enumerate(pairs(gtout)))
    if size(value,1)<6
        log2fc[idx] = missing
        species[idx] = missing
        accession_numbers[idx] = missing
        continue
    end
    log2fc[idx] = mean(value[value[!,:condition].=="A",:log2_abundance]) - mean(value[value[!,:condition].=="B",:log2_abundance])
    species[idx] = key[:species]
    accession_numbers[idx] = key[:protein]
end

tpquant = DataFrame(Dict(
    :log2fc => log2fc,
    :species => species,
    :protein => accession_numbers
))

gtpquant = groupby(tpquant,:species)

qbins = LinRange(-3, 3, 100)
stephist(gtpquant[(species = "HUMAN",)][!,:log2fc], normalize = :probability, bins = qbins)
stephist!(gtpquant[(species = "YEAST",)][!,:log2fc], normalize = :probability, bins = qbins)
stephist!(gtpquant[(species = "ECOLI",)][!,:log2fc], normalize =:probability, bins = qbins)
vline!([log2(1/3), log2(1), log2(2)])
=#