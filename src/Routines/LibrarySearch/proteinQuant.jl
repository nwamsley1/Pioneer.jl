#Filter certain modifications. (don't use these in protein quantitation)
sub_best_psms = best_psms[occursin.("M,Unimod:35)", best_psms[!,:structural_mods]).==false,:]
@time protein_quant = DataFrame(LFQ(sub_best_psms, :peak_area_normalized));
protein_quant[!,:n_peptides] = countPeptides(protein_quant[!,:peptides]);
filter!(x->(!ismissing(x.n_peptides))&(x.n_peptides>1), protein_quant);
protein_quant[!,:species] = [join(unique(split(species,';')),';') for species in protein_quant[!,:species]];
filter!(x->!occursin(';', x.species), protein_quant);
protein_quant[!,:q_value] = zeros(Union{Missing, Float32}, size(protein_quant, 1));
for i in ProgressBar(range(1, size(protein_quant,1)))
    key =         (
        UInt32(protein_quant[i,:experiments]),
        protein_quant[i,:protein]
    )
    if haskey(protein_to_q_value, key)
        protein_quant[i,:q_value] =     protein_to_q_value[
            key
        ] 
    else
        protein_quant[i,:q_value] = missing
    end
end
filter!(x->(!ismissing(x.q_value))&(x.q_value<=0.01),protein_quant);
filter!(x->x.target, protein_quant);
#gprotein_quant= groupby(protein_quant,[:target,:species,:protein])
