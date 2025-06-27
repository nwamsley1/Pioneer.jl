prec_results = DataFrame(
    Tables.columntable(
        Arrow.Table(
    "/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_test_1p/precursors_long.arrow"
)))

library_precursors = DataFrame(Arrow.Table(
    "/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_entrap.poin/precursors_table.arrow"
))

prec_results[!,:pair_id] = [library_precursors.pair_id[pid] for pid in prec_results[!,:precursor_idx]]

library_precursors[!, :partner_precursor_idx]

pid = UInt32(100)

random_precursors = rand(1:size(library_precursors, 1), 10)
for pid in random_precursors
    @show pid
    @show pid_partner = library_precursors[!,:partner_precursor_idx][pid]
    @show library_precursors[[pid, pid_partner], 
[:proteome_identifiers,:accession_numbers, :sequence, :base_pep_id, :pair_id, :partner_precursor_idx, :entrapment_group_id]]
end


library_precursors[library_precursors[!,:base_pep_id] .== 280365,
[:proteome_identifiers,:accession_numbers, :sequence, :base_pep_id, :pair_id, :partner_precursor_idx, :entrapment_group_id, :structural_mods, :prec_charge]]


base_peps = groupby(library_precursors, [:base_pep_id,:prec_charge,:is_decoy])

base_peps[20][!,
[:proteome_identifiers,:accession_numbers, :sequence, :base_pep_id, :pair_id, :partner_precursor_idx, :entrapment_group_id, :structural_mods, :prec_charge]]
