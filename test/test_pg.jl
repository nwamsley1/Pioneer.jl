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

pgtable = DataFrame(Arrow.Table("/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_test_1p/temp_data/passing_proteins/Rep1.arrow"))
psms_table = DataFrame(Arrow.Table("/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_test_1p/temp_data/passing_psms/Rep1.arrow"))

pg_passing = pgtable[(pgtable[!,:pg_qval].<=0.01).&(pgtable[!,:global_pg_qval].<=0.01).&(pgtable[!,:target].==true),:];
uniq_pgs_passing = unique(pg_passing[!,:protein_name]);

psms_pg_passing = psms_table[(psms_table[!,:pg_qval].<=0.01).&(psms_table[!,:global_qval_pg].<=0.01).&(psms_table[!,:target].==true),:];
uniq_pgs_passing_fr_psms = unique(psms_pg_passing[!,:inferred_protein_group]);

@info "Sets Equal? " length(uniq_pgs_passing) == length(uniq_pgs_passing_fr_psms) == length(intersect(uniq_pgs_passing, uniq_pgs_passing_fr_psms))

psms_pg_passing = psms_table[(psms_table[!,:pg_qval].<=0.01).&(psms_table[!,:global_qval_pg].<=0.01).&(psms_table[!,:target].==true).&(psms_table[!,:use_for_protein_quant].==true),:];
uniq_pgs_passing_fr_psms = unique(psms_pg_passing[!,:inferred_protein_group]);

@info "Sets Equal? " length(uniq_pgs_passing) == length(uniq_pgs_passing_fr_psms) == length(intersect(uniq_pgs_passing, uniq_pgs_passing_fr_psms))

psms_pg_passing = psms_table[(psms_table[!,:pg_qval].<=0.01
).&(psms_table[!,:global_qval_pg].<=0.01
).&(psms_table[!,:target].==true
).&(psms_table[!,:use_for_protein_quant].==true
).&(psms_table[!,:peak_area].>0
).&(psms_table[!,:global_qval].>0
).&(psms_table[!,:qval].>0)
,:];
uniq_pgs_passing_fr_psms = unique(psms_pg_passing[!,:inferred_protein_group]);

@info "Sets Equal? " length(uniq_pgs_passing) == length(uniq_pgs_passing_fr_psms) == length(intersect(uniq_pgs_passing, uniq_pgs_passing_fr_psms))



psms_long = DataFrame(Tables.columntable(Arrow.Table("/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_test_1p/precursors_long.arrow")));
psms_long = psms_long[psms_long[!,:file_name].=="Rep1", :];
psms_long_passing = psms_long[(psms_long[!,:pg_qval].<=0.01
).&(psms_long[!,:global_qval_pg].<=0.01
).&(psms_long[!,:target].==true
).&(psms_long[!,:use_for_protein_quant].==true
).&(psms_long[!,:peak_area].>0
).&(psms_long[!,:global_qval].>0
).&(psms_long[!,:qval].>0)
,:];

uniq_pgs_passing_fr_long_psms = unique(psms_long_passing[!,:inferred_protein_group]);
@info "Sets Equal? " length(uniq_pgs_passing) == length(uniq_pgs_passing_fr_long_psms ) == length(intersect(uniq_pgs_passing, uniq_pgs_passing_fr_long_psms ))

pgs_long = DataFrame(Tables.columntable(
    Arrow.Table("/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_test_1p/protein_groups_long.arrow")));
pgs_long = pgs_long[pgs_long[!,:file_name].=="Rep1", :];