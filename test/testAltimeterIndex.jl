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

precursors_arrow_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Altimeter111124_MixedSpecies_OlsenAstral_NoEntrapment_SplineTest.poin/precursors_for_chronologer.arrow"
precursors_table = Arrow.Table(precursors_arrow_path)
"MTNNDNIRTR"


length(unique(fragments_table[:precursor_idx]))
Int64(maximum(fragments_table[:precursor_idx]))
8961415