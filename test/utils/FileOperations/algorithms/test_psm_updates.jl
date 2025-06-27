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

# NOTE: The original tests in this file were disabled because they test functionality 
# that doesn't match the actual update_psms_with_scores function. 
#
# The function is designed to add protein group scores to PSMs based on protein 
# grouping (using protein_name, target, entrapment_group_id as keys), not for 
# generic PSM score updates by precursor_idx as the original tests assumed.
#
# The correct functionality is tested in test_psm_updates_simple.jl
#
# If you need to test generic PSM score updates by precursor_idx, that would
# require a different function to be implemented.

println("PSM Updates tests disabled - see test_psm_updates_simple.jl for working tests")