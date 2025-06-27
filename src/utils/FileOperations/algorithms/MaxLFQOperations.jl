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
MaxLFQ algorithm operations and validation.

This module is reserved for future MaxLFQ-specific file operations.
Currently, MaxLFQ functionality is implemented directly in MaxLFQSearch
which uses MSData and the existing LFQ algorithm from utils/maxLFQ.jl.
"""

using Arrow, DataFrames

# This module is currently empty as the MaxLFQ wrapper functions
# were not being used in the actual search pipeline.
# MaxLFQ functionality is instead handled by:
# - MaxLFQSearch search method which uses MSData directly
# - Existing LFQ algorithm in utils/maxLFQ.jl
# - Pipeline system for data preprocessing
# - Merge operations for combining files

# Future MaxLFQ-specific file operations would go here when needed.