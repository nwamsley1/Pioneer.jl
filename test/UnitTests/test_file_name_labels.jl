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

@testset "file names: full for tables, distinguishing for plot labels" begin
    paths = ["/d/20240101_lab_rep1.arrow", "/d/20240101_lab_rep2.arrow", "/d/20240101_lab_lowsignal.arrow"]
    full = Pioneer.parseFileNames(paths)
    @test full == ["20240101_lab_rep1", "20240101_lab_rep2", "20240101_lab_lowsignal"]
    @test Pioneer.distinguishingFileNames(full) == ["rep1", "rep2", "lowsignal"]

    # Differing token counts: nothing is stripped.
    @test Pioneer.distinguishingFileNames(["a_b_c", "a_b"]) == ["a_b_c", "a_b"]
    # Nothing distinguishes the runs: the full name stands in.
    @test Pioneer.distinguishingFileNames(["same_name", "same_name"]) == ["same_name", "same_name"]
    @test Pioneer.distinguishingFileNames(["only_one"]) == ["only_one"]
end
