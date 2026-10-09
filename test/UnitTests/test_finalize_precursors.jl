# finalize_precursor_table's index columns == the in-memory add_pair_indices! / add_entrapment_indices!.

using Test, Random, DataFrames, Arrow
using Pioneer

@testset "finalize_precursor_table indices" begin
    rng = Xoshiro(7)
    for trial in 1:20
        n = rand(rng, 1:400)
        pair_id = UInt32.(rand(rng, 1:max(1, n ÷ 2), n))               # groups of 1, 2 and more
        group = UInt8.(rand(rng, 0:2, n)); decoy = rand(rng, Bool, n)
        epair = Union{Missing, UInt32}[decoy[i] ? missing : UInt32(rand(rng, 1:max(1, n ÷ 3))) for i in 1:n]
        df = DataFrame(pair_id = pair_id, entrapment_pair_id = epair, entrapment_group_id = group, is_decoy = decoy)
        Pioneer.add_pair_indices!(df)
        Pioneer.add_entrapment_indices!(df)
        partner = Pioneer.pair_partners(pair_id)
        @test isequal(Union{Missing, Int64}[p == 0 ? missing : p for p in partner], df.partner_precursor_idx)
        etarget = Pioneer.entrapment_target_rows(epair, group, decoy)
        @test isequal(Union{Missing, UInt32}[t == 0 ? missing : t for t in etarget], df.entrapment_target_idx)
    end
    # whole table: names, order and index columns over several record batches
    dir = mktempdir(); in_path = joinpath(dir, "precursors.arrow"); out_path = joinpath(dir, "precursors_table.arrow")
    t = (accession_number = ["A", "B", "C", "D"], start_idx = [UInt32[1], UInt32[2, 3], UInt32[4], UInt32[5]],
         mods = Union{Missing, String}["", missing, "(1,M,Unimod:35)", ""], precursor_charge = UInt8[2, 2, 3, 3],
         decoy = [false, true, false, true], entrapment_group_id = UInt8[0, 0, 0, 0],
         pair_id = Union{Missing, UInt32}[1, 1, 2, 3], entrapment_pair_id = Union{Missing, UInt32}[1, missing, 2, missing],
         koina_sequence = ["A", "B", "M[UNIMOD:35]", "D"], collision_energy = fill(26f0, 4),
         isotopic_mods = Union{Missing, String}[missing, missing, missing, missing],
         isotope_mods = Union{Missing, String}[missing, missing, missing, missing])
    open(Arrow.Writer, in_path) do w
        Arrow.write(w, map(c -> c[1:2], t)); Arrow.write(w, map(c -> c[3:4], t))
    end
    @test Pioneer.finalize_precursor_table(in_path, out_path; entrapment_targets = true, isotope_mods = false) == (4, 2)
    out = Arrow.Table(out_path)
    @test collect(propertynames(out)) == [:accession_numbers, :start_idx, :structural_mods, :prec_charge, :is_decoy,
                                          :entrapment_group_id, :pair_id, :entrapment_pair_id,
                                          :partner_precursor_idx, :entrapment_target_idx]   # build-only columns dropped
    with_iso = joinpath(dir, "with_isotopes.arrow")
    Pioneer.finalize_precursor_table(in_path, with_iso; entrapment_targets = false, isotope_mods = true)
    @test collect(propertynames(Arrow.Table(with_iso))) == [:accession_numbers, :start_idx, :structural_mods, :prec_charge,
        :is_decoy, :entrapment_group_id, :pair_id, :entrapment_pair_id, :isotopic_mods, :isotope_mods, :partner_precursor_idx]
    @test isequal(collect(out.partner_precursor_idx), [2, 1, missing, missing])
    @test isequal(collect(out.entrapment_target_idx), Union{Missing, UInt32}[1, missing, 3, missing])
    @test [collect(s) for s in out.start_idx] == [UInt32[1], UInt32[2, 3], UInt32[4], UInt32[5]]
end
